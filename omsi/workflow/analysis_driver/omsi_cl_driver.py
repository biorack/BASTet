"""
Module used to help with driving the execution of omsi-based analyses.
"""
import sys
import argparse
from omsi.analysis.base import analysis_base
from omsi.dataformat.omsi_file.common import omsi_file_common
from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.workflow.base import analysis_driver_base
import omsi.shared.mpi_helper as mpi_helper
import numpy as np
import os
from omsi.shared.log import log_helper


class RawDescriptionDefaultHelpArgParseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                                 argparse.RawDescriptionHelpFormatter):
    """
    Simple derived formatter class for use with argparse used by the
    omsi_cl_driver class. This formatter combines the default
    argparse.ArgumentDefaultsHelpFormatter and
    argparse.RawDescriptionHelpFormatter
    for formatting arguments and help descriptions.

    """
    pass


class omsi_cl_driver(analysis_driver_base):
    """
    Command-line analysis driver.

    :cvar output_save_arg_name: Name of the optional keyword argument for specifying the name
        and target for the analysis

    :cvar analysis_class_arg_name: Name of the optional first positional argument to be used
        to define the analysis class to be used.

    :cvar profile_arg_name: Name of the keyword argument used to enable profiling of the analysis

    :cvar profile_mem_arg_name: Name of the keyword argument used to enable profiling of memory usage of an analysis


    :ivar analysis_class: The class (subclass of analysis_base) defining the analysis to be executed
    :ivar add_analysis_class_arg: Boolean indicating whether an optional positional command line argument
        should be used to determine the analysis class (or whether the analysis class will be set explicitly)
    :ivar add_output_arg: Boolean indicating whether an optional keyword argument should be added to define
        the output target for the analysis.
    :ivar add_profile_arg: Add the optional --profile keyword argument for profiling the analysis
    :ivar profile_analysis: Boolean indicating whether we should profile the analysis
    :ivar parser: The argparse.ArgumentParser instance used for defining command-line arguments
    :ivar required_argument_group: argparse.ArgumentParser argument group used to define required command line arguments
    :ivar custom_argument_groups: Dict of custom argparse.ArgumentParser argument groups specified by the analysis
    :ivar output_target: Specification of the output target where the analysis result should be stored
    :ivar analysis_arguments: Dictionary defining the input arguments to be used for the analysis

    """
    analysis_class_arg_name = '__analysis_class'
    """The name where the positional argument for defining the analysis class will be stored."""

    output_save_arg_name = 'save'
    """Name of the key-word argument used to define"""

    profile_arg_name = 'profile'
    """Name of the keyword argument used to enable profiling of the analysis"""

    profile_mem_arg_name = 'memprofile'
    """Name of the keyword argument used to enable profiling of memory usage of an analysis"""

    log_level_arg_name = 'loglevel'
    """Name of the keyword argument used to specify the level of logging to be used"""

    def __init__(self,
                 analysis_class,
                 add_analysis_class_arg=False,
                 add_output_arg=True,
                 add_log_level_arg=True,
                 add_profile_arg=False,
                 add_mem_profile_arg=False):
        """

        :param analysis_class: The analysis class for which we want to execute the analysis.
            The analysis class must derive from omsi.analysis.analysis_base. May be None
            in case that we use the command-line to define the analysis class via the optional
            positional argument for the command class (i.e., set add_analysis_class_arg to True).
        :type analysis_class: omsi.analysis.base.analysis_base
        :param add_analysis_class_arg: Boolean indicating whether we will use the positional
            command-line argument to determine the analysis class name
        :param add_output_arg: Boolean indicating whether we should add the optional keyword
            argument for defining the output target for the analysis.
        :param add_profile_arg: Boolean indicating whether we should add the optional keyword
            argument for enabling profiling of the analysis.
        :param add_mem_profile_arg: Boolean indicating whether we should add the optional
            keyword argument for enabling memory profiling of the analysis.

        :raises: A ValueError is raised in the case of conflicting inputs, i.e., if
            i) analysis_class==None and add_analysis_class_arg=False, i.e., the analysis class is not determined or
            ii) analysis_class!=None and add_analysis_class_arg=False, i.e, the analysis class is determined
            via two separate mechanisms.

        """
        # Check if the input is valid
        if analysis_class is None and not add_analysis_class_arg:
            raise ValueError('The analysis class must be either set explicitly or determined from the command line.')
        if analysis_class is not None and add_analysis_class_arg:
            raise ValueError('Conflicting inputs: analysis_class set and add_analysis_class_arg set to True.')
        super(omsi_cl_driver, self).__init__(analysis_class)
        # self.analysis_class = analysis_class  # Initialized by the super constructor call
        self.add_analysis_class_arg = add_analysis_class_arg
        self.add_output_arg = add_output_arg
        self.add_profile_arg = add_profile_arg
        self.add_mem_profile_arg = add_mem_profile_arg
        self.add_log_level_arg = add_log_level_arg
        self.parser = None
        self.required_argument_group = None
        self.output_target = None
        self.profile_analysis = False
        self.profile_analysis_mem = False
        self. __output_target_self = None  # Path to output target created by the driver if we need to remove it
        self.analysis_arguments = {}
        self.custom_argument_groups = {}
        self.mpi_root = 0

    def get_analysis_class_from_cl(self):
        """
        Internal helper function used to get the analysis class object based on the
        analysis_class_arg_name positional argument from the command line.

        *Side effects:* The function sets ``self.analysis_class`

        :raises: ImportError in case that the analysis module cannot be loaded
        :raises: AttributeError in case that the analysis class cannot be extracted from the module
        """
        if len(sys.argv) < 2 or sys.argv[1].startswith('--'):
            raise ValueError("Missing required input argument defining the analysis to be executed missing")

        # Get the analysis class we need to operate on as the first positional argument
        # Get the name of the analysis class and format the string to remove common formatting problems.
        analysis_script = sys.argv[1].replace('/', '.')
        if analysis_script.endswith('.py'):
            analysis_script = analysis_script.rstrip('.py')
        if analysis_script.startswith('.'):
            analysis_script = analysis_script.lstrip('.')

        # Determine the name of the module and name from the string
        analysis_class_name = analysis_script.split('.')[-1]
        analysis_module_name = analysis_script.rstrip(analysis_class_name)[:-1]
        if not analysis_module_name.startswith('omsi.analysis'):
            analysis_module_name = 'omsi.analysis.' + analysis_module_name

        # Import the module
        try:
            analysis_module_object = __import__(analysis_module_name, globals(), locals(), [analysis_class_name], -1)
        except ImportError as e:
            log_helper.error(__name__, e.message)
            log_helper.error(__name__, "Could not locate module " + analysis_module_name)
            log_helper.error(__name__, "Please check the name of the module. Maybe there is a spelling error.")
            raise

        # Determine the self.analysis parameter
        try:
            self.analysis_class = getattr(analysis_module_object, analysis_class_name)
        except AttributeError as e:
            log_helper.error(__name__, e.message)
            log_helper.error(__name__, "Could not locate " + analysis_class_name + " in " + analysis_module_name)
            log_helper.error(__name__, "Please check the name of the analysis. Maybe there is a spelling error.")
            raise

    def initialize_argument_parser(self):
        """
        Internal helper function used to initialize the argument parser.
        NOTE: self.analysis_class must be set before calling this function.

        *Side effects:* The function sets ``self.parser`` and ``self.required_argument_group``

        """
        # Initialize the analysis object to collect information about the arguments
        analysis_object = None if self.analysis_class is None else self.analysis_class()

        # Setup the argument parser
        if self.analysis_class is not None:
            parser_description = "class description: \n\n" + \
                                 self.analysis_class.__doc__ + " \n\n" + \
                                 "execution description: \n\n" + \
                                 self.analysis_class.execute_analysis.__doc__
        else:
            parser_description = "NO ANALYSIS CLASS GIVEN"
        parser_epilog = "how to specify ndarray data? \n" \
                        "---------------------------- \n" +\
                        "n-dimensional arrays stored in OpenMSI data files may be specified as \n" + \
                        "input parameters via the following syntax: \n" + \
                        "      -- MSI data: <filename>.h5:/entry_#/data_# \n" + \
                        "      -- Analysis data: <filename>.h5:/entry_#/analysis_#/<dataname> \n" + \
                        "      -- Arbitrary dataset: <filename>.h5:<object_path>\n" + \
                        "E.g. a valid definition may look like: 'test_brain_convert.h5:/entry_0/data_0'\n" + \
                        "In rear cases we may need to manually define an array (e.g., a mask)\n" + \
                        "Here we can use standard python syntax, e.g, '[1,2,3,4]' or '[[1, 3], [4, 5]]' \n" + \
                        "\n\n" + \
                        "This command-line tool has been auto-generated using BASTet (Berkeley Analysis & Storage Toolkit)"

        self.parser = argparse.ArgumentParser(description=parser_description,
                                              epilog=parser_epilog,
                                              formatter_class=RawDescriptionDefaultHelpArgParseFormatter,
                                              add_help=False)  # We'll add the help later in

        # Create the argument group for required arguments
        self.required_argument_group = self.parser.add_argument_group(title="required analysis arguments")

        # Create custom argument groups from the analysis
        if analysis_object is not None:
            arg_group_dict = {arg_param.get_group_name(): arg_param.get_group_description()
                              for arg_param in analysis_object.get_all_parameter_data()
                              if arg_param.get_group_name() is not None}
            for group_name, group_description in arg_group_dict.iteritems():
                self.custom_argument_groups[group_name] = self.parser.add_argument_group(title=group_name,
                                                                                         description=group_description)

        # Create the optional positional argument for the analysis class
        if self.add_analysis_class_arg:
            self.parser.add_argument(self.analysis_class_arg_name,
                                     help='The omsi.analysis class.')

        # Add the argument for defining where we should save the analysis
        if self.add_output_arg:
            output_arg_help = 'Define the file and experiment where the analysis should be stored. ' + \
                              'A new file will be created if the given file does not exists ' + \
                              'but the directory does. The filename is expected to be of the from:  ' + \
                              '<filename>:<entry_#> . If no experiment index is given, then ' + \
                              'experiment index 0 (i.e, entry_0) will be assumed by default. A valid' + \
                              'path may, e.g, be "test.h5:/entry_0" or jus "test.h5"'
            self.parser.add_argument("--"+self.output_save_arg_name,
                                     action='store',
                                     default=None,
                                     type=str,
                                     required=False,
                                     help=output_arg_help)

        # Add the optional keyword argument for enabling profiling of the analysis
        if self.add_profile_arg:
            profile_arg_help = 'Enable runtime profiling of the analysis. NOTE: This is intended for ' + \
                               'debugging and investigation of the runtime behavior of an analysis.' + \
                               'Enabling profiling entails certain overheads in performance'
            self.parser.add_argument("--"+self.profile_arg_name,
                                     action='store_true',
                                     default=False,
                                     required=False,
                                     help=profile_arg_help)
        # Add the optional keyword argument for enabling memory profiling of the analysis
        if self.add_mem_profile_arg:
            profile_mem_arg_help = 'Enable runtime profiling of the memory usage of analysis. ' + \
                                   'NOTE: This is intended for debugging and investigation of ' + \
                                   'the runtime behavior of an analysis. Enabling profiling ' + \
                                   'entails certain overheads in performance.'
            self.parser.add_argument("--"+self.profile_mem_arg_name,
                                     action='store_true',
                                     default=False,
                                     required=False,
                                     help=profile_mem_arg_help)
        # Add the optional logging argument
        if self.add_log_level_arg:
            self.parser.add_argument("--"+self.log_level_arg_name,
                                     action='store',
                                     default='INFO',
                                     required=False,
                                     help='Specify the level of logging to be used.',
                                     choices=log_helper.log_levels.keys())


    def parse_cl_arguments(self):
        """
        The function assumes that the command line parser has been setup using the initialize_argument_parser(..)

        This function parses all arguments that are specific to the command-line parser itself. Analysis
        arguments are added and parsed later by the add_and_parse_analysis_arguments(...) function.
        The reason for this is two-fold: i) to separate the parsing of analysis arguments and arguments of the
        command-line driver and ii) if the same HDF5 file is used as input and output target, then we need to
        open it first here in append mode before it gets opened in read mode later by the arguments.

        *Side effects:* The function sets ``self.output_target`` and ``self.profile_analysis``

        """
        # Parse the arguments and convert them to a dict using vars
        parsed_arguments = vars(self.parser.parse_known_args()[0])

        # Clean up the arguments to remove default arguments of the driver class
        # before we hand the arguments to the analysis class
        if self.analysis_class_arg_name in parsed_arguments:
            parsed_arguments.pop(self.analysis_class_arg_name)

        # Process the --save argument to determine where we should save the output
        if self.output_save_arg_name in parsed_arguments and  mpi_helper.get_rank() == self.mpi_root:
            # Determine the filename and experiment group from the path
            self.output_target = parsed_arguments.pop(self.output_save_arg_name)
            if self.output_target is not None:
                output_filename, output_object_path = omsi_file_common.parse_path_string(self.output_target)
                # Create the output file
                if output_filename is None:
                    raise ValueError("ERROR: Invalid save parameter specification " + self.output_target)
                elif os.path.exists(output_filename) and not os.path.isfile(output_filename):
                    raise ValueError("ERROR: Save parameter not specify a file.")
                if not os.path.exists(output_filename):
                    out_file = omsi_file(output_filename, mode='a')
                    self.output_target = out_file.create_experiment()
                    self. __output_target_self = output_filename
                else:
                    out_file = omsi_file(output_filename, mode='r+')
                    if output_object_path is not None:
                        self.output_target = omsi_file_common.get_omsi_object(out_file[output_object_path])
                    else:
                        if out_file.get_num_experiments() > 0:
                            self.output_target = out_file.get_experiment(0)
                        else:
                            self.output_target = out_file.create_experiment()
        else:
            self.output_target = parsed_arguments.pop(self.output_save_arg_name)

        # Process the --profile profiling argument
        if self.profile_arg_name in parsed_arguments:
            self.profile_analysis = parsed_arguments.pop(self.profile_arg_name)

        # Process the --memprofile argument
        if self.profile_mem_arg_name in parsed_arguments:
            self.profile_analysis_mem = parsed_arguments.pop(self.profile_mem_arg_name)
        # The --loglovel argument
        if self.log_level_arg_name in parsed_arguments:
            user_log_level = parsed_arguments.pop(self.log_level_arg_name)
            if user_log_level in log_helper.log_levels.keys():
                log_helper.set_log_level(level=log_helper.log_levels[user_log_level])
            else:
                log_helper.error(module_name=__name__, message="Invalid log level specified")

    def add_and_parse_analysis_arguments(self):
        """
        The function assumes that the command line parser has been setup using the initialize_argument_parser(..)

        This function is responsible for adding all command line arguments that are specific to the analysis and
        to then parse those argument and save the relevant data in the self.analysis_arguments dictionary.
        Command-line arguments that are specific to the command line driver are removed, so that only
        arguments that can be consumed by the analysis are handed to the analysis.

        *Side effects:* The function sets ``self.analysis_arguments``

        """
        # Add all arguments from our analysis object
        analysis_object = None if self.analysis_class is None else self.analysis_class()
        if analysis_object is not None:
            for arg_param in analysis_object.get_all_parameter_data():
                arg_name = "--" + arg_param['name']
                arg_action = 'store'
                arg_default = arg_param['default']
                arg_required = arg_param['required'] and arg_default is None
                arg_type = arg_param['dtype']
                arg_choices = arg_param['choices']
                arg_help = arg_param['help']
                arg_dest = arg_param['name']
                arg_group = arg_param.get_group_name()

                # Determine the group the argument belongs to
                argument_group = self.required_argument_group if arg_required else self.parser
                if arg_group in self.custom_argument_groups:
                    argument_group = self.custom_argument_groups[arg_group]

                # Add the argument to the proper group
                argument_group.add_argument(arg_name,               # <-- Required, user specified arg name
                                            action=arg_action,      #     Constant. We define this not the user.
                                            # nargs=1,                    Don't use. Leave as default
                                            # const=None,                 Don't use. We don't use this type of action
                                            default=arg_default,    # <-- Optional default value for the argument
                                            type=arg_type,          # <-- Optional dtype of the argument
                                            choices=arg_choices,    # <-- Optional Key may be missing.
                                            required=arg_required,  # <-- Optional
                                            help=arg_help,          # <-- Required
                                            # metavar               #     Don't use. Positional analysis arguments
                                            #                       #     are not allowed
                                            dest=arg_dest)          #     Automatically determined by the name
        # Add the help argument
        self.parser.add_argument('-h', '--help',
                                action='help',
                                default=argparse.SUPPRESS,
                                help='show this help message and exit')
        parsed_arguments = vars(self.parser.parse_args())
        parsed_arguments.pop(self.analysis_class_arg_name, None)
        parsed_arguments.pop(self.profile_arg_name, None)
        parsed_arguments.pop(self.output_save_arg_name, None)
        parsed_arguments.pop(self.profile_mem_arg_name, None)
        parsed_arguments.pop(self.log_level_arg_name, None)
        self.analysis_arguments = parsed_arguments

    def print_settings(self):
        """
        Print the analysis settings.
        """
        log_helper.info(__name__, "Inputs:")
        for key, value in self.analysis_arguments.iteritems():
            log_helper.info(__name__, "   " + unicode(key) + " = " + unicode(value))
        if self.output_target is not None:
            if isinstance(self.output_target, omsi_file_common):
                h5py_object = omsi_file_common.get_h5py_object(self.output_target)
                log_helper.info(__name__, "Save to: " + unicode(h5py_object.file.filename) + u":" + unicode(h5py_object.name))
            else:
                log_helper.info(__name__, "Save to: " + unicode(self.output_target))

    def remove_output_target(self):
        """
        This function is used to delete any output target files created by the
        command line driver. This is done in case that an error occurred and
        we do not want to leave garbage files left over.

        *Side effects* The function modifies ``self.output_target``

        :return: Boolean indicating whether we succesfully cleaned up the output
        """
        success = False
        if self.__output_target_self is not None:
            try:
                os.remove(self.__output_target_self)
                log_helper.info(__name__, "Successfully removed output target: " + unicode(self.__output_target_self))
                success = True
            except:
                log_helper.error(__name__, "Clean-up of output failed. File may be left on system: "
                                 + unicode(self.__output_target_self))
        elif self.output_target is not None:
            log_helper.info(__name__, "Output target not removed because it was not created " +
                                      "by the analysis but potentially modified by it")
        else:
            success = True
        return success

    def main(self):
        """
        Default main function for running an analysis from the command line.
        The default implementation exposes all specified analysis parameters as command
        line options to the user. The default implementation also provides means to
        print a help text for the function.

        :raises: ValueError is raised in case that the analysis class is unknown
        """
        # Get the analysis object if needed
        if self.add_analysis_class_arg:
            try:
                self.get_analysis_class_from_cl()
            except (ImportError, AttributeError, ValueError):
                pass

        # Initialize the argument parser
        if self.parser is None:
            self.initialize_argument_parser()

        # Check if we have a valid analysis class
        if self.analysis_class is None:
            print self.parser.print_help()
            raise ValueError('Could not determine the analysis class.')
        if not issubclass(self.analysis_class, analysis_base):
            print self.parser.print_help()
            raise ValueError('Analysis class is not a subclass of analysis_base.')

        try:
            # Parse the command line arguments to determine the command line driver settings
            self.parse_cl_arguments()
            # Add and parse the command line arguments specific to the analysis to determine the analysis settings
            self.add_and_parse_analysis_arguments()
        except:
            self.remove_output_target()
            raise

        # Print the analysis settings
        if mpi_helper.get_rank() == self.mpi_root:
            self.print_settings()

        # Call the execute function of the analysis
        try:
            # Create the analysis object
            analysis_object = self.analysis_class()
            # Enable time and usage profiling if requested
            if self.profile_analysis:
                try:
                    analysis_object.enable_time_and_usage_profiling(self.profile_analysis)
                except ImportError as e:
                    log_helper.warning(__name__, "Profiling of time and usage not available due to missing packages.")
                    log_helper.warning(__name__, e.message)
            # Enable memory profiling if requested
            if self.profile_analysis_mem:
                try:
                    analysis_object.enable_memory_profiling(self.profile_analysis_mem)
                except ImportError as e:
                    log_helper.warning(__name__, "Profiling of memory usage not available due to missing packages")
                    log_helper.warning(__name__, e.message)
            # Execute the analysis
            log_helper.debug(__name__, 'Analysis arguments: ' + str(self.analysis_arguments))
            analysis_object.execute(**self.analysis_arguments)
        except:
            if mpi_helper.get_rank() == self.mpi_root:
                self.remove_output_target()
            raise

        # Finalize the saving of results on rank our mpi root rank. NOTE: When running in serial
        # the condition of  mpi_helper.get_rank() ==  self.mpi_root evaluates to True because
        # our mpi_root is 0 and the mpi_helper returns 0 for the rank when running in serial.
        if mpi_helper.get_rank() == self.mpi_root:
            # Print the profiling results of time and usage
            if self.profile_analysis:
                print ""
                print "PROFILING DATA: TIME AND USAGE"
                print ""
                analysis_object.get_profile_stats_object(consolidate=True).print_stats()

            # Print the profiling results for memory usage
            if self.profile_analysis_mem:
                print ""
                print "PROFILING DATA: MEMORY"
                print ""
                print analysis_object.get_memory_profile_info()

            # Print the time it took to run the analysis
            try:
                # Parallel case: We need to compile/collect timing data from all cores
                if isinstance(analysis_object.run_info['execution_time'] , list):
                    # Time for each task to execute
                    log_helper.info(__name__, "Time in seconds for each analysis process: " +
                                     str(analysis_object.run_info['execution_time']))
                    # Start times of each task
                    log_helper.info(__name__, "Time when each of the processes started: " +
                                              str(analysis_object.run_info['start_time']))
                    # Stop times for each task

                    log_helper.info(__name__, "Time when each of the processes finished: " +
                                              str(analysis_object.run_info['end_time']))

                    # Compile the time to execute string
                    exec_time_array = np.asarray(analysis_object.run_info['execution_time'], dtype=float)
                    max_exec_time = str(exec_time_array.max())
                    min_exec_time = str(exec_time_array.min())
                    mean_exec_time = str(exec_time_array.mean())
                    exec_time_string = max_exec_time + " s " + \
                                       "    ( min = " + min_exec_time + " , mean = " + mean_exec_time + " )"
                # Serial case: We only have a single time to worry about
                else:
                    exec_time_string = str(analysis_object.run_info['execution_time']) + " s"
                log_helper.info(__name__, "Time to execute analysis: " + exec_time_string)
            except:
                raise

        # Save the analysis to file
        if self.output_target is not None:
            from omsi.dataformat.omsi_file.analysis import omsi_analysis_manager
            omsi_analysis_manager.create_analysis_static(analysis_parent=self.output_target,
                                                         analysis=analysis_object)

            #if mpi_helper.get_rank() == self.mpi_root:
            #    self.output_target.create_analysis(analysis_object)
            #else:
            #    try:
            #        analysis_object.write_analysis_data()
            #    except NotImplementedError:
            #        pass
            # TODO we should compute the minimum and maximum start time and compute the total runtime that way as well
            # TODO add MPI Barrier at the beginning to make sure everyone has started up before we do anything


if __name__ == "__main__":

    # Create an command-line driver and call the main function to run the analysis
    omsi_cl_driver(analysis_class=None,
                   add_analysis_class_arg=True,
                   add_output_arg=True,
                   add_profile_arg=True,
                   add_mem_profile_arg=True).main()
