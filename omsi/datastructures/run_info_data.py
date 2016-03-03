"""
Module with helper data structures for recording runtime provenance data
"""

import warnings
import platform
import datetime
import sys
import time
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.log import log_helper

try:
    from cProfile import Profile
    PROFILE_AVAILABLE = True
except ImportError:
    try:
        from profile import Profile
        PROFILE_AVAILABLE = True
    except ImportError:
        PROFILE_AVAILABLE = False
try:
    import pstats
except ImportError:
    PSTATS_AVAILABLE = False

try:
    import memory_profiler
    PROFILE_MEMORY_AVAILABLE = True
except ImportError:
    PROFILE_MEMORY_AVAILABLE = False

try:
    import StringIO
except ImportError:
    PROFILE_AVAILABLE = False
    PROFILE_MEMORY_AVAILABLE = False


# TODO Expand the data that is being recorded, e.g. psutil data etc.
# TODO Add options to record inputs and outputs of functions and the function itself


class run_info_dict(dict):
    """
    Simple dictionary class for collecting runtime information

    The typical use is as follows:

    >> my_run_info = run_info_dict()
    >> my_run_info(my_function)(my_parameters)

    With this, all runtime information is automatically collected in my_run_info.
    We can enable time-and-usage and memory profiling simply by calling
    enable_profile_time_and_usage(...) or  enable_profile_memory(...), respectively,
    before we run our function.

    We can also use the data structure directly and control the population ourselves,
    however, memory profiling is not supported by default in this case but we need to
    set and run the memory profiler ourselves, since memory_profiler expects that it
    can wrap the function

    """
    DEFAULT_TIME_FORMAT = '%Y-%m-%d %H:%M:%S.%f'

    def __init__(self, *args, **kwargs):
        super(run_info_dict, self).__init__(*args, **kwargs)
        self.__profile_time_and_usage = False
        self.__profile_memory = False
        self.__time_and_use_profiler = None
        self.__memory_profiler = None
        self.mpi_comm = mpi_helper.get_comm_world()
        self.mpi_root = 0
        self.gather_data = True

    def __call__(self, func):
        """

        :param func: The function to be wrapped for execution
        :return: A wrapped function for which we track the runtime information in self
        """
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            # Pre-execute recording
            self.clear()                # Clear all runtime information data and profilers
            self.record_preexecute()    # Record system provenance and pre-execution data
            start_time = time.time()    # Start the execution timer
            # Execute the function
            if not self.get_profile_memory():
                result = func(*args, **kwargs)  # Execute the function without memory profiling
            else:
                self.__memory_profiler = memory_profiler.LineProfiler()
                result = self.__memory_profiler(func)(*args, **kwargs)  # Execute the function with memory profiling
            # Post-execute recording
            execution_time = time.time() - start_time                   # Compute the execution time
            self.record_postexecute(execution_time=execution_time)      # Record post-execution data
            self.clean_up()                                             # Clean up empty data
            if self.gather_data:
                self.gather()                                           # Gather the data from all MPI ranks
            # Return the result
            return result

        # Return our wrapped function
        return wrapper

    def clear(self):
        """
        Clear the dictionary and other internal parameters

        Side Effects

            * Remove all key/value pairs from the dict
            * Set self.__time_and_use_profiler to None
            * Set self.__memory_profiler to None
            * Set self.__profile_memory to False if invalid (i.e, if set to True but memory profiling is unavailable)
            * Set self.__profile_time_and_usage to False if invalid (i.e., if set to True but profiling is unavailable)
        """
        # Make sure profiling settings are valid
        if self.get_profile_memory() and not PROFILE_MEMORY_AVAILABLE:
            self.enable_profile_time_and_usage(False)
        if self.get_profile_time_and_usage() and not PROFILE_AVAILABLE:
            self.enable_profile_memory(False)
        # Remove old profilers
        self.__time_and_use_profiler = None
        self.__memory_profiler = None
        # Clear all data from the dictionary
        return super(run_info_dict, self).clear()

    def enable_profile_memory(self, enable=True):
        """
        Enable/disable profiling of memory usage

        :param enable: boolean to enable (True) or disable (False) memory profiling

        """
        if PROFILE_MEMORY_AVAILABLE:
            if not enable and self.__profile_memory:
                log_helper.debug(__name__, "Disabled memory profiling. ",
                                 root=self.mpi_root, comm=self.mpi_comm)
            if enable and not self.__profile_memory:
                log_helper.debug(__name__, "Enabled memory profiling. ",
                                 root=self.mpi_root, comm=self.mpi_comm)
            self.__profile_memory = enable
        else:
            self.__profile_memory = False
            if enable:
                log_helper.warning(__name__, 'Profiling of memory usage not available.' +
                                   ' Missing memory_profiler or StringIO package')

    def enable_profile_time_and_usage(self, enable=True):
        """
        Enable/disable time and usage profiling

        :param enable: boolean to enable (True) or disable (False) time and usage profiling

        """
        if PROFILE_AVAILABLE:
            if not enable and self.__profile_time_and_usage:
                log_helper.debug(__name__, "Disabled time and usage profiling. ",
                                 root=self.mpi_root, comm=self.mpi_comm)
            if enable and not self.__profile_time_and_usage:
                log_helper.debug(__name__, "Enabled time and usage profiling. ",
                                 root=self.mpi_root, comm=self.mpi_comm)
            self.__profile_time_and_usage = enable
        else:
            self.__profile_time_and_usage = False
            if enable:
                log_helper.warning(__name__, 'Profiling of time and usage not available.' +
                                   ' Missing profile and/or pstats package')

    def get_profile_time_and_usage(self):
        """
        Check whether time and usage profiling is enabled

        :return: Boolean indicating whether time and usage profiling is enabled
        """
        return self.__profile_time_and_usage

    def get_profile_memory(self):
        """
        Check whether profiling of memory usage is enabled

        :return: Boolean indicating whether memory profiling is enabled
        """
        return self.__profile_memory

    def record_preexecute(self):
        """
        Record basic runtime information in this dict before the exeuction is started.


        Function used to record runtime information prior to executing the process we want to track, e.g.,
        the `execute_analysis(...)` of a standard analysis.

        The function may be overwritten in child classes to add recording of
        additional runtime information. All runtime data should be recorded in the
        main dict (i.e, self). This ensures in the case of standard analysis that
        the data is stored in the HDF5 file. Other data should be stored in separate
        variables that we may add to the object.

        When overwriting the function we should typically call super(...,self).runinfo_record_pretexecute()
        last in the custom version to ensure that the start_time is properly recorded right before
        the execution of the analysis.

        """
        log_helper.debug(__name__, 'Recording pre-execution runtime data', root=self.mpi_root, comm=self.mpi_comm)
        # Record basic runtime environment information using the platform module
        try:
            self['architecture'] = unicode(platform.architecture())
            self['java_ver'] = unicode(platform.java_ver())
            self['libc_ver'] = unicode(platform.libc_ver())
            self['linux_distribution'] = unicode(platform.linux_distribution())
            self['mac_ver'] = unicode(platform.mac_ver())
            self['machine'] = unicode(platform.machine())
            self['node'] = unicode(platform.node())
            self['platform'] = unicode(platform.platform())
            self['processor'] = unicode(platform.processor())
            self['python_branch'] = unicode(platform.python_branch())
            self['python_build'] = unicode(platform.python_build())
            self['python_compiler'] = unicode(platform.python_compiler())
            self['python_implementation'] = unicode(platform.python_implementation())
            self['python_revision'] = unicode(platform.python_revision())
            self['python_version'] = unicode(platform.python_version())
            self['release'] = unicode(platform.release())
            self['system'] = unicode(platform.system())
            self['uname'] = unicode(platform.uname())
            self['version'] = unicode(platform.version())
            self['win32_ver'] = unicode(platform.win32_ver())
        except:
            warnings.warn("WARNING: Recording of platform provenance failed: " + str(sys.exc_info()))

        # Attempt to record the svn version information
        try:
            import subprocess
            self['svn_ver'] = subprocess.check_output('svnversion').rstrip('\n')
        except ImportError:
            log_helper.warning(__name__, 'Recording of svn version not possible. subprocess not installed',
                               root=self.mpi_root, comm=self.mpi_comm)
        except:
            warnings.warn("Recording of svn version information failed: "+str(sys.exc_info()))

        # Attempt to record software library version
        try:
            import numpy as np
            self['numpy_version_full_version'] = unicode(np.version.full_version)
            self['numpy_version_release'] = unicode(np.version.release)
            self['numpy_version_git_revision'] = unicode(np.version.git_revision)
        except ImportError:
            log_helper.warning(__name__, 'Recording of numpy version not possible.',
                               root=self.mpi_root, comm=self.mpi_comm)

        # Attempt to record psutil data
        try:
            import psutil
            self['logical_cpu_count'] = unicode(psutil.cpu_count())
            self['cpu_count'] = unicode(psutil.cpu_count(logical=False))
            process = psutil.Process()
            self['open_files'] = unicode(process.open_files())
            self['memory_info_before'] = unicode(process.memory_info())
        except ImportError:
            log_helper.warning(__name__, 'psutil not installed. Recording of part of runtime information not possible',
                               root=self.mpi_root, comm=self.mpi_comm)
        except:
            warnings.warn("Recording of psutil-based runtime information failed: "+str(sys.exc_info()))

        # Record the start time for the analysis
        self['start_time'] = unicode(datetime.datetime.now())

        # Enable time and usage profiling if requested
        if self.__profile_time_and_usage:
            self.__time_and_use_profiler = Profile()
            self.__time_and_use_profiler.enable()

    def record_postexecute(self, execution_time=None):
        """
        Function used to record runtime information after the task we want to track is comleted, e.g.
        the `execute_analysis(...)` function of a standard analysis.

        The function may be overwritten in child classes to add recording of
        additional runtime information.

        When overwriting the function we should call super(...,self).runinfo_record_postexecute(execution_time)
        in the custom version to ensure that the execution and end_time are properly
        recorded.

        :param execution_time: The total time it took to execute the analysis. May be None, in which
            case the function will attempt to compute the execution time based on the start_time
            (if available) and the the current time.

        :param comm: Used for logging only. The MPI communicator to be used. Default value is None,
            in which case MPI.COMM_WORLD is used.

        """
        log_helper.debug(__name__, 'Recording post-execution runtime data', root=self.mpi_root, comm=self.mpi_comm)
        # Finalize recording of post execution provenance
        self['end_time'] = unicode(datetime.datetime.now())
        if execution_time is not None:
            self['execution_time'] = unicode(execution_time)
        elif 'start_time' in self:
            start_time = run_info_dict.string_to_time(self['start_time'])
            stop_time = run_info_dict.string_to_time(self['end_time'])
            self['execution_time'] = unicode(stop_time - start_time)    # TODO: This only gives execution time in full seconds right now
        else:
            self['execution_time'] = None
        # Attempt to record psutil data
        try:
            import psutil
            process = psutil.Process()
            self['memory_info_after'] = unicode(process.memory_info())
        except ImportError:
            log_helper.warning(__name__, 'psutil not installed. Recording of part of runtime information not possible',
                               root=self.mpi_root, comm=self.mpi_comm)
        except:
            warnings.warn("Recording of psutil-based runtime information failed: "+str(sys.exc_info()))

        # Record the time and use profiling data if possible
        if self.__time_and_use_profiler is not None:
            self.__time_and_use_profiler.disable()
            self.__time_and_use_profiler.create_stats()
            self['profile'] = unicode(self.__time_and_use_profiler.stats)
            # Save the summary statistics for the profiling data
            stats_io = StringIO.StringIO()
            profiler_stats = pstats.Stats(self.__time_and_use_profiler, stream=stats_io).sort_stats('cumulative')
            profiler_stats.print_stats()
            self['profile_stats'] = stats_io.getvalue()

        # Record the memory profiling data if possible
        if self.__memory_profiler is not None and self.get_profile_memory():
            log_helper.debug(__name__, 'Recording memory profiling data', root=self.mpi_root, comm=self.mpi_comm)
            mem_stats_io = StringIO.StringIO()
            memory_profiler.show_results(self.__memory_profiler, stream=mem_stats_io)
            self['profile_mem'] = unicode(self.__memory_profiler.code_map)
            self['profile_mem_stats'] = mem_stats_io.getvalue()

    def clean_up(self):
        """
        Clean up the runinfo object. In particular remove empty keys that
        either recorded None or recorded just an empty string.

        This function may be overwritten to also do clean-up needed
        due to additional custom runtime instrumentation.

        When overwriting this function we should call super(..., self).runinfo_clean_up()
        at the end of the function to ensure that the runinfo dictionary
        is clean, i.e., does not contain any empty entries.

        """
        log_helper.debug(__name__, 'Clean up runtime data', root=self.mpi_root, comm=self.mpi_comm)
        # Remove empty items from the run_info dict
        for ri_key, ri_value in self.items():
            try:
                if ri_value is None or len(ri_value) == 0:
                    self.pop(ri_key)
            except:
                pass

    def gather(self):
        """
        Simple helper function to gather the runtime information---that has been collected on
        multiple processes when running using MPI---on a single root process

        :return: If we have more than one processes then this function returns a
            dictionary with the same keys as usual for the run_info but the
            values are now lists with one entry per mpi processes. If we only have
            a single process, then the run_info object will be returned without
            changes. NOTE: Similar to mpi gather, the function only collects
            information on the root. All other processes will return just their
            own private runtime information.

        """
        if mpi_helper.MPI_AVAILABLE:
            if self.mpi_comm.Get_size() > 1:
                log_helper.debug(__name__, 'Gather runtime data from parallel tasks',
                                 root=self.mpi_root, comm=self.mpi_comm)
                self['mpi_rank'] = self.mpi_comm.Get_rank()
                run_data = self.mpi_comm.gather(self, self.mpi_root)
                if self.mpi_comm.Get_rank() == self.mpi_root:
                    merged_run_data = {}
                    for run_dict in run_data:
                        for key in run_dict:
                            try:
                                merged_run_data[key].append(run_dict[key])
                            except KeyError:
                                merged_run_data[key] = [run_dict[key]]
                    return merged_run_data
        return self

    def get_profile_stats_object(self, consolidate=True, stream=None):
        """
        Based on the execution profile of the execute_analysis(..) function get
        ``pstats.Stats`` object to help with the interpretation of the data.

        :param consolidate: Boolean flag indicating whether multiple stats (e.g., from multiple cores)
            should be consolidated into a single stats object. Default is True.
        :param stream: The optional stream parameter to be used fo the pstats.Stats object.

        :return: A single pstats.Stats object if consolidate is True. Otherwise the function
            returns a list of pstats.Stats objects, one per recorded statistic. None is returned
            in case that the stats objects cannot be created or no profiling data is available.
        """
        from ast import literal_eval
        if stream is None:
            import sys
            stream = sys.stdout

        if 'profile' in self:
            # Parse the profile data (which is stored as a string) or in the case of MPI we may
            # have a list of strings from each MPI processes
            if isinstance(self['profile'], list):
                # Convert the profile from each MPI process independently
                profile_data = [literal_eval(profile) for profile in self['profile']]
            else:
                # If we only have a single stat, then convert our data to a list, so that we can
                # handle the single and multiple statistics case in the same way in the remainder of this function
                profile_data = [literal_eval(self['profile']), ]

            # Create a list of profile objects that the pstats.Stats class understands
            profile_dummies = []
            for profile_i in profile_data:
                # Here we are creating for each statistic a dummy class on the fly that holds our
                # profile_data in the stats attributes and has an empty create_stats function.
                # This trick allows us to create a pstats.Stats object without having to write our
                # stats data to file or having to create a cProfile.Profile object first. Writing
                # the data to file involves overhead and is ugly and creating a profiler and
                # overwriting its stats is potentially problematic
                profile_dummies.append(type('Profile',
                                            (object,),
                                            {'stats': profile_i, 'create_stats': lambda x: None})())

            # Create the statistics object and return it
            if consolidate:
                profile_stats = pstats.Stats(*profile_dummies, stream=stream)
                return profile_stats
            else:
                profile_stats = [pstats.Stats(profile_i, stream=stream) for profile_i in profile_dummies]
                return profile_stats
        else:
            return None

    @staticmethod
    def string_to_structime(time_string, time_format=None):
        """
        Covert a time string to a time.struct_time using time.strptime

        :param time_string: String with the time, e.g, with the start time of a program.
        :param time_format: The time format to be used or None in which case run_info_dict.DEFAULT_TIME_FORMAT
            will be used.

        """
        return time.strptime(time_string,
                             time_format if time_format is not None else run_info_dict.DEFAULT_TIME_FORMAT)

    @staticmethod
    def string_to_time(time_string, time_format=None):
        """
        Convert a time string to local time object using time.mktime.

        :param time_string: String with the time, e.g, with the start time of a program.
        :param time_format: The time format to be used or None in which case run_info_dict.DEFAULT_TIME_FORMAT
            will be used.

        """
        return time.mktime(time.strptime(time_string,
                                         time_format if time_format is not None else run_info_dict.DEFAULT_TIME_FORMAT))

