"""
Module containing base classes for workflow executors.

Workflow executors control the execution of workflows. The setup of workflows is often controlled
either by a workflow driver or the user.
"""

from omsi.workflow.common import analysis_task_set
from omsi.shared.run_info_data import run_info_dict
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.analysis_data import parameter_data
from omsi.shared.log import log_helper


class workflow_executor_base(object):
    """
    Base class used to execute a workflow of one or many analyses. This is an object-based execution,
    i.e., the user defines a set of analyses to be executed.

    We are given a set of existing analysis objects for which we need to coordinate the execution.

    :ivar analysis_objects: Private set of analysis objects to be executed
    :type analysis_objects: analysis_task_set

    :cvar DEFAULT_EXECUTOR_CLASS: Define the derived workflow_executor_base class to be used as default executor
        The default value is None, in which case the greedy_workflow_executor is used. This variable is used
        by the get_default_executor function to instantiate a  default workflow executor on request. Using this
        variable we can change the default executor to our own preferred executor, e.g., to change the executor
        used by the omsi.analysis.base.analysis_base functions execute_all(...) and execute_recursive(...)

    To implement a derived workflow executor, we need to implement a derived class that implements the main()
    function.

    """

    DEFAULT_EXECUTOR_CLASS = None
    """
    The default executor class to be used
    """

    @classmethod
    def from_script_files(cls,
                          script_files):
        """
        Same as from_script, but the scripts are read from the given list of files.

        NOTE: This function executes scripts using exec(..), i.e., there are NO safeguards against malicious codes.

        :param script_files: List of strings with the paths to the script files. If only a single
            script is used, then a single string may be used as well.

        :return: Instance of the current workflow executor class for running the given workflow
        """
        analysis_objects = analysis_task_set.from_script_files(script_files)
        new_executor = cls(analysis_objects)
        return new_executor

    @classmethod
    def from_scripts(cls,
                     scripts):
        """
        Create and initalize a worklflow executor of the current class type to execute the workflow defined
        in the given set of scripts.

        This function using analysis_task_set.from_scripts to evaluate the workflow scripts to extract all
        analyses to be created.

        NOTE: This function executes scripts using exec(..), i.e., there are NO safeguards against malicious codes.

        :param scripts: The script with the setup of the workflow. This should only include
            the definition of analyses and their inputs.

        :return: Instance of the current workflow executor class for running the given workflow
        """
        analysis_objects = analysis_task_set.from_scripts(scripts)
        new_executor = cls(analysis_objects)
        return new_executor

    @classmethod
    def get_default_executor(cls, analysis_objects=None):
        """
        Create an instance of the default workflow executor to be used.

        :param analysis_objects: A set or unique list of analysis objects to be executed by the workflow

        :return: Instance of the default workflow executor
        """
        default_executor_class = cls.get_default_executor_class()
        executor = default_executor_class(analysis_objects=analysis_objects)
        return executor

    @classmethod
    def get_default_executor_class(cls):
        """
        Get the default executor class

        :return: Derived class of workflow_executor_base
        """
        if cls.DEFAULT_EXECUTOR_CLASS is None:
            from omsi.workflow.executor.greedy_executor import greedy_executor
            return greedy_executor
        else:
            return cls.DEFAULT_EXECUTOR_CLASS

    def __init__(self,
                 analysis_objects=None):
        """
        Initialize the workflow executor

        :param analysis_objects: A list of analysis objects to be executed
        """
        log_helper.debug(__name__, "Creating workflow executor")
        if analysis_objects is not None:
            if not isinstance(analysis_objects, list) and not isinstance(analysis_objects, set):
                analysis_objects = [analysis_objects, ]
        log_helper.log_var(__name__, analysis_objects=analysis_objects, level='DEBUG')
        self.run_info = run_info_dict()
        self.track_runinfo = True
        self.analysis_tasks = analysis_task_set(analysis_objects) if analysis_objects is not None else analysis_task_set()
        self.mpi_comm = mpi_helper.get_comm_world()
        self.mpi_root = 0
        self.parameters = []

    def __call__(self):
        """
        Same as execute

        :return: The output of main
        """
        return self.execute()

    def __iter__(self):
        """
        Convenience iterator function used to iterate over all analysis tasks.
        Same as self.analysis_tasks.__iter__
        :return:
        """
        return self.analysis_tasks.__iter__()

    def __getitem__(self, item):
        """
        Convenience function used to access analysis objects directly.
        Same as self.analysis_tasks.__getitem__
        :param item:
        :return: Output of self.analysis_tasks.__getitem__ implemented by omsi.workflow.common
        """
        if isinstance(item, basestring):
            for param in self.parameters:
                if param['name'] == item:
                    return param.get_data_or_default()
        try:
            super(workflow_executor_base, self)[item]
        except TypeError:
            raise KeyError('Invalid parameter key')

    def __setitem__(self, key, value):
        """
        Set worflow driver parameter options directly via slicing

        Overwrite this function in child classes to implement custom setting behavior, e.g., error
        checking for valid values before setting a non-standard parameter.

        :param key: name of the parameters
        :param value: new value

        :raise: ValueError if an invalid value is given
        :raise: KeyError if an invalid key is given
        """
        # Check if we have a valid key
        param_set = False
        if isinstance(key, basestring):
            for param in self.parameters:
                if param['name'] == key:
                    log_helper.debug(__name__, 'Setting parameter ' + key, root=self.mpi_root, comm=self.mpi_comm)
                    param['data'] = value
                    param_set = True
        if not param_set:
            try:
                super(workflow_executor_base, self)[key] = value
            except TypeError:
                raise KeyError('Invalid parameter key')

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError("Invalid attribute given")

    def get_all_parameter_data(self):
        """
        Overwrite this function to describe input parameters that we want the executor to expose to the
        outside, e.g., workflow dirvers where a users needs to be able ot set parameters of the executor
        via the command line or other form of UI.

        :return: List of omsi.shared.analysis_data.parameter_data objects describing the options. Options
            are set using slicing via the __setitem__ function
        """
        return self.parameters

    def get_parameter_data(self,
                           index):
        """
        Given the index return the associated dataset to be written to the HDF5 file

        :param index : Return the index entry of the private member parameters.

        :raises: IndexError is raised when the index is out of bounds
        """
        return self.get_all_parameter_data()[index]

    def main(self):
        """
        Implement the execution of the workflow. We should always call execute(..) or __call__(..) to run
        the workflow. This function is intended to implementd the executor-specific exeuction behavior and
        must be implemented in child classes.
        """
        raise NotImplementedError("Child classes must implement the main function")

    def execute(self):
        """
        Execute the workflow. This uses the main() function to run the actual workflow.
        """
        import time
        log_helper.debug(__name__, "Execute", root=self.mpi_root, comm=self.mpi_comm)
        if self.track_runinfo:
            # 1) Record basic execution provenance information prior to running the analysis
            self.run_info.clear()
            self.run_info.record_preexecute()
            start_time = time.time()

        # 2) Execute the workflow
        re = self.main()

        # 3) Record post-execution information, e.g., the execution time
        if self.track_runinfo:
            execution_time = time.time() - start_time
            self.run_info.record_postexecute(execution_time=execution_time)
            self.run_info.clean_up()
            self.run_info = self.run_info.gather(root=self.mpi_root,
                                                 comm=self.mpi_comm)
            try:
                log_helper.debug(__name__, 'Execution time: ' + str(self.run_info['execution_time']) + "s",
                                 root=self.mpi_root, comm=self.mpi_comm)
            except (KeyError, ValueError):
                pass
        # 4) Return the result of the execution execution
        return re

    def add_analysis(self,
                     analysis_object):
        """
        Add a given analysis to the set of object to be executed by the workflow

        Shorthand for: self.analysis_tasks.add(analysis_object)
        """
        self.analysis_tasks.add(analysis_object)

    def add_analysis_from_scripts(self,
                                  script_files):
        """
        Evaluate the list of scripts and add all (i.e., zero, one, or multiple) analyses to this workflow

        NOTE: This function executes scripts using exec(..), i.e., there are NO safeguards against malicious codes.

        :param script_files: List of strings with the paths to the script files. If only a single
            script is used, then a single string may be used as well.

        """
        new_analysis_objects = analysis_task_set.from_script_files(script_files)
        if new_analysis_objects is not None and len(new_analysis_objects) > 0:
            log_helper.debug("Adding %i new analyses to the workflow from scripts" % len(new_analysis_objects),
                             root=self.mpi_root, comm=self.mpi_comm)
            self.analysis_tasks = self.analysis_tasks.union(new_analysis_objects)
        else:
            log_helper.debug("No analysis found in scripts", root=self.mpi_root, comm=self.mpi_comm)


    def add_analysis_all(self):
        """
        Add all known analyses to the workflow.

        Shorthand for: self.analysis_tasks.add_analysis_all()
        """
        self.analysis_tasks.add_all()

    def add_analysis_dependencies(self):
        """
        Add the dependencies of all analyses to the workflow in case they are missing.

        This function is recursive, step-by-step adding all dependencies of the workflow to the list of
        tasks to be executed, until no more dependencies are found.

        Usually this function is called by the workflow executor itself before running the analysis and
        should not need to be called by the user.

        :returns: Integer indicating the number of dependencies added to the list of tasks
        """
        return self.analysis_tasks.add_analysis_dependencies()

    def get_analyses(self):
        """
        Get the list of analyses to be run.

        Shorthand for: self.analysis_tasks
        """
        return self.analysis_tasks

    def get_analysis(self, index):
        """
        Get the analysis with the given index

         Shorthand for: self.analysis_tasks[index]

        :param index: Integer index of the analysis

        :return: omsi.analysis.base.analysis_base object

        :raises: IndexError in case that the index is invalid
        """
        return self.analysis_tasks[index]

    def clear(self):
        """
        Remove all analyses from the workflow.

        Shorthand for: self.analysis_tasks.clear()
        """
        log_helper.debug(__name__, "Clearing the workflow", root=self.mpi_root, comm=self.mpi_comm)
        self.analysis_tasks.clear()

    def add_parameter(self,
                      name,
                      help,
                      dtype=unicode,
                      required=False,
                      default=None,
                      choices=None,
                      data=None,
                      group=None):
        """
        Add a new parameter for the workflow executor. This function is typically used in the constructor
        of a derived analysis to specify the parameters of the workflow executor.

        :param name: The name of the parameter
        :param help: Help string describing the parameter
        :param type: Optional type. Default is string.
        :param required: Boolean indicating whether the parameter is required (True) or optional (False). Default False.
        :param default: Optional default value for the parameter. Default None.
        :param choices: Optional list of choices with allowed data values. Default None, indicating no choices set.
        :param data: The data assigned to the parameter. None by default.
        :param group: Optional group string used to organize parameters. Default None, indicating that
            parameters are automatically organized by driver class (e.g. in required and optional parameters)

        :raises: ValueError is raised if the parameter with the given name already exists.
        """
        log_helper.debug(__name__, "Adding workflow executor parameter. " + str(name),
                         root=self.mpi_root, comm=self.mpi_comm)
        for param in self.parameters:
            if param['name'] == name:
                log_helper.debug(__name__, "Parameter already exists. " + str(name), root=self.mpi_root, comm=self.mpi_comm)
                raise ValueError('A parameter with the name ' + unicode(name) + " already exists.")
        self.parameters.append(parameter_data(name=name,
                                              help=help,
                                              dtype=dtype,
                                              required=required,
                                              default=default,
                                              choices=choices,
                                              data=data,
                                              group=group))
