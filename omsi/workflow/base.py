"""
Base module for defining analysis drivers
"""
from omsi.shared.log import log_helper
from omsi.shared.run_info_data import run_info_dict
import omsi.shared.mpi_helper as mpi_helper


class analysis_driver_base(object):
    """
    Base class used to drive the execution of a specific type of analysis. This a class-based
    execution, i.e, the user defines only the type of analysis and inputs but does not actually
    create the analysis.

    This is a class-based execution paradigm, i.e., we are given only the class of analysis
    and we need to handle the inputs, outputs, execution etc. of the whole analysis.
    """
    def __init__(self,
                 analysis_class):
        """
        Initialize the analysis driver

        :param analysis_class: The analysis class for which we want to execute the analysis.
            The analysis class must derive from omsi.analysis.analysis_base. May be None
            in case that we use the command-line to define the analysis class via the optional
            positional argument for the command class (i.e., set add_analysis_class_arg to True).
        :type analysis_class: omsi.analysis.base
        """
        super(analysis_driver_base, self).__init__()
        self.analysis_class = analysis_class

    def main(self):
        """
        The main function for running the analysis.
        """
        raise NotImplementedError("Child classes must implement the main function")

    def execute(self):
        """
        Same as main
        """
        self.main()


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
        if cls.DEFAULT_EXECUTOR_CLASS is None:
            from omsi.workflow.executor.greedy_executor import greedy_workflow_executor
            executor = greedy_workflow_executor(analysis_objects=analysis_objects)
        else:
            executor = cls.DEFAULT_EXECUTOR_CLASS(analysis_objects=analysis_objects)
        return executor

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
        self.analysis_tasks = analysis_task_set(analysis_objects) if analysis_objects is not None else analysis_task_set()
        self.mpi_comm = mpi_helper.get_comm_world()
        self.mpi_root = 0

    def __call__(self):
        """
        Same as main

        :return: The output of main
        """
        return self.execute()

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
        log_helper.debug(__name__, "Execute workflow")
        # 1) Record basic execution provenance information prior to running the analysis
        self.run_info.clear()
        self.run_info.record_preexecute()

        # 2) Execute the workflow
        start_time = time.time()
        re = self.main()
        execution_time = time.time() - start_time

        # 3) Record post-execution information, e.g., the execution time
        self.run_info.record_postexecute(execution_time=execution_time)
        self.run_info.clean_up()
        self.run_info = self.run_info.gather(root=self.mpi_root,
                                             comm=self.mpi_comm)
        try:
            log_helper.info(__name__, 'Execution time: ' + str(self.run_info['execution_time']) + "s")
        except (KeyError, ValueError):
            pass
        # 4) Return the result of the workflow execution
        return re

    def add_analysis(self,
                     analysis_object):
        """
        Add a given analysis to the set of object to be executed by the workflow

        Shorthand for: self.analysis_tasks.add(analysis_object)
        """
        self.analysis_tasks.add(analysis_object)

    def add_all(self):
        """
        Add all known analyses to the workflow.

        Shorthand for: self.analysis_tasks.add_all()
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
        log_helper.debug(__name__, "Clearing the workflow")
        self.analysis_tasks.clear()


class analysis_task_set(set):
    """
    Define a python set of analyses to be executed as a workflow.
    """
    @classmethod
    def from_script_files(cls,
                          script_files):
        """
        Same as from_script, but the script is read from the given files.

        :param script_files: List of strings with the paths to the script files. If only a single
            script is used, then a single string may be used as well.

        :return: An instance of workflow_base with the specification of the workflow to be executed
        """
        scripts = []
        inscripts = [script_files, ] if isinstance(script_files, basestring) else script_files
        for script_path in inscripts:
            infile = open(script_path, 'r')
            scripttext = infile.read()
            print scripttext
            scripts.append(scripttext)
            infile.close()
        return cls.from_scripts(scripts=scripts)

    @classmethod
    def from_scripts(cls,
                     scripts):
        """
        Evaluate the workflow script to extract all analyses to be run.

        NOTE: This function executes using eval(..), i.e., there are NO safeguards against malicious codes.

        :param scripts: The script with the setup of the workflow. This should only include
            the definition of analyses and their inputs.

        :return: An instance of workflow_base with the specification of the workflow to be executed
        """
        from omsi.analysis.base import analysis_base
        # Check which analyses exist already as we only want to include the analyses created by the script
        current_analyses = set([ana_obj for ana_obj in analysis_base.get_analysis_instances()])
        # Evaluate the workflows script
        if isinstance(scripts, basestring):
            exec(scripts)
        else:
            for script in scripts:
                exec(script)
        # Create the list of new analyses that have been created
        all_analyses = set([ana_obj for ana_obj in analysis_base.get_analysis_instances()])
        new_analyses = set(all_analyses - current_analyses)
        new_workflow = cls(new_analyses)
        return new_workflow

    def __init__(self,
                 analysis_objects=None):
        """
        Initialize the set of analysis tasks to be performed as part of the workflow

        :param analysis_objects: Set or list of unique analysis objects to be added.
            Duplicates will be removed.
        """
        super(analysis_task_set, self).__init__()
        self.update(analysis_objects)

    def __getitem__(self, item):
        """
        Get the analysis object with the given integer index of analysis_identifier
        :param item: Integer index of the analysis or string with the analysis identifier to be retrieved.
            NOTE: If analyses identifiers are user-defined, and if they are not unique, then the subset
            with matching identifiers will be retrieved.

        :return: Instance of omsi.analysis.base.analysis_base with the matching analysis object. In case that
            multiple analyses match the key, then a analysis_task_set object will be returned instead
        """
        if isinstance(item, int):
            return list(self)[item]
        elif isinstance(item, basestring):
            matching_analyses = []
            for ana_obj in self:
                if ana_obj.analysis_identifier == item:
                    matching_analyses.append(ana_obj)
            if len(matching_analyses) == 1:
                return matching_analyses[0]
            if len(matching_analyses) > 1:
                return analysis_task_set(matching_analyses)
        else:
            raise KeyError('Invalid key type')
        raise KeyError('Invalid key')

    def update(self, analysis_objects):
        """
        Return the set with elements added from the given set of analysis_objects.

        This is the same as set.update() but we ensure that only analysis_base objects
        are added.

        :param analysis_objects: List or set of analysis_base objects to be added to workflow set

        :raise: ValueError is raised in case that objects that are not instances of analysis_base are to be added.

        :return: self with elements added to self.
        """
        from omsi.analysis.base import analysis_base
        # Nothing to be added
        if analysis_objects is None:
            return self
        # Check that all elements to be added are instances of analysis_base
        for ana_obj in analysis_objects:
            if not isinstance(ana_obj, analysis_base):
                raise ValueError('Object is not instance of analysis_base')
        # Update the workflow set
        return super(analysis_task_set, self).update(analysis_objects)

    def add(self, analysis_object):
        """
        Add a given analysis to the set of object to be executed by the workflow

        This is the same as set.add() but we ensure that only analysis_base objects
        are added.

        :param analysis_object: Analysis object to be added to the execution.
            All dependencies of the analysis will also be executed as part of the
            execution.
        :type analysis_object: omsi.analysis.base.analysis_base

        :raises: ValueError is raised if the given analysis_object is invalid
        """
        from omsi.analysis.base import analysis_base
        if isinstance(analysis_object, analysis_base):
            if analysis_object in self:
                log_helper.debug(__name__, "Analysis already in the list of tasks")
                return
            log_helper.info(__name__, "Adding analysis object to the workflow set. " + str(analysis_object))
            super(analysis_task_set, self).add(analysis_object)
        else:
            raise ValueError('Analysis is not of type omsi.analysis.base.analysis_base')

    def add_analysis_dependencies(self):
        """
        Add the dependencies of all analyses to the workflow set in case they are missing.

        This function is recursive, step-by-step adding all dependencies of the workflow to the list of
        tasks to be executed, until no more dependencies are found.

        Usually this function is called by the workflow executor itself before running the analysis and
        should not need to be called by the user.

        :returns: Integer indicating the number of dependencies added to the list of tasks
        """
        # Add all direct dependencies
        add_objects = self.get_additional_analysis_dependencies()
        for new_analysis in add_objects:
            self.add(new_analysis)

        # Recursively add all indirect dependencies if add new objects
        added_dependencies = len(add_objects)
        if added_dependencies > 0:
            added_dependencies += self.add_analysis_dependencies()
        # Return the number of dependencies we added
        return added_dependencies

    def get_additional_analysis_dependencies(self):
        """
        Compute a set of all dependencies of the current list of analyses (excluding analyses that
        are already in the the list of tasks.

        :return: Python set of all analysis dependencies
        """
        from omsi.dataformat.omsi_file.common import omsi_file_common
        from omsi.analysis.base import analysis_base

        missing_dependencies = set()
        for analysis_obj in self:
            for dependency_param_obj in analysis_obj.get_all_dependency_data():
                dependency_analysis = dependency_param_obj['data']['omsi_object']
                if isinstance(dependency_analysis, analysis_base):
                    if dependency_analysis not in self:
                        missing_dependencies.add(dependency_analysis)
                elif isinstance(dependency_analysis, omsi_file_common):
                    pass  # Ignore dependencies on data files. We do not need to execute those
                else:
                    log_helper.warning(__name__, 'Unknown dependency object type that cannot be processed by workflow.'
                                       + str(dependency_analysis))
        return missing_dependencies

    def get_all_dependency_data(self):
        """
        Get the complete list of all direct and indirect dependencies of all analysis tasks.

        NOTE: These are only the direct dependencies as specified by the analysis itself.
        Use  get_all_dependency_data_recursive(..) to also get the indirect dependencies of
        the analysis due to dependencies of the dependencies themselves.

        :returns: List of parameter_data objects that define dependencies.

        """
        all_analyses = self.get_additional_analysis_dependencies()
        all_analyses.update(self)
        dependencies = [ana.get_all_dependency_data() for ana in all_analyses]
        return dependencies

    def get_all_parameter_data(self,
                               exclude_dependencies=False):
        """
        Get the complete list of all parameters

        :param exclude_dependencies: Boolean indicating whether we should exclude parameters
            that define dependencies from the list

        :return: List of omsi.analysis.analysis_data.parameter_data objects with the description of the parameters
        """
        parameters = [ana.get_all_parameters(exclude_dependencies) for ana in self]
        return parameters

    def add_all(self):
        """
        Add all known analyses to the workflow set

        :return: The updated analysis_task_set with all analyses added
        """
        from omsi.analysis.base import analysis_base
        for ana in analysis_base.get_analysis_instances():
            self.add(ana)
        return self











