"""
Base module for defining analysis drivers
"""
from omsi.shared.log import log_helper

class analysis_driver_base(object):
    """
    Base class used to drive the execution of a specific type of analysis.

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


class workflow_driver_base(object):
    """
    Base class used to drive a workflow of existing analyses.

    We are given a set of existing analysis objects for which we need to coordinate the execution.

    :ivar __analysis_objects: Private set of analysis objects to be executed

    """

    DEFAULT_DRIVER_CLASS = None

    @classmethod
    def get_default_driver(cls, analysis_objects=None):
        """
        Create an instance of the default workflow driver to be used.

        :param analysis_objects: A list of analysis objects to be executed.

        :return: Instance of the default workflow driver
        """
        if cls.DEFAULT_DRIVER_CLASS is None:
            from omsi.workflow.analysis_driver.greedy_workflow_driver import greedy_workflow_driver
            driver = greedy_workflow_driver(analysis_objects=analysis_objects)
        else:
            driver = cls.DEFAULT_DRIVER_CLASS(analysis_objects=analysis_objects)
        return driver


    def __init__(self,
                 analysis_objects=None):
        """
        Initalize the workflow driver

        :param analysis_objects: A list of analysis objects to be executed
        """
        log_helper.debug(__name__, "Creating workflow driver")
        if analysis_objects is not None:
            if not isinstance(analysis_objects, list) and not isinstance(analysis_objects, set):
                analysis_objects = [analysis_objects, ]
        self.__analyses_objects = set(analysis_objects) if analysis_objects is not None else set()


    def main(self):
        """
        The main function for running the analysis workflow
        """
        raise NotImplementedError("Child classes must implement the main function")

    def execute(self):
        """
        Same as main
        """
        self.main()

    def add_analysis(self,
                     analysis_object):
        """
        Add a given analysis to the set of object to be executed by the workflow

        :param analysis_object: Analysis object to be added to the execution.
            All dependencies of the analysis will also be executed as part of the
            execution.
        :type analysis_object: omsi.analysis.base.analysis_base

        :raises: ValueError is raised if the given analysis_object is invalid
        """
        from omsi.analysis.base import analysis_base

        if isinstance(analysis_object, analysis_base):
            if analysis_object in self.__analyses_objects:
                log_helper.debug(__name__, "Analysis already in the list of tasks")
                return
            log_helper.info(__name__, "Adding analysis object to the workflow. " + str(analysis_object))
            self.__analyses_objects.add(analysis_object)
        else:
            raise ValueError('Analysis is not of type omsi.analysis.base.analysis_base')

    def add_all(self):
        """
        Add all known analyses to the workflow
        """
        from omsi.analysis.base import analysis_base
        for ana in analysis_base.get_analysis_instances:
            self.add_analysis(ana)

    def add_analysis_dependencies(self):
        """
        Add the dependencies of all analyses to the workflow in case they are missing.

        This function is recursive, step-by-step adding all dependencies of the workflow to the list of
        tasks to be executed, until no more dependencies are found.

        Usually this function is called by the workflow driver itself before running the analysis and
        should not need to be called by the user.

        :returns: Integer indicating the number of dependencies added to the list of tasks
        """
        from omsi.analysis.base import analysis_base
        from omsi.dataformat.omsi_file.common import omsi_file_common

        # Add all direct dependencies
        add_objects = set()
        for analysis_obj in self.__analyses_objects:
            for dependency_param_obj in analysis_obj.get_all_dependency_data():
                dependency_analysis = dependency_param_obj['data']['omsi_object']
                if isinstance(dependency_analysis, analysis_base):
                    if dependency_analysis not in self.__analyses_objects:
                        add_objects.add(dependency_analysis)
                        log_helper.debug(__name__, "Added dependency to the workflow. " + str(dependency_analysis))
                elif isinstance(dependency_analysis, omsi_file_common):
                    pass  # Ignore dependencies on data files. We do not need to execute those
                else:
                    log_helper.warning(__name__, 'Unknown dependency object type that cannot be processed by workflow.'
                                       + str(dependency_analysis))
        for new_analysis in add_objects:
            self.add_analysis(new_analysis)

        # Recursively add all indiret dependencies if add new objects
        added_dependencies = len(add_objects)
        if added_dependencies > 0:
            added_dependencies += self.add_analysis_dependencies()
        # Return the number of dependencies we added
        return added_dependencies

    def get_analyses(self):
        """
        Get the list of analyses to be run.
        """
        return self.__analyses_objects

    def clear(self):
        """
        Remove all analyses from the workflow.
        """
        log_helper.debug(__name__, "Clearing the workflow")
        self.__analyses_objects.clear()


