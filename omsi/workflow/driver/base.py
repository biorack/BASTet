"""
Module with base classes for workflow drivers.

Workflow drivers are responsible for the creation and initialization of workflows. The execution
of workflows is then controlled by the workflow executor.
"""


class driver_base(object):
    """
    Primitve base class for driving the execution of an object
    """
    def __init__(self):
        super(driver_base, self).__init__()

    def main(self):
        """
        The main function for running the analysis.
        """
        raise NotImplementedError("Child classes must implement the main function")

    def __call__(self):
        """
        Same as execute

        :return: The output of main
        """
        return self.execute()

    def execute(self):
        """
        Same as main
        """
        self.main()


class analysis_driver_base(driver_base):
    """
    Base class defining the minimal interface for drivers of a single analysis based on the type/class of the analysis.

    This is a class-based execution, i.e, the user defines only the type of analysis and inputs but does not actually
    create the analysis.

    Derived classes must implement the main(...) function where the analysis is created and executed.

    :ivar analysis_class: The analysis class for which we want to execute the analysis.
            The analysis class must derive from omsi.analysis.analysis_base. May be None
            in case that we use other means to set the analysis_class, e.g., via the
            command-line.
    :type analysis_class: omsi.analysis.base
    """
    def __init__(self,
                 analysis_class):
        """
        Initialize the analysis driver

        :ivar analysis_class: The analysis class for which we want to execute the analysis.
            The analysis class must derive from omsi.analysis.analysis_base. May be None
            in case that we use other means to set the analysis_class, e.g., via the
            command-line.
        :type analysis_class: omsi.analysis.base
        """
        super(analysis_driver_base, self).__init__()
        self.analysis_class = analysis_class

    def main(self):
        """
        The main function for running the analysis.
        """
        raise NotImplementedError("Child classes must implement the main function")


class workflow_driver_base(object):
    """
    Base class defining the minimal interface for drivers of complex analysis workflows.

    Workflows may be specified via scripts or given via a set of analysis objects.

    Derived classes must implement the main(...) function where the analysis is created and executed.

    :ivar workflow_executor: The executor of the workflow.
    :type workflow_executor: omsi.workflow.executor.base.workflow_executor_base

    """
    def __init__(self,
                 workflow_executor):
        """
        Initalize the workflow driver

        :param workflow_executor: The executor of the workflow we want to drive.
        """
        self.workflow_executor = workflow_executor

    def main(self):
        """
        The main function for running the analysis.
        """
        raise NotImplementedError("Child classes must implement the main function")
