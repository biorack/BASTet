"""
Module with base classes for workflow drivers.

Workflow drivers are responsible for the creation and initialization of workflows. The execution
of workflows is then controlled by the workflow executor.
"""


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