"""
Base module for defining analysis drivers
"""

class omsi_driver_base(object):
    """
    Based class used to drive omsi-based analyses.
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
        super(omsi_driver_base, self).__init__()
        self.analysis_class = analysis_class

    def main(self):
        """
        The main function for running the analysis.
        """
        raise NotImplementedError("Child classes must implement the main function")