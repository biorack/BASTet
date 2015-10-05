"""
Module used to help with the execution of complex analyses workflows
"""
from omsi.workflow.base import workflow_driver_base
from omsi.shared.run_info_data import run_info_dict
from omsi.shared.log import log_helper

class greedy_workflow_driver(workflow_driver_base):
    """
    Execute a set of analysis objects and their dependencies
    """
    def __init__(self, analysis_objects=None):
        """
        Initalize the workflow driver

        :param analysis_objects: A list of analysis objects to be executed
        """
        super(greedy_workflow_driver, self).__init__(analysis_objects)
        self.run_info = run_info_dict()

    def main(self):
        """Execute the analysis workflow"""
        if len(self.get_analyses()) == 0:
            log_helper.info(__name__, "The workflow is empty")
            return

        # Add all dependencies to the workflow
        log_helper.debug(__name__, "Executing the workflow")
        log_helper.info(__name__, "Adding all dependencies")
        self.add_analysis_dependencies()

        # Record the runtime information
        log_helper.debug(__name__, "Recording runtime information")
        self.run_info.clear()
        self.run_info.record_preexecute()

        # Execute the workflow in a greedy fashion (i.e., execute whichever analysis is ready and has not be run yet)
        log_helper.debug(__name__, "Running the analysis workflow")
        all_analyses = self.get_analyses()
        iterations = 0
        while True:
            # Run all analyses that are ready
            for analysis in all_analyses:
                if analysis.update_analysis and len(analysis.check_ready_to_execute()) == 0:
                    log_helper.debug(__name__, "Execute analysis: " + str(analysis))
                    analysis.execute()
            # Check if there is any other tasks that we need to execte now
            num_tasks = 0
            num_tasks_ready = 0
            for analysis in all_analyses:
                if analysis.update_analysis:
                    num_tasks += 1
                    if len(analysis.check_ready_to_execute()) == 0:
                        num_tasks_ready += 1
            if num_tasks == 0:
                log_helper.info(__name__, "Completed executing the workflow.")
                break
            if num_tasks > 0 and num_tasks_ready == 0:
                log_helper.warning(__name__, "Workflow could not be fully executed. " + str(num_tasks) +
                                   " remain in the queue but cannot be completed due to unresolved dependencies.")
            iterations += 1

        log_helper.log_var(__name__, iterations=iterations, level='DEBUG')

        # Record the runtime information after we are done with the workflow
        self.run_info.record_postexecute()
        self.run_info.gather()


