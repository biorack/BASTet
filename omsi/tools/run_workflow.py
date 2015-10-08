"""
Simple helper tool to run an analysis. This is essentially just a short-cut
to the omsi/workflow/analysis_driver/omsi_cl_diver module
"""
from omsi.workflow.driver.cl_workflow_driver import cl_workflow_driver

if __name__ == "__main__":

    # Create a command-line driver and call main to run the analysis
    cl_workflow_driver(workflow_executor=None,
                       add_script_arg=True,
                       add_output_arg=True,
                       add_log_level_arg=True,
                       add_profile_arg=True,
                       add_mem_profile_arg=True).main()
