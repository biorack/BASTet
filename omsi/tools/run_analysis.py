"""
Simple helper tool to run an analysis. This is essentially just a short-cut
to the omsi/workflow/analysis_driver/omsi_cl_diver module
"""
from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver

if __name__ == "__main__":

    # Create a command-line driver and call main to run the analysis
    cl_analysis_driver(analysis_class=None, add_analysis_class_arg=True, add_output_arg=True, add_log_level_arg=True).main()
