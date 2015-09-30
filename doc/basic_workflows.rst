Defining and Executing Analysis Workflows
=========================================

Defining a workflow:
--------------------

Defining a workflow is simple, we just need to create the analyses we want to execute and define the input parameters. We can define input parameters prior to execution simply using standard dict-like assignement.

.. code-block:: python
    :linenos:

    from omsi.dataformat.omsi_file import *
    from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
    from omsi.analysis.multivariate_stats.omsi_nmf import omsi_nmf

    # Open a file to get some MSI data
    f = omsi_file('/Users/oruebel/Devel/openmsi-data/msidata/20120711_Brain.h5' , 'r')
    d = f.get_experiment(0).get_msidata(0)

    # Specify the analysis workflow
    # Create a global peak finding analysis
    a1 = omsi_findpeaks_global()      # Create the analysis
    a1['msidata'] = d                 # Set the input msidata
    a1['mzdata'] = d.mz               # Set the input mz data
    # Create an NMF that processes our peak cube
    a2 = omsi_nmf()                   # Create the NMF analysis
    a2['msidata'] = a1['peak_cube']   # Set the input data to the peak cube
    a2['numIter'] = 2                 # Set input to perform 2 iterations only


The above example creates a simple analysis workflow in which we compute a peak-cube from a raw MSI dataset and then compute an NMF from the peak cube. NOTE: The above code only specifies our workflow, but none of the analyses we specified are actually executed yet.

Executing a single analysis
---------------------------

To execute a single analysis, we can simply call the ``execute()`` function of our analysis. Note, the execute may raise and ``AnalysisReadyError`` in case that the inputs of the analsis are not ready. E.g.:

.. code-block:: python
    :linenos:

    a2.execute()   # Will fail with ``AnalysisReadyError``

.. code-block:: python
    :linenos:

    a1.execute()   # Will successfully execute a1

Executing a partial workflow
----------------------------

To execute a single analysis including any missing dependencies, we can simply call the ``execute_recursive()`` function. E.g.:

.. code-block:: python
    :linenos:

    a2.execute_recursive()   # Will successfully execute a1

The above will execute ``a1`` as well as ``a2`` since ``a2`` depends on ``a1``.

**NOTE:** Recursive execution will only execute other analyses that are actually needed to complete our analysis and analysis results of dependent analyses that have been executed before will be reused. E.g., if we would call ``a2.execute_recursive()`` again, then only ``a2`` would be executed again.

**NOTE:** When executing multiple dependent analyses, then the execution is typically controlled by a workflow driver py:meth:`omsi.workflow.analysis_driver.base.workflow_driver_base`. By default, ``execute_recursive(..)`` will automatically create a default driver. If we want to customize the driver to be used then we can simply assign a driver to the analysis before-hand by setting the py:var:`omsi.analysis.base.analysis_base.driver`` instance variable.

Executing all analyses
----------------------

To run all analyses that have been created---independent of whether they depend on each other or not---we can simply call :py:meth:`omsi.analysis.base.analysis_base.execute_all`.

.. code-block:: python
    :linenos:

    a1.execute_all()   # Execute all analyses

The above will execute any analysis that have not up-to-date. NOTE: In contrast to py:meth:`omsi.analysis.base.analysis_base.execute` and py:meth:`omsi.analysis.base.analysis_base.execute_recursive`, this is a class-level method and not an object-method. Again, the function uses a workflow driver, which we can customize by providing as driver as input to the function.

Explicitly executing workflows
------------------------------

To explicitly execute a subset of analyses (and all their dependencies) we can explicitly define a driver for the workflow we want to execute:

.. code-block:: python
    :linenos:

    from omsi.workflow.analysis_driver.greedy_workflow_driver import greedy_workflow_driver
    driver = greedy_workflow_driver()  # Create a driver
    driver.add_analysis(a1)            # Add one ore more analyses
    driver.add_analysis(a2)
    driver.execute()                   # Execute the workflow and its dependencies

.. code-block:: python
    :linenos:

    driver2 = greedy_workflow_driver()
    driver2.add_all()  # Add all analyses
    driver2.execute()  # Execute all analyses


