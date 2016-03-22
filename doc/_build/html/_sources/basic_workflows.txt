Defining and Executing Analysis Workflows
=========================================

Figure :ref:`workflow_illustration`, illustrates the basic steps of using analysis workflows, i.e,.:

1) Create the analysis tasks
2) Define the analysis inputs
3) Execute

.. _workflow_illustration:

.. figure:: _static/workflow_illustration.png
   :scale: 100 %
   :alt: Figure showing and example workflow

   Illustration of an example workflow for image normalization


In the following we will use a simple analysis---workflow in which we compute a peak-cube from a raw MSI dataset and then compute an NMF from the peak cube---to illustrate the main steps involved for performing complex analysis workflows.


Step 1: Create the analysis tasks:
----------------------------------

First we need to create our main analysis objects.

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
    # Create an NMF that processes our peak cube
    a2['numIter'] = 2                 # Set input to perform 2 iterations only

Step 2: Define analysis inputs:
-------------------------------

We can define the input parameters of analysis simply using standard dict-like assignment. Any dependencies between analysis tasks or OpenMSI files are created automatically for us.

.. code-block:: python
    :linenos:

    # Define the inputs of the global peak finder
    a1['msidata'] = d                 # Set the input msidata
    a1['mzdata'] = d.mz               # Set the input mz data
    # Define the inputs of the NMF
    a2['msidata'] = a1['peak_cube']   # Set the input data to the peak cube
    a2['numIter'] = 2                 # Set input to perform 2 iterations only



NOTE: So far we have only specified our workflow. We have not executed any analysis yet, nor have we loaded any actual data yet.


Step 3: Execute
---------------

Finally we need to execute our analyses. For this we have various options, depending on which parts of our workflow we want to execute.


Executing a single analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To execute a single analysis, we can simply call the ``execute()`` function of our analysis. Note, the execute may raise and ``AnalysisReadyError`` in case that the inputs of the analsis are not ready. E.g.:

.. code-block:: python
    :linenos:

    a2.execute()   # Will fail with ``AnalysisReadyError``

.. code-block:: python
    :linenos:

    a1.execute()   # Will successfully execute a1

Executing a single sub-workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To execute a single analysis including any missing dependencies, we can simply call the ``execute_recursive()`` function. E.g.:

.. code-block:: python
    :linenos:

    a2.execute_recursive()   # Will successfully execute a1

The above will execute ``a1`` as well as ``a2`` since ``a2`` depends on ``a1``.

**NOTE:** Recursive execution will only execute other analyses that are actually needed to complete our analysis and analysis results of dependent analyses that have been executed before will be reused. E.g., if we would call ``a2.execute_recursive()`` again, then only ``a2`` would be executed again.

**NOTE:** When executing multiple dependent analyses, then the execution is typically controlled by a workflow executor py:meth:`omsi.workflow.executor`. By default, ``execute_recursive(..)`` will automatically create a default driver. If we want to customize the driver to be used then we can simply assign a driver to the analysis before-hand by setting the py:var:`omsi.analysis.base.analysis_base.driver`` instance variable.

Executing all analyses
^^^^^^^^^^^^^^^^^^^^^^

To run all analyses that have been created---independent of whether they depend on each other or not---we can simply call :py:meth:`omsi.analysis.base.analysis_base.execute_all`.

.. code-block:: python
    :linenos:

    a1.execute_all()   # Execute all analyses

The above will execute any analysis that have not up-to-date. NOTE: In contrast to py:meth:`omsi.analysis.base.analysis_base.execute` and py:meth:`omsi.analysis.base.analysis_base.execute_recursive`, this is a class-level method and not an object-method. Again, the function uses a workflow driver, which we can customize by providing as driver as input to the function.

Executing multiple sub-workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To explicitly execute a subset of analyses (and all their dependencies) we can explicitly define a driver for the workflow we want to execute:

.. code-block:: python
    :linenos:

    from omsi.workflow.driver.greedy_executor import greedy_executor
    driver = greedy_executor()  # Create a driver
    driver.add_analysis(a1)            # Add one ore more analyses
    driver.add_analysis(a2)
    driver.execute()                   # Execute the workflow and its dependencies

.. code-block:: python
    :linenos:

    driver2 = greedy_executor()
    driver2.add_analysis_all()  # Add all analyses
    driver2.execute()  # Execute all analyses


Example: Normalizing an image
-----------------------------

The goal of this example is to 1) illustrate the general concepts of how we can define analysis workflows and 2) illustrate the use of simple wrapped functions in combination with integrated analytics to create complex analysis workflows. The example shown below defines a basic image normalization workflow in which we:

1. Compute a reduced peak cube from an MSI image using the global peak finding analysis provided by BASTet
2. Use a simple wrapped function to compute the total intensity image for the peak cube dataset computed in step 1
3. use a simple wrapped function to normalize the peak cube computed in step 1 using the total intensity image computed in step 2

This is the same workflow as shown in Figure :ref:`workflow_illustration`.


.. code-block:: bash
    :linenos:

    # Illustration of the basic image normalization workflow defined below:
    #
    # +-----------a1------------+       +-------------a2----------------+       +-----------a3--------------+
    # +---global-peak-finder----+       +------total_intensities--------+       +---normalize_intensities---+
    # |                         |       |                               |       |                           |
    # | msidata       peak_cube +---+---> msidata      total_intensities+-------> norm_factors     output_0 |
    # |                         |   |   |                               |       |                           |
    # | mzdata                  |   |   | axis=2                        |   +---> msidata                   |
    # +-------------------------+   |   +-------------------------------+   |   +---------------------------+
    #                               |                                       |
    #                               |                                       |
    #                               +---------------------------------------+


.. code-block:: python
    :linenos:
    :emphasize-lines: 21,22,23,26,27,28,31,32,33,43

    import numpy as np
    from omsi.shared.log import log_helper
    log_helper.set_log_level('DEBUG')
    from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
    from omsi.dataformat.omsi_file.main_file import omsi_file
    from omsi.analysis.generic import analysis_generic

    # Define a simple function to compute the total intensity image
    def total_intensity(msidata, axis=2):
        import numpy as np
        return np.sum(msidata, axis=axis)

    # Define a simple function to normalize an MSI data cube by per-spectrum normalization factors
    def normalize_intensities(msidata, normfactors):
        import numpy as np
        return msidata / normfactors[:,:,np.newaxis]

    # Get an ezample MSI image
    f = omsi_file('/Users/oruebel/Devel/openmsi-data/msidata/20120711_Brain.h5' , 'r')
    d = f.get_experiment(0).get_msidata(0)

    # Define the global peak finder
    a1 = omsi_findpeaks_global()
    a1['msidata'] = d
    a1['mzdata'] = d.mz

    # Define compute of total intensity image
    a2 = analysis_generic.from_function(analysis_function=total_intensity,
                                        output_names=['total_intensities'])
    a2['msidata'] = a1['peak_cube']

    # Define the normalization of the peak cube
    a3 = analysis_generic.from_function(normalize_intensities)
    a3['msidata'] = a1['peak_cube']
    a3['normfactors'] = a2['total_intensities']

    # To run the workflow we now have several basic options
    #
    # 1) a3.execute_recursive()  : Recursively execute the last analysis and all its dependencies (i.e., a1, a2)
    # 2) a1.execute_all() : Tell any analysis to execute all available analyses (i.e., a1,a2,a3)
    # 3) Create our own workflow driver to control the execution of the analyses
    # 4) Manually call execute on a1, a2, and a3 in order of their dependencies

    # Execute the workflow
    a3.execute_recursive()


Workflow Tools
==============

Similar to the :py:mod:`omsi.workflow.driver.cl_analysis_driver` (and the corresponding tool :py:mod:`omsi.tools.run_analysis`) for running single analysis tasks, BASTet provides basic tools for executing complete workflows via the concept of workflow dirvers. Users may implement their own drivers using the approbriate base classes :py:mod:`omsi.workflow.driver.base`.

Some basic drivers and tools are already available with BASTet, e.g., the :py:mod:`omsi.workflow.driver.cl_workflow_driver` module (and the corresponding tool :py:mod:`omsi.tools.run_workflow`) defines a driver for driving and executing one or multiple workflows defined via workflow scripts, directly from the command-line.

Workflow Scripts
----------------

Workflow scripts are regular python scripts that include the i) creation of the analusis objects, and ii) full or partial definition of analysis parameters but usually **NOT** the actual execution of any of the analyses. Following our example from earlier, we may simply save the following code in python source file, e.g, `normalize_image.py`.


.. code-block:: python
    :linenos:

    import numpy as np
    from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
    from omsi.dataformat.omsi_file.main_file import omsi_file
    from omsi.analysis.generic import analysis_generic

    # Define a simple function to compute the total intensity image
    def total_intensity(msidata, axis=2):
        import numpy as np
        return np.sum(msidata, axis=axis)

    # Define a simple function to normalize an MSI data cube by per-spectrum normalization factors
    def normalize_intensities(msidata, normfactors):
        import numpy as np
        return msidata / normfactors[:,:,np.newaxis]

    # Define the global peak finder
    a1 = omsi_findpeaks_global()

    # Define compute of total intensity image
    a2 = analysis_generic.from_function(analysis_function=total_intensity,
                                        output_names=['total_intensities'])
    a2['msidata'] = a1['peak_cube']

    # Define the normalization of the peak cube
    a3 = analysis_generic.from_function(normalize_intensities)
    a3['msidata'] = a1['peak_cube']
    a3['normfactors'] = a2['total_intensities']


When using our command-line tool, all parameters that are not defined for any of the analyses are automatically exposed via command-line options. In contrast to our previous example, we here, e.g., do not set the input msidata and mzdata parameters for our global peak finder (a1). In this way, we can now easily set the input file we want to process directly via the command line. In cases where we want to expose a parameter via the command line but still want to provide a good default setting for the user, we can set the default value of a parameter via, e.g, ``a1.get_parameter_data_by_name('peakheight')['default'] = 3``.

To execute our above example from the command line we can now simply do the following:

.. code-block:: bash

    python run_workflow.py --script normalize_image.py
                           --ana_0:msidata $HOME/20120711_Brain.h5:/entry_0/data_0
                           --ana_0:mzdata  $HOME/20120711_Brain.h5:/entry_0/data_0/mz

In order to avoid collisions between parameters with the same name for different analyses, the tool prepends the unique ``analysis_identifier`` to each parameter. Since we did not set any explicit ``analysis_identifier` (e.g, via ``a1.analysis_identifier='a1'``), the tool automatically generated unique identifiers (i.e, ``ana_0``, ``ana_1``, and ``ana_3`` for our 3 analyses). To view all available command line option we can simply call the script with ``--help``. If one or more workflow scipts are given (here via seperate ``--script`` parameters), then all unfilled options of those workflows and the corresponding analyses will be listed as. E.g.


.. code-block:: python
    :linenos:

    newlappy:tools oruebel$ python run_workflow.py --script normalize_image.py --help
    usage: run_workflow.py --script SCRIPT [--save SAVE] [--profile]
                           [--memprofile]
                           [--loglevel {INFO,WARNING,CRITICAL,ERROR,DEBUG,NOTSET}]
                           --ana_0:msidata ANA_0:MSIDATA --ana_0:mzdata
                           ANA_0:MZDATA
                           [--ana_0:integration_width ANA_0:INTEGRATION_WIDTH]
                           [--ana_0:peakheight ANA_0:PEAKHEIGHT]
                           [--ana_0:slwindow ANA_0:SLWINDOW]
                           [--ana_0:smoothwidth ANA_0:SMOOTHWIDTH]
                           [--ana_1:axis ANA_1:AXIS]
                           [--reduce_memory_usage REDUCE_MEMORY_USAGE]
                           [--synchronize SYNCHRONIZE] [-h]

    Execute analysis workflow(s) based on a given set of scripts

    required arguments:
      --script SCRIPT       The workflow script to be executed. Multiple scripts
                            may be added via separate --script arguments (default:
                            None)

    optional arguments:
      --save SAVE           Define the file and experiment where all analysis
                            results should be stored. A new file will be created
                            if the given file does not exists but the directory
                            does. The filename is expected to be of the from:
                            <filename>:<entry_#> . If no experiment index is
                            given, then experiment index 0 (i.e, entry_0) will be
                            assumed by default. A validpath may, e.g, be
                            "test.h5:/entry_0" or jus "test.h5" (default: None)
      --profile             Enable runtime profiling of the analysis. NOTE: This
                            is intended for debugging and investigation of the
                            runtime behavior of an analysis.Enabling profiling
                            entails certain overheads in performance (default:
                            False)
      --memprofile          Enable runtime profiling of the memory usage of
                            analysis. NOTE: This is intended for debugging and
                            investigation of the runtime behavior of an analysis.
                            Enabling profiling entails certain overheads in
                            performance. (default: False)
      --loglevel {INFO,WARNING,CRITICAL,ERROR,DEBUG,NOTSET}
                            Specify the level of logging to be used. (default:
                            INFO)
      -h, --help            show this help message and exit

    ana_0:omsi.analysis.findpeaks.omsi_findpeaks_global:analysis settings:
      Analysis settings

      --ana_0:integration_width ANA_0:INTEGRATION_WIDTH
                            The window over which peaks should be integrated
                            (default: 0.1)
      --ana_0:peakheight ANA_0:PEAKHEIGHT
                            Peak height parameter (default: 2)
      --ana_0:slwindow ANA_0:SLWINDOW
                            Sliding window parameter (default: 100)
      --ana_0:smoothwidth ANA_0:SMOOTHWIDTH
                            Smooth width parameter (default: 3)

    ana_0:omsi.analysis.findpeaks.omsi_findpeaks_global:input data:
      Input data to be analyzed

      --ana_0:msidata ANA_0:MSIDATA
                            The MSI dataset to be analyzed (default: None)
      --ana_0:mzdata ANA_0:MZDATA
                            The m/z values for the spectra of the MSI dataset
                            (default: None)

    ana_1 : generic:
      --ana_1:axis ANA_1:AXIS

    optional workflow executor options:
      Additional, optional settings for the workflow execution controls

      --reduce_memory_usage REDUCE_MEMORY_USAGE
                            Reduce memory usage by pushing analyses to file each
                            time they complete, processing dependencies out-of-
                            core. (default: False)
      --synchronize SYNCHRONIZE
                            Place an MPI-barrier at the beginning of the exection
                            of the workflow. This can be useful when we require
                            that all MPI ranks are fully initalized. (default:
                            False)

    how to specify ndarray data?
    ----------------------------
    n-dimensional arrays stored in OpenMSI data files may be specified as
    input parameters via the following syntax:
          -- MSI data: <filename>.h5:/entry_#/data_#
          -- Analysis data: <filename>.h5:/entry_#/analysis_#/<dataname>
          -- Arbitrary dataset: <filename>.h5:<object_path>
    E.g. a valid definition may look like: 'test_brain_convert.h5:/entry_0/data_0'
    In rear cases we may need to manually define an array (e.g., a mask)
    Here we can use standard python syntax, e.g, '[1,2,3,4]' or '[[1, 3], [4, 5]]'

    This command-line tool has been auto-generated by BASTet (Berkeley Analysis & Storage Toolkit)




