omsi Package
============

Subpackages
-----------

The omsi package defines the main BASTet software stack used, e.g., by the OpenMSI project.
It contains a large range of functionality for interacting analyses, OpenMSI HDF5 data files,
analysis workflows and other related infrastructure. The following is a rough overview
of the various packages and modules

* **Analysis** :py:mod:`omsi.analysis` Package containing the base classes that facilitate the integration of new analysis with the BASTet software stack (e.g, the file format) and collection of specific analysis functionality.

.. autosummary::

   omsi.analysis
   omsi.analysis.compound_stats
   omsi.analysis.findpeaks
   omsi.analysis.msi_filtering
   omsi.analysis.multivariate_stats

* **Data Format** :py:mod:`omsi.dataformat` Package for implementation and specification of file formats. In particular this package contains the base API for interacting with OpenMSI HDF5 datasets.

.. autosummary::

   omsi.dataformat
   omsi.dataformat.omsi_file

* **Workflow** :py:mod:`omsi.workflow` Package with modules for specification and execution of analysis tasks and complex analysis workflows.

.. autosummary::

   omsi.workflow
   omsi.worflow.dirver
   omsi.workflow.executor

* **Data Structures** :py:mod:`omsi.datastructure` Package with a collection of various data structures and related classes used throughout the software stack, e.g., for metadata, analysis parameter data, runtime information data etc.

.. autosummary::

   omsi.datastructures
   omsi.datastructures.metadata

* **Shared** :py:mod:`omsi.shared` Package used to implement shared functionality and helper functions.

.. autosummary::

   omsi.shared
   omsi.shared.thirdparty

* **Tools** :py:mod:`omsi.tools` Package for collecting tools (e.g,. command-line programs) to help with particular tasks. This includes, e.g, tools for data conversion, document generation, etc.

.. autosummary::

   omsi.tools
   omsi.tools.misc

* **Templates** :py:mod:`omsi.templates` This package provides a collection of code templates to ease the development of additional components, e.g., analysis modules. As such, this package is NOT intended for direct usage but is rather just a library of code templates.

.. autosummary::

   omsi.templates

* **Examples** :py:mod:`omsi.examples` Package with a collection of various misc. example scripts.

.. autosummary::

   omsi.examples





.. toctree::

    omsi.analysis
    omsi.dataformat
    omsi.datastructures
    omsi.shared
    omsi.workflow
    omsi.tools
    omsi.examples
    omsi.templates

