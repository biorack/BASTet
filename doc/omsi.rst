omsi Package
============

Subpackages
-----------

The omsi package defines the main OpenMSI python software stack. It contains a large range
of functionality for interacting with OpenMSI HDF5 files. The package consists of the
folliwing subpackages designed for:

    * **Analysis** (``omsi.analysis``) This package contains the base classes that facilitate the integration of new analysis with the software stack and in particular the OpenMSI HDF5 data format.
    * **Data Format** (``omsi.dataformat``) This package contains a series of modules for interacting with different mass spectrometry imaging formats. In particular this package contains the base API for interacting with OpenMSI HDF5 datasets.
    * **Workflow** (``omsi.workflow``) This pacakge contains a series of modules to define and execute analyses and define and execute complex analysis workflows.
    * **Tools** (``omsi.tools``) This module contains various helpful tools based on the omsi software stack. This includes, e.g., tools for converting img data to HDF5 (imgToHDF5).
    * **Shared** (``omsi.shared``) This module contains functionality that is required by different codes (e.g., the omsi_server and omsi) that does not quite fit into the other modules but that should not be replicated in multiple places.
    * **Examples** (``omsi.examples``) Various example scripts.

.. toctree::

    omsi.analysis
    omsi.dataformat
    omsi.shared
    omsi.workflow
    omsi.tools
    omsi.examples

