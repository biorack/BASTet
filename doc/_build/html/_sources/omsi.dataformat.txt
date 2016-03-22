dataformat Package
==================

Main module for specification of data formats. This module also contains the
`omsi_file` module which specifies the OpenMSI HDF5 data format. In addition
it defines the base class for third-party file readers (i.e, `file_reader_base`)
and implements various basic file readers for third-party formats, e.g,
`img_file` and `mzml_file` for IMG and MZML data files respectively (among others).


.. autosummary::

    omsi.dataformat.omsi_file
    omsi.dataformat.omsi_file.analysis
    omsi.dataformat.omsi_file.common
    omsi.dataformat.omsi_file.dependencies
    omsi.dataformat.omsi_file.experiment
    omsi.dataformat.omsi_file.format
    omsi.dataformat.omsi_file.instrument
    omsi.dataformat.omsi_file.main_file
    omsi.dataformat.omsi_file.metadata_collection
    omsi.dataformat.omsi_file.methods
    omsi.dataformat.omsi_file.msidata
    omsi.dataformat.file_reader_base
    omsi.dataformat.bruckerflex_file
    omsi.dataformat.img_file
    omsi.dataformat.imzml_file
    omsi.dataformat.mzml_file

.. toctree::

    omsi.dataformat.omsi_file


:mod:`file_reader_base` Module
------------------------------

.. automodule:: omsi.dataformat.file_reader_base
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`img_file` Module
----------------------

.. automodule:: omsi.dataformat.img_file
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`mzml_file` Module
-----------------------

.. automodule:: omsi.dataformat.mzml_file
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`bruckerflex_file` Module
------------------------------

.. automodule:: omsi.dataformat.bruckerflex_file
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:



