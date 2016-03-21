OMSI Dataformat Package
========================

Main module for specification of the OpenMSI HDF5-based data format.
The module contains various sub-modules, with the main goal to
organize different categories of data.

Naming conventions for objects inside the HDF5 file are defined
in the `omsi_file.format` module. These are used by the manager
API classes to then implement the specific format.

The basic idea behind the design of the OpenMSI file format is the
concept of managed objects. Managed objects are objects in an HDF5 file
(usually HDF5 Groups --similar to directories just within a file)
that have a corresponding interface class in the API. These
classes in the API always start with the prefix `omsi_file_` and
inherit from `omsi_file.common.omsi_file_common`.

To make it easy to nest different objects, we also have the concept
of manager helper classes, which encapsulate common functionality
for creation and interaction with the objects when they are
contained in another object. Manager helper classes follow the
following naming convection `omsi_<objectname>_manager`, where
objectname is name of the object to be managed. E.g,
`omsi_instrument_manager` is used to help place instrument
groups inside another managed object. This is done by inheriting
from the given manager helper class.

This means, multiple inheritance is used in order to nest other
managed modules with other interfaces. This allows us to easily
encapsulate common interaction features in centralized locations
and construct more complex containers simply via inheritance.

The user of multiple inheritance and `super` can be tricky in python.
To simplify the use and ensure stability we use the following conventions:

  * All `omsi_file_*` manager classes must except the h5py.Group object
    they manage as input and call `super(..)__init__(managed_group)`
    with the managed group as parameter in their `__init__`.
  * All `omsi_<objectname>_manager` manager helper classes must except
    the h5py.Group that contains the object(s) that should be managed
    using the helper class as input and call `super(..)__init__(managed_group)`
    with the managed group as parameter in their in `__init__`.
  * All `omsi_file_*` manager classes must inherit from `omsi_file.common.omsi_file_common`
  * All `omsi_<objectname>_manager` manager helper classes must inherit from `object`
    (i.e., we use new-style classes).


:mod:`omsi_file` Package
------------------------

.. automodule:: omsi.dataformat.omsi_file
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`format` Module
--------------------

.. automodule:: omsi.dataformat.omsi_file.format
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`common` Module
--------------------

.. automodule:: omsi.dataformat.omsi_file.common
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`main_file` Module
-----------------------

.. automodule:: omsi.dataformat.omsi_file.main_file
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`experiment` Module
------------------------

.. automodule:: omsi.dataformat.omsi_file.experiment
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`metadata_collection` Module
---------------------------------

.. automodule:: omsi.dataformat.omsi_file.metadata_collection
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`instrument` Module
------------------------

.. automodule:: omsi.dataformat.omsi_file.instrument
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`methods` Module
---------------------

.. automodule:: omsi.dataformat.omsi_file.methods
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`msidata` Module
---------------------

.. automodule:: omsi.dataformat.omsi_file.msidata
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`analysis` Module
----------------------

.. automodule:: omsi.dataformat.omsi_file.analysis
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

:mod:`dependencies` Module
--------------------------

.. automodule:: omsi.dataformat.omsi_file.dependencies
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
