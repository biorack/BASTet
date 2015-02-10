"""
Module for common data format classes and functionality.
"""

from omsi.dataformat.omsi_file.format import *

class omsi_file_object_manager(object):
    """
    Base class used to define manager helper classes  used to manage contained managed objects.
    Managed objects are HDF5.Groups (or Datasets) with a corresponding manager API class
    and may be nested within other Managed objects.

    **What is a manager helper class?**

    Manager classes are used in the design of `omsi.dataformat.omsi_file` to encapsulate
    functionality needed for management of other manager objects. The expected use
    of this class, hence, is through multiple inheritance where the main base
    class is `omsi.dataformat.omsi_file.common.omsi_file_common`. This is important
    due to the use of super to accomodate multiple inheritance to allow object
    to manage an arbitrary number of other object and inherit from other object as well.

    **Use of super()**

    This class inherits only from object but calls super in the __init__(manager_group)
    with the manager_group as only input parameter, in the expectation that this
    class is used using multiple inheritance with `omsi_file_common` as main base class .

    Multiple inheritance is used in `omsi.dataformat.omsi_file module` when a class contains
    other managed objects and uses the manager classes (such as this one)
    to get all the features needed to manage those objects.

    All child classes of omsi_file_common call super(..).__init__(manager_group)
    and all manager helper classes (such as this one) use a single input parameter
    indicating the manager h5py.Group object that contains the given object.

    """
    def __init__(self, *args, **kwargs):
        """
        """
        super(omsi_file_object_manager, self).__init__(*args, **kwargs)


class omsi_file_common(object):
    """
    Base class for definition of file format modules for the OpenMSI data format.

    **Use of super()**

    This class inherits only from object and calls super in the __init__ without
    parameters. In the standard design pattern of the `omsi.dataformat.omsi_file module`,
    it is, therefore, the last class we inherit from in the case of multiple inheritance.

    Multiple inheritance is used in `omsi.dataformat.omsi_file module` when a class contains
    other managed objects and uses the manager classes, e.g, omsi_instrument_mangager etc.
    to get all the features needed to manage those objects.

    All child classes of omsi_file_common also call super(..).__init__(manager_group) but
    using a single input parameter indicating the manager h5py.Group object that
    contains the given object.

    """

    def __init__(self, managed_group):
        super(omsi_file_common, self).__init__()
        if managed_group is None:
            raise ValueError("Managed group must be initialized with an h5py.Group object. None given.")
        self.managed_group = managed_group
        self.name = self.managed_group.name

    @classmethod
    def is_managed(cls, in_object):
        """
        Check whether the given object is managed by any omsi API class.

        :param in_object: The object to be checked
        :type in_object: Any omsi_file API object or h5py.Dataset or h5py.Group or h5py.File object.
        """
        from omsi.dataformat.omsi_file.experiment import omsi_file_experiment
        from omsi.dataformat.omsi_file.main_file import omsi_file
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies, omsi_file_dependencydata
        from omsi.dataformat.omsi_file.instrument import omsi_file_instrument
        from omsi.dataformat.omsi_file.methods import omsi_file_methods
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        from omsi.dataformat.omsi_file.msidata import omsi_file_msidata

        managed_object = omsi_file_common.get_omsi_object(in_object)
        if isinstance(managed_object, omsi_file) or \
           isinstance(managed_object, omsi_file_experiment) or \
           isinstance(managed_object, omsi_file_methods) or\
           isinstance(managed_object, omsi_file_instrument) or \
           isinstance(managed_object, omsi_file_analysis) or \
           isinstance(managed_object, omsi_file_msidata) or \
           isinstance(managed_object, omsi_file_dependencies) or \
           isinstance(managed_object, omsi_file_dependencydata):
            return True
        else:
            return False

    @classmethod
    def get_h5py_object(cls, omsi_object, resolve_dependencies=False):
        """
        This static method is a convenience function used to retrieve the corresponding h5py
        interface object for any omsi file API object.

        :param omsi_object: omsi file API input object for which the corresponding h5py.Group, h5py.File, or
                        h5py.Dataset object should be retrieved. The omsi_object input may itself also be
                        a h5py.Group, h5py.File, or h5py.Dataset, in which case omsi_object itself is returned
                        by the function.
        :param resolve_dependencies: Set to True if omsi_file_dependencydata objects should be resolved to retrieve
                        the dependent object the dependency is pointing to. Dependencies are resolved recursively,
                        i.e., if a dependency points to another dependency then that one will be resolved as well.
                        Default value is False, i.e., the omis_file_dependency object itself is returned.
        :returns:  h5py.Group, h5py.File, or h5py.Dataset corresponding to the given omsi_object.

        :raises ValueError: A ValueError is raised in case that an unsupported omsi_object object is given, i.e.,
                            the input object is not a omsi_file API object nor a h5py Group, File, or
                            Dataset object.
        """
        from omsi.dataformat.omsi_file.main_file import omsi_file
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencydata

        if isinstance(omsi_object, omsi_file):
            h5pyobject = omsi_object.get_h5py_file()
        elif isinstance(omsi_object, omsi_file_dependencydata):
            if resolve_dependencies:
                h5pyobject = omsi_file_common.get_h5py_object(omsi_object.get_dependency_omsiobject(),
                                                              resolve_dependencies)
            else:
                h5pyobject = omsi_object.get_managed_group()
        elif isinstance(omsi_object, omsi_file_common):
            # This case covers:
            #  omsi_file_dependencies,  omsi_file_instrument, omsi_file_methods, omsi_file_analysis,
            #  omsi_file_msidata, omsi_file_experiment
            h5pyobject = omsi_object.get_managed_group()
        elif isinstance(omsi_object, h5py.File) or \
                isinstance(omsi_object, h5py.Group) or \
                isinstance(omsi_object, h5py.Dataset):
            h5pyobject = omsi_object
        else:
            raise ValueError("Unsupported input object given to omsi_file.get_h5py_object.")
        return h5pyobject

    @classmethod
    def get_omsi_object(cls, h5py_object, resolve_dependencies=False):
        """
        This static method is convenience function used to retrieve the corresponding interface class for a
        given h5py group object.

        :param h5py_object: h5py object for which the corresponding omsi_file API object should be generated.
        :param resolve_dependencies: Set to True if omsi_file_dependencydata objects should be resolved to retrieve
                        the dependent object the dependency is pointing to. Dependencies are resolved recursively,
                        i.e., if a dependency points to another dependency then that one will be resolved as well.
                        Default value is False, i.e., the omis_file_dependency object itself is returned.

        :returns: None in case no corresponding object was found. Otherwise an instance of:

            * omsi_file : If the given object is a h5py.File object
            * omsi_file_experiment : If the given object is an experiment groupt
            * omsi_file_methods : If the given object is a method group
            * omsi_file_instrument : If the given object is an instrument group
            * omsi_file_analysis : If the given object is an analysis group
            * omsi_file_msidata : If the given object is a MSI data group
            * omsi_file_dependencydata : If the fiven object is a dependency group
            * The input h5py_object: If the given objet is a h5py.Dataset
            * None: In case that an unknown type is given. E.g. a user may have \
                    given an unmanaged Group object that does not have a \
                    corresponding omsi file API object.
        """
        from omsi.dataformat.omsi_file.experiment import omsi_file_experiment
        from omsi.dataformat.omsi_file.main_file import omsi_file
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies, omsi_file_dependencydata
        from omsi.dataformat.omsi_file.instrument import omsi_file_instrument
        from omsi.dataformat.omsi_file.methods import omsi_file_methods
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        from omsi.dataformat.omsi_file.msidata import omsi_file_msidata

        # If the input object is already an omsi API object then return it as
        # is
        if isinstance(h5py_object, omsi_file) or \
           isinstance(h5py_object, omsi_file_experiment) or \
           isinstance(h5py_object, omsi_file_methods) or\
           isinstance(h5py_object, omsi_file_instrument) or \
           isinstance(h5py_object, omsi_file_analysis) or \
           isinstance(h5py_object, omsi_file_msidata) or \
           isinstance(h5py_object, omsi_file_dependencies) or \
           isinstance(h5py_object, omsi_file_dependencydata):
            return h5py_object

        if isinstance(h5py_object, h5py.File):
            return omsi_file(h5py_object)
        elif isinstance(h5py_object, h5py.Group):
            # Check if the group has an explicit type attribute
            try:
                # Try to determine the type of the group based on the
                # attributes
                type_attribute = h5py_object.attrs[omsi_format_common.type_attribute]
                if type_attribute == "omsi_file_experiment":
                    return omsi_file_experiment(h5py_object)
                elif type_attribute == "omsi_file_methods" or type_attribute == "omsi_file_methods":
                    return omsi_file_methods(h5py_object)
                elif type_attribute == "omsi_file_instrument":
                    return omsi_file_instrument(h5py_object)
                elif type_attribute == "omsi_file_analysis":
                    return omsi_file_analysis(h5py_object)
                elif type_attribute == "omsi_file_msidata":
                    return omsi_file_msidata(h5py_object)
                elif type_attribute == "omsi_file":
                    return omsi_file(h5py_object)
                elif type_attribute == "omsi_file_dependencydata":
                    omsiobject = omsi_file_dependencydata(h5py_object)
                    if resolve_dependencies:
                        return omsi_file_common.get_omsi_object(
                            omsi_file_common.get_h5py_object(omsiobject.get_dependency_omsiobject(),
                                                             resolve_dependencies))
                    else:
                        return omsiobject
                elif type_attribute == "omsi_file_dependencies":
                    return omsi_file_dependencies(h5py_object)
                else:
                    return None
            except:
                # If the attribute is missing, then try to determin the type
                # based on he name of group
                groupname = h5py_object.name.split("/")[-1]
                parentgroupname = h5py_object.parent.name.split("/")[-1]
                if groupname.startswith(omsi_format_experiment.exp_groupname):
                    return omsi_file_experiment(h5py_object)
                elif groupname.startswith(omsi_format_methods.methods_groupname) or \
                        groupname.startswith(omsi_format_methods.methods_old_groupname):
                    return omsi_file_methods(h5py_object)
                elif groupname.startswith(omsi_format_instrument.instrument_groupname):
                    return omsi_file_instrument(h5py_object)
                elif groupname.startswith(omsi_format_analysis.analysis_groupname):
                    return omsi_file_analysis(h5py_object)
                elif groupname.startswith(omsi_format_data.data_groupname):
                    return omsi_file_msidata(h5py_object)
                elif groupname.startswith(omsi_format_dependencies.dependencies_groupname):
                    return omsi_file_dependencies(h5py_object)
                elif parentgroupname.startswith(omsi_format_dependencies.dependencies_groupname):
                    omsiobject = omsi_file_dependencydata(h5py_object)
                    if resolve_dependencies:
                        return omsi_file_common.get_omsi_object(
                            omsi_file_common.get_h5py_object(omsiobject.get_dependency_omsiobject(),
                                                             resolve_dependencies))
                    else:
                        return omsiobject
                elif groupname == "":  # We are at the root group
                    return omsi_file(h5py_object.file)
                else:
                    return None
        elif isinstance(h5py_object, h5py.Dataset):
            return h5py_object
        else:
            return None

    @classmethod
    def get_num_items(cls, file_group, basename=""):
        """
        Get the number of object with the given basename at the given path

        :param file_group: The h5py object to be examined
        :param basename: The name that should be searched for.

        :returns: Number of objexts with the given basename at the given path
        """

        numitems = 0
        # Iterate through all groups of the root folder
        for item_obj in file_group.items():
            if item_obj[0].startswith(basename):
                numitems += 1
        return numitems

    def __getitem__(self, key):
        """
        Support direct read interaction with the method h5py group
        """
        return self.managed_group[key]

    def __setitem__(self, key, value):
        """
        Support direct write interaction with the method h5py group.
        """
        # If the version of h5py does not support automatic storage of strings and we have a string then do-it-yourself
        if (isinstance(value, unicode) or isinstance(value, str)) and not omsi_format_common.str_type_unicode:
            key_dataset = self.managed_group.require_dataset(name=unicode(key),
                                                             shape=(1,),
                                                             dtype=omsi_format_common.str_type)
            key_dataset[0] = str(value)
        else:
            self.managed_group[key] = value

    def __eq__(self, value):
        """
        Check whether the two objects have the same h5py name.
        """
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """
        Check whether the two objects have different h5py names
        """
        return not self.__eq__(value)

    def get_managed_group(self):
        """
        Return the h5py object with the analysis data.

        The returned object can be used to read data directly from the HDF5 file.
        Write operations to the analysis group can be performed only if the
        associated omsi_file was opened with write permissions.

        :returns: h5py object for the analysis group.
        """
        return self.managed_group

    def get_version(self):
        """
        Get the omsi version for the representation of this object in the HDF5 file
        """
        try:
            return self.managed_group.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """
        Get the timestamp when the analysis group was created in the HDF5 file.

        :returns: Python timestamp string generated using time.ctime().
                  None may be returned in case that the timestamp does not exists
                  or cannot be retrieved from the file for some reason.

        """
        try:
            return self.managed_group.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def items(self):
        """
        Get the list of items associdated with the h5py.Group object managed by this object
        """
        return self.managed_group.items()
