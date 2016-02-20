"""
Module for common data format classes and functionality.
"""
from omsi.dataformat.omsi_file.format import *

# TODO Extend parse_path_string(..) and create_path_string(...) to allow for the inclusing of data selections


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

    :ivar managed_group: The h5py.Group object managed by the class
    :ivar name: The path to the object in the hdf5 file. Same as managed_group.name
    :ivar file: The h5py.File object the managed_group is associated with. Same as managed_group.file

    """

    def __init__(self, managed_group):
        super(omsi_file_common, self).__init__()
        if managed_group is None:
            raise ValueError("Managed group must be initialized with an h5py.Group object. None given.")
        self.managed_group = managed_group
        self.name = self.managed_group.name
        self.file = self.managed_group.file

    @classmethod
    def is_managed(cls, in_object):
        """
        Check whether the given object is managed by any omsi API class.

        :param in_object: The object to be checked
        :type in_object: Any omsi_file API object or h5py.Dataset or h5py.Group or h5py.File object.
        """
        try:
            managed_object = omsi_file_common.get_omsi_object(in_object)
        except ValueError:
            managed_object = None
        if isinstance(managed_object, omsi_file_common):
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
                        This may also be a string describing the requested object based on a combination of the
                        path to the file and a path ot the object <filename.h5>:<object_path>
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
            * The input h5py_object: If the given object is a h5py.Dataset or h5py.Group
            * None: In case that an unknown type is given.
        """
        from omsi.dataformat.omsi_file.experiment import omsi_file_experiment
        from omsi.dataformat.omsi_file.main_file import omsi_file
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies, omsi_file_dependencydata
        from omsi.dataformat.omsi_file.instrument import omsi_file_instrument
        from omsi.dataformat.omsi_file.methods import omsi_file_methods
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        from omsi.dataformat.omsi_file.msidata import omsi_file_msidata
        from omsi.dataformat.omsi_file.metadata_collection import omsi_file_metadata_collection

        # If the input object is already an omsi API object then return it as is
        if isinstance(h5py_object, omsi_file) or \
                isinstance(h5py_object, omsi_file_experiment) or \
                isinstance(h5py_object, omsi_file_methods) or\
                isinstance(h5py_object, omsi_file_instrument) or \
                isinstance(h5py_object, omsi_file_analysis) or \
                isinstance(h5py_object, omsi_file_msidata) or \
                isinstance(h5py_object, omsi_file_dependencies) or \
                isinstance(h5py_object, omsi_file_dependencydata) or \
                isinstance(h5py_object, omsi_file_metadata_collection):
            return h5py_object
        # IF we have an h5py.File then create and omsi_file
        if isinstance(h5py_object, h5py.File):
            return omsi_file(h5py_object)
        # If we have an h5py.Group, then try to create a corresponding omsi API object
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
                elif type_attribute == 'omsi_file_metadata_collection':
                    return omsi_file_metadata_collection(h5py_object)
                else:
                    return h5py_object
            except:
                # If the attribute is missing, then try to determine the type
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
                elif parentgroupname.startswith(omsi_format_metadata_collection.metadata_collection_groupname_default):
                    return omsi_file_metadata_collection(h5py_object)
                elif groupname == "":  # We are at the root group
                    return omsi_file(h5py_object.file)
                else:
                    return h5py_object
        # If we have an hpy.Dataset then we don't have a corresponding API object. Return as is
        elif isinstance(h5py_object, h5py.Dataset):
            return h5py_object
        elif isinstance(h5py_object, basestring):
            import os
            filename, object_path = cls.parse_path_string(h5py_object)
            if filename is not None and os.path.exists(filename) and os.path.isfile(filename):
                try:
                    curr_omsi_file = omsi_file(filename, 'r')
                except:
                    return None
            else:
                return None

            if object_path is not None:
                try:
                    file_object = curr_omsi_file.managed_group[object_path]
                except KeyError:
                    return None
                return omsi_file_common.get_omsi_object(file_object, resolve_dependencies=True)
            else:
                if isinstance(curr_omsi_file, omsi_file):
                    return curr_omsi_file
                else:
                    raise ValueError('omsi_file_common.Invalid path or file')
        else:
            return None

    @staticmethod
    def parse_path_string(path):
        """
        Given a string of the form <filename.h5>:<object_path> retrieve
        the name of the file and the object path.

        :param path: The string defining the file and object path.

        :return: Tuple with the filename and the object path. Both may
            be None depending on whether an object_path is given and
            whether the path string is valid.

        :raises: ValueError in case that an invalid string is given
        """
        if not isinstance(path, basestring):
            raise ValueError('The given path is not a string.')
        file_path = None
        object_path = None
        # Case 1: We have a string of the form <filename>:<objectname>
        if ".h5" not in path and ":" in path:
            split_string = path.split(':')
            if len(split_string) == 2:
                file_path = split_string[0]
                object_path = split_string[1]
            elif len(split_string) == 1:
                file_path = split_string[0]
                object_path= None
            else:
                raise ValueError("Invalid path string given")
        # Case 2: We have a string of the form filename.h5
        elif path.endswith(".h5") and ".h5:" not in path:
            file_path = path
            object_path = None
        # Case 3: We have a string of the form <filename>.h5:<objectname>
        elif ".h5:" in path:
            split_string = path.split('.h5:')
            file_path = split_string[0] + ".h5"
            if len(split_string) == 2:
                object_path = split_string[1]
            elif len(split_string) > 2:
                raise ValueError("Invalid path string given")
        # Case 4: We have a string that contains only the path to the object itself
        else:
            file_path = None
            object_path = path
        return file_path, object_path

    @staticmethod
    def create_path_string(filename, objectname):
        """
        Given the name of the file and the object path within the file,
        create a string describing the external reference to the datra

        :param filename: The full or relative path to the file
        :param objectname: The object path in the HDF5 file

        :return: String describing the path to the object
        """
        import os
        import warnings
        if not filename.endswith('.h5'):
            warnings.warn('Filename does not end with .h5 as expected.')
        objectname_string = objectname if objectname is not None else ""
        filename_string = filename if filename is not None else ""
        return filename_string + ":" + objectname_string

    @staticmethod
    def same_file(filename1, filename2):
        """
        Check whether two files are the same.

        This function uses the os.path.samefile(...) method to compare files and
        falls back to comparing the absolute paths of files if samefile should
        fail or cannot be imported.

        :param filename1: The name of the first file
        :param filename2: The name of the second file
        :return:
        """
        try:
            from os.path import samefile, abspath
            try:
                return samefile(filename1, filename2)
            except:
                pass
        except ImportError:  # E.g. we are on Windows
            pass
        return abspath(filename1) == abspath(filename2)

    @classmethod
    def get_num_items(cls, file_group, basename=""):
        """
        Get the number of object with the given basename at the given path

        :param file_group: The h5py object to be examined
        :param basename: The name that should be searched for.

        :returns: Number of objects with the given basename at the given path
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
            equal &= self.same_file(value.file.filename, self.file.filename)
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
