"""
Module for management of method specific data in OMSI data files
"""
from omsi.dataformat.omsi_file.format import omsi_format_common, omsi_format_methods
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager


class omsi_methods_manager(omsi_file_object_manager):
    """
    This is a file format manager helper class
    used to define common functionality needed for methods-related data.
    Usually, a class that defines a format that contains an omsi_file_methods object
    will inherit from this class (in addition to omsi_file_common) to acquire the common
    features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar method_parent: The parent h5py.Group object containing the method object to be managed

    """
    def __init__(self, methods_parent=None):
        super(omsi_methods_manager, self).__init__(methods_parent)
        self.method_parent = methods_parent

    def create_method_info(self, method_name=None, flush_io=True):
        """
        Add information about the method imaged to the experiment.
        Note, if a methods group already exists, then that group
        will be used. If method_name is not None, then the existing
        name will be overwritten by the new value.

        :param method_name: Optional name of the method
        :type method_name: str, None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 buffers are
                         flushed so that all data has been written to file

        :returns: h5py object of the newly created method group.
        """
        return omsi_file_methods.___create___(self.method_parent,
                                              method_name=method_name,
                                              flush_io=flush_io)

    def get_method_info(self, check_parent=True):
        """
        Get the omsi_file_methods object with the method information.

        :param check_parent: If no method group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the method. (default=True)

        :returns: omsi_file_methods object for the requested method info. The function returns
                  None in case no method information was found for the experiment

        """
        if self.method_parent is not None:
            try:
                return omsi_file_methods(self.method_parent[unicode(omsi_format_methods.methods_groupname)])
            except:
                try:
                    return omsi_file_methods(self.method_parent[unicode(omsi_format_methods.methods_old_groupname)])
                except:
                    if check_parent:
                        try:
                            return omsi_file_common.get_omsi_object(self.method_parent.parent).get_method_info()
                        except:
                            pass
        return None

    def has_method_info(self, check_parent=False):
        """
        Check whether custom method information is available for this dataset.

        :param check_parent: If no method group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the method. (default=False)
        :returns: Boolean indicating whether method info is available.
        """
        return self.get_method_info(check_parent=check_parent) is not None


class omsi_file_methods(omsi_file_common):
    """
    Class for managing method specific data.

    **Use of super():**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common`.
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the __init__ function calls
    super(...).__init__(manager_group) with a single  parameter indicating the
    parent group.

    **Inherited Instance Variables**

    :ivar managed_group: The group that is managed by this object
    :ivar name: Name of the managed group

    """
    @classmethod
    def ___create___(cls, parent_group, method_name=None, flush_io=True):
        """
        Add information about the method imaged to the experiment.
        If method_name is not None, then the existing
        name will be overwritten by the new value.

        :param parent_group: The h5py.Group where the method object should be created
        :type parent_group: h5py.Group
        :param method_name: Optional name of the method
        :type method_name: str, None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 buffers are
                         flushed so that all data has been written to file

        :returns: h5py object of the newly created method group.
        """
        if omsi_format_methods.methods_groupname in parent_group.items():
            out_method = omsi_file_methods(methods_group=parent_group[omsi_format_methods.methods_groupname])
            if method_name is not None:
                out_method.set_method_name(method_name)
        else:
            import time
            method_group = parent_group.require_group(omsi_format_methods.methods_groupname)
            method_group.attrs[omsi_format_common.type_attribute] = "omsi_file_methods"
            method_group.attrs[omsi_format_common.version_attribute] = omsi_format_methods.current_version
            method_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
            out_method = omsi_file_methods.__create_method_info__(method_group=method_group, method_name=method_name)
        if flush_io:
            parent_group.file.flush()
        return out_method

    @classmethod
    def __create_method_info__(cls, method_group, method_name=None):
        """
        Add information about the method  imaged to the experiment

        NOTE: This is a private helper function used to populate the given group for the \
        method data. Use the corresponding omsi_file_experiment.create_method_info(...) \
        to create a new method group and populate it with data.

        :param method_group: h5py group object that should be populated with the method data.
        :param method_name: Optional name of the method
        :type method_name: string, None

        :returns: h5py object of the newly created method group.
        """
        methodnamedata = method_group.require_dataset(name=unicode(
            omsi_format_methods.methods_name), shape=(1,), dtype=omsi_format_common.str_type)
        if method_name is None:
            if len(methodnamedata[0]) == 0:
                methodnamedata[0] = "undefined"
        else:
            if omsi_format_common.str_type_unicode:
                methodnamedata[0] = method_name
            else:
                methodnamedata[0] = str(method_name)
        return omsi_file_methods(method_group)

    def __init__(self, methods_group):
        """
            Initialize the method object given the h5py object of the method group

            :param methods_group: The h5py object with the method group of the omsi hdf5 file.
            """
        super(omsi_file_methods, self).__init__(methods_group)
        # Initialized by super call
        # self.managed_group = methods_group
        # self.name = self.managed_group.name

    def get_method_name(self):
        """
        Get the HDF5 dataset with the name of the method.

        To retrieve the name string use get_method_name()[...]

        :returns: h5py object where the method name is stored.
                 Returns None in case no method name is found.
        """
        if self.managed_group is None:
            return None
        try:
            method_name = self.managed_group[unicode(omsi_format_methods.methods_name)]
        except:
            method_name = None
        return method_name

    def set_method_name(self, name_string):
        """
        Overwrite the name string for the method with the given name string

        :param name_string: The new method name.
        :type name_string: string
        """
        # Get the name of the intrument
        name = self.get_method_name()
        # Create the dataset for the instrument name if it does not exist
        if name is None:
            name = self.managed_group.require_dataset(name=unicode(
                omsi_format_methods.methods_name), shape=(1,), dtype=omsi_format_common.str_type)

        if omsi_format_common.str_type_unicode:
            name[0] = name_string
        else:
            name[0] = str(name_string)

