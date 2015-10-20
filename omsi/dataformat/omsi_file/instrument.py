"""
Module for managing instrument related data in OMSI files.
"""
from omsi.dataformat.omsi_file.format import omsi_format_instrument
from omsi.dataformat.omsi_file.common import omsi_file_common
from omsi.dataformat.omsi_file.metadata_collection import omsi_metadata_collection_manager
from omsi.dataformat.omsi_file.metadata_collection import omsi_file_metadata_collection
from omsi.datastructures.metadata.metadata_data import metadata_value, metadata_dict


class omsi_instrument_manager(omsi_metadata_collection_manager):
    """
    Instrument manager helper class
    used to define common functionality needed for instrument-related data.
    Usually, a class that defines a format that contains an omsi_file_methods object
    will inherit from this class (in addition to omsi_file_common) to acquire the common
    features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar instrument_parent: The h5py.Group parent object containing the instrument object to be managed.

    """
    def __init__(self, instrument_parent):
        super(omsi_instrument_manager, self).__init__(instrument_parent)
        self.instrument_parent = instrument_parent

    def create_instrument_info(self,
                               instrument_name=None,
                               mzdata=None,
                               flush_io=True):
        """
        Add information about the instrument used for creating the images for this experiment.

        :param instrument_name: The name of the instrument
        :type instrument_name: string, None
        :param mzdata: Numpy array of the mz data values of the instrument
        :type mzdata: numpy array or None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all
                         data has been written to file

        :returns: The function returns the h5py HDF5 handler to the instrument info group created for the experiment.

        """
        return omsi_file_instrument.__create__(parent_group=self.instrument_parent,
                                               instrument_name=instrument_name,
                                               mzdata=mzdata,
                                               flush_io=flush_io)

    def get_instrument_info(self,
                            check_parent=True):
        """
        Get the HDF5 group opbject with the instrument information.

        :param check_parent: If no method group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the method. (default=True)

        :returns:  omsi_file_instrument object for the requested instrument info. The function returns \
                   None in case no instrument information was found for the experiment
        """
        try:
            return omsi_file_instrument(self.instrument_parent[unicode(omsi_format_instrument.instrument_groupname)])
        except KeyError:
            # Check whether the parent group has information about the instrument
            if check_parent:
                return omsi_file_common.get_omsi_object(self.instrument_parent.parent).get_instrument_info()
        except:
            pass
        return None

    def has_instrument_info(self,
                            check_parent=False):
        """
        Check whether custom instrument information is available for this dataset.

        :param check_parent: If no instrument group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the instrument. (default=False)

        :returns: Boolean indicating whether instrument info is available.
        """
        return self.get_instrument_info(check_parent=check_parent) is not None


class omsi_file_instrument(omsi_file_metadata_collection):
    """
    Class for managing instrument specific data

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
    def __create__(cls,
                   parent_group,
                   instrument_name=None,
                   mzdata=None,
                   flush_io=True):
        """
        Create an instrument group and populate it with the given data.

        :param parent_group: The parent h5py group where the instrument group should be created in.
        :type parent_group. h5py.Group
        :param instrument_name: The name of the instrument
        :type instrument_name: string, None
        :param mzdata: Numpy array of the mz data values of the instrument
        :type mzdata: numpy array or None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all
                         data has been written to file

        :returns: The function returns the h5py HDF5 handler to the instrument info group created for the experiment.

        """
        if instrument_name is not None or mzdata is not None:
            all_meta = metadata_dict()
            if instrument_name is not None:
                all_meta[omsi_format_instrument.instrument_name] = \
                    metadata_value(value=instrument_name,
                                   name=omsi_format_instrument.instrument_name,
                                   description='Name of the instrument')
            if mzdata is not None:
                all_meta[omsi_format_instrument.instrument_name] = \
                    metadata_value(value=mzdata,
                                   name=omsi_format_instrument.instrument_mz_name,
                                   description='The global m/z axis for the recordings')

        else:
            all_meta = None

        # Initialize the group and populate the data using the create method of the parent class
        metadata_obj = omsi_file_metadata_collection.___create___(
            parent_group=parent_group,
            group_name=omsi_format_instrument.instrument_groupname,
            metadata=all_meta,
            type_attr_value="omsi_file_instrument",
            version_attr_value=omsi_format_instrument.current_version,
            flush_io=flush_io)

        if flush_io:
            parent_group.file.flush()
        return omsi_file_instrument.__create_instrument_info___(instrument_group=metadata_obj.managed_group)

    @classmethod
    def __create_instrument_info___(cls,
                                    instrument_group):
        """
        Populate the empty instrument group with the given data.

        NOTE: This is a private helper function used to populate the instrument group with data. \
        Use the corresponding omsi_file_experiment.create_instrument_info(...) function to  \
        generate a new instrument information data in HDF5

        :return: omsi_file_instrument object that manages the given instrument_group
        """
        return omsi_file_instrument(instrument_group)

    def __init__(self,
                 instrument_group):
        """
            Initalize the instrument object given the h5py object of the instrument group

            :param instrument_group: The h5py object with the instrument group of the omsi hdf5 file.
            """
        super(omsi_file_instrument, self).__init__(instrument_group)
        # Initialized by super call
        # self.managed_group = instrument_group
        # self.name = self.managed_group.name

    def has_instrument_name(self):
        """
        Check whether a name has been saved for the instrument

        :return: bool
        """
        return self.get_instrument_name() is not None

    def get_instrument_name(self):
        """
        Get the HDF5 dataset with the name of the instrument.

        To get the string of the instrument name use:
        get_instrument_name()[...]

        :returns: h5py object to the dataset with the instrument name.
                  Returns None in case no method name is found.
        """
        if self.managed_group is None:
            return None
        try:
            instrument_name = self.managed_group[unicode(omsi_format_instrument.instrument_name)]
        except KeyError:
            instrument_name = None
        return instrument_name

    def get_instrument_mz(self):
        """
        Get the HDF5 dataset with the mz data for the instrument.

        To get the numpy array of the full mz data use:
        get_instrument_mz()[:]

        :returns: Returns the h5py object with the instrument mz data.
                 Returns None in case no mz data was found for the instrument.
        """
        if self.managed_group is None:
            return None
        try:
            return self.managed_group[unicode(omsi_format_instrument.instrument_mz_name)]
        except KeyError:
            return None

    def set_instrument_name(self,
                            name):
        """
        Overwrite the current identifier string for the experiment with the given string.

        :param name: The new instrument name.
        :type name: string.
        """
        self.add_metadata(metadata=metadata_value(value=name,
                                                  name=omsi_format_instrument.instrument_name,
                                                  description='Name of the instrument',
                                                  unit=None,
                                                  ontology=None))
