"""
Module for managing instrument related data in OMSI files.
"""
from omsi.dataformat.omsi_file.format import omsi_format_common, omsi_format_instrument
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager


class omsi_instrument_manager(omsi_file_object_manager):
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
        return omsi_file_instrument.__create__(instrument_parent=self.instrument_parent,
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
        except:
            # Check whether the parent group has information about the instrument
            if check_parent:
                try:
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


class omsi_file_instrument(omsi_file_common):
    """Class for managing instrument specific data

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
                   instrument_parent,
                   instrument_name=None,
                   mzdata=None,
                   flush_io=True):
        """
        Create an instrument group and populate it with the given data.

        :param instrument_parent: The parent h5py group where the instrument group should be created in.
        :type instrument_parent. h5py.Group
        :param instrument_name: The name of the instrument
        :type instrument_name: string, None
        :param mzdata: Numpy array of the mz data values of the instrument
        :type mzdata: numpy array or None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all
                         data has been written to file

        :returns: The function returns the h5py HDF5 handler to the instrument info group created for the experiment.

        """
        import time
        # Create the group for instrument specific data
        instrument_group = instrument_parent.require_group(omsi_format_instrument.instrument_groupname)
        instrument_group.attrs[omsi_format_common.type_attribute] = "omsi_file_instrument"
        instrument_group.attrs[omsi_format_common.version_attribute] = omsi_format_instrument.current_version
        instrument_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())

        instrument_object = omsi_file_instrument.__create_instrument_info___(instrument_group=instrument_group,
                                                                             instrument_name=instrument_name,
                                                                             mzdata=mzdata)
        if flush_io:
            instrument_parent.file.flush()
        return instrument_object

    @classmethod
    def __create_instrument_info___(cls,
                                    instrument_group,
                                    instrument_name=None,
                                    mzdata=None):
        """
        Populate the empty instrument group with the given data.

        NOTE: This is a private helper function used to populate the instrument group with data. \
        Use the corresponding omsi_file_experiment.create_instrument_info(...) function to  \
        generate a new instrument information data in HDF5

        :param instrument_group: The h5py group to which the instrument group data should be written to.
        :param instrument_name: The name of the instrument used.
        :param mzdata: The mz data for the instrument.
        """
        # Name of the instrument
        instrumentnamedata = instrument_group.require_dataset(name=unicode(
            omsi_format_instrument.instrument_name), shape=(1,), dtype=omsi_format_common.str_type)
        if instrument_name is None:
            if len(instrumentnamedata[0]) == 0:
                instrumentnamedata[0] = "undefined"
        else:
            if omsi_format_common.str_type_unicode:
                instrumentnamedata[0] = instrument_name
            else:
                instrumentnamedata[0] = str(instrument_name)
        # MZ data for the instrument
        if mzdata is not None:
            instrumentmzdata = instrument_group.require_dataset(
                name=omsi_format_instrument.instrument_mz_name, shape=mzdata.shape, dtype=mzdata.dtype)
            instrumentmzdata[:] = mzdata[:]

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
        except:
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
            return self.managed_group[
                unicode(omsi_format_instrument.instrument_mz_name)]
        except:
            return None

    def set_instrument_name(self,
                            name):
        """
        Overwrite the current identifier string for the experiment with the given string.

        :param name: The new instrument name.
        :type name: string.
        """

        # Get the name of the intrument
        namedataset = self.get_instrument_name()
        # Create the dataset for the id name if it does not exist
        if namedataset is None:
            namedataset = self.managed_group.require_dataset(name=unicode(
                omsi_format_instrument.instrument_name), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            namedataset[0] = name
        else:
            namedataset[0] = str(name)
