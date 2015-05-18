"""
Module for managing MSI data in OMSI data files
"""

import math
from omsi.dataformat.omsi_file.format import omsi_format_msidata, \
    omsi_format_msidata_partial_cube, \
    omsi_format_msidata_partial_spectra, \
    omsi_format_common
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager
from omsi.dataformat.omsi_file.dependencies import omsi_dependencies_manager
from omsi.dataformat.omsi_file.methods import omsi_methods_manager
from omsi.dataformat.omsi_file.instrument import omsi_instrument_manager
from omsi.dataformat.omsi_file.metadata_collection import omsi_metadata_collection_manager
import numpy as np


class omsi_msidata_manager(omsi_file_object_manager):
    """
    MSI-data  manager helper class used to define common functionality needed for
    msidata-related data. Usually, a class that defines a format that contains an
    omsi_file_msidata object will inherit from this class (in addition to omsi_file_common)
    to acquire the common features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar msidata_parent: The h5py.Group parent object containing the instrument object to be managed.

    """
    def __init__(self, msidata_parent):
        """
        Initatize the manger object.

        :param msidata_parent: The h5py.Group parent object for the msi data.
        """
        super(omsi_msidata_manager, self).__init__(msidata_parent)
        self.msidata_parent = msidata_parent

    def get_num_msidata(self):
        """
        Get the number of raw mass spectrometry images stored for a given experiment

        :returns: Integer indicating the number of msi datasets available for the experiment.
        """

        return omsi_file_common.get_num_items(self.msidata_parent,
                                              omsi_format_msidata.data_groupname)

    def get_msidata(self,
                    data_index,
                    fill_space=True,
                    fill_spectra=True,
                    preload_mz=True,
                    preload_xy_index=True):
        """
        Get the dataset with the given index for the given experiment.

        For more detailed information about the use of the fill_space and fill_spectra and preload_mz and
        preload_xy_index options, see the init function of omsi.dataformat.omsi_file_msidata.

        :param data_index: Index of the dataset.
        :type data_index: unsigned int
        :param fill_space: Define whether the data should be padded in space (filled with 0's)
                        when accessing the data using [..] operator so that the data behaves like a 3D cube.
        :param fill_spectra: Define whether the spectra should completed by adding 0's so that
                        all spectra retrived via the [..] opeator so that always spectra of the full
                        length are returned.
        :param preload_mz: Should the data for the mz axis be kept in memory or loaded on the fly when needed.
        :param preload_xy_index: Should the xy index (if available) be preloaderd into memory or should the
                        required data be loaded on the fly when needed.
        :returns: omsi_file_msidata object for the given data_index or None in case the data
                 with given index does not exist or the access failed for any
                 other reason.
        """
        msidata_name = unicode(omsi_format_msidata.data_groupname + str(data_index))
        try:
            return omsi_file_msidata(self.msidata_parent[msidata_name],
                                     fill_space=fill_space,
                                     fill_spectra=fill_spectra,
                                     preload_mz=preload_mz,
                                     preload_xy_index=preload_xy_index)
        except KeyError:
            return None

    def get_msidata_by_name(self, data_name):
        """
        Get the h5py data object for the the msidata with the given name.

        :param data_name: The name of the dataset
        :type data_name: string

        :returns: h5py object of the dataset or None in case the dataset is not found.
        """
        # Iterate through all items of the experiment
        for item_obj in self.msidata_parent.items():
            if item_obj[0] == data_name:
                return self.msidata_parent[item_obj[0]]
        return None

    def create_msidata_full_cube(self,
                                 data_shape,
                                 data_type='f',
                                 mzdata_type='f',
                                 chunks=None,
                                 compression=None,
                                 compression_opts=None,
                                 flush_io=True):
        """
        Create a new mass spectrometry imaging dataset for the given experiment written as a full 3D cube.

        :param data_shape: Shape of the dataset. Eg. shape=(10,10,10) creates a 3D dataset with
                         10 entries per dimension
        :param data_type:  numpy style datatype to be used for the dataset.
        :param mzdata_type: numpy style datatype to be used for the mz data array.
        :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk
                       sizes to be used in x,y, and m/z explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has
                        been written to file

        :returns: The following two empty (but approbriately sized) h5py datasets are returned in order
                 to be filled with data:

            * ``data_dataset`` : Primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * ``mz_dataset`` : h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the group
        data_group = omsi_file_msidata.__create_msidata_group__(parent_group=self.msidata_parent,
                                                                msidata_index=None,
                                                                flush_io=False)
        data_dataset, mz_dataset = omsi_file_msidata.__create_msidata_full_cube__(data_group=data_group,
                                                                                  data_shape=data_shape,
                                                                                  data_type=data_type,
                                                                                  mzdata_type=mzdata_type,
                                                                                  chunks=chunks,
                                                                                  compression=compression,
                                                                                  compression_opts=compression_opts)
        if flush_io:
            self.msidata_parent.file.flush()
        return data_dataset, mz_dataset, data_group

    def create_msidata_partial_cube(self,
                                    data_shape,
                                    mask,
                                    data_type='f',
                                    mzdata_type='f',
                                    chunks=None,
                                    compression=None,
                                    compression_opts=None,
                                    flush_io=True):
        """
        Create a new mass spectrometry imaging dataset for the given experiment written as a partial 3D cube
        of complete spectra.

        :param data_shape: Shape of the dataset. Eg. shape=(10,10,10) creates a 3D dataset with
                          10 entries per dimension
        :param mask: 2D boolean NumPy array used as mask to indicate which (x,y) locations have
                    spectra associated with them.
        :param data_type:  numpy style datatype to be used for the dataset.
        :param mzdata_type: numpy style datatype to be used for the mz data array.
        :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk sizes
                       to be used in x,y, and m/z explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has
                       been written to file

        :returns: The following two empty (but approbriately sized) h5py datasets are returned in order to
                 be filled with data:

            * ``data_dataset`` : Primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * ``mz_dataset`` : h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        :returns: The following already complete dataset

            * ``xy_index_dataset`` : This dataset indicates for each xy location to which index in \
                                     ``data_dataset`` the location corresponds to. This dataset is \
                                     needed to identify where spectra need to be written to.

        :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the group
        data_group = omsi_file_msidata.__create_msidata_group__(parent_group=self.msidata_parent,
                                                                msidata_index=None,
                                                                flush_io=False)
        data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset = \
            omsi_file_msidata.__create_msidata_partial_cube__(data_group=data_group,
                                                              data_shape=data_shape,
                                                              mask=mask,
                                                              data_type=data_type,
                                                              mzdata_type=mzdata_type,
                                                              chunks=chunks,
                                                              compression=compression,
                                                              compression_opts=compression_opts)
        if flush_io:
            self.msidata_parent.file.flush()
        return data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset, data_group

    def create_msidata_partial_spectra(self,
                                       spectra_length,
                                       len_global_mz,
                                       data_type='f',
                                       mzdata_type='f',
                                       chunks=None,
                                       compression=None,
                                       compression_opts=None,
                                       flush_io=True):
        """
        Create a new mass spectrometry imaging dataset for the given experiment written as a partial 3D
        cube of partial spectra.

        :param spectra_length: 2D boolean NumPy array used indicating for each (x,y) locations the length
                of the corresponding partial spectrum.
        :param len_global_mz: The total number of m/z values in the global m/z axis for the full 3D cube
        :param data_type: The dtype for the MSI dataset
        :param mzdata_type: The dtype for the mz dataset
        :param mzdata_type: numpy style datatype to be used for the mz data array.
        :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk sizes to
                       be used in x,y, and m/z explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data
                        has been written to file

        :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to be
            filled with data:

            * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * ``mz_index_dataset`` : h5py dataset with the mz_index values
            * ``mz_dataset`` : h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        :returns: The following already complete dataset

            * ``xy_index_dataset`` : This dataset indicates for each xy location at which index in ``data_dataset``
                            the corresponding spectrum starts. This dataset is needed to identify where spectra
                            need to be written to.
            * ``xy_index_end_dataset`` : This dataset indicates for each xy location at which index in
                            ``data_dataset`` the corresponding spectrum ends (exclusing the given value).
                            This dataset is needed to identify where spectra need to be written to.

        :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the
        # group
        data_group = omsi_file_msidata.__create_msidata_group__(parent_group=self.msidata_parent,
                                                                msidata_index=None,
                                                                flush_io=False)
        data_dataset, mz_index_dataset, mz_dataset, xy_index_start_dataset, xy_index_end_dataset = \
            omsi_file_msidata.__create_msidata_partial_spectra__(data_group=data_group,
                                                                 spectra_length=spectra_length,
                                                                 len_global_mz=len_global_mz,
                                                                 data_type=data_type,
                                                                 mzdata_type=mzdata_type,
                                                                 chunks=chunks,
                                                                 compression=compression,
                                                                 compression_opts=compression_opts)
        if flush_io:
            self.msidata_parent.file.flush()
        return data_dataset, mz_index_dataset, mz_dataset, xy_index_start_dataset, xy_index_end_dataset, data_group


class omsi_file_msidata(omsi_dependencies_manager,
                        omsi_methods_manager,
                        omsi_instrument_manager,
                        omsi_metadata_collection_manager,
                        omsi_file_common):
    """
    Interface for interacting with mass spectrometry imaging datasets stored in omis HDF5 files.
    The interface allows users to interact with the data as if it where a 3D cube even if data
    is missing. Full spectra may be missing in cases where only a region of interest in space
    has been imaged. Spectra may further be pre-processed so that each spectrum has only information
    about its peaks so that each spectrum has it's own mz-axis.

    To load data ue standard array syntax, e.g., [1,1,:] can be used to retrieve the spectrums at
    location (1,1).

    **Use of super():**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common`.
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the __init__ function calls
    super(...).__init__(manager_group) with a single  parameter indicating the
    parent group.


    **Current limitations:**

    * The estimates in def __best_dataset__(self,keys) are fairly crude at this point
    * The __getitem__ function for the partial_spectra case is not implemented yet.
    * The __setitem__ function for the partial spectra case is not implemented yet (Note, it \
      should also support dynamic expansion of the cube by adding previously missing spectra).
    * For the partial cube case, assignement using __setitem__ function is only supported to \
      valid spectra, i.e., spectra that were specified as occupied during the intital creation process.

    **Public object variables:**

    :ivar shape: Define the full 3D shape of the dataset (i.e., even if the data is stored in sparse manner)
    :ivar dtype: The numpy datatyp of the main MSI data. This is the same as dataset.dtype
    :ivar name: The name of the corresponding groupt in the HDF5 file. Used to generate hard-links to the group.
    :ivar format_type: Define according to which standard the data is stored in the file
    :ivar datasets: List of h5py objects containing possibly multiple different version of the same MSI data \
                    (spectra). There may be multiple versions stored with different layouts in order to \
                    optimize the selection process.
    :ivar mz: dataset with the global mz axis information. If prelaod_mz is set in the constructor, then this \
              is a numpy dataset with the preloaded data. Otherwise, this is the h5py dataset pointing to \
              the data on disk.
    :ivar xy_index: None if format_type is 'full_cube'. Otherwise, this is the 2D array indicating for \
                    each x/y location the index of the spectrum in dataset. If prelaod_xy_index is set \
                    in the constructor, then this is a numpy dataset with the preloaded data. Otherwise, \
                    this is the h5py dataset pointing to the data on disk. Negative (-1) entries indicate \
                    that no spectrum has been recored for the given pixel.
    :ivar inv_xy_index:  2D dataset with n rows and 2 columns indicating for each spectrum i the (x,y) \
                         pixel index the spectrum belongs to. This index is stored for convenience purposes but \
                         is not actually needed for data access.
    :ivar mz_index: None if format_type is not 'partial_spectra'. Otherwise this is a dataset of the same size \
                    as the spectra data stored in dataset. Each entry indicates an index into the mz dataset \
                    to determine the mz_data value for a spectrum. This means mz[ mx_index ] gives the true \
                    mz value.
    :ivar xy_index_end: None if format_type is not 'partial_spectra'. Otherwise this is a 2D array indicating \
                   for each x/y location the index where the given spectrum ends in the dataset. If \
                   prelaod_xy_index is set in the constructor, then this is a numpy dataset with the \
                   preloaded data. Otherwise, this is the h5py dataset pointing to the data on disk. \
                   Negative (-1) entries indicate that no spectrum has been recored for the given pixel.

    **Private object variables:**

    :ivar _data_group: Store the pointer to the HDF5 group with all the data
    :ivar _fill_xy: Define whether the data should be reconstructed as a full image cube Set using the \
                    set_fill_space function(..)
    :ivar _fill_mz: Define whether spectra should be remapped onto a global m/z axis. Set using the \
                    set_fill_spectra function(..)

    """
    @classmethod
    def __create_msidata_group__(cls,
                                 parent_group,
                                 msidata_index=None,
                                 flush_io=True):
        """
        Create a new empty group for storing MSI data.

        :param parent_group: The h5py.Group where the method object should be created
        :type parent_group: h5py.Group
        :param msidata_index: The integer index of the msidata to be created. Set to None if the
                         function should determine the index automatically.
        :type msidata_index: unsigned int or None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                       so that all data has been written to file
        :type flush_io: True

        :return: h5py group object for storing msi data
        """
        import time
        if msidata_index is None:
            msidata_index = omsi_file_common.get_num_items(parent_group,
                                                           omsi_format_msidata.data_groupname)
        data_group = parent_group.require_group(omsi_format_msidata.data_groupname + str(msidata_index))
        data_group.attrs[omsi_format_common.type_attribute] = "omsi_file_msidata"
        data_group.attrs[omsi_format_common.version_attribute] = omsi_format_msidata.current_version
        data_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            parent_group.file.flush()
        return data_group

    def __init__(self,
                 data_group,
                 fill_space=True,
                 fill_spectra=True,
                 preload_mz=False,
                 preload_xy_index=False):
        """
        Initialize the omsi_msidata object.

        The fill options are provided to enable a more convenient access to the data independent of how the
        data is stored in the file. If the fill options are enabled, then the user can interact with the
        data as if it where a 3D cube while missing is data is filled in by the given fill value.

        The prelaod options provided here refer to generally smaller parts of the data for which it may be
        more efficient to load the data and keep it around rather than doing repeated reads. If the object is
        used only for a single read and destroyed afterwards, then disabling the preload options may give a
        slight advantage but in most cases enabling the preload should be Ok (default).

        :param data_group: The h5py object for the group with the omsi_msidata.

        :param fill_space: Define whether the data should be padded in space (filled with 0's) when
                           accessing the data using [..] operator so that the data behaves like a 3D cube.

        :param fill_spectra: Define whether the spectra should completed by adding 0's so that all
                             spectra retrieved via the [..] operator so that always spectra of the full
                             length are returned. This option is provided to ease extension of the class
                             to cases where only partial spectra are stored in the file but is not used at this point.

        :param preload_mz: Should the data for the mz axis be kept in memory or loaded on the fly when needed.

        :param preload_xy_index: Should the xy index (if available) be preloaderd into memory or should the
                                 required data be loaded on the fly when needed.
        """
        super(omsi_file_msidata, self).__init__(data_group)
        # The following initialization are performed by the super call
        # self.managed_group = data_group
        # self.name = self.managed_group.name
        # self.methods_parent = self.managed_group
        # self.instrument_parent = self.managed_group
        # self.dependencies = ...
        self.shape = []
        self.dtype = None  # Which datatype does the main MSI data have
        self.format_type = omsi_format_msidata.format_types[
            str(self.managed_group .get(unicode(omsi_format_msidata.format_name))[0])]
        self.datasets = []
        self.mz = None
        self.xy_index = None
        self.inv_xy_index = None
        self.xy_index_end = None
        self.mz_index = None
        self._fill_xy = fill_space
        self._fill_mz = fill_spectra
        self.is_valid = False
        if self.format_type is not None:
            # Initalize the dataset
            self.datasets = [self.managed_group[unicode(x[0])] for x in self.managed_group.items()
                             if x[0].startswith(omsi_format_msidata.dataset_name)]
            self.dtype = self.datasets[0].dtype
            # Initalize the mz data
            if preload_mz:
                self.mz = self.managed_group[unicode(omsi_format_msidata.mzdata_name)][:]
            else:
                self.mz = self.managed_group[
                    unicode(omsi_format_msidata.mzdata_name)]
            # Initalize the data shape
            if self.format_type == omsi_format_msidata.format_types['full_cube']:
                self.shape = self.datasets[0].shape
                self.is_valid = True
            else:
                self.shape = self.managed_group[unicode(omsi_format_msidata_partial_cube.shape_name)][:]

            # We are done initializing all variables for the 'full_cube' format
            # Initialize the xy index, inv_xy_index and mz-index for the other
            # formats
            if self.format_type != omsi_format_msidata.format_types['full_cube']:
                if preload_xy_index:
                    self.xy_index = self.managed_group[unicode(omsi_format_msidata_partial_cube.xy_index_name)][:]
                    self.inv_xy_index = \
                        self.managed_group[unicode(omsi_format_msidata_partial_cube.inv_xy_index_name)][:]
                    if self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                        self.xy_index_end = \
                            self.managed_group[unicode(omsi_format_msidata_partial_spectra.xy_index_end_name)][:]
                        self.mz_index = self.managed_group[
                            unicode(omsi_format_msidata_partial_spectra.mz_index_name)]
                        # We have successfully loaded all data for the
                        # partial_spectra case
                        self.is_valid = True
                    elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
                        # We have successfully loaded all data for the
                        # partial_cube case
                        self.is_valid = True
                else:
                    self.xy_index = self.managed_group[unicode(omsi_format_msidata_partial_cube.xy_index_name)]
                    self.inv_xy_index = self.managed_group[unicode(omsi_format_msidata_partial_cube.inv_xy_index_name)]
                    if self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                        self.xy_index_end = \
                            self.managed_group[unicode(omsi_format_msidata_partial_spectra.xy_index_end_name)]
                        self.mz_index = self.managed_group[unicode(omsi_format_msidata_partial_spectra.mz_index_name)]
                        # We have successfully loaded all data for the
                        # partial_spectra case
                        self.is_valid = True
                    elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
                        # We have successfully loaded all data for the
                        # partial_cube case
                        self.is_valid = True

    def get_h5py_datasets(self, index=0):
        """
        Get the h5py dataset object for the given dataset.

        :param index: The index of the dataset.

        :returns: h5py object for the requested dataset.

        :raises: and Index error is generated in case an invalid index is given.
        """
        return self.datasets[index]

    def get_h5py_mzdata(self):
        """
        Get the h5py object for the mz datasets.

        :returns: h5py object of the requested mz dataset.
        """
        return self.managed_group[unicode(omsi_format_msidata.mzdata_name)]

    def __setitem__(self, key, value):
        """
        The __getitem__ function is used in python to implement the [..] operator for setting data values.
        This function allows the user to write the data as if it were a 3D cube using array notation [:,:,:]
        even if the data may be stored in a different fashion in the file. If less than three
        selections are specified then the remaining dimensions are assumed to be ":", i.e., all.

        :param key: Three elements of type index, list or slice defining the data selection
        """
        # The object is not fully initalized
        if self.managed_group is None:
            raise ValueError("The msidata object has not been initialized.")

        # Complete the input selection if it is only partially specified. In this way we can
        # assume in the following code that we always have three key selection
        # parameters
        if not isinstance(key, tuple):
            key = (key, slice(None), slice(None))
        elif len(key) == 1:
            key = (key, slice(None), slice(None))
        elif len(key) == 2:
            key = (key[0], key[1], slice(None))
        elif len(key) != 3:
            raise ValueError("Invalid selection")

        # Check the data format and call the approbriate getitem function
        if self.format_type == omsi_format_msidata.format_types['full_cube']:
            return self.__setitem_fullcube__(key, value)
        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
            return self.__setitem_partialcube__(key, value)
        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
            return self.__setitem_partialspectra__(key, value)

    def __setitem_fullcube__(self, key, value):
        """
        Private helper function used in case that the data is stored as a full cube

        :param key: Three elements of type index, list or slice defining the data selection
        """
        # Update the data in all available version of the dataset
        for dset in self.datasets:
            dset[key] = value

    def __setitem_partialcube__(self, key, value):
        """
        Private helper function used in case that the data is stored as a partial cube

        :param key: Three elements of type index, list or slice defining the data selection
        """
        # Check which elements need to be updated
        indexlist = self.xy_index[key[0], key[1]]
        # Check whether only valid spectra are selected
        validselectsize = (indexlist >= 0).sum()
        try:
            totalsize = indexlist.size
            indexlist = indexlist.reshape(indexlist.size)
        except:
            totalsize = 1
        if validselectsize != totalsize:
            raise KeyError(
                "Assignment to uninitalized spectra is currently not supported")

        # Update the data in all available version of the dataset
        for dset in self.datasets:
            dset[indexlist.tolist(), key[2]] = value.reshape(totalsize, self.__num_elements__(key[2]))

    def __setitem_partialspectra__(self, key, value):
        """
        Private helper function used in case that the data is stored as a partial cube of partial spectra.

        :param key: Three elements of type index, list or slice defining the data selection
        """
        raise NotImplementedError("Assignement of partial spectra datasets is currently not implemented")

    def __getitem__(self, key):
        """
        The __getitem__ function is used in python to implement the [..] operator. This function
        allows the user to access the data as if it were a 3D cube using array notation [:,:,:]
        even if the data may be stored in a different fashion in the file. If less than three
        selections are specified then the remaining dimensions are assumed to be ":", i.e., all.

        :param key: Three elements of type index, list or slice defining the data selection
        """

        # The object is not fully initalized
        if self.managed_group is None:
            raise ValueError("The msidata object has not been initialized.")

        # Complete the input selection if it is only partially specified. In this way we can
        # assume in the following code that we always have three key selection parameters
        if isinstance(key, basestring):
            return super(omsi_file_msidata, self).__getitem__(key)
        elif not isinstance(key, tuple):
            key = (key, slice(None), slice(None))
        else:
            if len(key) == 1:
                key = (key, slice(None), slice(None))
            elif len(key) == 2:
                key = (key[0], key[1], slice(None))
            elif len(key) != 3:
                raise ValueError("Invalid selection")

        # Check the data format and call the appropriate getitem function
        if self.format_type == omsi_format_msidata.format_types['full_cube']:
            return self.__getitem_fullcube__(key)
        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
            return self.__getitem_partialcube__(key)
        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
            return self.__getitem_partialspectra__(key)

    def __getitem_fullcube__(self, key):
        """
        Private helper function used in case that the data is stored as a full cube

        :param key: Three elements of type index, list or slice defining the data selection
        """
        # Get the dataset that is best suited for the selection
        dset = self.__best_dataset__(key)
        # Access data using h5py
        return dset[key]

    def __getitem_partialcube__(self, key):
        """
        Private helper function used in case that the data is stored as a partial cube.

        :param key: Three elements of type index, list or slice defining the data selection
        """
        # Get the dataset that is best suited for the selection
        dset = self.__best_dataset__(key)
        # Compute the list elements to be loaded in (x,y)
        indexlist = self.xy_index[key[0], key[1]]
        if indexlist.size > 1:
            cleanindexlist = indexlist[indexlist >= 0]
        else:
            indexlist = np.asarray([indexlist])
            if indexlist[0] < 0:
                cleanindexlist = np.empty(0, dtype=int)
            else:
                cleanindexlist = np.asarray([indexlist], dtype=int)
        cleanindexlist = cleanindexlist.reshape(cleanindexlist.size)
        # Load the data
        if cleanindexlist.size > 0:
            data = dset[cleanindexlist, key[2]]
        else:
            data = np.empty(0, dtype=dset.dtype)

        # Check if we need to complete the data with 0's. If the selection defines a range in
        # x and y (i.e we select a full rectangular region in space), then reconstruct the
        # full region by filling any missing data with 0's.
        if self._fill_xy:
            # Compute how many elements we need to retrieve in the mz-axis
            mzsize = self.__num_elements__(key[2])
            renshape = indexlist.shape + (mzsize, )
            filldata = np.zeros(renshape, dtype=dset.dtype)
            filldata[indexlist >= 0] = data.reshape(
                (cleanindexlist.size, self.__num_elements__(key[2])))
            return filldata
            # if len( indexlist.shape )==2 :
            # xsize  = indexlist.shape[0]
            # ysize  = indexlist.shape[1]
            # filldata = np.zeros( shape=(xsize ,  ysize , mzsize) , dtype=d.dtype )
            # if not isinstance( key[0] , list ) and not isinstance( key[1] , list  ) :
            #    filldata[ self.inv_xy_index[cleanindexlist,0]-self.__offset__(key[0]) , self.inv_xy_index[cleanindexlist,1]-self.__offset__(key[1])  ] = data.reshape((cleanindexlist.size, self.__num_elements__(key[2]) ))
            #    return filldata
            # if isinstance( key[0] , list ) :
        else:
            # Otherwise return a list of numpy arrays with the requested data
            # completed with None objects for missing data
            filldata = [None] * indexlist.size
            for i in xrange(0, indexlist.size):
                if indexlist[i] >= 0:
                    filldata[i] = data[cleanindexlist[i], :]
            return filldata

    def __getitem_partialspectra__(self, key):
        """
        Private helper function used in case that the data is stored as full or partial cube with partial spectra.

        :param key: Three elements of type index, list or slice defining the data selection
        """
        raise NotImplementedError(
            "Selection for partial spectra datasets is currently not implemented. See __getitem_partialspectra__(...)")

        # Get the dataset that is best suited for the selection
        # d = self.__best_dataset__(key)
        # Determine the data values to be loaded in xy
        # indexList_start = self.xy_index[ key[0] , key[1] ]
        # indexList_end   = self.xy_index_end[ key[0] , key[1] ]
        # index_select = (indexList_start>=0)
        # if not isinstance( indexList_start , int ) :
        #   cleanIndexList_start = indexList_start[ index_select ]
        #   cleanIndexList_end = indexList_end[ index_select ]
        # else :
        #    if indexList_start < 0 :
        #        cleanIndexList_start = np.empty(0, dtype=int)
        #        cleanIndexList_end   = np.empty(0, dtype=int)
        #    else :
        #       cleanIndexList_start = np.asarray( [ indexList_start ] , dtype=int )
        #       cleanIndexList_end   = np.asarray( [ indexList_end ]   , dtype=int )

        # Load the full specta for each x/y location
        # mz_data = [[]] * cleanIndexList.shape[0]
        # for i in xrange(0,cleanIndexList.shape[0]) :
        #     mz_data[i] = self.mz_index[ cleanIndexList_start[i]:cleanIndexList_end[i] ]

        # Determine which parts of the m/z data we actually need
        # if isinstance(key[2] , int ) :
        #   mz_select = [ i==key[2] for i in mz_data  ]
        # if isinstance(key[2] , slice ) :
        #   if key[2].step>1:
        #       raise NotImplementedError("Selections with a stepping are currently not supported by the interface for the m/z dimension and particle_spectra data")
        #   Select full specta
        #   if (key[2].start is None or key[2].start == 0) and \
        #       (key[2].end is None or key[2].end == mz_indexList.shape[1]) and \
        #       (key[2].step is None or key[2].step == 1 ) :
        #           mz_select = [ np.ones(dtype='bool', shape=i.shape ) for i in mz_data ]
        #   Select an mz-range of the spectra
        #   else :
        #       mz_select = [ np.logical_and( i>=key[2].start , i<key[2].end ) for i in mz_data  ]
        #       #mz_select = np.where( np.logical_and( mz_indexList>=key[2].start , mz_indexList<key[2].end )  , mz_indexList , -1 )
        # elif isinstance(key[2] , list ) :
        #   Treat the list the same way as an index selection if we only have one entry in the list
        #   if len( key[2] ) == 1 :
        #       mz_select = [ i==key[2][0] for i in mz_data  ]
        #   Check if the list defines a continues selection
        #   for i in xrange(1,len(key[2])) :
        #       if key[2][i] != (key[2][i-1]+1) :
        #           raise NotImplementedError("Selections using discontinoues lists are currently not supported by the interface for the m/z dimension and particle_spectra data.")
        #           Implementing discontinoues lists for selection will require a different treatment also when filling the data
        # Treat the list as a continues slice-based selection.
        # mz_select = [ np.logical_and( i>=key[2][0] , i<key[2][-1] ) for i in mz_data  ]

        # Load the requested spectra
        # if self._fill_mz and self._fill_xy and (len( indexList_start.shape )==2)  :
        #   Remap both the spatial and m/z coordinates so that the data behaves completely like a 3D cube
        #   Compute how many elements we need to retrieve in the mz-axis
        #   mzsize = self.__num_elements__(key[2])
        #   ysize  = indexList_start.shape[0]
        #   xsize  = indexList_start.shape[1]
        #   filldata = np.zeros( shape=(xsize, ysize, mzsize) , dtype = d.dtype )
        #   mzi = 0
        #   for xi in xrange(0,xsize) :
        #       for yi in xrange(0,ysize) :
        #           if indexList_start[xi,yi] >= 0 :
        #               filldata[xi,yi,mz_select[mzi]] = d[ indexList_start[i]:indexList_end[i] ][ mz_select[mzi] ]
        #               mzi = mzi+1
        # elif self._fill_mz: #This includes the case where both fill_mz and fill_xy are set but fill_xy is not needed.
        #   mzsize = self.__num_elements__(key[2])
        #   filldata = [np.zeros( mzsize) ] * len(mz_select)
        #    for i in xrange(0, len(mz_select) ):
        #       filldata[i][mz_select[i]] =  d[ cleanIndexList_start[i]:cleanIndexList_end[i] ][ mz_select[i] ]
        #   return filldata
        # elif self._fill_xy: #Only fill_xy is set but not the fill_mz
        #    Remap only the spatial coordinates
        #   raise NotImplementedError("fill_xy without fill_mz set is currently not supported by the API")
        # else : #This is the case where neither fill_mz nor fill_xy are set
        #    Just load the data and return both the data and mz indicies
        #   filldata     =  [[]] * cleanIndexList.shape[0]
        #   fillmzdata   =  [[]] * cleanIndexList.shape[0]
        #   for i in xrange(0,cleanIndexList.shape[0]) :
        #       filldata[i]   = d[ cleanIndexList_start[i]:cleanIndexList_end[i] ][mz_select[i]]
        #       fillmzdata[i] = mz_data[i][mz_select[i]]
        #   return filldata, fillmzdata

    def set_fill_space(self, fill_space):
        """
        Define whether spatial selection should be filled with 0's to retrieve full image slices

        :param fill_space: Boolean indicating whether images should be filled with 0's
        """
        self._fill_xy = fill_space

    def set_fill_spectra(self, fill_spectra):
        """
        Define whether spectra should be filled with 0's to map them to the global mz axis when retrieved.

        :param fill_spectra: Define whether m/z values should be filled with 0's.
        """
        self._fill_mz = fill_spectra

    @classmethod
    def __create_msidata_full_cube__(cls,
                                     data_group,
                                     data_shape,
                                     data_type='f',
                                     mzdata_type='f',
                                     chunks=None,
                                     compression=None,
                                     compression_opts=None,
                                     dependencies_data_list=None):
        """
        Create a new msi data group with all necessary datasets for storing a full 3D cube.

        NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
        corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

        Required input parameters

        :param data_group: The h5py.Group for which the dataset should be initalized.
        :param data_shape: The 3D shape of the MSI dataset in x,y, and m/z
        :param data_type: The dtype for the MSI dataset
        :param mzdata_type: The dtype for the mz dataset

        Optional layout optimization parameters (these refer to the main MSI dataset only)

        :param chunks:  Specify whether chunking should be used (True,False), or specify the chunk sizes
                        to be used explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).

        Other optional parameters

        :param dependencies_data_list: List of dependency_dict objects to be stored as dependencies.
                                       Default is empty list []

        :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to
                  be filled with data:

            * data_dataset : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * mz_dataset : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        """
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(omsi_format_msidata.format_name), shape=(1,),
                                                    dtype=omsi_format_common.str_type)
        format_dataset[0] = 'full_cube'
        # mz data
        mz_dataset = data_group.require_dataset(name=omsi_format_msidata.mzdata_name,
                                                shape=(data_shape[2], ),
                                                dtype=mzdata_type)
        # Main MSI dataset
        if compression is None:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=data_shape, dtype=data_type,
                                                      chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=data_shape,
                                                      dtype=data_type,
                                                      chunks=chunks,
                                                      compression=compression)
        else:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=data_shape,
                                                      dtype=data_type,
                                                      chunks=chunks,
                                                      compression=compression,
                                                      compression_opts=compression_opts)
        # Create any dependencies
        _ = omsi_file_dependencies.__create__(parent_group=data_group,
                                              dependencies_data_list=dependencies_data_list)
        # Return the datasets that need to be written
        return data_dataset, mz_dataset

    @classmethod
    def __create_msidata_partial_cube__(cls,
                                        data_group,
                                        data_shape,
                                        mask,
                                        data_type='f',
                                        mzdata_type='f',
                                        chunks=None,
                                        compression=None,
                                        compression_opts=None,
                                        dependencies_data_list=None):
        """
        Create a new msi data group with all necessay datasets for storing a full 3D cube

        NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
        corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

        Required input parameters

        :param data_group: The hHDF5 group for which the dataset should be initalized.
        :param data_shape: The 3D shape of the MSI dataset in x,y, and m/z
        :param mask: 2D boolean NumPy array used as mask to indicate which (x,y) locations have spectra
                     associated with them.
        :param data_type: The dtype for the MSI dataset
        :param mzdata_type: The dtype for the mz dataset

        Optional layout optimization parameters (these refer to the main MSI dataset only)

        :param chunks:  Specify whether chunking should be used (True,False), or specify the
                        chunk sizes to be used explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).

        Other optional parameters

        :param dependencies_data_list: List of dependency_dict objects to be stored as dependencies.
                                       Default is empty None which is mapped to an empty list []

        :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to
                  be filled with data:

            * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        :returns: The following already complete dataset

            * ``xy_index_dataset`` : This dataset indicates for each xy location to which index in ``data_dataset`` \
                                     the location corresponds to. This dataset is needed to identify where spectra \
                                     need to be written to.
            * ``inv_xy_index_dataset`` : This datasets indicates for each spectrum its xy locating. The dataset is \
                                     of the from n x 2 with inv_xy_index_dataset[:,0] being the x indices and \
                                     inv_xy_index_dataset[:,1] being the y indicies.

        """
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(
            omsi_format_msidata.format_name), shape=(1,), dtype=omsi_format_common.str_type)
        format_dataset[0] = 'partial_cube'
        # mz data
        mz_dataset = data_group.require_dataset(
            name=omsi_format_msidata.mzdata_name, shape=(data_shape[2], ), dtype=mzdata_type)

        # Main MSI dataset
        xyshape = data_shape[0] * data_shape[1]
        numspectra = mask.sum()
        mzshape = data_shape[2]
        if isinstance(chunks, tuple):
            chunks = (chunks[0] * chunks[1], chunks[2])
        if compression is None:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=(numspectra, mzshape),
                                                      dtype=data_type,
                                                      chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=(numspectra, mzshape),
                                                      dtype=data_type, chunks=chunks,
                                                      compression=compression)
        else:
            data_dataset = data_group.require_dataset(name=(omsi_format_msidata.dataset_name + "0"),
                                                      shape=(numspectra, mzshape),
                                                      dtype=data_type,
                                                      chunks=chunks,
                                                      compression=compression,
                                                      compression_opts=compression_opts)

        # Determine the minimal dtype for the xy index
        xy_index_dtype = 'int16'
        if xyshape < np.iinfo('int16').max:
            xy_index_dtype = 'int16'
        elif xyshape < np.iinfo('int32').max:
            xy_index_dtype = 'int32'
        elif xyshape < np.iinfo('int64').max:
            xy_index_dtype = 'int64'

        # Create the xy_index_dataset
        xy_index_dataset = data_group.require_dataset(name=omsi_format_msidata_partial_cube.xy_index_name,
                                                      shape=(data_shape[0], data_shape[1]),
                                                      dtype=xy_index_dtype)
        xy_index_dataset[:, :] = -1
        xy_index_dataset[mask] = np.arange(0, numspectra)

        # Create the inverse inv_xy_index dataset
        inv_xy_index_dataset = data_group.require_dataset(name=omsi_format_msidata_partial_cube.inv_xy_index_name,
                                                          shape=(numspectra, 2),
                                                          dtype=xy_index_dtype)
        for x_index in xrange(0, data_shape[0]):
            for y_index in xrange(0, data_shape[1]):
                if xy_index_dataset[x_index, y_index] >= 0:
                    inv_xy_index_dataset[xy_index_dataset[x_index, y_index], :] = np.asarray([x_index, y_index])

        # Create the shape dataset
        shape_dataset = data_group.require_dataset(name=omsi_format_msidata_partial_cube.shape_name,
                                                   shape=(3, ),
                                                   dtype='uint32')
        shape_dataset[:] = np.asarray(data_shape, dtype='uint32')

        # Create any dependencies
        _ = omsi_file_dependencies.__create__(parent_group=data_group,
                                              dependencies_data_list=dependencies_data_list)

        # Return the datasets that need to be written
        return data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset

    @classmethod
    def __create_msidata_partial_spectra__(cls,
                                           data_group,
                                           spectra_length,
                                           len_global_mz,
                                           data_type='f',
                                           mzdata_type='f',
                                           chunks=None,
                                           compression=None,
                                           compression_opts=None,
                                           dependencies_data_list=None):
        """
        Create a new msi data group with all necessay datasets for storing a full or partial cube of partial spectra

        NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
        corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

        Required input parameters

        :param data_group: The hDF5 group for which the dataset should be initalized.
        :param spectra_length: 2D boolean NumPy array used indicating for each (x,y) locations the length of
                               the corresponding partial spectrum.
        :param len_global_mz: The total number of m/z values in the global m/z axis for the full 3D cube
        :param data_type: The dtype for the MSI dataset
        :param mzdata_type: The dtype for the mz dataset

        Optional layout optimization parameters (these refer to the main MSI dataset only)

        :param chunks:  Specify whether chunking should be used (True,False), or specify the
                        chunk sizes to be used explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                            Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                            For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                            a number between zero and nine (inclusive) to indicate the tradeoff between speed
                            and compression ratio (zero is fastest, nine is best ratio).

        Other optional parameters

        :param dependencies_data_list: List of dependency_dict objects to be stored as dependencies.
                                       Default is empty list []

        :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order
                  to be filled with data:

            * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
            * ``mz_index_dataset`` : The
            * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        :returns: The following already complete dataset

            * ``xy_index_dataset`` : This dataset indicates for each xy location at which index in ``data_dataset``
                                    the corresponding spectrum starts. This dataset is needed to identify where
                                    spectra need to be written to.
            * ``xy_index_end_dataset`` : This dataset indicates for each xy location at which index in
                            ``data_dataset`` the corresponding spectrum ends (exclusing the given value).
                            This dataset is needed to identify where spectra need to be written to.

        """
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(
            omsi_format_msidata.format_name), shape=(1,), dtype=omsi_format_common.str_type)
        format_dataset[0] = 'partial_spectra'
        # mz data
        mz_dataset = data_group.require_dataset(
            name=omsi_format_msidata.mzdata_name, shape=(len_global_mz, ), dtype=mzdata_type)

        # Main MSI dataset
        mask = (spectra_length > 0)
        numspectra = mask.sum()
        total_len_spectra = spectra_length[mask].sum()
        if isinstance(chunks, tuple):
            chunks = (chunks[0] * chunks[1] * chunks[2], )
        if compression is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks,
                compression=compression)
        else:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        # Determine the minimal dtype for the mz index
        mz_index_dtype = 'uint16'
        if len_global_mz < np.iinfo('uint16').max:
            mz_index_dtype = 'uint16'
        elif len_global_mz < np.iinfo('uint32').max:
            mz_index_dtype = 'uint32'
        elif len_global_mz < np.iinfo('uint64').max:
            mz_index_dtype = 'uint64'

        # Create the mz-index dataset
        if compression is None:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks)
        elif compression_opts is None:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks,
                compression=compression)
        else:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        # Determine the minimal dtype for the xy index
        xy_index_dtype = 'uint16'
        xyshape = spectra_length.shape[0] * spectra_length.shape[1]
        if xyshape < np.iinfo('uint16').max:
            xy_index_dtype = 'uint16'
        elif xyshape < np.iinfo('uint32').max:
            xy_index_dtype = 'uint32'
        elif xyshape < np.iinfo('uint64').max:
            xy_index_dtype = 'uint64'

        # Create the xy index start dataset
        xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.xy_index_name,
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        xy_index_dataset[:] = -1
        xy_index_dataset[mask] = np.insert(
            np.cumsum(spectra_length[mask]), 0, 0)[0:-1]

        # Create the xy index end dataset
        xy_index_end_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_spectra.xy_index_end_name,
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        xy_index_end_dataset[:] = -1
        xy_index_end_dataset[mask] = np.insert(
            np.cumsum(spectra_length[mask]), 0, 0)[1:]

        # Create the inverse inv_xy_index dataset
        inv_xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.inv_xy_index_name,
            shape=(numspectra, 2),
            dtype=xy_index_dtype)
        temp_count_xy_index = np.zeros(
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        temp_count_xy_index[:] = -1
        temp_count_xy_index[mask] = np.insert(
            np.cumsum(np.ones(shape=spectra_length.shape, dtype='uint16')[mask]), 0, 0)[0:-1]
        for x_index in xrange(0, spectra_length.shape[0]):
            for y_index in xrange(0, spectra_length.shape[1]):
                if temp_count_xy_index[x_index, y_index] >= 0:
                    inv_xy_index_dataset[temp_count_xy_index[x_index, y_index], :] = np.asarray([x_index, y_index])

        # Create the shape dataset
        shape_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.shape_name,
            shape=(3, ),
            dtype='uint32')
        shape_dataset[:] = np.asarray(mask.shape, dtype='uint32')

        # Create any dependencies
        _ = omsi_file_dependencies.__create__(parent_group=data_group,
                                              dependencies_data_list=dependencies_data_list)

        # Return the datasets that need to be written
        return data_dataset, mz_index_dataset, mz_dataset, xy_index_dataset, xy_index_end_dataset

    def create_optimized_chunking(self,
                                  chunks=None,
                                  compression=None,
                                  compression_opts=None,
                                  copy_data=True,
                                  print_status=False,
                                  flush_io=True):
        """
        Helper function to allow one to create optimized copies of the dataset with different internal data
        layouts to speed up selections. The function expects that the original data has already been written
        to the data group. The function takes

        :param chunks:  Specify whether chunking should be used (True,False), or specify the chunk sizes
                       to be used explicitly.
        :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                          Can also use an integer in range(10) indicating gzip.
        :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                          For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                          a number between zero and nine (inclusive) to indicate the tradeoff between speed
                          and compression ratio (zero is fastest, nine is best ratio).
        :param copy_data: Should the MSI data be copied by this function to the new dataset or not. If False, then
                         it is up to the user of the function to copy the appropriate data into the returned h5py
                         dataset (not recommended but may be useful for performance optimization).
        :param print_status: Should the function print the status of the conversion process to the command line?
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has
                       been written to file

        :returns: h5py dataset with the new copy of the data
        """
        if len(self.datasets) == 0:
            raise ValueError("No datasets are currently stored for the dataset group that could be replicated")

        # Get the donor dataset. Try tp get the main dataset first. If it is
        # not found take the first one from the list.
        try:
            dset = self.managed_group .get(unicode(omsi_format_msidata.dataset_name + "0"))
        except:
            dset = self.datasets[0]

        # Check and adjust the chunking
        if isinstance(chunks, tuple):
            if self.format_type == omsi_format_msidata.format_types['partial_cube']:
                chunks = (chunks[0] * chunks[1], chunks[2])
                # Correct the chunking if it is larger than the actual data
                if chunks[0] > dset.shape[0]:
                    chunks = (dset.shape[0], chunks[1])
                if chunks[1] > dset.shape[1]:
                    chunks = (chunks[0], dset.shape[1])
            elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                chunks = (chunks[0] * chunks[1] * chunks[2], )
                # Correct the chunking if the chunk is larger than the actual
                # data
                if chunks[0] > dset.shape[0]:
                    chunks = (dset.shape[0], )

        # Generate the new dataset in the HDF5 file
        new_data_name = omsi_format_msidata.dataset_name + str(len(self.datasets))
        if compression is None:
            data_dataset = self.managed_group.require_dataset(
                name=new_data_name, shape=dset.shape, dtype=self.dtype, chunks=chunks)
        elif compression_opts is None:
            data_dataset = self.managed_group.require_dataset(
                name=new_data_name, shape=dset.shape, dtype=self.dtype, chunks=chunks, compression=compression)
        else:
            data_dataset = self.managed_group.require_dataset(
                name=new_data_name,
                shape=dset.shape,
                dtype=self.dtype,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        if copy_data:
            self.copy_dataset(source=dset,
                              destination=data_dataset,
                              print_status=print_status)

        if flush_io:
            self.managed_group.file.flush()

        # Add the dataset to the list of datasets
        self.datasets.append(data_dataset)

        # Return the dataset. This is needed in case the caller decided to set
        # copy_data to False and copy the data themselfs
        return data_dataset

    def copy_dataset(self,
                     source,
                     destination,
                     print_status=False):
        """
        Helper function used to copy a source msi dataset one chunk at a time to the destination dataset.
        The data copy is done one destination chunk at a time to achieve chunk-aligned write.

        :param source: The source h5py dataset
        :param destination: The h5py desitnation h5py dataset.
        :param print_status: Should the function print the status of the conversion process to the command line?

        """
        # Write the data from the donor to the target dataset
        if print_status:
            import sys
        if self.format_type == omsi_format_msidata.format_types['full_cube']:

            chunks = destination.chunks
            numchunksx = int(math.ceil(float(source.shape[0]) / float(chunks[0])))
            numchunksy = int(math.ceil(float(source.shape[1]) / float(chunks[1])))
            numchunksz = int(math.ceil(float(source.shape[2]) / float(chunks[2])))
            numchunks = numchunksx * numchunksy * numchunksz
            itertest = 0
            for x_chunk_index in xrange(0, numchunksx):
                xstart = x_chunk_index * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                for y_chunk_index in xrange(0, numchunksy):
                    ystart = y_chunk_index * chunks[1]
                    yend = min(ystart + chunks[1], source.shape[1])
                    for z_chunk_index in xrange(0, numchunksz):
                        zstart = z_chunk_index * chunks[2]
                        zend = min(zstart + chunks[2], source.shape[2])
                        destination[xstart:xend, ystart:yend, zstart:zend] = \
                            source[xstart:xend, ystart:yend, zstart:zend]
                        itertest += 1
                        if print_status:
                            sys.stdout.write("[" +
                                             str(int(100. * float(itertest) / float(numchunks))) +
                                             "%]" + "\r")
                            sys.stdout.flush()

        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:

            chunks = destination.chunks
            numchunksx = int(math.ceil(float(source.shape[0]) / float(chunks[0])))
            numchunksy = int(math.ceil(float(source.shape[1]) / float(chunks[1])))
            numchunks = numchunksx * numchunksy
            itertest = 0
            for x_chunk_index in xrange(0, numchunksx):
                xstart = x_chunk_index * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                for y_chunk_index in xrange(0, numchunksy):
                    ystart = y_chunk_index * chunks[1]
                    yend = min(ystart + chunks[1], source.shape[1])
                    destination[xstart:xend, ystart:yend] = source[xstart:xend, ystart:yend]
                    itertest += 1
                    if print_status:
                        sys.stdout.write("[" +
                                         str(int(100. * float(itertest) / float(numchunks))) +
                                         "%]" + "\r")
                        sys.stdout.flush()

        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:

            chunks = destination.chunks
            numchunksx = int(
                math.ceil(float(source.shape[0]) / float(destination.chunks[0])))
            numchunks = numchunksx
            itertest = 0
            for x_chunk_index in xrange(0, numchunksx):
                xstart = x_chunk_index * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                destination[xstart:xend] = source[xstart:xend]
                itertest += 1
                if print_status:
                    sys.stdout.write("[" +
                                     str(int(100. * float(itertest) / float(numchunks))) +
                                     "%]" + "\r")
                    sys.stdout.flush()

    def __best_dataset__(self, keys, print_info=False):
        """
        Compute the index of the dataset that is best suited for executing the given selection

        :param keys: List of three keys indicting which elements should be selected from the MSI dataset.

        :returns: Integer indicating the index of the dataset to be used for the given selection.
        """
        # If there is only one version of the dataset then use that one
        if len(self.datasets) == 1:
            return self.datasets[0]
        # Get the basic properties of the selection
        sizex = self.__num_elements__(keys[0])
        sizey = self.__num_elements__(keys[1])
        sizez = self.__num_elements__(keys[2])
        offsetx = self.__offset__(keys[0])
        offsety = self.__offset__(keys[1])
        offsetz = self.__offset__(keys[2])

        # Try to determine how many chunks need to be touched in order to fullfill the given selection
        # Currently this is just an estimate assuming that we have spatially consecutive selections,
        # i.e., if we have a stepping >1 in a slice or lists of non-neighboring elements, then this
        # estimate will be wrong. This may result in a suggestion of a sub-optimal dataset with
        # respect to performance.
        # For the particle_cube case the estimate is also not accurate for spatial selections since
        # both x and y are linearized in this case. However, the estimate should still distinguish
        # properly between selection in space vs. spectra.
        # For the partial spectra case, currently no chunking optimization is
        # currently available.
        suggestion = self.datasets[0]
        if self.format_type == omsi_format_msidata.format_types['full_cube']:

            currenttouch = -1
            currentload = -1
            for dset in self.datasets:
                chunking = dset.chunks
                chunkload = chunking[0] * chunking[1] * chunking[2]
                # Account for the misalignment of the selection with chunk
                # boundaries
                adjusted_offset_x = offsetx - (math.floor(float(offsetx) / chunking[0]) * chunking[0])
                adjusted_offset_y = offsety - (math.floor(float(offsety) / chunking[1]) * chunking[1])
                adjusted_offset_z = offsetz - (math.floor(float(offsetz) / chunking[2]) * chunking[2])
                # Determine how many chunks need to be touched, assuming that
                # we have continues selection.
                touchx = math.ceil(float(sizex + adjusted_offset_x) / chunking[0])
                touchy = math.ceil(float(sizey + adjusted_offset_y) / chunking[1])
                touchz = math.ceil(float(sizez + adjusted_offset_z) / chunking[2])
                touchtotal = touchx * touchy * touchz
                # Determine the amount of data that need to be loaded
                loadtotal = touchtotal * chunkload
                # If less chunks have to be touched for this dataset then use
                # that one instead
                if loadtotal < currentload or currentload < 0:
                    suggestion = dset
                    currenttouch = touchtotal
                    currentload = loadtotal

        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:

            currenttouch = -1
            currentload = -1
            for dset in self.datasets:
                chunking = dset.chunks
                chunkload = chunking[0] * chunking[1]
                # Account for the misalignment of the selection with chunk boundaries
                adjusted_offset_x = offsetx - (math.floor(float(offsetx) / chunking[0]) * chunking[0])
                adjusted_offset_y = offsety - (math.floor(float(offsety) / chunking[0]) * chunking[0])
                adjusted_offset_z = offsetz - (math.floor(float(offsetz) / chunking[1]) * chunking[1])
                # Determine how many chunks need to be touched, assuming that
                # we have continues selection.
                touchx = math.ceil(float(sizex + adjusted_offset_x) / chunking[0])
                touchy = math.ceil(float(sizey + adjusted_offset_y) / chunking[0])
                touchz = math.ceil(float(sizez + adjusted_offset_z) / chunking[1])
                touchtotal = touchx * touchy * touchz
                # Determine the amount of data that need to be loaded
                loadtotal = touchtotal * chunkload
                # If less chunks have to be touched for this dataset then use that one instead
                if loadtotal < currentload or currentload < 0:
                    suggestion = dset
                    currenttouch = touchtotal
                    currentload = loadtotal

        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:

            # ToDo currently no optimization is implemented for the partial
            # spectra case
            pass

        if print_info:
            print "Selection: " + str(keys)
            print "Suggest: " + str(suggestion.chunks)
            print "Alternatives: " + str([dset.chunks for dset in self.datasets])

        return suggestion

    @staticmethod
    def __offset__(key):
        """
        Determine the start offset of the given single key

        :param key: List, slice or interger indicating a single selection

        :returns: Integer indicating the lower bound of the selection defined by the key.
        """
        if isinstance(key, int):
            return key
        elif isinstance(key, list):
            return min(key)
        elif isinstance(key, slice):
            start = key.start
            if start is None:
                start = 0
            return start
        else:
            return 0

    def __num_elements__(self, key):
        """
        Compute the number of elements selected by a given single selection key

        :param key: List, slice or integer indicating a single selection

        :returns: The number of elements selected by the given key.
        """
        if isinstance(key, int):
            return 1
        elif isinstance(key, list) or isinstance(key, np.ndarray):
            return len(key)
        elif isinstance(key, slice):
            start = key.start
            if start is None:
                start = 0
            end = key.stop
            if end is None:
                end = self.shape[2]
            step = key.step
            if key.step is None:
                step = 1
            return math.ceil(float(end - start) / float(step))
        else:
            raise ValueError("Unexpected key value " + str(key) + " " + str(type(key)))

