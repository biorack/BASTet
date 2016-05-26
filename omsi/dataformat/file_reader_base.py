"""
Module for base classes for implementation and integration of third-party file readers.

ToDo:

* get_number_of_regions(...) should be updated to return a list of regions, one per dataset
* Need to add base class for multi dataset formats
* Need to add base class for multi dataset+region formats
* Need to implement new file format for combined raw data file (ie., multiple raw files in one folder).

"""
from omsi.datastructures.metadata.metadata_data import metadata_dict

class file_reader_base(object):
    """
    Base-class used to define the basic interface that file-readers
    for a new format need to implement.

    **__init__ interface:**

    To avoid the need for custom code subclasses should be able to be constructed \
    by providing just the basename parameter and optional requires_slicing parameter. \
    If additional inputs are needed, then file conversion and management scripts \
    may need to be modified to account for the custom requirements. \

    **Required Attributes:**

    :ivar data_type: String indicating the data type to be used (e.g., uint16)
    :ivar shape: Tuple indicating the shape of the data
    :ivar mz: Numpy array with the m/z axis data. In the case of multi-data this is a
        list of numpy arrays, one per dataset.
    :ivar basename: The basename provided for opening the file.

    **Required Interface Functions:**

    * ``close_file`` : Close any opened files
    * ``is_valid_dataset`` : Check whether a given dir/file is valid under the current format
    * ``spectrum_iter`` : Generator function that iterates over all the spectra in the current data cube and yield the \
        numpy array with the intensity and integer x, y position of the spectrum

    **Optional Interface Functions:**

    * ``__getitem__`` : Implement array slicing for files. Reuqired if the requires_slicing parameter
            should be supported.
    * ``supports_regions`` : Specify whether the format supports multiple regions (default=False)
    * ``supports_multidata``: Specify whether the format supports multiple datasets (default=False)
    * ``supports_multiexperiment``: Specify whether the format supports multiple experiments (default=False)

    """

    def __init__(self, basename, requires_slicing=True):
        """
        Construct the base class and define required attributes.

        :param basename: The name of the file or directory with the file to be opened by the reader
        :param requires_slicing: Boolean indicating whether the user requires array slicing via
            the __getitem__ function to work or not. This is an optimization, because many MSI
            data formats do not easily support arbitrary slicing of data but rather only
            iteration over spectra.
        """
        self.data_type = ''
        self.shape = ()
        self.mz = None
        self.basename = basename
        self.requires_slicing = requires_slicing

    def __getitem__(self, key):
        """Enable slicing of files"""
        raise NotImplementedError('Slicing not supported by the file reader.')

    def spectrum_iter(self):
        """
        Enable iteration of the spectra of the current data cube.

        Iterate over all the spectra in the current data cube and yield the
        numpy array with the intensity and integer x, y position of the spectrum.

        NOTE: As this is a generator one needs to use yield.

        :returns: The function yields for each spectrum the following information:

            * tuple of (x,y) or (x,y,z) position of the spectrum
            * Numpy array with the spectrum
        """
        raise NotImplementedError('Iteration over spectra not supported by the file reader.')

    def close_file(self):
        """
        Close the file.
        """
        raise NotImplementedError('close_file function not implemented.')

    def get_dataset_metadata(self):
        """
        Get dict of additional metadata associated with the current dataset

        NOTE: In the case that multiple regions and/or datasets are supported,
        this function should return the metadata of the currently selected
        dataset only. If no particular dataset is selected, then all should
        be returned.

        :return: Instance of omsi.shared.metadata_data.metadata_dict
        """
        return metadata_dict()

    def get_instrument_metadata(self):
        """
        Get dict of additional metadata associated with the current instrument

        NOTE: In the case that multiple regions and/or datasets are supported,
        this function should return the metadata of the currently selected
        dataset only. If no particular dataset is selected, then all should
        be returned.

        :return: Instance of omsi.shared.metadata_data.metadata_dict
        """
        return metadata_dict()

    def get_method_metadata(self):
        """
        Get dict of additional metadata associated with the current instrument

        NOTE: In the case that multiple regions and/or datasets are supported,
        this function should return the metadata of the currently selected
        dataset only. If no particular dataset is selected, then all should
        be returned.

        :return: Instance of omsi.shared.metadata_data.metadata_dict
        """
        return metadata_dict()

    @classmethod
    def is_valid_dataset(cls, name):
        """
        Classmethod used to check whether a given directory (or file)
        defines as valid data file of the file format specified by
        the current child class

        :param name: Name of the dir or file.
        :type name: String

        :returns: Boolean indicating whether the given file or folder is a valid file.

        """
        raise NotImplementedError('is_valid_dataset function not implemented')

    @classmethod
    def size(cls, name):
        """
        Classmethod used to check the estimated size for the given file/folder.

        :param name: Name of the dir or file.
        :type name: String

        :returns: Integer indicating the size in byte or None if unknown.
        """
        return None

    def get_number_of_regions(self):
        """
        File readers with multi region support must overwrite this function
        to retrieve the true number of regions in the file.
        Default implementation returns 1.
        """
        return 1

    def get_number_of_datasets(self):
        """
        File readers with multi dataset support must overwrite this function
        to retrieve the true number of raw datasets in the file.
        Default implementation returns 1.
        """
        return 1

    @classmethod
    def supports_regions(cls):
        """
        Define whether the file format support multiple regions.
        """
        return False

    @classmethod
    def supports_multidata(cls):
        """
        Define whether the file format support multiple independent datasets.
        """
        return False

    @classmethod
    def supports_multiexperiment(cls):
        """
        Define whether the file format supports multiple independent experiments,
        each of which may contain multiple datasets.
        """
        return False

    def __del__(self):
        """Close the file before garbage collection"""
        self.close_file()

    @classmethod
    def format_name(cls):
        """
        Define the name of the format.

        :returns: String indicating the name of the format.
        """
        return cls.__name__

    @staticmethod
    def available_formats():
        """
        Get dictionary of all available file formats that implement the file_format_base API.

        :returns: Dictionary where the keys are the names of the formats and the values are
                  the corresponding file reader classes.
        """
        # Create list of all derived formats
        formatreaders = {sub.format_name(): sub for sub in file_reader_base.__subclasses__()}
        formatreaders.update({sub.format_name(): sub for sub in file_reader_base_with_regions.__subclasses__()})
        formatreaders.update({sub.format_name(): sub for sub in file_reader_base_multidata.__subclasses__()})
        # Remove extended interface classes which do not implement an actual format
        formatreaders.pop('file_reader_base_with_regions')
        formatreaders.pop('file_reader_base_multidata')
        # Retrun the list of formats
        return formatreaders


class file_reader_base_with_regions(file_reader_base):
    """
    Base-class used to define the basic interface for file-readers
    used to implement new file formats with support for multiple
    imaging regions per file. This class extends file_reader_base,
    and accordingly all required attributes and functions of
    file_reader_base must be implemented by subclasses.

    **Additional required attributes:**

    * ``select_region`` : Integer indicating which region should be selected. \
                    If set to None, indicates that the data should be treated \
                    as a whole. If set to a region index, then the data should \
                    be treated by the reader as if it only pertains to that \
                    region, ie., the shape of the data should be set accordingly \
                    and __getitem__ should behave as such as well.
    * ``region_dicts`` : List of dictionaries, where each dictionary describes a \
                    given region (e.g,. the origin and extend for rectangular \
                    regions.

    """
    def __init__(self, basename, requires_slicing):
        """
        Construct the base class and define required attributes.
        """
        super(file_reader_base_with_regions, self).__init__(basename, requires_slicing)
        self.select_region = None
        self.region_dicts = []

    def set_region_selection(self, region_index=None):
        """
        Define which region should be selected for local data reads.

        :param region_index: The index of the region that should be read. The shape of the
            data will be adjusted accordingly. Set to None to select all regions and treat
            the data as a single full 3D image.
        """
        raise NotImplementedError('set_region_selection function not implemented')

    def get_region_selection(self):
        """Get the index of the selected region"""
        return self.select_region

    def get_number_of_regions(self):
        """Get the number of available regions"""
        return len(self.region_dicts)

    def get_regions(self):
        """
        Get list of all region dictionaries defining for each region the origin
        and extend of the region. See also self.region_dicts.
        """
        return self.region_dicts

    def get_dataset_dependencies(self):
        """
        Get the dependencies between the current region and any of the
        other region datasets stored in the current file. If self.select_region
        is not set, then the function is expected to return a list of lists
        with all dependencies for all datasets.

        :return: List of dependencies (or list of lists of dependencies if self.select_dataset is None)
            where each dependency is a dict of the following form:

            * 'omsi_object': None,         # The omsi file API object where the data is stored. Often None.
            * 'link_name': ms2_link_name,  # Name for the dependency link to be used
            * 'basename': basename,        # Basename of the file
            * 'region': None,              # Index of the region in the dataset or None
            * 'dataset': ind2,             # Index of the dataset withing the file or None
            * 'help':scan_types[ms1scan],  # Help describing the depdency
            *  dependency_type': ... }     # Type of dependency see dependency_dict.dependency_type for available types
        """
        raise NotImplementedError('Determine the dependencies to other data blocks for the the current block')

    @classmethod
    def supports_regions(cls):
        """
        Define whether the file format support multiple regions.
        """
        return True


class file_reader_base_multidata(file_reader_base):
    """
    Base-class used to define the basic interface for file-readers
    used to implement new file formats with support for multiple
    dataset (e.g, MSI dataset with multiple spectrum types).
    This class extends file_reader_base, and accordingly all
    required attributes and functions of
    file_reader_base must be implemented by subclasses.

    In addition to the file_reader_base functions we need to implement
    the get_number_of_datasets(...) and get_dataset_dependencies(...)
    functions.

    :ivar select_dataset: Unsigned integer indicating the currently selected dataset
    """

    def __init__(self, basename, requires_slicing):
        """
        Construct the base class and define required attributes.
        """
        super(file_reader_base_multidata, self).__init__(basename, requires_slicing)
        self.select_dataset = None

    def set_dataset_selection(self, dataset_index):
        """
        Define the current dataset to be read.
        """
        try:
            num_datasets = self.get_number_of_datasets()
        except NotImplementedError:
            num_datasets = -1
        if num_datasets >= 0 and dataset_index >= self.get_number_of_datasets():
            raise ValueError("Given dataset_index is out of range.")
        self.select_dataset = dataset_index

    def get_number_of_datasets(self):
        """
        Get the number of available datasets.
        """
        raise NotImplementedError('Determine the number of different data blocks')

    def get_dataset_dependencies(self):
        """
        Get the dependencies between the current dataset and any of the
        other datasets stored in the current file. If self.select_dataset
        is not set, then the function is expected to return a list of lists
        with all dependencies for all datasets.

        :return: List of dependencies (or list of lists of dependencies if self.select_dataset is None)
            where each dependency is a dict of the following form:

            * 'omsi_object': None,         # The omsi file API object where the data is stored. Often None.
            * 'link_name': ms2_link_name,  # Name for the dependency link to be used
            * 'basename': basename,        # Basename of the file
            * 'region': None,              # Index of the region in the dataset or None
            * 'dataset': ind2,             # Index of the dataset withing the file or None
            * 'help':scan_types[ms1scan],  # Help describing the depdency
            *  dependency_type': ... }     # Type of dependency see dependency_dict.dependency_type for available types

        """
        raise NotImplementedError('Determine the dependencies to other data blocks for the the current block')

    @classmethod
    def supports_multidata(cls):
        """
        Define whether the file format supports multiple data blocks.
        """
        return True

