class file_reader_base(object):
    """Base-class used to define the basic interface that file-readers
       for a new format need to implement.

       __init__ interface:
       ^^^^^^^^^^^^^^^^^^^
       To avoid the need for custom code subclasses should be able to be constructed \
       by providing just the basename parameter and optional readall parameter. \
       If additional inputs are needed, then file conversion and management scripts \
       may need to be modified to account for the custom requirements. \

       * ``basename`` : Name of the directory with the files
       * ``readdata`` : Optional parameter used for optimization. When set to false \
                       indicates that no data will be read using the instance but  \
                       the instance may only be used for basic metadata checking, \
                       e.g., shape. \



       Required Attributes:
       ^^^^^^^^^^^^^^^^^^^^

       :ivar data_type: String indicating the data type to be used (e.g., uint16)
       :ivar shape: Tuple indicating the shape of the data
       :ivar mz: Numpy array with the m/z axis data
       :ivar basename: The basename provided for opening the file.

       Required Interface Functions:
       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       * ``__getitem__`` : Implement array slicing for files
       * ``close_file`` : Close any opened files
       * ``is_valid_dataset`` : Check whether a given dir/file is valid under the current format

       Optional Interface Functions:
       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       * ``supports_regions`` : Specify whether the format supports multiple regions (default=False)

    """

    def __init__(self, basename, readdata=True):
        """
        Construct the base class and define required attributes.
        """
        self.data_type = ''
        self.shape = ()
        self.mz = None
        self.basename = basename
        self.readdata = readdata

    def __getitem__(self, key):
        """Enable slicing of files"""
        raise NotImplementedError('Slicing not supported by the file reader.')

    def close_file(self):
        """
        Close the file.
        """
        raise NotImplementedError('close_file function not implemented.')

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
    def supports_regions(cls):
        """
        Define whether the file format support multiple regions.
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
        formatReaders =  {sub.format_name(): sub for sub in file_reader_base.__subclasses__()}
        formatReaders.update({sub.format_name(): sub for sub in file_reader_base_with_regions.__subclasses__()})
        # Remove extended interface classes which do not implement an actual format
        formatReaders.pop('file_reader_base_with_regions')
        # Retrun the list of formats
        return  formatReaders


class file_reader_base_with_regions(file_reader_base):
    """Base-class used to define the basic interface for file-readers
       used to implement new file formats with support for multiple
       imaging regions per file. This class extends file_reader_base,
       and accordingly all required attributes and functions of
       file_reader_base must be implemented by subclasses.

       Additional required attributes:
       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    def __init__(self, basename, readdata):
        """
        Construct the base class and define required attributes.
        """
        super(file_reader_base_with_regions, self).__init__(basename, readdata)
        self.select_region = None
        self.region_dicts = []

    def set_region_selection(self, region_index=None):
        """Define which region should be selected for local data reads.

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
        """ Get list of all region dictionaries defining for each region the origin
            and extend of the region. See also self.region_dicts.
        """
        return self.region_dicts

    @classmethod
    def supports_regions(cls):
        """
        Define whether the file format support multiple regions.
        """
        return True




