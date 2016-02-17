"""
This module provides functionality for reading imzML mass spectrometry image files.
"""

# import basic packages
import numpy as np
import os

# import xml parser
from pyimzml.ImzMLParser import ImzMLParser
from omsi.dataformat.file_reader_base import *
from omsi.shared.log import log_helper

class imzml_file(file_reader_base):
    """
    Interface for reading a single 2D imzml file with a single distinct scan types.
    :ivar available_imzml_types: Dict of available mzml flavors.
    """
    # TODO: Rewrite the entire class to call the pyimzml.ImzMLParser only once for efficiency.
    # TODO:   this code is derived from a pymzml-based parser for .mzML files which didn't read the whole file on import

    available_imzml_types = {'unknown': 'unknown', 'thermo': 'thermo', 'bruker': 'bruker'}

    def __init__(self, basename, requires_slicing=True):
        """
        Open an imzml file for data reading.

        :param  basename:   The name of the mzml file. If basename is a directory, then the first mzML file found
                             in the directory will be used instead.
        :type   basename:   string

        :param  requires_slicing:   Should the complete data be read into memory
                             (this makes slicing easier). (default is True)
        :type   requires_slicing:   bool
        """
        # Determine the correct base
        if os.path.isdir(basename):
            filelist = self.get_files_from_dir(basename)
            if len(filelist) > 0:
                basename = filelist[0]
            else:
                raise ValueError("No valid imzML file found in the given directory.")

        # Call super constructor. This sets self.basename and self.readall
        super(imzml_file, self).__init__(basename=basename, requires_slicing=requires_slicing)
        self.mzml_type = self.__compute_filetype(filename=self.basename)
        self.data_type = 'uint32'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans(filename=self.basename)
        log_helper.info(__name__, 'Read %s scans from imzML file.' % self.num_scans)
        self.scan_types = self.__compute_scan_types(filename=self.basename)
        self.scan_dependencies = self.__compute_scan_dependencies(scan_types=self.scan_types)

        # Read x, y coordinates from file
        self.coordinates = self.__compute_coordinates(filename=basename)

        # Determine scan parameters ## TODO: after solving imzML generation prob, fix this for multicube data
        self.scan_params = self.__parse_scan_parameters(self)

        # Compute step size
        coordarray = np.asarray(self.coordinates)
        self.x_pos = np.unique(coordarray[:, 0])
        self.y_pos = np.unique(coordarray[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        # Compute the mz axis
        self.mz_all = self.__compute_mz_axis(filename=self.basename,
                                             mzml_filetype=self.mzml_type,
                                             scan_types=self.scan_types)

        # Determine the shape of the dataset ## TODO: after solving imzML generation prob, fix this for multicube data

        self.shape_all_data = [len(np.unique(self.x_pos)), len(np.unique(self.y_pos)), len(self.mz_all)]
        self.shape = None
        self.mz = None
        self.select_dataset = None  # Used only for *.img files?

        # Read the data into memory
        self.data = None
        if requires_slicing:
            self.__read_all(filename=basename)

        self.mz = self.mz

    def __read_all(self, filename):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.  Data is now a list of datacubes.
        """

        self.data = np.zeros(shape=self.shape_all_data, dtype=self.data_type)
        log_helper.info(__name__, 'Datacube shape is %s' % [self.data.shape])

        reader = ImzMLParser(filename)
        log_helper.debug(__name__,'READING ALL DATA!! GIVE ME RAM (please)!')
        for ind in xrange(0, len(reader.coordinates)):
            xidx, yidx = reader.coordinates[ind]
            mz, intens = reader.getspectrum(ind)
            try:
                self.data[xidx-1, yidx-1, :] = intens
            except:
                print xidx, yidx, len(mz), len(intens)

    def spectrum_iter(self):
        """
        Generator function that yields a position and associated spectrum for a selected datacube type.
        :yield: (xidx, yidx) a tuple of ints representing x and y position in the image
        :yield: yi,          a numpy 1D-array of floats containing spectral intensities at the given position
                                and for the selected datacube type
        """
        reader = ImzMLParser(self.basename)
        if self.select_dataset is None:
            # multiple datatypes are not supported in mzML files
            self.select_dataset = 0

        for idx in xrange(0, len(reader.coordinates)):
            xidx, yidx = reader.coordinates[idx]
            mz, intens = reader.getspectrum(idx)

            yield (xidx-1, yidx-1), intens  # -1 because pyimzML coordinates are 1-based

    @classmethod
    def __compute_mz_axis(cls, filename, mzml_filetype, scan_types):
        ## TODO completely refactor this to make it smartly handle profile or centroid datasets
        ## TODO: centroid datasets should take in a user parameter "Resolution" and resample data at that resolution
        ## TODO: profile datasets should work as is
        ## TODO: checks for profile data vs. centroid data on the variation in length of ['m/z array']
        """
        Internal helper function used to compute the mz axis of each scantype
        Returns a list of numpy arrays
        """
        reader = ImzMLParser(filename)
        mz_axes, intens = reader.getspectrum(0)

        for ind, loc in enumerate(reader.coordinates):
            mz, intens = reader.getspectrum(ind)
            if mz == mz_axes:
                pass
            else:
                log_helper.error(__name__, 'Inconsistent m/z axis from scan to scan. ' +
                                           'Currently OpenMSI supports only continuous-mode imzML.')
                raise AttributeError

        return mz_axes

    @classmethod
    def __compute_filetype(cls, filename):
        """
        Internal helper function used to compute the filetype.
        TODO: figure out if this will ever be useful for imzML
        """
        return cls.available_imzml_types['unknown']

    @staticmethod
    def __compute_coordinates(filename):
        """
        Internal helper function used to compute the coordinates for each scan.

        :returns: 2D numpy integer array of shape (numScans,2) indicating for each scan its x and y coordinate
        """
        reader = ImzMLParser(filename)
        return reader.coordinates

    @classmethod
    def test(cls):
        """
        Test method
        """
        pass

    @staticmethod
    def __compute_num_scans(filename=None):
        """
        Internal helper function used to compute the number of scans in the imzml file.
        """
        reader = ImzMLParser(filename)
        return len(reader.coordinates)

    @staticmethod
    def __compute_scan_types(filename=None):
        """
        Internal helper function used to compute a list of unique scan types in the imzml file.
        """
        # TODO: FIGURE OUT HOW TO WRITE THIS FUNCTION.  I cannot figure out how to use any publicly available
        # TODO:               imzML parser to convert thermo .raw files to .imzml or to convert .mzML files to .imzml
        # TODO:               http://www.cs.bham.ac.uk/~ibs/imzMLConverter/ is based on Java and does not run on my
        # TODO:               macbook pro / Yosemite OSX  with java successfully installed.
        # TODO:               'The Java JAR file "imzMLConverter.jar" could not be launched.  Check the Console for
        # TODO:                possible error messages.'
        # TODO:               http://www.imzml.org/index.php?option=com_content&view=article&id=312&Itemid=99
        # TODO:               is a Windows .exe and did not run on the Windows machine I tried it on:
        # TODO:               "Exception EOIeSysError in Modul imzMLConvert.exe bei 00078531.  Class not registered."

        reader = ImzMLParser(filename)
        scantypes = []
        #for _ in reader:
        #    try:
        #        scanfilter = reader.next()['scanList']['scan'][0]['filter string']
        #        if scanfilter not in scantypes:
        #            scantypes.append(scanfilter)
        #    except:
        #        pass
        #return scantypes
        scantypes = ['1']
        return scantypes

    @staticmethod
    def __parse_scan_parameters(self):
        """
        Internal helper function used to parse out scan parameters from the scan filter string
        """
        # TODO: implement this function for imzML files.  Right now inability to generate imzML from
        # TODO:    arbitrary files is blocking

        scan_params = []
        return scan_params

    @staticmethod
    def __compute_scan_dependencies(scan_types=None):
        """
        Takes a scan_types list and returns a list of tuples (x, y) indicating that scan_type[y] depends on scan_type[x]
        Internal helper function used to parse out scan parameters from the scan filter string
        """
        # TODO: implement this function for imzML files.  Right now inability to generate imzML from
        # TODO:   arbitrary files is blocking
        dependencies = []
        return dependencies

    def __getitem__(self, key):
        """Enable slicing of img files"""
        if self.data is not None:
            if self.select_dataset is None:
                return self.data[key]
            else:
                return self.data[self.select_dataset][key]
        else:
            raise ValueError("Slicing is currently only supported when the object has been initialized with readall")

    def close_file(self):
        """Close the mzml file"""
        pass

    @classmethod
    def is_valid_dataset(cls, name):
        """Check whether the given file or directory points to a img file.

           :param name: Name of the dir or file.
           :type name: String

           :returns: Boolean indicating whether the given file or folder is a valid img file.
        """
        if os.path.isdir(name):  # If we point to a directory, check if the dir contains an imzML file
            filelist = cls.get_files_from_dir(name)
            return len(filelist) > 0
        else:
            try:
                # Try to open the file
                temp = ImzMLParser(name)
                del temp
                return True
            except (RuntimeError, IOError) as e:
                # print ('pyimzml could not parse your file.')
                return False


    @classmethod
    def size(cls, name, max_num_reads=1000):
        """
        Classmethod used to check the estimated size for the given file/folder.
        For mzml this is an estimate of the final size of the full 3D datacube.
        For efficiency the number of scans is estimated based on the size of
        the first 1000 scans.

        :param name: Name of the dir or file.
        :type name: unicode
        :param max_num_reads: The maximum number of spectrum reads to be performed to estimate the file size
        :type max_num_reads: int

        :returns: None
        """
        # TODO: implement sensible size estimation scheme for imzML, perhaps based on associated *.ibd file size

        basename = None
        if os.path.isdir(name):  # If we point to a directory, check if the dir contains an mzML file
            filelist = cls.get_files_from_dir(name)
            if len(filelist) > 0:
                basename = filelist[0]
        else:
            basename = name
        if basename is not None:
           return None

    @classmethod
    def get_files_from_dir(cls, dirname):
        """
        Get a list of all basenames of all img files in a given directory.
        Note: The basenames include the dirname.
        """
        filelist = []
        for l in os.listdir(dirname):
            currname = os.path.join(dirname, l)
            filename_only, extension = os.path.splitext(currname)
            if os.path.isfile(currname) and currname.endswith(".imzML"):
                if os.path.isfile(filename_only + '.ibd'):
                    filelist.append(currname)
                else:
                    log_helper.warning(__name__, 'Could not find binary .ibd file for file %s . Skipping conversion of this file.' % currname)

        return filelist

    def get_number_of_datasets(self):
        """
        Get the number of available datasets.
        """
        return len(self.mz_all)

    def set_dataset_selection(self, dataset_index):
        """
        Define the current dataset to be read.
        """
        super(imzml_file, self).set_dataset_selection(dataset_index)
        self.shape = self.shape_all_data[self.select_dataset]
        self.mz = self.mz_all[self.select_dataset]

    def get_dataset_dependencies(self):
        """
        Get the dependencies between the current dataset and any of the
        other datasets stored in the current file.
        """
        # TODO Implement dependencies between current dataset given by self.select_dataset and all other datasets
        return []