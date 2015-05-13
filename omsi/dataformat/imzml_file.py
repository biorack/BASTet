"""
This module provides functionality for reading imzML mass spectrometry image files.
"""

import numpy as np
import re
import os

from pyteomics import mzml
from pyimzml.ImzMLParser import ImzMLParser  ## TODO: install pyimzml at NERSC

from omsi.dataformat.file_reader_base import file_reader_base_multidata


class imzml_file(file_reader_base_multidata):
    """
    Interface for reading a single 2D imzml file with a single distinct scan types.
    :ivar available_imzml_types: Dict of available mzml flavors.
    """

    available_imzml_types = {'unknown': 'unknown', 'thermo': 'thermo', 'bruker': 'bruker'}

    def __init__(self, basename, readdata=True):
        """
        Open an imzml file for data reading.

        :param basename: The name of the mzml file. If basename is a directory, then the first mzML file found
                             in the directory will be used instead.
        :type basename: string

        :param readdata: Should the complete data be read into memory
                             (this makes slicing easier). (default is True)
        :type readdata: bool

        :param resolution: For profile data only, the minimum m/z spacing to use for creating the "full" reprofiled
                            data cube
        :type resolution: float
        """
        # Determine the correct base
        if os.path.isdir(basename):
            filelist = self.get_files_from_dir(basename)
            if len(filelist) > 0:
                basename = filelist[0]
            else:
                raise ValueError("No valid mzML file found in the given directory.")

        # Call super constructor. This sets self.basename and self.readall
        super(imzml_file, self).__init__(basename=basename, readdata=readdata)
        self.mzml_type = self.__compute_filetype(filename=self.basename)
        self.data_type = 'uint32'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans(filename=self.basename)
        print 'Read %s scans from mzML file.' % self.num_scans
        self.scan_types = self.__compute_scan_types(filename=self.basename)
        self.scan_dependencies = self.__compute_scan_dependencies(scan_types=self.scan_types)
        print 'Found %s different scan types in mzML file.' % len(self.scan_types)
        self.coordinates = self.__compute_coordinates(filename=basename)
        self.scan_params = self.__parse_scan_parameters(self)

        # TODO: redo this block for imzML / Compute the spatial configuration of the matrix

        coordarray = np.asarray(self.coordinates)

        self.x_pos = np.unique(coordarray[:, 0])
        self.y_pos = np.unique(coordarray[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        # Compute the mz axis
        self.mz_all = self.__compute_mz_axis(filename=self.basename, mzml_filetype=self.mzml_type,
                                         scan_types=self.scan_types)

        # Determine the shape of the dataset, result is a list of shapes for each datacube

        self.shape_all_data = [len(np.unique(self.x_pos)), len(np.unique(self.y_pos)), len(self.mz_all)]
        self.shape = None
        self.mz = None

        # Read the data into memory
        self.data = None
        if readdata:
            self.__read_all(filename=basename)

    def __read_all(self, filename):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.  Data is now a list of datacubes.
        """

        self.data = np.zeros(shape=self.shape_all_data, dtype=self.data_type)
        print self.data.shape

        reader = ImzMLParser(filename)

        for ind in xrange(0, len(reader.coordinates)):
            xidx, yidx = reader.coordinates[ind]
            mz, intens = reader.getspectrum(ind)
            try:
                print len(intens)
                print type(intens)
                self.data[xidx-1, yidx-1, :] = intens
            except:
                print xidx, yidx, len(mz), len(intens)


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
                print 'Inconsistent m/z axis from scan to scan'
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

    # TODO: add methods for classifying scans by type: MS1 vs. MS2_precursorA vs. MS2_precursorB ?

    @classmethod
    def test(cls):
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
        TODO: FIGURE OUT HOW TO WRITE THIS FUNCTION.  I cannot figure out how to use any publicly available
                        imzML parser to convert thermo .raw files to .imzml or to convert .mzML files to .imzml

                        http://www.cs.bham.ac.uk/~ibs/imzMLConverter/ is based on Java and does not run on my
                        macbook pro / Yosemite OSX  with java successfully installed.
                        'The Java JAR file "imzMLConverter.jar" could not be launched.  Check the Console for possible error messages.'

                        http://www.imzml.org/index.php?option=com_content&view=article&id=312&Itemid=99 is a Windows .exe
                        and did not run on the Windows machine I tried it on.  "Exception EOIeSysError in Modul imzMLConvert.exe bei 00078531.  Class not registered."

        """
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
        TODO: implement this function for imzML files.  Right now inability to generate imzML from
        arbitrary files is blocking
        """
        # precursor m/z
        # example scan filter: MS2: ITMS + p MALDI Z ms2 907.73@cid60.00 [500.00-700.00]
        ## example scan filter: MS1: FTMS + p MALDI Full ms [850.00-1000.00]
        ## imzML files that I can access do not contain filter strings

        scan_params = []
        return scan_params

    @staticmethod
    def __compute_scan_dependencies(scan_types=None):  #Why is the "self" arg needed here & below?
        """
        Takes a scan_types list and returns a list of tuples (x, y) indicating that scan_type[y] depends on scan_type[x]
        Internal helper function used to parse out scan parameters from the scan filter string
        TODO: implement this function for imzML files.  Right now inability to generate imzML from
        arbitrary files is blocking
        """
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
        if os.path.isdir(name):  # If we point to a directory, check if the dir contains an mzML file
            filelist = cls.get_files_from_dir(name)
            return len(filelist) > 0
        else:
            try:
                # Try to open the file and iterate over it
                reader = mzml.read(name)
                for _ in reader:
                    pass
                return True
            except:
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

        :returns: Integer indicating the size in byte or None if unknown.
        """
        basename = None
        if os.path.isdir(name):  # If we point to a directory, check if the dir contains an mzML file
            filelist = cls.get_files_from_dir(name)
            if len(filelist) > 0:
                basename = filelist[0]
        else:
            basename = name
        if basename is not None:
            num_scans = -1
            # Try to compute the number of scans by looking at the spectrumList count entry in the file
            try:
                size_line = os.popen('head -n 120 "' + basename + '" | grep "spectrumList count="').read()
                if len(size_line) > 0:
                    size_text = size_line.split('spectrumList count=')[1].split('"')[1]
                    if size_text.isdigit():
                        num_scans = int(size_text)
            except:
                pass
            if num_scans < 0:
                # Estimate the number of scans by reading the first 1000 spectra
                index = 0
                prev_tell = 0
                sizes = []
                reader = mzml.read(basename)
                for _ in reader:
                    if index >= max_num_reads:
                        break
                    current_tell = reader.file.file.tell()
                    sizes.append(current_tell - prev_tell)
                    prev_tell = current_tell
                    index += 1
                npsizes = np.asarray(sizes)
                filesize = os.stat(basename).st_size
                scansize = (npsizes.max() - npsizes.min()) / 2.
                num_scans = int(filesize/scansize)
            mz_axis_len = cls.__compute_mz_axis(filename=basename,
                                                mzml_filetype=cls.__compute_filetype(filename=basename)).shape[0]
            return num_scans*mz_axis_len

            # temp_mzml_file = cls(basename=basename, readdata=False)
            # itemsize = np.dtype(temp_mzml_file.data_type).itemsize
            # size = np.asarray(temp_mzml_file.shape).prod() * itemsize
            # print ('MZML size', size)
            # return size
        else:
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
            if os.path.isfile(currname) and currname.endswith(".mzML"):
                filelist.append(currname)
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
        super(mzml_file, self).set_dataset_selection(dataset_index)
        self.shape = self.shape_all_data[self.select_dataset]
        self.mz = self.mz_all[self.select_dataset]

    def get_dataset_dependencies(self):
        """
        Get the dependencies between the current dataset and any of the
        other datasets stored in the current file.
        """
        # TODO We need to implement the creation of dependencies between the current dataset given by self.select_dataset and all other datasets
        return []