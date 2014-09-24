"""
This module provides functionality for reading mzml mass spectrometry image files.


filename = '/Users/oruebel/Devel/openmsi-data/mzML_Data/N2A2_Serratia_spots_extract_TI.mzML'
"""
import numpy as np
from pyteomics import mzml
import re
import os
from omsi.dataformat.file_reader_base import file_reader_base


class mzml_file(file_reader_base):
    """
    Interface for reading a single 2D mzml file.

    :ivar available_mzml_types: Dict of available mzml flavors.
    """

    available_mzml_types = {'unknown': 'unknown',
                            'bruker': 'bruker',
                            'thermo': 'thermo'}

    def __init__(self, basename, readdata=True):
        """
        Open an img file for data reading.

        :param basename: The name of the mzml file. If basename is a directory, then the first mzML file found
                             in the directory will be used instead.
        :type basename: string

        :param readdata: Should the complete data be read into memory
                             (this makes slicing easier). (default is True)
        :type readdata: bool

        """
        # Determine the correct base
        if os.path.isdir(basename):
            filelist = self.get_files_from_dir(basename)
            if len(filelist) > 0:
                basename = filelist[0]
            else:
                raise ValueError("No valid mzML file found in the given directory.")

        # Call super constructor. This sets self.basename and self.readall
        super(mzml_file, self).__init__(basename=basename, readdata=readdata)
        self.mzml_type = self.__compute_filetype(filename=self.basename)
        self.data_type = 'uint32'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans(filename=self.basename)
        self.coordinates = self.__compute_coordinates()

        #Compute the spatial configuration of the matrix
        self.x_pos = np.unique(self.coordinates[:, 0])
        self.y_pos = np.unique(self.coordinates[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        #Compute the mz axis
        self.mz = self.__compute_mz_axis(filename=self.basename, mzml_filetype=self.mzml_type)

        #Determine the shape of the dataset
        self.shape = (self.x_pos.shape[0], self.y_pos.shape[0], self.mz.shape[0])

        #Read the data into memory
        self.data = None
        if readdata:
            self.__read_all()

    def __read_all(self):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.
        """
        reader = mzml.read(self.basename)
        self.data = np.zeros(shape=self.shape, dtype=self.data_type)
        spectrumid = 0
        for spectrum in reader:
            x = spectrum['m/z array']
            try:
                y = spectrum['intensity array']
            except KeyError:
                # TODO Do we need user switched for reading the different spectra MS1 vs. profile spectrum?
                #The ['MS1 spectrum'] stores the calibrated data
                y = spectrum['MS1 spectrum']
                #['profile spectrum'] stores raw data
            yi = np.interp(self.mz, x, y, 0, 0)
            xidx = np.nonzero(self.x_pos == self.coordinates[spectrumid, 0])[0]
            yidx = np.nonzero(self.y_pos == self.coordinates[spectrumid, 1])[0]
            self.data[xidx, yidx, :] = yi  # TODO Note if the data is expected to be of float precision then self.data_type needs to be set accordingly.
            spectrumid += 1

    @classmethod
    def __compute_mz_axis(cls, filename, mzml_filetype):
        """
        Internal helper function used to compute the mz axis
        """
        reader = mzml.read(filename)
        spectrum = reader.next()
        if mzml_filetype == cls.available_mzml_types['thermo']:
            mzmin = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
            mzmax = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']
            # TODO Should this be a user-definable parameter?
            ppm = 5
            f = np.ceil(1e6*np.log(mzmax/mzmin)/ppm)
            return np.logspace(np.log10(mzmin), np.log10(mzmax), f)
        elif mzml_filetype == cls.available_mzml_types['bruker']:
            return spectrum['m/z array']
        else:
            raise ValueError('Unknown mzml format')

    @classmethod
    def __compute_filetype(cls, filename):
        """
        Internal helper function used to compute the filetype.
        """
        spectrum = next(mzml.read(filename))
        if 'spotID' in spectrum:
            return cls.available_mzml_types['thermo']
        elif 'id' in spectrum:
            return cls.available_mzml_types['bruker']
        else:
            return cls.available_mzml_types['unknown']

    def __compute_coordinates(self):
        """
        Internal helper function used to compute the coordinates for each scan.

        :returns: 2D numpy integer array of shape (numScans,2) indicating for each scan its x and y coordinate
        """
        reader = mzml.read(self.basename)
        coords = np.zeros(shape=(self.num_scans, 2), dtype='uint32')
        if self.mzml_type == self.available_mzml_types['thermo']:
            spectrumid = 0
            for spectrum in reader:
                spotid = spectrum['spotID']
                coords[spectrumid, :] = map(int, spotid.split(',')[-1].split('x'))
                spectrumid += 1
        elif self.mzml_type == self.available_mzml_types['bruker']:
            spectrumid = 0
            for spectrum in reader:
                spotdesc = spectrum['id'].split('_x002f_')[1]
                matchobj = re.findall('\d+', spotdesc)
                coords[spectrumid, 0] = int(matchobj[2])
                coords[spectrumid, 1] = int(matchobj[3])
                spectrumid += 1
        return coords

    @staticmethod
    def __compute_num_scans(filename=None):
        """
        Internal helper function used to compute the number of scans in the mzml file.
        """
        reader = mzml.read(filename)
        return sum(1 for _ in reader)

    def __getitem__(self, key):
        """Enable slicing of img files"""
        if self.data is not None:
            return self.data[key]
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
        if os.path.isdir(name):  # If we point to a director, check if the dir contains an mzML file
            filelist = cls.get_files_from_dir(name)
            return len(filelist) > 0
        else:
            try:
                #Try to open the file and iterate over it
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
        if os.path.isdir(name):  # If we point to a director, check if the dir contains an mzML file
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

            #temp_mzml_file = cls(basename=basename, readdata=False)
            #itemsize = np.dtype(temp_mzml_file.data_type).itemsize
            #size = np.asarray(temp_mzml_file.shape).prod() * itemsize
            #print ('MZML size', size)
            #return size
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
