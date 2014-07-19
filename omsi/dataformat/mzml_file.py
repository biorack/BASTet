"""
This module provides functionality for reading mzml mass spectrometry image files.


filename = '/Users/oruebel/Devel/openmsi-data/mzML_Data/N2A2_Serratia_spots_extract_TI.mzML'
"""
import numpy as np
from pyteomics import mzml
import re
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
        """Open an img file for data reading.

            :param basename: The name of the mzml file
            :type basename: string

             :param readdata: Should the complete data be read into memory
                              (this makes slicing easier). (default is True)
            :type readdata: bool

        """
        # Call super constructor. This sets self.basename and self.readall
        super(mzml_file, self).__init__(basename=basename, readdata=readdata)
        self.mzml_type = self.__compute_filetype()
        self.mz = self.__compute_mz_axis()
        self.data_type = 'uint16'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans()
        self.coordinates = self.__compute_coordinates()

        #Compute the spatial configuration of the matrix
        self.x_pos = np.unique(self.coordinates[:, 0])
        self.y_pos = np.unique(self.coordinates[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        #Compute the mz axis
        self.mz = self.__compute_mz_axis()

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

    def __compute_mz_axis(self):
        """
        Internal helper function used to compute the mz axis
        """
        reader = mzml.read(self.basename)
        spectrum = reader.next()
        if self.mzml_type == self.available_mzml_types['thermo']:
            mzmin = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
            mzmax = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']
            # TODO Should this be a user-definable parameter?
            ppm = 5
            f = np.ceil(1e6*np.log(mzmax/mzmin)/ppm)
            return np.logspace(np.log10(mzmin), np.log10(mzmax), f)
        elif self.mzml_type == self.available_mzml_types['bruker']:
            return spectrum['m/z array']
        else:
            raise ValueError('Unknown mzml format')

    def __compute_filetype(self):
        """
        Internal helper function used to compute the filetype.
        """
        spectrum = next(mzml.read(self.basename))
        if 'spotID' in spectrum:
            return self.available_mzml_types['thermo']
        elif 'id' in spectrum:
            return self.available_mzml_types['bruker']
        else:
            return self.available_mzml_types['unknown']

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

    def __compute_num_scans(self):
        """
        Internal helper function used to compute the number of scans in the mzml file.
        """
        reader = mzml.read(self.basename)
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
        try:
            #Try to open the file and iterate over it
            reader = mzml.read(name)
            for _ in reader:
                pass
            return True
        except:
            return False

