"""
This module provides functionality for reading imzML mass spectrometry image files.
"""

# import basic packages
import numpy as np
import os
from scipy import interpolate
# USE xmltodict to expand support for metadata from the imzml files

# import xml parser
from pyimzml.ImzMLParser import ImzMLParser
from omsi.dataformat.file_reader_base import *
from omsi.datastructures.metadata.metadata_data import metadata_dict, metadata_value
from omsi.shared.log import log_helper
try:
    import xmltodict
except ImportError:
    import omsi.shared.third_party.xmltodict

class imzml_file(file_reader_base):
    """
    Interface for reading a single 2D imzml file with a single distinct scan types.
    """
    # TODO: Rewrite the entire class to call the pyimzml.ImzMLParser only once for efficiency.

    available_imzml_types = {'unknown': 'unknown',
                            'continuous': 'continuous',
                            'processed': 'processed'}



    def __init__(self, basename, requires_slicing=True, resolution=15):
        """
        Open an imzml file for data reading.

        :param  basename:   The name of the mzml file. If basename is a directory, then the first mzML file found
                             in the directory will be used instead.
        :type   basename:   string

        :param  requires_slicing:   Should the complete data be read into memory
                             (this makes slicing easier). (default is True)
        :type   requires_slicing:   bool

        :param resolution: For processed data only, the minimum m/z spacing to use for creating the "full" reprofiled
                            data cube
        :type resolution: float
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
        self.resolution=resolution

        # Compute the mz axis, pixel coordinates data type etc.
        self.coordinates, self.mz, self.data_type, self.imzml_type, self.dataset_metadata, self.instrument_metadata, \
        self.method_metadata = self.__compute_file_info(filename=self.basename, resolution=self.resolution)

        self.num_scans = self.coordinates.size
        log_helper.info(__name__, 'Read %s scans from imzML file.' % self.num_scans)

        # Compute step size
        self.x_pos = np.unique(self.coordinates[:, 0])
        self.y_pos = np.unique(self.coordinates[:, 1])
        self.x_pos_min = self.x_pos.min()
        self.y_pos_min = self.y_pos.min()

        # self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        # Determine the shape of the dataset ## TODO: after solving imzML generation prob, fix this for multicube data
        num_x = self.x_pos.max() - self.x_pos.min() + 1
        num_y = self.y_pos.max() - self.y_pos.min() + 1

        self.shape = (num_x,num_y, self.mz.size)

        # Read the data into memory
        self.data = None
        if requires_slicing:
            self.__read_all(filename=basename)

        log_helper.info(__name__, "IMZML file type: " + str(self.imzml_type))
        log_helper.info(__name__, "IMZML data type: " + str(self.data_type))



    def __read_all(self, filename):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.  Data is now a list of datacubes.
        """

        self.data = np.zeros(shape=self.shape, dtype=self.data_type)
        log_helper.info(__name__, 'Datacube shape is %s' % [self.data.shape])
        reader = ImzMLParser(filename)
        log_helper.debug(__name__,'READING ALL DATA!! GIVE ME RAM (please)!')

        # Compute the bin edges for reinterpolation if needed
        if self.imzml_type == self.available_imzml_types['processed']:
            shift = np.diff(self.mz).mean()
            bin_edges = np.append(self.mz, self.mz[-1]+ shift)
        else:
            bin_edges = None
        for ind in xrange(0, len(reader.coordinates)):
            xidx, yidx = reader.coordinates[ind]
            # Coordinates may start at arbitrary locations, hence, we need to substract the minimum to recenter at (0,0)
            xidx -= self.x_pos_min
            yidx -= self.y_pos_min
            # Read the spectrum
            mz, intens = reader.getspectrum(ind)
            # Reinterpolate intensities if we are in processed mode
            if bin_edges is not None:
                f = interpolate.interp1d(mz,intens,fill_value=0,bounds_error=False)
                intens = f(self.mz)
                #intens, bin_edges_new = np.histogram(mz, bins=bin_edges, weights=intens)
            # Save the intensity values in our data cube
            self.data[xidx, yidx, :] = intens

    def spectrum_iter(self):
        """
        Generator function that yields a position and associated spectrum for a selected datacube type.
        :yield: (xidx, yidx) a tuple of ints representing x and y position in the image
        :yield: yi,          a numpy 1D-array of floats containing spectral intensities at the given position
                                and for the selected datacube type
        """
        reader = ImzMLParser(self.basename)
        for idx in xrange(0, len(reader.coordinates)):
            xidx, yidx, zidx = reader.coordinates[idx]
            # Coordinates may start at arbitrary locations, hence, we need to substract the minimum to recenter at (0,0)
            xidx -= self.x_pos_min
            yidx -= self.y_pos_min
            mz, intens = reader.getspectrum(idx)
            # Rehistogram the data if we are in procesed mode
            if self.imzml_type == self.available_imzml_types['processed']:
                # shift = np.diff(self.mz).mean()
                # bin_edges = np.append(self.mz, self.mz[-1]+ shift)
                f = interpolate.interp1d(mz,intens,fill_value=0,bounds_error=False)
                intens = f(self.mz)
                # intens, bin_edges_new = np.histogram(mz, bins=bin_edges, weights=intens)

            yield (xidx, yidx), np.asarray(intens)

    @classmethod
    def __compute_file_info(cls, filename, resolution):
        ## TODO completely refactor this to make it smartly handle profile or centroid datasets
        ## TODO: centroid datasets should take in a user parameter "Resolution" and resample data at that resolution
        ## TODO: profile datasets should work as is
        ## TODO: checks for profile data vs. centroid data on the variation in length of ['m/z array']
        """
        Internal helper function used to compute the mz axis, data type for the intensities, format type

        :return: Numpy array with mz axis
        :return: string with data type
        :return: imzml file type
        :return:
        """
        reader = ImzMLParser(filename)
        # Read the first spectrum
        mz_axes, intens = reader.getspectrum(0)   # NOTE: mz_axes is a tuple
        # Read the coordinates
        coordinates = np.asarray(reader.coordinates)

        # #Start the data at [0,0,0]
        # coordinates[:,0] = coordinates[:,0] - np.amin(coordinates,axis=0)[0]
        # coordinates[:,1] = coordinates[:,1] - np.amin(coordinates,axis=0)[1]
        # coordinates[:,2] = coordinates[:,2] - np.amin(coordinates,axis=0)[2]

        # Determine the data type for the internsity values
        dtype = np.asarray(intens).dtype.str

        # Compute the mz axis and file type
        file_type = cls.available_imzml_types['continuous']
        min_mz, max_mz = np.amin(mz_axes), np.amax(mz_axes)
        for ind in range(coordinates.shape[0]):      #for ind, loc in enumerate(reader.coordinates):
            mz, intens = reader.getspectrum(ind)
            if mz == mz_axes:
                pass
            else:
                file_type = cls.available_imzml_types['processed']
                if min_mz > np.amin(mz):
                    min_mz = np.amin(mz)
                if max_mz < np.amax(mz):
                    max_mz = np.amax(mz)
        # Reinterpolate the mz-axis if we have a processed mode imzml file
        if file_type == cls.available_imzml_types['processed']:
            f = np.ceil(1e6 * np.log(max_mz/min_mz)/resolution)
            mz_axes = np.logspace(np.log10(min_mz), np.log10(max_mz), f)
            log_helper.info(__name__, "Reinterpolated m/z axis for processed imzML file")

        # Construct the imzml metadata information
        dataset_metadata = metadata_dict()
        instrument_metadata = metadata_dict()
        method_metadata = metadata_dict()
        for k, v in reader.imzmldict.iteritems():
            dataset_metadata[k] = metadata_value(name=k,
                                                 value=v,
                                                 unit=None,
                                                 description=k,
                                                 ontology=None)

        # Delete the parser and read the metadata
        del reader

        # Parse the metadata for the file. We try to parse only the header and ignore the
        # <run > group in the XML file to avoid going throught the whole file again
        # while extracting the majority of the relevant metadata
        try:
            with open(filename, 'r') as ins:
                metdata_header = ''
                for line in ins:
                    if '<run' in line:
                        break
                    else:
                        metdata_header += line
                metdata_header += '</mzML>'
                metdata_header_dict = xmltodict.parse(metdata_header)['mzML']
                for k, v in metdata_header_dict.iteritems():
                    store_value = metadata_value(name=k,
                                                 value=v,
                                                 unit=None,
                                                 description=str(k) + " extracted from imzML XML header.",
                                                 ontology=None)
                    if k == 'instrumentConfigurationList':
                        instrument_metadata[k] = store_value
                    elif k == 'dataProcessingList':
                        method_metadata[k] = store_value
                    elif k == 'scanSettingsList':
                        dataset_metadata[k] = store_value
                    elif k == 'softwareList':
                        method_metadata[k] = store_value
                    elif k =='sampleList':
                        method_metadata[k] = store_value
                    else:
                        dataset_metadata[k] = store_value
                dataset_metadata['imzml_xml_metadata_header'] = metadata_value(name='imzml_xml_metadata_header',
                                                                               value=metdata_header,
                                                                               unit=None,
                                                                               description='XML imzML header',
                                                                               ontology=None)
        except:
            log_helper.warning(__name__, "Extraction of additional imzML metadata failed")

        return coordinates, np.asarray(mz_axes), dtype, file_type, dataset_metadata, instrument_metadata, method_metadata

    @classmethod
    def test(cls):
        """
        Test method
        """
        pass

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
        """
        Check whether the given file or directory points to a img file.

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
            if os.path.isfile(currname) and currname.lower().endswith(".imzml"):
                if os.path.isfile(filename_only + '.ibd'):
                    filelist.append(currname)
                else:
                    log_helper.warning(__name__, 'Could not find binary .ibd file for file %s . Skipping conversion of this file.' % currname)

        return filelist

    def get_dataset_metadata(self):
        """
        Get dict of additional metadata associated with the current dataset.

        Inherited from file_reader_base

        :return: Instance of omsi.shared.metadata_data.metadata_dict

        """
        return self.dataset_metadata

    def get_instrument_metadata(self):
        """
        Get dict of additional metadata associated with the current instrument

        Inherited from file_reader_base

        :return: Instance of omsi.shared.metadata_data.metadata_dict
        """
        return self.instrument_metadata

    def get_method_metadata(self):
        """
        Get dict of additional metadata associated with the current instrument

        Inherited from file_reader_base

        :return: Instance of omsi.shared.metadata_data.metadata_dict
        """
        return self.method_metadata

