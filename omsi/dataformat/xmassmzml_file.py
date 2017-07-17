"""
This module provides functionality for reading mzml mass spectrometry image files.

filename = '/Users/oruebel/Devel/openmsi-data/mzML_Data/N2A2_Serratia_spots_extract_TI.mzML'
"""
import re
import os

import numpy as np

from pyteomics import mzml

from omsi.dataformat.file_reader_base import file_reader_base_multidata
from omsi.datastructures.dependency_data import dependency_dict
from omsi.datastructures.metadata.metadata_data import metadata_dict, metadata_value
from omsi.datastructures.metadata.metadata_ontologies import METADATA_ONTOLOGIES
from omsi.shared.log import log_helper


class xmassmzml_file(file_reader_base_multidata):
    """
    Interface for reading a single 2D mzml file with several distinct scan types.

    :ivar available_mzml_types: Dict of available mzml flavors.
    """

    def __init__(self, basename, requires_slicing=True, resolution=5):
        """
        Open an img file for data reading.

        :param basename: The name of the mzml file. If basename is a directory, then the first mzML file found
                             in the directory will be used instead.
        :type basename: string

        :param requires_slicing: Should the complete data be read into memory
                             (this makes slicing easier). (default is True)
        :type requires_slicing: bool

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
        # self.basename = basename
        # self.requires_slicing = requires_slicing
        
        # Call super constructor. This sets self.basename and self.readall
        super(xmassmzml_file, self).__init__(basename=basename, requires_slicing=requires_slicing)
        self.resolution = resolution
        self.data_type = 'uint32'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans(filename=self.basename)
        log_helper.info(__name__, 'Read %s scans from mzML file.' % self.num_scans)
        log_helper.debug(__name__, 'Compute coordinates')
        self.coordinates = self.__compute_coordinates(filename=self.basename,num_scans=self.num_scans)
        # Compute the spatial configuration of the matrix
        self.x_pos = np.unique(self.coordinates[:, 0])
        self.y_pos = np.unique(self.coordinates[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        # Compute the mz axis
        log_helper.debug(__name__, 'Compute mz axes')
        self.mz = self.__compute_mz_axis(filename=self.basename)
        log_helper.debug(__name__, 'mz axes computed')

        # Determine the shape of the dataset, result is a list of shapes for each datacube
        # self.shape_all_data = [(self.x_pos.shape[0], self.y_pos.shape[0], mz.shape[0]) for mz in self.mz_all

        log_helper.debug(__name__, 'Compute shape')
        self.shape = (self.x_pos.shape[0], self.y_pos.shape[0], len(self.mz))#self.shape[0])
        # self.shape = None
        # self.mz = None

        # Read the data into memory
        # self.data = None
        log_helper.debug(__name__, 'read all')
        if requires_slicing:
            self.data = self.__read_all()
        log_helper.debug(__name__, 'Finished with init')


    def __read_all(self):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.  Data is now a list of datacubes
        """

        # self.data = np.zeros(shape=self.shape_all_data[scan_idx], dtype=self.data_type) for scan_idx, scantype in enumerate(self.scan_types)
        data = np.zeros(self.shape)
        reader = mzml.read(self.basename)
        spectrumid = 0
        # if not self.scan_profiled[scan_idx]:
        #     shift = np.diff(self.mz_all[scan_idx]).mean()
        #     bin_edges = np.append(self.mz_all[scan_idx], self.mz_all[scan_idx][-1]+ shift)
        # else:
        #     bin_edges = None

        for spectrum in reader:
            # if spectrum['scanList']['scan'][0]['filter string'] == scantype:
            x = spectrum['m/z array']
                # try:
            y = spectrum['intensity array']
                # except KeyError:
                    # raise KeyError
                # if bin_edges is None:
            yi = np.interp(self.mz_all[scan_idx], x, y, 0, 0)  # Re-interpolate the data in profiled mode
                # else:
                     # yi, _ = np.histogram(x, bins=bin_edges, weights=y)   # Re-histogram the data in centroided mode
                # xidx = np.nonzero(self.x_pos == self.coordinates[spectrumid, 0])[0]
                # yidx = np.nonzero(self.y_pos == self.coordinates[spectrumid, 1])[0]
                # try:
            data[self.coordinates[spectrumid, 0], self.coordinates[spectrumid, 1], :] = yi
                # except:
                    # log_helper.debug(__name__, spectrumid, scan_idx, scantype, self.mz_all[scan_idx].shape)
        # TODO Note if the data is expected to be of float precision then self.data_type needs to be set accordingly
            # if spectrumid%1000 == 0:
                # log_helper.info(__name__, 'Processed data for %s spectra to datacube for scan type %s' % (spectrumid, scantype))
            spectrumid += 1
        return data

    def spectrum_iter(self):
        """
        Generator function that yields a position and associated spectrum for a selected datacube type.

        :yield: (xidx, yidx) a tuple of ints representing x and y position in the image
        :yield: yi,          a numpy 1D-array of floats containing spectral intensities at the given position \
                             and for the selected datacube type

        """
        reader = mzml.read(self.basename)
        if self.select_dataset is None:
           raise ValueError('Select a dataset to continue!')
        dataset_index = self.select_dataset
        for idx, spectrum in enumerate(reader):
            mz = self.mz_all[0]
            x = spectrum['m/z array']
            try:
                y = spectrum['intensity array']
            except KeyError:
                raise KeyError('Key "intensity array" not found in this mzml file')

            yi = np.interp(mz, x, y, 0, 0)      # Interpolate the data onto the new axes in profiles mode
            # else:
            #     shift = np.diff(mz).mean()
            #     bin_edges = np.append(mz, mz[-1]+ shift)
            #     yi, _ = np.histogram(x, bins=bin_edges, weights=y)   # Re-histogram the data in centroided mode
            xidx = np.nonzero(self.x_pos == self.coordinates[idx, 0])[0][0]
            yidx = np.nonzero(self.y_pos == self.coordinates[idx, 1])[0][0]

            yield (xidx, yidx), yi


    @classmethod
    def __compute_mz_axis(cls, filename):
        ## TODO completely refactor this to make it smartly handle profile or centroid datasets
        ## TODO: centroid datasets should take in a user parameter "Resolution" and resample data at that resolution
        ## TODO: profile datasets should work as is
        ## TODO: checks for profile data vs. centroid data on the variation in length of ['m/z array']
        """
        Internal helper function used to compute the mz axis of each scantype
        Returns a list of numpy arrays
        """
        reader = mzml.read(filename)
        mz_list = []
        counter = 0
        for spectrum in reader:
            mz_list.append(np.asarray(spectrum['m/z array']))
            counter += 1
        mzdiff = 10000000.0
        mz_min = 10000000.0
        mz_max = 0.0
        for mz in mz_list:
            d  = np.diff(mz).min()
            if d < mzdiff:
                mzdiff = d
            m = mz.min()
            if m < mz_min:
                mz_min = m
            m = mz.max()
            if m > mz_max:
                mz_max = m
            
        mz_axes = np.arange(start=mz_min, stop=mz_max, step=mzdiff)
        return mz_axes

    @classmethod
    def __compute_coordinates(self,filename,num_scans):
        """
        Internal helper function used to compute the coordinates for each scan.

        :returns: 2D numpy integer array of shape (numScans,2) indicating for each scan its x and y coordinate
        """
        spectrumid = 0
        reader = mzml.read(filename)
        coords = np.zeros(shape=(num_scans, 2), dtype='uint32')
        with open(filename,'r') as origin_file:
            for line in origin_file:
                s = re.findall(r'location="', line)
                if s:
                    m = re.search(r'_[0-9]+x_[0-9]+y_', line,)
                    if m:
                        coord_str = m.group()
                        coord_str = coord_str.strip('_').split('_')
                        coord_str = [int(c[:-1]) for c in coord_str]
                        coords[spectrumid, 0] = coord_str[0]
                        coords[spectrumid, 1] = coord_str[1]
                        spectrumid += 1
        return coords


    # @classmethod
    # def test(cls):
    #     """
    #     Test method
    #     """
    #     pass

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
            if self.select_dataset is None:
                return self.data[key]
            else:
                return self.data[self.select_dataset][key]
        else:
            raise ValueError("Slicing is currently only supported when the object has been initialized with readall")

    def close_file(self):
        """Close the mzml file"""
        pass

    # @classmethod
    # def is_valid_dataset(cls, name):
    #     """Check whether the given file or directory points to a img file.

    #        :param name: Name of the dir or file.
    #        :type name: String

    #        :returns: Boolean indicating whether the given file or folder is a valid img file.
    #     """
    #     if os.path.isdir(name):  # If we point to a directory, check if the dir contains an mzML file
    #         filelist = cls.get_files_from_dir(name)
    #         return len(filelist) > 0
    #     else:
    #         try:
    #             # Try to open the file and iterate over it
    #             reader = mzml.read(name)
    #             for _ in reader:
    #                 pass
    #             del reader
    #             return True
    #         except:
    #             return False

    # @classmethod
    # def size(cls, name, max_num_reads=1000):
    #     """
    #     Classmethod used to check the estimated size for the given file/folder.
    #     For mzml this is an estimate of the final size of the full 3D datacube.
    #     For efficiency the number of scans is estimated based on the size of
    #     the first 1000 scans.

    #     :param name: Name of the dir or file.
    #     :type name: unicode
    #     :param max_num_reads: The maximum number of spectrum reads to be performed to estimate the file size
    #     :type max_num_reads: int

    #     :returns: Integer indicating the size in byte or None if unknown.
    #     """
    #     basename = None
    #     if os.path.isdir(name):  # If we point to a directory, check if the dir contains an mzML file
    #         filelist = cls.get_files_from_dir(name)
    #         if len(filelist) > 0:
    #             basename = filelist[0]
    #     else:
    #         basename = name
    #     if basename is not None:
    #         num_scans = -1
    #         # Try to compute the number of scans by looking at the spectrumList count entry in the file
    #         try:
    #             size_line = os.popen('head -n 120 "' + basename + '" | grep "spectrumList count="').read()
    #             if len(size_line) > 0:
    #                 size_text = size_line.split('spectrumList count=')[1].split('"')[1]
    #                 if size_text.isdigit():
    #                     num_scans = int(size_text)
    #         except:
    #             pass
    #         if num_scans < 0:
    #             # Estimate the number of scans by reading the first 1000 spectra
    #             index = 0
    #             prev_tell = 0
    #             sizes = []
    #             reader = mzml.read(basename)
    #             for _ in reader:
    #                 if index >= max_num_reads:
    #                     break
    #                 current_tell = reader.file.file.tell()
    #                 sizes.append(current_tell - prev_tell)
    #                 prev_tell = current_tell
    #                 index += 1
    #             npsizes = np.asarray(sizes)
    #             filesize = os.stat(basename).st_size
    #             scansize = (npsizes.max() - npsizes.min()) / 2.
    #             num_scans = int(filesize/scansize)
    #         mz_axis_len = cls.__compute_mz_axis(filename=basename,
    #                                             mzml_filetype=cls.__compute_filetype(filename=basename),
    #                                             scan_types=cls.__compute_scan_types(filename=basename)).shape[0]
    #         return num_scans*mz_axis_len

    #         # temp_mzml_file = cls(basename=basename, requires_slicing=False)
    #         # itemsize = np.dtype(temp_mzml_file.data_type).itemsize
    #         # size = np.asarray(temp_mzml_file.shape).prod() * itemsize
    #         # print ('MZML size', size)
    #         # return size
    #     else:
    #         return None

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

