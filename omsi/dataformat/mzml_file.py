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


class mzml_file(file_reader_base_multidata):
    """
    Interface for reading a single 2D mzml file with several distinct scan types.

    :ivar available_mzml_types: Dict of available mzml flavors.
    """

    available_mzml_types = {'unknown': 'unknown',
                            'bruker': 'bruker',
                            'thermo': 'thermo'}

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

        # Call super constructor. This sets self.basename and self.readall
        super(mzml_file, self).__init__(basename=basename, requires_slicing=requires_slicing)
        self.resolution = resolution
        self.mzml_type = self.__compute_filetype(filename=self.basename)
        self.data_type = 'uint32'  # TODO What data type should we use for the interpolated data?
        self.num_scans = self.__compute_num_scans(filename=self.basename)
        log_helper.info(__name__, 'Read %s scans from mzML file.' % self.num_scans)
        log_helper.debug(__name__, 'Compute scan types and indicies')
        self.scan_types, self.scan_indices, self.scan_profiled = self.__compute_scan_types_and_indices(filename=self.basename)
        log_helper.debug(__name__, 'Compute scan dependencies')
        self.scan_dependencies = self.__compute_scan_dependencies(scan_types=self.scan_types,
                                                                  basename=basename)
        log_helper.info(__name__, 'Found %s different scan types in mzML file.' % len(self.scan_types))
        log_helper.debug(__name__, 'Compute coordinates')
        self.coordinates = self.__compute_coordinates()
        log_helper.debug(__name__, 'Parse scan parameters')
        self.scan_params = self.__parse_scan_parameters(self)

        # Compute the spatial configuration of the matrix
        self.x_pos = np.unique(self.coordinates[:, 0])
        self.y_pos = np.unique(self.coordinates[:, 1])
        self.step_size = min([min(np.diff(self.x_pos)), min(np.diff(self.y_pos))])

        # Compute the mz axis
        log_helper.debug(__name__, 'Compute mz axes')
        self.mz_all = self.__compute_mz_axis(filename=self.basename,
                                             mzml_filetype=self.mzml_type,
                                             scan_types=self.scan_types,
                                             resolution=self.resolution)

        # Determine the shape of the dataset, result is a list of shapes for each datacube
        self.shape_all_data = [(self.x_pos.shape[0], self.y_pos.shape[0], mz.shape[0]) for mz in self.mz_all]
        self.shape = None
        self.mz = None

        # Read the data into memory
        self.data = None
        if requires_slicing:
            self.__read_all()

    def __read_all(self):
        """
        Internal helper function used to read all data. The
        function directly modifies the self.data entry.  Data is now a list of datacubes
        """

        self.data = [np.zeros(shape=self.shape_all_data[scan_idx], dtype=self.data_type) for scan_idx, scantype in enumerate(self.scan_types)]

        for scan_idx, scantype in enumerate(self.scan_types):
            reader = mzml.read(self.basename)
            spectrumid = 0
            if not self.scan_profiled[scan_idx]:
                shift = np.diff(self.mz_all[scan_idx]).mean()
                bin_edges = np.append(self.mz_all[scan_idx], self.mz_all[scan_idx][-1]+ shift)
            else:
                bin_edges = None

            for spectrum in reader:
                if spectrum['scanList']['scan'][0]['filter string'] == scantype:
                    x = spectrum['m/z array']
                    try:
                        y = spectrum['intensity array']
                    except KeyError:
                        raise KeyError
                    if bin_edges is None:
                        yi = np.interp(self.mz_all[scan_idx], x, y, 0, 0)  # Re-interpolate the data in profiled mode
                    else:
                         yi, _ = np.histogram(x, bins=bin_edges, weights=y)   # Re-histogram the data in centroided mode
                    xidx = np.nonzero(self.x_pos == self.coordinates[spectrumid, 0])[0]
                    yidx = np.nonzero(self.y_pos == self.coordinates[spectrumid, 1])[0]
                    try:
                        self.data[scan_idx][xidx, yidx, :] = yi
                    except:
                        log_helper.debug(__name__, spectrumid, scan_idx, scantype, self.mz_all[scan_idx].shape)
            # TODO Note if the data is expected to be of float precision then self.data_type needs to be set accordingly
                if spectrumid%1000 == 0:
                    log_helper.info(__name__, 'Processed data for %s spectra to datacube for scan type %s' % (spectrumid, scantype))
                spectrumid += 1

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
            current_scan_index = self.scan_indices[idx]

            if current_scan_index == dataset_index:
                # get the mz axis of the datacube the scan belongs to
                mz = self.mz_all[self.scan_indices[idx]]
                x = spectrum['m/z array']
                try:
                    y = spectrum['intensity array']
                except KeyError:
                    raise KeyError('Key "intensity array" not found in this mzml file')

                if self.scan_profiled[current_scan_index]:
                    yi = np.interp(mz, x, y, 0, 0)      # Interpolate the data onto the new axes in profiles mode
                else:
                    shift = np.diff(mz).mean()
                    bin_edges = np.append(mz, mz[-1]+ shift)
                    yi, _ = np.histogram(x, bins=bin_edges, weights=y)   # Re-histogram the data in centroided mode
                xidx = np.nonzero(self.x_pos == self.coordinates[idx, 0])[0][0]
                yidx = np.nonzero(self.y_pos == self.coordinates[idx, 1])[0][0]

                yield (xidx, yidx), yi


    @classmethod
    def __compute_mz_axis(cls, filename, mzml_filetype, scan_types, resolution):
        ## TODO completely refactor this to make it smartly handle profile or centroid datasets
        ## TODO: centroid datasets should take in a user parameter "Resolution" and resample data at that resolution
        ## TODO: profile datasets should work as is
        ## TODO: checks for profile data vs. centroid data on the variation in length of ['m/z array']
        """
        Internal helper function used to compute the mz axis of each scantype
        Returns a list of numpy arrays
        """
        reader = mzml.read(filename)

        if mzml_filetype == cls.available_mzml_types['thermo']:

            mz_axes = [np.array([]) for _ in scan_types]
            # all_centroid = True
            for spectrum in reader:
                scanfilt = spectrum['scanList']['scan'][0]['filter string']
                scantype_idx = scan_types.index(scanfilt)
                mz = spectrum['m/z array']
                try:
                    len_axes = len(mz_axes[scantype_idx])
                except TypeError:
                    len_axes = 1
                if spectrum.has_key('profile spectrum'):
                    # all_centroid = False
                    if len(mz) > len_axes:
                        mzdiff = np.diff(mz).min()
                        mzmin = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
                        mzmax = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']
                        mz_axes[scantype_idx] = np.arange(start=mzmin, stop=mzmax, step=mzdiff)
                        mz_axes[scantype_idx] = np.append(arr=mz_axes[scantype_idx], values=mzmax)
                else:
                    if len(mz) > len_axes:
                        mzmin = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
                        mzmax = spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']
                        f = np.ceil(1e6 * np.log(mzmax/mzmin)/resolution)
                        mz_axes[scantype_idx] = np.logspace(np.log10(mzmin), np.log10(mzmax), f)
                        # ['count', 'index', 'highest observed m/z', 'm/z array', 'total ion current', 'ms level', 'spotID', 'lowest observed m/z', 'defaultArrayLength', 'intensity array', 'centroid spectrum', 'positive scan', 'MS1 spectrum', 'spectrum title', 'base peak intensity', 'scanList', 'id', 'base peak m/z']

            return mz_axes

        # assume bruker instruments have constant m/z axis from scan to scan
        elif mzml_filetype == cls.available_mzml_types['bruker']:

            mz_axes = [np.array([]) for _ in scan_types]

            for spectrum in reader:
                scanfilt = spectrum['scanList']['scan'][0]['filter string']
                scantype_idx = scan_types.index(scanfilt)
                # grossly inefficient reassignment of m/z array at each scan
                mz_axes[scantype_idx] = spectrum['m/z array']

            return mz_axes

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


    @classmethod
    def test(cls):
        """
        Test method
        """
        pass

    @staticmethod
    def __compute_num_scans(filename=None):
        """
        Internal helper function used to compute the number of scans in the mzml file.
        """
        reader = mzml.read(filename)
        return sum(1 for _ in reader)

    def __compute_scan_types_and_indices(self, filename=None):
        """
        Internal helper function used to compute a list of unique scan types in the mzml file.
        Also computes a numpy 1d array of ints which index every scan to relevant datacube.
        """
        reader = mzml.read(filename)
        scantypes = []
        scan_indices = []
        scan_profiled = []
        for idx, spectrum in enumerate(reader):
            try:
                scanfilter = spectrum['scanList']['scan'][0]['filter string']
                if scanfilter not in scantypes:
                    scantypes.append(scanfilter)
                    scan_profiled.append(spectrum.has_key('profile spectrum'))
                scan_indices.append(scantypes.index(scanfilter))
            except:
                log_helper.debug(__name__, idx)

        assert len(scan_indices) == self.num_scans
        return scantypes, scan_indices, scan_profiled

    @staticmethod
    def __parse_scan_parameters(self):
        """
        Internal helper function used to parse out scan parameters from the scan filter string
        """
        # precursor m/z
        # example scan filter: MS2: ITMS + p MALDI Z ms2 907.73@cid60.00 [500.00-700.00]
        ## example scan filter: MS1: FTMS + p MALDI Full ms [850.00-1000.00]

        scan_params = []

        for scan_idx, scantype in enumerate(self.scan_types):

            #MSnValueOfN
            n = filter(None, re.findall('(?<=ms)\d*', scantype))
            if n:
                MSnValueOfN = int(n[0])
            else:
                MSnValueOfN = 1

            #precursor
            ms2pre = filter(None, re.findall('[\d.]+(?=@)', scantype))
            if ms2pre:
                ms2_precursor = float(ms2pre[0])
            else:
                ms2_precursor = None
            #dissociation type
            dissot = filter(None, re.findall('(?<=\d@)[A-z]*', scantype))
            if dissot:
                dissociationtype = dissot[0]
            else:
                dissociationtype = 'None'
            #dissociation energy
            dissoe = filter(None, re.findall('(?<='+dissociationtype+')'+'[\d.]+', scantype))
            if dissoe:
                dissociationenergy = float(dissoe[0])
            else:
                dissociationenergy = None
            #polarity
            pol = filter(None, re.findall('([+][ ]p)', scantype))
            if pol:
                polarity = 'pos'
            else:
                pol = filter(None, re.findall('([-][ ]p)', scantype))
                if pol:
                    polarity = 'neg'
                else:
                    polarity = 'unk'
            #put all params in dictionary
            paramdict = metadata_dict()
            msn_von_ontology = METADATA_ONTOLOGIES['msn_value_of_n']
            paramdict['msn_value_of_n'] = metadata_value(name='msn_value_of_n',
                                                         value=MSnValueOfN,
                                                         unit=msn_von_ontology['unit'],
                                                         description=msn_von_ontology['description'],
                                                         ontology=msn_von_ontology)
            if dissociationenergy:
                paramdict['dissociation_energy'] = metadata_value(name='dissociation_energy',
                                                                  value=dissociationenergy,
                                                                  unit='V',
                                                                  description='Dissociation energy')
            if ms2_precursor is not None:
                paramdict['msn_precursor_mz'] = metadata_value(name='msn_precursor_mz',
                                                               value=ms2_precursor,
                                                               unit='m/z',
                                                               description='The precursor m/z value')
            paramdict['dissociation_type'] = metadata_value(name='dissociation_type',
                                                            value=dissociationtype,
                                                            unit=None,
                                                            description='Dissociation type')
            polarity_ontology = METADATA_ONTOLOGIES['polarity']
            paramdict['polarity'] = metadata_value(name='polarity',
                                                   value=polarity,
                                                   unit=polarity_ontology['unit'],
                                                   description=polarity_ontology['description'],
                                                   ontology=polarity_ontology)

            scan_params.append(paramdict)

        return scan_params

    @staticmethod
    def __compute_scan_dependencies(scan_types=None,
                                    basename=None):
        """
        Takes a scan_types list and returns a list of tuples (x, y) indicating that scan_type[y] depends on scan_type[x]
        """
        dependencies = [[] for i in range(len(scan_types))]

        #find MS1 scans   # TODO this algorithm could be much smarter; now will work only on single MS scans

        ms1scanlist = []

        for ind, stype in enumerate(scan_types):
            if stype.find('Full ms') != -1:    #Regular
                ms1scanlist.append(ind)

        #make MS2 scans dependent on MS1 scans
        for ms1scan in ms1scanlist:
            for ind2, stype in enumerate(scan_types):
                if stype.find('ms2') != -1:     #MS2 scan filter strings from thermo contain the string 'ms2'
                    ms2_to_ms1_dependency = {
                        'omsi_object': None,
                        'link_name': 'MS1',
                        'basename': basename,
                        'region': None,
                        'dataset': ms1scan,
                        'help': scan_types[ms1scan],
                        'dependency_type': dependency_dict.dependency_types['co_modality']}
                    ms2_link_name = 'MS2_' + scan_types[ind2].split('ms2')[-1].lstrip(' ').rstrip(' ').replace(' ', '_')
                    ms1_to_ms2_dependency = {
                        'omsi_object': None,
                        'link_name': ms2_link_name,
                        'basename': basename,
                        'region': None,
                        'dataset': ind2,
                        'help':scan_types[ms1scan],
                        'dependency_type': dependency_dict.dependency_types['co_modality']}
                    dependencies[ind2].append(ms2_to_ms1_dependency)
                    dependencies[ms1scan].append(ms1_to_ms2_dependency)

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
                del reader
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
                                                mzml_filetype=cls.__compute_filetype(filename=basename),
                                                scan_types=cls.__compute_scan_types(filename=basename)).shape[0]
            return num_scans*mz_axis_len

            # temp_mzml_file = cls(basename=basename, requires_slicing=False)
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

        Inherited from
        """
        if self.select_dataset is not None:
            return self.scan_dependencies[self.select_dataset]
        else:
            return self.scan_dependencies

    def get_dataset_metadata(self):
        """
        Get dict of additional metadata associated with the current dataset.

        Inherited from file_reader_base.file_reader_base_multidata.

        :return: Dict where keys are strings and associated values to be stored as
            metadata with the dataset.

        """
        if self.select_dataset is not None:
            return self.scan_params[self.select_dataset]
        else:
            return self.scan_params
