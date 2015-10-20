"""This module provides functionality for reading bruker flex mass spectrometry image files

Limitations:

1) Currently the reader assumes a single global m/z axis for all spectra.
2) The read of acqu files does not convert <..> entries to python values but leaves them as strings.
3) __read_spotlist__ converts the regions to start with a 0 index. This is somewhat inconsistent in the spot list file.\
  The spotname seems to number regions starting with 0 while the region list numbers them starting with 1.
4) __read_spotlist__ computes the folder where the spots are located based on the filename of the spotlist. \
   The question is whether this is always the case??? The advantage is that we do not rely on the regions.xml \
   file which contains absolute paths which are in most cases invalid as the data has been copied between different \
   machines in many cases and is stored in different locations on each of the machines.
5) __read_spotlist__ currenly assumes that there is only one fid file per spot
6) __read_spotlist__ currenlty only looks at where the acqu and fid file is located for the spot. Other files \
  are currently ignored.
7) __read_spotlist__ (and hence the reader at large) currently assumes that we have 2D images only.
8) __read_spotlist__ currently generates maps for the image that assume that x and y pixel indices start at \
   0. Not all images may record data until the border, so that this strategy may add empty spectra rather than \
   generating a new bounding box for the image.
9) __read_spotlist__ assumes in the variable spotname_encoding a maximum of 24 characters in the spotname \
   R01X080Y013. This should in general be more than sufficient as this allows for 7 characters for each R, X, \
   Y entry, however, if this is not enough then this behaviour needs ot be changed.
10) __getitem__ currently only works if we have read the full data into memory. An on-demand load should be supported \
    as well.
11) We can currently only selected either a single region or the full data but we cannot selected multiple regions \
    at once. E.g. if a dataset contains 3 regions then we can either select all regions at once or region 1,2, or 3 \
    but one cannot selected region 1+2, 1+3, or 2+3.

import bruckerflex_file
spotlist = "/Users/oruebel/Devel/msidata/Bruker_Data/UNC IMS Data/20130417 Bmycoides Paenibacillus Early SN03130/" + \
           "2013 Bmyc Paeni Early LP/2013 Bmyc Paeni Early LP Spot List.txt"
exppath ="/Users/oruebel/Devel/msidata/Bruker_Data/UNC IMS Data/20130417 Bmycoides Paenibacillus Early SN03130/" + \
        "2013 Bmyc Paeni Early LP/2013 Bmyc Paeni Early LP/0_R00X012Y006/1/1SLin"
f = bruckerflex_file.bruckerflex_file( spotlist_filename = spotlist)
f.s_read_fid( exppath+"/fid" , f.data_type )
testacqu = f.s_read_acqu( exppath+"/acqu" )
testmz   = f.s_mz_from_acqu( testacqu )
testspotlist = f.s_read_spotlist(spotlist)


.....
from bruckerflex_file import *
dirname = "/Users/oruebel/Devel/msidata/Bruker_Data/UNC IMS Data/20130417 Bmycoides Paenibacillus Early SN03130/" + \
          "20130417 Bmyc Paeni Early LN/20130417 Bmyc Paeni Early LN"
a = bruckerflex_file( dirname )

"""

import os
import math

import numpy as np

from omsi.datastructures.dependency_data import dependency_dict
from omsi.dataformat.file_reader_base import file_reader_base_with_regions
from omsi.shared.log import log_helper


class bruckerflex_file(file_reader_base_with_regions):
    """Interface for reading a single bruker flex image file.

       The reader supports standard array slicing for data read. I.e., to read a
       spectrum use [x,y,:] to read an ion image using [:,:,mz].

       The reader supports multiple regions, i.e., reading of different independent
       regions that were imaged as part of the same dataset. Using the get_regions, get_number_of_regions,
       set_region_selection and get_region_selection the user can interact with the
       region settings. Using the set_region_selection the user can define whether
       the data of the complete image should be read (set_region(None), default) or
       whether data from a single region should be read. When a region is selected,
       then the reader acts as if the region were the complete dataset, i.e., the
       shape variable is addjusted to fit the selected region and the __get_item__
       method (which is used to implement array-like slicing (e.g., [1,1,:])) behaves
       as if the selected region where the full data.

    """

    def __init__(self, basename, fid_encoding='int32', requires_slicing=True):
        """Open an img file for data reading.

            :param basename: Name of the textfile with the spotlist. Alternatively this may also be the \
                                      folder with the spots.
            :type basename: string
            :param requires_slicing: Should the complete data be read into memory (this makes slicing easier). (default is True)
            :type requires_slicing: bool
            :param fid_encoding: String indicating in which binary format the intensity values (fid files) are stored. \
                                 (default value is 'int32')
            :type fid_encoding: string

            :var self.basename: Name of the file with the spotlist
            :var self.pixel_dict: dictionary with the pixel array metadata (see also s_read_spotlist(...)). Some of \
                 the main keys of the dictionary are, e.g. (see also s_read_spotlist(...)) :

                * 'spotfolder' : String indicating the folder where all the spot-data is located
                * 'fid' : 2D numpy masked array of strings indicating for (x,y) pixel the fid file with  \
                          the intensity values.
                * 'acqu' : 2D numpy masked array of strings indicating for each (x,y) pixel the acqu file with the \
                           metadata for that pixel.
                * 'regions' : 2D numpy masked array of integers indicating for each pixels the index of the region \
                           it belongs to.
                * 'xpos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'ypos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'spotname' : 2D masked numpy array of strings with the names of the spot corresponding to a pixel.

            :var self.data_type: the encoding used for intensity values
            :var self.shape: The 3D shape of the MSI data volume for the currently selected region.
            :var self.full_shape: Shape of the full 3D MSI dataset including all regions imaged.
            :var self.metadata: Dictionary with metadata from the acqu file
            :var self.mz: The 1D numpy array with the m/z axis information
            :var self.data: If requires_slicing is set to true then this 3D array includes the complete data of the MSI data \
                  cube. Missing data values (e.g., from regions not imaged during the aquistion processes) are \
                  completed with zero values.
            :var self.region_dicts: Dictionary with description of the imaging regions

            :raises ValueError: In case that no valid data is found.
        """
        super(bruckerflex_file, self).__init__(basename, requires_slicing)
        self.spotlist_filename = basename
        if os.path.isdir(self.spotlist_filename):
            self.spotlist_filename = None
            self.pixel_dict = self.s_spot_from_dir(basename)
            if not self.pixel_dict:
                raise ValueError(
                    "Invalid bruckerflex file. No valid spots found in the given directory")
        else:
            self.pixel_dict = self.s_read_spotlist(self.spotlist_filename)
        self.data_type = fid_encoding  # Data type of the fid files

        # ToDo this reads a single spectrum (ie., fid file) to figure out the
        # length of the spectra. This obviously means that we assume a single
        # m/z axis for all spectra
        temp_spectrum = self.s_read_fid(self.pixel_dict['fid'].compressed()[0],
                                        self.data_type)
        self.full_shape = [self.pixel_dict['regions'].shape[0],
                           self.pixel_dict['regions'].shape[1], temp_spectrum.shape[0]]
        self.shape = [self.pixel_dict['regions'].shape[0],
                      self.pixel_dict['regions'].shape[1],
                      temp_spectrum.shape[0]]  # Number of pixels in x,y, and z

        # ToDo this reads a single acqu metadata file to figure out the m/z
        # axis. This obviously means that we assume a single m/z axis for all
        # spectra
        self.metadata = self.s_read_acqu(self.pixel_dict['acqu'].compressed()[0])
        self.mz = self.s_mz_from_acqu(self.metadata)

        # Compute the region bouding boxes
        self.region_dicts = self.__compute_regions_extends___(self.pixel_dict)
        self.select_region = None

        # Should we read all the intensity values into memory. If so, allocate
        # space.
        if requires_slicing:
            self.data = np.zeros(self.shape)
            dmask = self.pixel_dict['fid'].mask
            for xindex in range(0, self.shape[0]):
                for yindex in range(0, self.shape[1]):
                    if not dmask[xindex, yindex]:
                        self.data[xindex, yindex, :] = self.s_read_fid(
                            self.pixel_dict['fid'][xindex, yindex], self.data_type)
        else:
            self.data = None

        # A numpy vector with the m/z values of the instrument
        self.instrument_mz = 0

    @classmethod
    def is_valid_dataset(cls, name):
        """Determine whether the given file or name specifies a bruckerflex file

           :param name: name of the file or dir
           :type name: string

           :return: Boolean indicating whether the name is a valid bruckerflex
        """
        if os.path.isfile(name) and name.endswith("Spot List.txt"):
            return True
        elif os.path.isdir(name):
            # Check if there are any spots in the directory
            spotlist = cls.s_spot_from_dir(name, spot_folder_only=True)
            if spotlist is not None:
                if len(spotlist) > 0:
                    return True
        return False

    def set_region_selection(self, region_index=None):
        """Define which region should be selected for local data reads.

           :param region_index: The index of the region that should be read. The shape of the
                    data will be adjusted accordingly. Set to None to select all regions and treat
                    the data as a single full 3D image.
        """
        if region_index is None:
            self.shape = self.full_shape
            self.select_region = None
        elif region_index < self.get_number_of_regions():
            self.select_region = region_index
            self.shape = [self.region_dicts[self.select_region]["extend"][0],
                          self.region_dicts[self.select_region]["extend"][1], self.full_shape[2]]

    def __getitem__(self, key):
        """Enable slicing of bruker files"""
        # If the full data has been loaded
        if self.data is not None:
            # Select the region of interest, i.e., the reader acts as if the
            # selected region where the complete image of interest.
            if self.select_region is not None:
                rd_origin = self.region_dicts[self.select_region]["origin"]
                rd_extend = self.region_dicts[self.select_region]["extend"]
                low_x = rd_origin[0]
                low_y = rd_origin[1]
                high_x = low_x + rd_extend[0] + 1
                high_y = low_y + rd_extend[1] + 1
                # Select the region and then selection the user data within
                # that region
                return self.data[low_x:high_x, low_y:high_y, :][key]
            else:
                # Return the data from the full image data.
                return self.data[key]
        else:
            raise ValueError(
                "Slicing is only supported when the object has been initialized with readall")

    def close_file(self):
        """Close the img file"""
        pass

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
        dependency_list = []
        if self.select_region is None:
            for region_index, region in enumerate(self.region_dicts):
                main_to_region = {'omsi_object': None,
                                  'link_name': 'msi_region_' + str(region_index),
                                  'basename': self.basename,
                                  'region': region_index,
                                  'dataset': None,
                                  'help': 'Dependency pointing to a subregion of the MSI dataset',
                                  'dependency_type': dependency_dict.dependency_types['contains']}
                dependency_list.append(main_to_region)
        else:
            region_to_main = {'omsi_object': None,
                              'link_name': 'msi_full',
                              'basename': self.basename,
                              'region': None,
                              'dataset': None,
                              'help': 'Dependency pointing to the full MSI dataset this dataset is a sub-region of',
                              'dependency_type': dependency_dict.dependency_types['subset']}
            dependency_list.append(region_to_main)
        return dependency_list

    def __read_spectrum__(self, filename, selection):
        """Read data from a single fid file with the intensities.

           :param filename: String indicating the name+path to the fid file.
           :type filename: string
           :param selection: This may be a python slice or a list of indecies to be read.
           :type selection: slice or list, i.e., a selection that numpy understands

           :returns:
            -- 1D numpy array of intensity values read from the file.
            -- 1D numpy array of the mz values for the part of the spectrum read.
        """
        # Read the requested intensity values from the fid file
        intensity = self.s_read_fid(filename, self.data_type, selection)

        # Construct the m/z axis for an acqu file
        # Get the directory where the file is located
        dirname = os.path.dirname(os.path.abspath(filename))
        acqu_filename = dirname + "/acqu"
        acqu_dict = self.s_read_acqu(acqu_filename)
        curr_mz = self.s_mz_from_acqu(acqu_dict)
        return intensity, curr_mz[selection]

    @staticmethod
    def __compute_regions_extends___(pixel_maps):
        """Compute the region extends from the given region map.

            :param pixel_maps: This is the self.pixel_dict which is computed in the __init__ function using
                               the s_read_spotlist function.

            :returns:
        """
        region_map = pixel_maps['regions']
        minr = np.min(region_map)
        maxr = np.max(region_map)
        regions = []
        for region_index in range(minr, maxr + 1):
            region_location = np.argwhere(region_map == region_index)
            min_x = np.min(region_location[:, 0])
            min_y = np.min(region_location[:, 1])
            max_x = np.max(region_location[:, 0])
            max_y = np.max(region_location[:, 1])
            extend_x = max_x - min_x
            extend_y = max_y - min_y
            bounding_box = {"origin": np.array([min_x, min_y]), "extend": np.array([extend_x, extend_y])}
            regions.append(bounding_box)

        return regions

    @staticmethod
    def s_spot_from_dir(in_dir, spot_folder_only=False):
        """Similar to  s_read_spotlist but instead of using a spotlist file the structure of the data is parsed
           directly from the structure of the direcory containint all spots.

            :param in_dir: Name of the directory with all spots
            :type in_dir: string
            :param spot_folder_only: If set to True, then the function only constructs the spot folders but
                                     does not check for acqu files etc. If set to True, only the spotfolder
                                     list will be returned.

            :returns: The function returns None in case that no valid spots were found. Returns a list of strings
                      with the spotfolders if spot_folder_only is set to True. Otherwise, the function returns a
                      number of different items in from of a python dictionary. Most data is stored as 2D spatial
                      maps, indicting for each (x,y) location the corresponding data. Most data is stored as 2D
                      masked numpy arrays. The masked of the array indicated whether data has been recorded for
                      a given pixel or not. The dict contains the following keys:

                * 'spotfolder' : String indicating the folder where all the spot-data is located
                * 'fid' : 2D numpy masked array of strings with fid file name for each (x,y) pixel (intensities)
                * 'acqu' : 2D numpy masked array of strings with acqu file name  for each (x,y) pixel (metadata)
                * 'regions' : 2D numpy masked array of integers with the index of the region for each (x,y) pixel.
                * 'xpos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'ypos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'spotname' : 2D masked numpy array of strings with the name of the spot corresponding to a pixel.

        """

        # Compute list of spots from folder structure
        # The base folder where the spots are supposed to be located
        spotfolder = os.path.abspath(in_dir)
        spotfolder_list = []  # List of folders with the spots
        spotname_list = []  # List of names for the spots
        pixel_x = []  # X location of each spot
        pixel_y = []  # Y location of each spot
        pixel_r = []  # Region index of each spot
        for l in os.listdir(in_dir):
            ap = os.path.abspath(os.path.join(in_dir, l))
            try:
                if os.path.isdir(ap):
                    rxsplit = l.split("R")[1].split("X")
                    pr = int(rxsplit[0])
                    sp = rxsplit[1].split("Y")
                    px = int(sp[0])
                    py = int(sp[1])
                    pixel_r.append(pr)
                    pixel_x.append(px)
                    pixel_y.append(py)
                    spotfolder_list.append(ap)
                    spotname_list.append(l)
            except:
                pass

        # Check if we need to generate anything else
        if spot_folder_only:
            return spotfolder_list

        # Check if we have any files at all
        if len(spotfolder_list) == 0:
            return None

        # Convert compute lists to numpy
        # numSpots = len(spotfolder_list)
        xpos = np.asarray(pixel_x)
        ypos = np.asarray(pixel_y)
        region = np.asarray(pixel_r)
        spotname = np.asarray(spotname_list)

        # Create maps for all entries
        xi_min = 0
        xi_max = np.max(xpos)
        yi_min = 0
        yi_max = np.max(ypos)
        xsize = xi_max - xi_min + 1
        ysize = yi_max - yi_min + 1

        # Create 2D image maps of the data (this is to make the reading process
        # easier)
        region_map = np.empty(shape=(xsize, ysize), dtype=region.dtype)
        region_map.fill(-1)
        region_map[xpos, ypos] = region[:]
        mask = (region_map < 0)
        region_map = np.ma.masked_array(region_map, mask)
        xpos_map = np.empty(shape=(xsize, ysize), dtype=xpos.dtype)
        xpos_map[xpos, ypos] = xpos[:]
        xpos_map = np.ma.masked_array(xpos_map, mask)
        ypos_map = np.empty(shape=(xsize, ysize), dtype=ypos.dtype)
        ypos_map[xpos, ypos] = ypos[:]
        ypos_map = np.ma.masked_array(ypos_map, mask)
        spotname_map = np.empty(shape=(xsize, ysize), dtype=spotname.dtype)
        spotname_map[xpos, ypos] = spotname[:]
        spotname_map = np.ma.masked_array(spotname_map, mask=mask)

        # Get a list of all files and compute and image map to indicate where the fid and acqu files
        # for each pixel are located
        # filelist = bruckerflex_file.get_all_files(spotfolder)
        # Compute the longest filename. This is only needed to determine the
        # encoding for the fix-size numpy array.
        longest_name_length = 0
        for path, subdirs, files in os.walk(spotfolder):
            for name in files:
                fn = os.path.join(path, name)
                if len(fn) > longest_name_length:
                    longest_name_length = len(fn)
        fltype = 'a' + str(longest_name_length)
        # ToDo: If needed we need to change this to allow for multiple fid
        # files per spot
        fidfiles = np.empty(shape=(xsize, ysize), dtype=fltype)
        # ToDo: Add additional storage for other files that may be needed
        acqufiles = np.empty(shape=(xsize, ysize), dtype=fltype)
        for path, subdirs, files in os.walk(spotfolder):
            for f in files:
                # ToDo: We are only interested in fid and acqu files right now
                # but that may need to change
                if f.endswith("fid") or f.endswith("acqu") and \
                   not (f.endswith("_fid") or f.endswith("_acqu")):
                    # Compute the pixelindex for the file
                    pindex = path.lstrip(spotfolder).lstrip(
                        "/").split("/")[0].split("X")[1].split("Y")
                    xi = int(pindex[0])
                    yi = int(pindex[1])
                    # ToDo: Again if we are interested in more than fid and
                    # acqu files then we need to extend this check
                    if f.endswith("fid") and not f.endswith("_fid"):
                        fidfiles[xi, yi] = os.path.join(path, f)
                    elif f.endswith("acqu") and not f.endswith("_acqu"):
                        acqufiles[xi, yi] = os.path.join(path, f)

        fidfiles = np.ma.masked_array(fidfiles, mask)
        acqufiles = np.ma.masked_array(acqufiles, mask)

        # Return all spatial maps of the different data we collected
        re = {'spotfolder': spotfolder,
              'fid': fidfiles,
              'acqu': acqufiles,
              'regions': region_map,
              'xpos': xpos_map,
              'ypos': ypos_map,
              'spotname': spotname_map}

        return re

    @staticmethod
    def s_read_spotlist(spotlist_filename):
        """Parse the given spotlist file.

            :param spotlist_filename: Name of the textfile with the spotlist
            :type spotlist_filename: string

            :returns: The function returns a number of different items in from of a python dictionary.
                      Most data is stored as 2D spatial maps, indicting for each (x,y) location the corresponding data.
                      Most data is stored as 2D masked numpy arrays. The masked of the array indicated whether data
                      has been recorded for a given pixel or not. The dict contains the following keys:

                * 'spotfolder' : String indicating the folder where all the spot-data is located
                * 'fid' : 2D numpy masked array of strings with fid file name for each (x,y) pixel (intensities).
                * 'acqu' : 2D numpy masked array of strings with acqu file name for each (x,y) (metadata).
                * 'regions' : 2D numpy masked array of integers with region index for each (x,y) pixel.
                * 'xpos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'ypos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'spotname' : 2D masked numpy array of strings with spot name for each (x,y) pixel.

        """
        # Read all lines
        f = open(spotlist_filename, 'r')
        lines = f.readlines()
        f.close()
        # Remove empty lines and comments from the file so that only lines with spot data remain
        # Clean up the text to remove extra spaces, endline and tab symols at
        # the end.
        li = 0
        while li < len(lines):
            # Clean up line
            lines[li] = lines[li].rstrip("\n").rstrip("\r").rstrip(" ")
            if len(lines[li]) == 0 or lines[li].startswith("#"):
                lines.pop(li)
            else:
                li += 1
        del li
        num_spots = len(lines)

        # Initalize the different arrays
        xpos = np.zeros(shape=(num_spots,), dtype='int32')
        ypos = np.zeros(shape=(num_spots,), dtype='int32')
        # ToDo: We here assume that we only have 24 characters per R02X071Y048
        # spotname.
        spotname_encoding = 'a24'
        spotname = np.zeros(shape=(num_spots,), dtype=spotname_encoding)
        region = np.zeros(shape=(num_spots,), dtype='int32')
        pixelindex = np.zeros(shape=(num_spots, 2), dtype='int32')
        # Parse all the pixels
        for i in range(0, len(lines)):
            sl = lines[i].split(" ")
            xpos[i] = int(sl[0])
            ypos[i] = int(sl[1])
            spotname[i] = sl[2]
            region[i] = int(sl[3])
            # Extract the pixel locations from the spotname string, e.g.,
            # R00X024Y005
            sp = spotname[i].split("X")[1].split("Y")
            pixelindex[i, 0] = int(sp[0])
            pixelindex[i, 1] = int(sp[1])

        # Create maps for all entries
        xi_min = 0  # np.min(pixelindex[:,0])
        xi_max = np.max(pixelindex[:, 0])
        yi_min = 0  # np.min(pixelindex[:,1])
        yi_max = np.max(pixelindex[:, 1])
        xsize = xi_max - xi_min + 1
        ysize = yi_max - yi_min + 1

        # Create 2D image maps of the data (this is to make the reading process
        # easier)
        region_map = np.empty(shape=(xsize, ysize), dtype=region.dtype)
        region_map.fill(-1)
        region_map[pixelindex[:, 0], pixelindex[:, 1]] = region[:]
        mask = (region_map < 0)
        region_map = np.ma.masked_array(region_map, mask)
        xpos_map = np.empty(shape=(xsize, ysize), dtype=xpos.dtype)
        xpos_map[pixelindex[:, 0], pixelindex[:, 1]] = xpos[:]
        xpos_map = np.ma.masked_array(xpos_map, mask)
        ypos_map = np.empty(shape=(xsize, ysize), dtype=ypos.dtype)
        ypos_map[pixelindex[:, 0], pixelindex[:, 1]] = ypos[:]
        ypos_map = np.ma.masked_array(ypos_map, mask)
        spotname_map = np.empty(shape=(xsize, ysize), dtype=spotname_encoding)
        spotname_map[pixelindex[:, 0], pixelindex[:, 1]] = spotname[:]
        spotname_map = np.ma.masked_array(spotname_map, mask=mask)

        # Compute the folder where the spot files (fid,acqu etc. are located).
        spotlist_abspath = os.path.abspath(spotlist_filename)
        spotlist_dir, spotlist_name = os.path.split(spotlist_abspath)
        spotfolder = spotlist_dir + "/" + \
            spotlist_name.rstrip(" Spot List.txt")
        spotfolder = os.path.abspath(spotfolder)
        if not os.path.isdir(spotfolder):
            raise ValueError("Could not find spotfolder: "+spotfolder)
        # Get a list of all files and compute and image map to indicate where the fid and acqu files
        # for each pixel are located
        # filelist = bruckerflex_file.get_all_files(spotfolder)
        # Compute the longest filename. This is only needed to determine the
        # encoding for the fix-size numpy array.
        longest_name_length = 0
        for path, subdirs, files in os.walk(spotfolder):
            for name in files:
                fn = os.path.join(path, name)
                if len(fn) > longest_name_length:
                    longest_name_length = len(fn)
        fltype = 'a' + str(longest_name_length)
        # ToDo: If needed we need to change this to allow for multiple fid
        # files per spot
        fidfiles = np.empty(shape=(xsize, ysize), dtype=fltype)
        # ToDo: Add additional storage for other files that may be needed
        acqufiles = np.empty(shape=(xsize, ysize), dtype=fltype)
        for path, subdirs, files in os.walk(spotfolder):
            for f in files:
                # ToDo: We are only interested in fid and acqu files right now
                # but that may need to change
                if f.endswith("fid") or f.endswith("acqu") and \
                   not (f.endswith("_fid") or f.endswith("_acqu")):
                    # Compute the pixelindex for the file
                    pindex = path.lstrip(spotfolder).lstrip(
                        "/").split("/")[0].split("X")[1].split("Y")
                    xi = int(pindex[0])
                    yi = int(pindex[1])
                    # ToDo: Again if we are interested in more than fid and
                    # acqu files then we need to extend this check
                    if f.endswith("fid"):
                        fidfiles[xi, yi] = os.path.join(path, f)
                    elif f.endswith("acqu"):
                        acqufiles[xi, yi] = os.path.join(path, f)

        fidfiles = np.ma.masked_array(fidfiles, mask)
        acqufiles = np.ma.masked_array(acqufiles, mask)

        # Return all spatial maps of the different data we collected
        re = {'spotfolder': spotfolder,
              'fid': fidfiles,
              'acqu': acqufiles,
              'regions': region_map,
              'xpos': xpos_map,
              'ypos': ypos_map,
              'spotname': spotname_map}

        return re

        # Return all the different values
        # return xpos, ypos, spotname, region, pixelindex, spotfolder

        # Create a python dictionary for each region
        # Make sure that regions are numbered starting with 0.
        # unique_regions = np.unique(regions)
        # minRegionIndex = np.min( unique_regions )
        # unique_regions = unique_regions - minRegionIndex
        # region = region - np.min(region) - minRegionIndex
        # region_dicts = {}
        # for i in unique_regions :
        #    regionSelect = (region == i)
        #    rd = {}
        #    rd['xpos'] = xpos[regionSelect]
        #    rd['ypos'] = ypos[regionSelect]
        #    rd['spotname'] = spotname[regionSelect]
        #    rd['pixelindex'] = pixelindex[regionSelect,:]

    # @staticmethod
    # def get_all_files(directory) :
    #    """Simple method to get all subdirectories and files located in a folder"""
    #    fl = []
    #    for path, subdirs, files in os.walk(directory):
    #        for name in files:
    #            fl.append(os.path.join(path, name))
    @staticmethod
    def s_read_fid(filename, data_type='int32', selection=slice(None, None, None)):
        """ Read data from an fid file

            :param filename: String indicating the name+path to the fid file.
            :type filename: string
            :param data_type: The numpy datatype encoding used for the fid files. (default is 'int32').
                              In the instance of this class this is encoded in the data_type variable associated
                              with the instance.
            :param selection: This may be a python slice or a list of indecies to be read. Default value is to
                              read all (i.e., slice(None,None,None))
            :type selection: slice or list, i.e., a selection that numpy understands

            :returns: 1D numpy array of intensity values read from the file.
        """
        fidtype = np.dtype(data_type)
        num_values = os.path.getsize(filename) / fidtype.itemsize
        fid_file = np.memmap(filename=filename, dtype=fidtype,
                             shape=(num_values,), mode='r', order='C')
        values = np.array(fid_file[selection])
        del fid_file
        return values

    @staticmethod
    def s_read_acqu(filename):
        """Construct an m/z axis for the given acqu file.

           :param filename: String with the name+path for the acqu file.
           :type filename: string

           :returns: Return dictonary with the parsed metadata information

        """
        # Read the complete acqu file
        acqu = open(filename, 'r')
        lines = acqu.readlines()  # read all lines of the file into a list
        acqu.close()
        #
        # Parse the acqu file and store all data in a python dictonary
        #
        acqu_dict = {}
        curr_line = 0
        while curr_line < len(lines):
            # Skip lines with no data
            if len(lines[curr_line].rstrip("\n").rstrip("\r").rstrip(" ")) == 0:
                curr_line += 1
                continue
            # All variables should start with ## in the acqu file
            if not lines[curr_line].startswith("##"):
                log_helper.warning(__name__,  "WARNING: Error while reading line" + str(curr_line) + \
                      " of the acqu file. The error may have occured on the previous line")
                if curr_line > 0:
                    log_helper.debug(__name__, str(curr_line - 1) + ": " + lines[curr_line - 1])
                log_helper.debug(str(curr_line) + ": " + lines[curr_line])
                curr_line += 1
                continue

            sl = lines[curr_line].split("=")
            key = sl[0]
            # Remove beginning spaces and endline and tabs at the end from the
            # value
            value = sl[1].lstrip(' ').rstrip("\n").rstrip("\r").rstrip(" ")

            # Try to convert the value to a number
            is_number = False
            try:
                value = int(value)
                is_number = True
            except ValueError:
                try:
                    value = float(value)
                    is_number = True
                except ValueError:
                    pass

            # Check whether the entry defines a vector of numbers
            unicode_value = unicode(value)
            if not is_number and unicode_value.startswith("(") and unicode_value.endswith(")"):
                # How many values and lines do we need to read?
                sv = unicode_value.lstrip("(").rstrip(")").split("..")
                # low = int(sv[0])
                high = int(sv[1])
                num_vals = high + 1
                vals_per_line = 8
                num_lines = int(math.ceil(num_vals / float(vals_per_line)))
                # Read all the values into a list and convert the numbers if
                # possible
                value = []
                curr_line += 1
                # print str(num_lines) + " : " + sv[0] + "..." + sv[1]
                for _ in range(0, num_lines):

                    sl = lines[curr_line].rstrip("\n").rstrip("\r").rstrip(" ").split(" ")
                    try:
                        sconv = [int(ix) for ix in sl]
                    except ValueError:
                        try:
                            sconv = [float(ix) for ix in sl]
                        except ValueError:
                            sconv = sl
                    value = value + sconv
                    curr_line += 1

                acqu_dict[key] = value
            else:
                acqu_dict[key] = value
                curr_line += 1

        return acqu_dict

    @staticmethod
    def s_mz_from_acqu(acqu_dict):
        """Construct the m/z axis from the data stored in the acqu_dict dictionary.
           See also s_read_acqu

           :param acqu_dict: Python dictionary with the complete information from the acqu file. See s_read_acqu(...).

            :returns: 1D Numpy array of floats with the mz axis data.
        """
        # Get some of the numbers we need from the acqu dictionary
        time_number = acqu_dict["##$TD"]
        time_delay = acqu_dict["##$DELAY"]
        time_delta = acqu_dict["##$DW"]
        calibration_constants = np.array([acqu_dict["##$ML1"], acqu_dict["##$ML2"], acqu_dict["##$ML3"]])
        calibration_constants[0] = math.pow(
            (1e12 / calibration_constants[0]), 0.5)
        # Compute the tof axis

        tof = np.arange(start=0, stop=time_number) * time_delta + time_delay
        tof = calibration_constants[1] - tof

        # Compute the mz axis
        if calibration_constants[2] == 0:
            mz = np.power(tof, 2.0) / \
                (calibration_constants[0] * calibration_constants[0])
        else:
            mz = np.power(((-1 * calibration_constants[0]) +
                           (np.power(math.pow(calibration_constants[0], 2) -
                                     (4 * calibration_constants[2] * tof), 0.5))) / (2 * calibration_constants[2]), 2.)

        # Return the acqu_dict dictionary as well as the the computed m/z axis
        return mz
