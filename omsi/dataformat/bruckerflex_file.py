"""This module provides functionality for reading bruker flex mass spectrometry image files

      
   Limitations: 
   1) Currently the reader assumes a single global m/z axis for all spectra.
   2) The read of acqu files does not convert <..> entries to python values but leaves them as strings.
   3) __read_spotlist__ converts the regions to start with a 0 index. This is somewhat inconsisten in the
      spot list file. The spotname seems to number regions starting with 0 while the region list numbers
      them starting with 1.
   4) __read_spotlist__ computes the folder where the spots are located based on the filename of the spotlist. 
      The question is whether this is always the case??? The advantage is that we do not rely on the 
    regions.xml file which contains absolute paths which are in most cases invalid as the data has been
    copied between different machines in many cases and is stored in different locations on each of the 
    machines. 
   5) __read_spotlist__ currenly assumes that there is only one fid file per spot
   6) __read_spotlist__ currenlty only looks at where the acqu and fid file is located for the spot. Other 
      files are currently ignored.
   7) __read_spotlist__ (and hence the reader at large) currently assumes that we have 2D images only.
   8) __read_spotlist__ currently generates maps for the image that assume that x and y pixle indicies start at 0. Not all images may record data until the border, so that this strategy may add empty spectra rather than generating a new bounding box for the image. 
   9) __read_spotlist__ assumes in the variable spotname_encoding a maximum of 24 characters in the spotname R01X080Y013. This should in general be more than sufficient as this allows for 7 characters for each R, X, Y entry, however, if this is not enough then this behaviour needs ot be changed.
   10) __getitem__ currently only works if we have read the full data into memory. An on-demand load should be supported as well.
   11) We can currently only selected either a single region or the full data but we cannot selected multiple regions at once. E.g. if a dataset contains 3 regions then we can either select all regions at once or region 1,2, or 3 but one cannot selected region 1+2, 1+3, or 2+3.
   
   
    
import bruckerflex_file
spotlist = "/Users/oruebel/Devel/msidata/Bruker_Data/UNC IMS Data/20130417 Bmycoides Paenibacillus Early SN03130/2013 Bmyc Paeni Early LP/2013 Bmyc Paeni Early LP Spot List.txt"    
exppath ="/Users/oruebel/Devel/msidata/Bruker_Data/UNC IMS Data/20130417 Bmycoides Paenibacillus Early SN03130/2013 Bmyc Paeni Early LP/2013 Bmyc Paeni Early LP/0_R00X012Y006/1/1SLin"
f = bruckerflex_file.bruckerflex_file( spotlist_filename = spotlist)
f.s_read_fid( exppath+"/fid" , f.data_type )
testacqu = f.s_read_acqu( exppath+"/acqu" )
testmz   = f.s_mz_from_acqu( testacqu )
testspotlist = f.s_read_spotlist(spotlist) 
   
"""

import sys 
import os
import numpy as np
import math

class bruckerflex_file :
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

    def __init__( self , spotlist_filename, fid_encoding='int32', readall = True) :
                
        """Open an img file for data reading.

            :param spotlist_filename: Name of the textfile with the spotlist
            :type spotlist_filename: string
            :param readall: Should the complete data be read into memory (this makes slicing easier). (default is True)
            :type readall: bool
            :param data_type: String indicating in which binary format the intensity values (fid files) are stored. (default value is 'int32')
            :type data_type: string
            
            :var spotlist_filename: Name of the file with the spotlist
            :var pixel_dict: dictionary with the pixel array metadata (see also s_read_spotlist(...)). Some of the main keys of the dictionary are, e.g. (see also s_read_spotlist(...)) : 
                
                * 'spotfolder' : String indicating the folder where all the spot-data is located
                * 'fid' : 2D numpy masked array of strings indicating for (x,y) pixel the fid file with the intensity values.
                * 'acqu' : 2D numpy masked array of strings indicating for each (x,y) pixel the acqu file with the metadata for that pixel.
                * 'regions' : 2D numpy masked array of integers indicating for each pixels the index of the region it belongs to.
                * 'xpos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'ypos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'spotname' : 2D masked numpy array of strings indicating the name of the spot corresponding to a pixel.
                
            :var data_type: the encoding used for intensity values
            :var shape: The 3D shape of the MSI data volume for the currently selected region.
            :var full_shape: Shape of the full 3D MSI dataset including all regions imaged. 
            :var metadata: Dictionary with metadata from the acqu file
            :var mz: The 1D numpy array with the m/z axis information
            :var data: If readall is set to true then this 3D array includes the complete data of the MSI data cube. Missing data values (e.g., from regions not imaged during the aquistion processes) are completed with zero values. 
            :var regions_dicts: 

            
        """
        self.spotlist_filename = spotlist_filename
        self.pixel_dict = self.s_read_spotlist( self.spotlist_filename )
        self.data_type = fid_encoding #Data type of the fid files
        
        #ToDo this reads a single spectrum (ie., fid file) to figure out the length of the spectra. This obviously means that we assume a single m/z axis for all spectra
        tempSpec = self.s_read_fid( self.pixel_dict['fid'].compressed()[0] , self.data_type )
        self.full_shape = [ self.pixel_dict['regions'].shape[0]  , self.pixel_dict['regions'].shape[1], tempSpec.shape[0] ]
        self.shape      = [ self.pixel_dict['regions'].shape[0]  , self.pixel_dict['regions'].shape[1], tempSpec.shape[0] ]#Number of pixels in x,y, and z
        
        #ToDo this reads a single acqu metadata file to figure out the m/z axis. This obviously means that we assume a single m/z axis for all spectra
        self.metadata = self.s_read_acqu( self.pixel_dict['acqu'].compressed()[0] )
        self.mz = self.s_mz_from_acqu( self.metadata )
        
        #Compute the region bouding boxes
        self.regions_dicts = self.__compute_regions_extends___(self.pixel_dict)
        self.select_region = None
        
        #Should we read all the intensity values into memory. If so, allocate space.
        if readall :
            self.data = np.zeros( self.shape )
            dmask = self.pixel_dict['fid'].mask
            fidFiles = self.pixel_dict['fid']
            for xi in range( 0 , self.shape[0]) :
                for yi in range( 0 , self.shape[1] ) :
                    if not dmask[xi,yi] :
                        self.data[xi,yi,:] = self.s_read_fid( self.pixel_dict['fid'][xi,yi] ,  self.data_type )
        else :
            self.data = None
        
        self.instrument_mz=0 #A numpy vector with the m/z values of the instrument        
  
    def get_regions(self) :
        """ Get list of all region dictionaries defining for each region the origin 
            and extend of the region. See also self.regions_dicts.
        """
        return self.regions_dicts
            
    def get_number_of_regions(self) :
        """Get the number of available regions"""
        return len(self.regions_dicts)

    def set_region_selection(self, region_index=None) :
        """Define which region should be selected for local data reads. 
        
           :param region_index: The index of the region that should be read. The shape of the
                    data will be adjusted accordingly. Set to None to select all regions and treat
                    the data as a single full 3D image.
        """
        if region_index is None :
            self.shape = self.full_shape
            self.select_region = None
        elif region_index<self.get_number_of_regions() :
            self.select_region = region_index
            self.shape = [ self.regions_dicts[self.select_region ]["extend"][0] , self.regions_dicts[self.select_region ]["extend"][1] , self.full_shape[2] ]

    def get_region_selection(self) :
        """Get the index of the selected region"""
        return self.select_region

            
    def __getitem__(self, key) :
        """Enable slicing of bruker files"""
        #If the full data has been loaded 
        if self.data is not None :
            #Select the region of interest, i.e., the reader acts as if the selected region where the complete image of interest. 
            if self.select_region is not None :
                rd_ori = self.regions_dicts[ self.select_region ]["origin"]
                rd_ext = self.regions_dicts[ self.select_region ]["extend"]
                lx = rd_ori[0]
                ly = rd_ori[1]
                hx = lx+rd_ext[0]+1
                hy = ly+rd_ext[1]+1
                #Select the region and then selection the user data within that region
                return self.data[ lx:hx, ly:hy  , : ][key]
            else:
                #Return the data from the full image data.
                return self.data[key]
        else :
            raise ValueError("Slicing is only supported when the object has been initialized with readall")
        pass

    def close_file(self) :
        """Close the img file"""
        pass

    #def __del__(self) :
    #    """Close the file before garbage collection"""
    #    self.close_file()
    

    def __read_spectrum__( self , filename, selection) :
        """Read data from a single fid file with the intensities.
        
           :param filename: String indicating the name+path to the fid file. 
           :type filename: string
           :param selection: This may be a python slice or a list of indecies to be read. 
           :type selection: slice or list, i.e., a selection that numpy understands
            
           :returns:
            -- 1D numpy array of intensity values read from the file.
            -- 1D numpy array of the mz values for the part of the spectrum read. 
        """
        #Read the requested intensity values from the fid file
        intensity = s_read_fid(filename , self.data_type, selection)
    
        #Construct the m/z axis for an acqu file 
        dirname = os.path.dirname(os.path.abspath(filename)) #Get the directory where the file is located
        acqu_filename = dirname+"/acqu"
        acqu_dict  = s_read_acqu(acqu_filename)
        mz = s_mz_from_acqu(acqu_dict)
        return intensity, mz[selection]
    
    
    def __compute_regions_extends___(self, pixel_maps) :
        """Compute the region extends from the given region map.
    
            :param pixel_maps: This is the self.pixel_dict which is computed in the __init__ function using the s_read_spotlist function. 
       
            :returns:
        """
        region_map = pixel_maps['regions']
        minr = np.min( region_map )
        maxr = np.max( region_map )
        regions = []
        for r in range( minr, maxr+1) :
            rr = np.argwhere( region_map==r )
            minx = np.min( rr[:,0] )
            miny = np.min( rr[:,1] )
            maxx = np.max( rr[:,0] )
            maxy = np.max( rr[:,1] )
            extendX = maxx - minx
            extendY = maxy - miny
            boundingBox = {}
            boundingBox["origin"] = np.array( [minx,miny] )
            boundingBox["extend"] = np.array( [extendX , extendY] )
            regions.append( boundingBox )
        
        return regions
    

    @staticmethod
    def s_read_spotlist(spotlist_filename) :
        """Parse the given spotlist file.
            
            :param spotlist_filename: Name of the textfile with the spotlist
            :type spotlist_filename: string
            
            :returns: The function returns a number of different items in from of a python dictionary. Most data is stored as 2D spatial maps, indicting for each (x,y) location the corresponding data. Most data is stored as 2D masked numpy arrays. The masked of the array indicated whether data has been recorded for a given pixel or not. The dict contains the following keys: 
            
                * 'spotfolder' : String indicating the folder where all the spot-data is located
                * 'fid' : 2D numpy masked array of strings indicating for (x,y) pixel the fid file with the intensity values.
                * 'acqu' : 2D numpy masked array of strings indicating for each (x,y) pixel the acqu file with the metadata for that pixel.
                * 'regions' : 2D numpy masked array of integers indicating for each pixels the index of the region it belongs to.
                * 'xpos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'ypos' : 2D numpy masked array of integers indicated for each pixel its x position.
                * 'spotname' : 2D masked numpy array of strings indicating the name of the spot corresponding to a pixel.
                
        """
        #Read all lines
        f = open( spotlist_filename, 'r' )
        lines = f.readlines()
        f.close()
        #Remove empty lines and comments from the file so that only lines with spot data remain
        #Clean up the text to remove extra spaces, endline and tab symols at the end.
        li = 0
        while li < len(lines) :
            lines[li] = lines[li].rstrip("\n").rstrip("\r").rstrip(" ") #Clean up line
            if len(lines[li]) == 0 or lines[li].startswith("#") :
                lines.pop(li)
            else :
                li = li+1
        del(li)
        numSpots = len(lines)
        
        #Initalize the different arrays
        xpos = np.zeros( shape=(numSpots,) , dtype='int32' )
        ypos = np.zeros( shape=(numSpots,) , dtype='int32' )
        spotname_encoding = 'a24' #ToDo: We here assume that we only have 24 characters per R02X071Y048 spotname.
        spotname =  np.zeros( shape=(numSpots,) , dtype=spotname_encoding)
        region = np.zeros( shape=(numSpots,) , dtype='int32' )
        pixelindex = np.zeros( shape=(numSpots,2) , dtype='int32' )
        #Parse all the pixels
        for i in range( 0, len(lines) ) :
            sl = lines[i].split(" ")
            xpos[i] = int(sl[0])
            ypos[i] = int(sl[1])
            spotname[i] = sl[2]
            region[i] = int(sl[3])
            #Extract the pixel locations from the spotname string, e.g., R00X024Y005
            sp = spotname[i].split("X")[1].split("Y")
            pixelindex[i,0] = int(sp[0])
            pixelindex[i,1] = int(sp[1])
    
        #Create maps for all entries
        xi_min = 0 #np.min(pixelindex[:,0])
        xi_max = np.max(pixelindex[:,0])
        yi_min = 0 #np.min(pixelindex[:,1])
        yi_max = np.max(pixelindex[:,1])
        xsize = xi_max - xi_min +1
        ysize = yi_max - yi_min +1
        
        #Create 2D image maps of the data (this is to make the reading process easier)
        region_map = np.empty( shape=(xsize,ysize) , dtype=region.dtype)
        region_map.fill(-1)
        region_map[ pixelindex[:,0] , pixelindex[:,1] ] = region[:]
        mask = (region_map<0)
        region_map = np.ma.masked_array( region_map , mask ) 
        xpos_map =   np.empty( shape=(xsize,ysize) , dtype=xpos.dtype)
        xpos_map[ pixelindex[:,0] , pixelindex[:,1] ] = xpos[:]
        xpos_map = np.ma.masked_array( xpos_map , mask )
        ypos_map = np.empty( shape=(xsize,ysize) , dtype=ypos.dtype)
        ypos_map[ pixelindex[:,0] , pixelindex[:,1] ] = ypos[:]
        ypos_map = np.ma.masked_array( ypos_map , mask )
        spotname_map = np.empty( shape=(xsize,ysize) , dtype=spotname_encoding)
        spotname_map [ pixelindex[:,0] , pixelindex[:,1] ] =spotname[:]
        spotname_map = np.ma.masked_array( spotname_map , mask=mask )
        
        #Compute the folder where the spot files (fid,acqu etc. are located).
        spotlist_abspath = os.path.abspath( spotlist_filename)
        spotlist_dir, spotlist_name = os.path.split( spotlist_abspath )
        spotfolder = spotlist_dir+"/"+spotlist_name.rstrip(" Spot List.txt")
        spotfolder = os.path.abspath( spotfolder )
        #Get a list of all files and compute and image map to indicate where the fid and acqu files
        #for each pixel are located 
        #filelist = bruckerflex_file.get_all_files(spotfolder)
        #Compute the longest filename. This is only needed to determine the encoding for the fix-size numpy array. 
        longestNameLength = 0 
        for path, subdirs, files in os.walk(spotfolder):
           for name in files:
                fn = os.path.join(path, name)
                if len(fn) > longestNameLength :
                    longestNameLength = len(fn)
        fltype = 'a'+str(longestNameLength)
        fidfiles = np.empty( shape=(xsize,ysize) , dtype=fltype) #ToDo: If needed we need to change this to allow for multiple fid files per spot
        acqufiles = np.empty( shape=(xsize,ysize) , dtype=fltype) #ToDo: Add additional storage for other files that may be needed
        for path, subdirs, files in os.walk(spotfolder) :
            for f in files :
                #ToDo: We are only interested in fid and acqu files right now but that may need to change
                if f.endswith("fid") or f.endswith("acqu") :
                    #Compute the pixelindex for the file
                    pindex = path.lstrip(spotfolder).lstrip("/").split("/")[0].split("X")[1].split("Y")
                    xi = int(pindex[0])
                    yi = int(pindex[1])
                    #ToDo: Again if we are interested in more than fid and acqu files then we need to extend this check
                    if f.endswith("fid") :
                        fidfiles[xi,yi] = os.path.join( path , f)
                    elif f.endswith("acqu") :
                        acqufiles[xi,yi] = os.path.join( path , f)
    
        fidfiles  = np.ma.masked_array( fidfiles , mask )
        acqufiles = np.ma.masked_array( acqufiles , mask )
            
        #Return all spatial maps of the different data we collected
        re = {}
        re['spotfolder'] = spotfolder
        re['fid'] = fidfiles
        re['acqu'] = acqufiles
        re['regions'] = region_map
        re['xpos'] = xpos_map
        re['ypos'] = ypos_map
        re['spotname'] = spotname_map
        
        return re
                
        #Return all the different values 
        #return xpos, ypos, spotname, region, pixelindex, spotfolder

        #Create a python dictionary for each region
        #Make sure that regions are numbered starting with 0. 
        #unique_regions = np.unique(regions)
        #minRegionIndex = np.min( unique_regions )
        #unique_regions = unique_regions - minRegionIndex
        #region = region - np.min(region) - minRegionIndex
        #regions_dicts = {}
        #for i in unique_regions :
        #    regionSelect = (region == i)
        #    rd = {}
        #    rd['xpos'] = xpos[regionSelect]
        #    rd['ypos'] = ypos[regionSelect]
        #    rd['spotname'] = spotname[regionSelect]
        #    rd['pixelindex'] = pixelindex[regionSelect,:]


    #@staticmethod
    #def get_all_files(directory) :
    #    """Simple method to get all subdirectories and files located in a folder"""
    #    fl = []
    #    for path, subdirs, files in os.walk(directory):
    #        for name in files:
    #            fl.append(os.path.join(path, name))
    
    @staticmethod
    def s_read_fid(filename, data_type='int32' , selection=slice(None,None,None) ) :
        """ Read data from an fid file
        
            :param filename: String indicating the name+path to the fid file.
            :type filename: string
            :param data_type: The numpy datatype encoding used for the fid files. (default is 'int32'). In the instance of this class this is encoded in the data_type variable associated with the instance.
            :param selection: This may be a python slice or a list of indecies to be read. Default value is to read all (i.e., slice(None,None,None))
            :type selection: slice or list, i.e., a selection that numpy understands
        
            :returns: 1D numpy array of intensity values read from the file.
        """
        fidtype = np.dtype( data_type )
        numValues = os.path.getsize( filename ) / fidtype.itemsize
        fid_file = np.memmap( filename=filename , dtype=fidtype , shape=( numValues ,) , mode='r', order='C')
        values = np.array(fid_file[selection])
        del(fid_file)
        return values
    
    @staticmethod
    def s_read_acqu(filename ) :
        """Construct an m/z axis for the given acqu file.
        
           :param filename: String with the name+path for the acqu file.
           :type filename: string
            
           :returns: Return dictonary with the parsed metadata information
        
        """
        #Read the complete acqu file
        acqu = open( filename, 'r')
        lines = acqu.readlines() #read all lines of the file into a list
        acqu.close()
        #
        #Parse the acqu file and store all data in a python dictonary
        #
        acqu_dict = {}
        currLine = 0
        while currLine < len(lines) :
            if len(lines[currLine].rstrip("\n").rstrip("\r").rstrip(" ")) ==0 : #Skip lines with no data
                currLine = currLine+1
                continue
            if not lines[currLine].startswith("##") : #All variables should start with ## in the acqu file
                print "WARNING: Error while reading line"+str(currLine)+" of the acqu file. The error may have occured on the previous line"
                if currLine>0 :
                    print str(currLine-1) + ": " + lines[currLine-1]
                print str(currLine) + ": " + lines[currLine]
                currLine = currLine+1
                continue
            
            sl = lines[currLine].split("=")
            key = sl[0]
            value = sl[1].lstrip(' ').rstrip("\n").rstrip("\r").rstrip(" ") #Remove beginning spaces and endline and tabs at the end from the value

            #Try to convert the value to a number
            isNumber = False
            try:
                value = int(value)
                isNumber = True
            except ValueError:
                try :
                    value = float(value)
                    isNumber = True
                except:
                    isNumber = False

            #Check whether the entry defines a vector of numbers 
            if not isNumber and value.startswith("(") and value.endswith(")") :
                #How many values and lines do we need to read?
                sv = value.lstrip("(").rstrip(")").split("..")
                low = int(sv[0])
                high = int(sv[1])
                numVals = high+1
                valsPerLine = 8
                numLines = int( math.ceil(numVals / float(valsPerLine)) )
                #Read all the values into a list and convert the numbers if possible 
                value = []
                currLine = currLine+1
                #print str(numLines) + " : " + sv[0] + "..." + sv[1]
                for i in range( 0 , numLines ) :
                    
                    sl = lines[currLine].rstrip("\n").rstrip("\r").rstrip(" ").split(" ")
                    try :
                        sconv = [ int(ix) for ix in sl ]
                    except :
                        try :
                            sconv = [ float(ix) for ix in sl ]
                        except :
                            sconv = sl
                    value = value + sconv 
                    currLine = currLine+1
                
                acqu_dict[key] = value
            else :
                acqu_dict[key] = value
                currLine = currLine+1
    
        return acqu_dict
    
    @staticmethod
    def s_mz_from_acqu(acqu_dict):
        """Construct the m/z axis from the data stored in the acqu_dict dictionary.
           See also s_read_acqu
            
           :param acqu_dict: Python dictonary with the complete information from the acqu file. See also s_read_acqu(...)
            
            :returns: 1D Numpy array of floats with the mz axis data. 
        """
        #Get some of the numbers we need from the acqu dictionary
        timeNumber = acqu_dict["##$TD"]
        timeDelay  = acqu_dict["##$DELAY"]
        timeDelta  = acqu_dict["##$DW"]
        calibrationConstants = np.array( [ acqu_dict["##$ML1"] , acqu_dict["##$ML2"] , acqu_dict["##$ML3"] ] )
        calibrationConstants[0] = math.pow( (1e12 / calibrationConstants[0]) , 0.5 )
        #Compute the tof axis
        
        
        tof = np.arange(start=0, stop=timeNumber)*timeDelta + timeDelay
        tof = calibrationConstants[1] - tof
        
        
        #Compute the mz axis
        if calibrationConstants[2] == 0 :
            mz = np.power( tof , 2.0 ) / ( calibrationConstants[0] * calibrationConstants[0])
        else :
            mz = np.power( ( (-1*calibrationConstants[0]) + ( np.power(  math.pow(calibrationConstants[0],2) - ( 4*calibrationConstants[2]*tof) , 0.5)   )  ) / ( 2* calibrationConstants[2] ) , 2.)


        #Return the acqu_dict dictionary as well as the the computed m/z axis
        return mz 
