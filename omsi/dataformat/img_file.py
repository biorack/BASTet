"""This module provides functionality for reading img mass spectrometry image files"""

import sys 
import os
import numpy as np

class img_file : 
    """Interface for reading a single 2D img file

       The img format consists of three different files:
       i) hdr header file, ii) t2m which contains the m/z data, 
       iii) img data file.
    """

    def __init__( self , hdrFile=None , t2mFile=None, imgFile=None, basename=None ) :
                
        """Open an img file for data reading.

            :param hdrFile: The name of the hdr header file
            :type hdrFile: string

            :param t2mFile: The name of the t2mFile
            :type t2mFile: string

            :param imgFile: The name of the img data file
            :type imgFile: string
            
            :param basename: Instead of imgFile, t2mFile, and hdrFile one may also supply just a basname to be completed. 
            :type basename: string
            
            :raises ValueError: In case that basename and hdrFile, t2mFile, and imgFile are specified. 
        """
     
        self.data_type = 'uint16'
        self.shape     = [0,0,0]#Number of pixels in x,y, and z
        self.mz=0 #A numpy vector with the m/z values of the instrument        
  
        if basename and hdrFile and t2mFile and imgFile :
            raise ValueError( "Conflicting input. Provider either basename or explicit hdrFile,t2mFile,imgFile parameters but not both.")
        if basename :
            basefile = basename
            if os.path.isdir( basename ) :
                filelist = img_file.get_files_from_dir( basename )
                if len(filelist) > 0 :
                    basefile = os.path.join( basename , filelist[0] )
                else :
                    raise ValueError("No valid img file found in the given directory.")
            if os.path.exists( basefile+".hdr") and os.path.exists( basefile+".t2m") and os.path.exists( basefile+".img") :
                hdrFile=basefile+".hdr"
                t2mFile=basefile+".t2m"
                imgFile=basefile+".img"
            else :
                raise ValueError("No valid img file found for the given basename.")
        elif  hdrFile and t2mFile and imgFile :
            pass #Nothing to be done
        else :
            raise ValueError("Missing input parameter. Either provide: i) basename or ii) hdrFile, t2mFile, imgFile")
            
  
        #Initalize the z length
        try:
            t2m = open( t2mFile , 'rb' )
            self.mz  = np.fromfile( file=t2m , dtype = 'float32' , count=-1 )
            self.shape[2] = self.mz.shape[0]
            t2m.close()
        except:
           self.mz = np.empty(0)

        #Initalize th x and y length
        try:
           hdr = open( hdrFile , 'rb' )
           hdrData = np.fromfile( file=hdrFile , dtype='int16' , count=-1 )
           self.shape[0] = hdrData[23]
           self.shape[1] = hdrData[22]
           hdr.close() 
        except:
           pass
        
        self.shape = tuple(self.shape)
            
        #Open the img file with the spectrum data 
        self.img_filename = imgFile
        self.fileOpened = False
        try:
           self.m_img_file = np.memmap( filename= self.img_filename , dtype=self.data_type , shape=self.shape , mode='r', order='C') #open( imgFile , 'rb' )
           self.fileOpened = True
        except : 
           print "Error while opening the img file: "+imgFile
           raise
        
    
    def __getitem__(self, key) :
        """Enable slicing of img files"""
        if self.m_img_file is None :
            self.m_img_file = np.memmap( filename= self.img_filename , dtype=self.data_type , shape=self.shape , mode='r', order='C')
        return self.m_img_file[key]
        

    def close_file(self) :
        """Close the img file"""
        if self.fileOpened :
            del(self.m_img_file)
            self.m_img_file = None
            self.fileOpened = False

    @classmethod
    def is_img(cls , name) :
        """Check whether the given file or directory points to a img file.
        
           :param name: Name of the dir or file.
           :type name: String
           
           :returns: Boolean indicating whether the given file or folder is a valid img file.
        """
        #Check if this points to a basname for the img file
        if os.path.exists( name+".hdr") and os.path.exists( name+".t2m") and os.path.exists( name+".img") :
            return True
        #If we point to a director, check if the dir contains an img file
        elif os.path.isdir( name ) :
            filelist = cls.get_files_from_dir( name )
            print filelist
            if len(filelist)>0 :
                return True
                
        return False
   
    @classmethod
    def get_files_from_dir(cls, dirname) :
        """Get a list of all basenames of all img files in a given directory"""
        filelist = []
        for l in os.listdir(dirname)  :
            currName = os.path.join( dirname , l )
            if os.path.isfile( currName ) and \
                currName.endswith(".img") :
                basename = currName.rstrip(".img")
                if os.path.exists( basename+".hdr") and os.path.exists( basename+".t2m") and os.path.exists( basename+".img") :
                    filelist.append( basename )
        return filelist 

    def __del__(self) :
        """Close the file before garbage collection"""
        self.close_file()


	


