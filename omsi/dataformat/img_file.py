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

    def __init__( self , hdrFile , t2mFile, imgFile ) :
                
        """Open an img file for data reading.

            :param hdrfile: The name of the hdr header file
            :type hdrfile: string

            :param t2mfile: The name of the t2mFile
            :type t2mfile: string

            :param imgFile: The name of the img data file
            :type imgFile: string
        """
     
        self.data_type = 'uint16'
        self.shape     = [0,0,0]#Number of pixels in x,y, and z
        self.mz=0 #A numpy vector with the m/z values of the instrument        
  
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
   

    def __del__(self) :
        """Close the file before garbage collection"""
        self.close_file()


	


