"""Simple script to generate thumbnail images"""

import h5py
from omsi.dataformat.omsi_file import *
from PIL import Image
import numpy as np

def main(argv=None):
    '''Then main function'''

    import sys
    from sys import argv,exit
    if argv is None:
        argv = sys.argv
   
    #Check for correct usage
    if len(argv) <5 :
        print  "Usage: python makeThumb.py <HDF5-File> <NMF-Index1> <NMF-Index2> <NMF-Index3>"
        exit(0)
    
  
    #Settings
    nmfIndex = 1
    expIndex = 0
    datafile = argv[1] #"/project/projectdirs/openmsi/omsi_data/11042008_NIMS.h5" 
    n1 = int(argv[2])
    n2 = int(argv[3])
    n3 = int(argv[4])
    applyLogScale = True
    thumbnailFilename = datafile+"_0.png"

    #Open the file and required analysis dataset
    omsiFile =  omsi_file( datafile , 'r' )
    exp = omsiFile.get_exp(expIndex)
    ana = exp.get_analysis(nmfIndex)
    ho = ana['ho']

    #Load the required data
    Nx = ho.shape[0]
    Ny = ho.shape[1]
    d1 = ho[:,:,n1].reshape( (Nx,Ny) )
    d2 = ho[:,:,n2].reshape( (Nx,Ny) )
    d3 = ho[:,:,n3].reshape( (Nx,Ny) )

    #Scale the data
    if applyLogScale :
       d1 = np.log( d1+1 )
       d2 = np.log( d2+1 )
       d3 = np.log( d3+1 )
    
    #Normalize the data values
    d1 = d1 / float( np.max( d1 ) )
    d2 = d2 / float( np.max( d2 ) )
    d3 = d3 / float( np.max( d3 ) )

    #Generate the grayscale images
    im1 = Image.fromarray( d1.astype('float')*255 ).convert('L')
    im2 = Image.fromarray( d2.astype('float')*255 ).convert('L')
    im3 = Image.fromarray( d3.astype('float')*255 ).convert('L')
    
    #Generate the RGB image and save the file
    thumbnail = Image.merge( 'RGB' , (im1,im2,im3) )
    print "Save image"
    thumbnail.save( thumbnailFilename , 'PNG' ) 
    print thumbnailFilename


if __name__ == "__main__":
    main()

