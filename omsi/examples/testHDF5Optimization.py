"""Simple test script used to test the performance of different HDF5 
   optimizations (using chunking) to improve the performance of 
   hyperslap selections"""

from omsi.dataformat.img_file import img_file
from omsi.dataformat.omsi_file import *
import numpy as np
import os
import time
import random
import sys
import math

def main(argv=None):
    """Then main function"""

    import sys
    from sys import argv,exit
    
    if argv is None:
        argv = sys.argv
        
    #Check for correct usage
    if len(argv) !=4 :
        printHelp()
        exit(0)
        
    zval = int(argv[1] )
    omsiOutFile  = argv[2]
    resultsFile  = argv[3]
    
    xdim = 100
    ydim = 100
    zdim = 100000
    
    xtest = [1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10   ]  #range(10,50,10)
    ytest = [1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10   ]  #range(10,50,10)
    ztest = [zval, zval, zval, zval, zval, zval, zval, zval, zval ,zval ]  #range(100, 10000, 100)
    
    allRes = np.zeros( len(xtest) , dtype=[ ('x','int32') , ('y','int32') , ('z','int32') , ('write','f') , ('z_min','f') , ('z_avg','f') , ('z_median','f') , ('z_max','f') , ('xy_min','f') , ('xy_avg','f') , ('xy_median','f') , ('xy_max','f') , ('xyz_min','f') , ('xyz_avg','f') , ('xyz_median','f') , ('xyz_max','f') , ('filesize' , 'f') ]  )
    
    numTests = 3
    repeats = 50

    
    for tr in xrange( 0 , len(xtest) ) :
        
        allRes[tr]['x'] = xtest[tr]
        allRes[tr]['y'] = ytest[tr]
        allRes[tr]['z'] = ztest[tr]

    #XYZ selection 25,000 elements
    sliceWidthX=5
    sliceWidthY=5
    sliceWidthZ = 1000

    for tr in xrange( 0 , len(xtest) ) :
        
        allRes[tr]['x'] = xtest[tr]
        allRes[tr]['y'] = ytest[tr]
        allRes[tr]['z'] = ztest[tr]
        
        xchunk = xtest[tr]
        ychunk = ytest[tr]
        zchunk = ztest[tr]
        print "Chunking:"+str(xchunk)+":"+str(ychunk)+":"+str(zchunk)
        start = time.time()
        generateTestFile( omsiOutFile , xdim , ydim, zdim, xchunk, ychunk , zchunk )
        allRes[tr]['write'] =  allRes[tr]['write']+(time.time() - start)
        allRes[tr]['filesize'] = os.stat( omsiOutFile ).st_size
        omsiFile = omsi_file( omsiOutFile )
        data = omsiFile.get_experiment(0).get_msidata(0)
        print "Time for data write:"+str(allRes[tr]['write'] )
        
        #Select xyz slize
        resXYZ = np.zeros( repeats , dtype='f' )
        for te in xrange( 0 , repeats ) :
            
            valX = random.randint(0, xdim-sliceWidthX-1 )
            valY = random.randint(0, ydim-sliceWidthY-1 )
            valZ = random.randint(0, zdim-sliceWidthZ-1 )
            start = time.time()
            d=np.sum( data[valX:(valX+sliceWidthX),valY:(valY+sliceWidthY),valZ:(valZ+sliceWidthZ)] )
            resXYZ[te] = (time.time() - start)
        
        np.savetxt( "resXYZ_"+str(xchunk)+"_"+str(ychunk)+"_"+str(zchunk)+".txt" , resXYZ )
        allRes[tr]['xyz_min'] = np.min( resXYZ )
        allRes[tr]['xyz_avg'] = np.average( resXYZ )
        allRes[tr]['xyz_median'] = np.median( resXYZ )
        allRes[tr]['xyz_max'] = np.max( resXYZ )
        print "XYZ-Slicing: "+str(np.max(resXYZ)) +":"+str(np.average( resXYZ ))+":"+str(np.median( resXYZ ))+":"+str(np.min(resXYZ))
        
        omsiFile.close_file()
        os.remove( omsiOutFile  )
        
        
    #mz-slice selection 250,000 elements
    sliceWidthZ = 25 #xdim=100 , ydim=100
    for tr in xrange( 0 , len(xtest) ) :
        
        allRes[tr]['x'] = xtest[tr]
        allRes[tr]['y'] = ytest[tr]
        allRes[tr]['z'] = ztest[tr]
        
        xchunk = xtest[tr]
        ychunk = ytest[tr]
        zchunk = ztest[tr]
        print "Chunking:"+str(xchunk)+":"+str(ychunk)+":"+str(zchunk)
        start = time.time()
        generateTestFile( omsiOutFile , xdim , ydim, zdim, xchunk, ychunk , zchunk )
        allRes[tr]['write'] =  allRes[tr]['write']+(time.time() - start)
        omsiFile = omsi_file( omsiOutFile )
        data = omsiFile.get_experiment(0).get_msidata(0)
        print "Time for data write:"+str(allRes[tr]['write'] )

        #Perform different types of slicing opterations and keep track of the times
        #Select 10 slizes in z
        resZ = np.zeros( repeats , dtype='f' )
        for te in xrange( 0 , repeats ) :
            
            val = random.randint(0, zdim-sliceWidthZ-1 )
            start = time.time()
            d=np.sum( data[:,:,val:(val+sliceWidthZ)] )
            resZ[te] = (time.time() - start)
        
        np.savetxt( "resZ_"+str(xchunk)+"_"+str(ychunk)+"_"+str(zchunk)+".txt" , resZ )
        allRes[tr]['z_min'] = np.min( resZ )
        allRes[tr]['z_avg'] = np.average( resZ )
        allRes[tr]['z_median'] = np.median( resZ )
        allRes[tr]['z_max'] = np.max( resZ )
        print "Z-Slicing: "+str(np.max(resZ)) +":"+str(np.average( resZ ))+":"+str(np.median( resZ ))+":"+str(np.min(resZ))
        
        omsiFile.close_file()
        os.remove( omsiOutFile  )
    
    #XY selection of multiple spectra. Selection size= 2,500,000 elements
    for tr in xrange( 0 , len(xtest) ) :
        
        xchunk = xtest[tr]
        ychunk = ytest[tr]
        zchunk = ztest[tr]
        print "Chunking:"+str(xchunk)+":"+str(ychunk)+":"+str(zchunk)
        start = time.time()
        generateTestFile( omsiOutFile , xdim , ydim, zdim, xchunk, ychunk , zchunk )
        allRes[tr]['write'] =  allRes[tr]['write']+(time.time() - start)
        omsiFile = omsi_file( omsiOutFile )
        data = omsiFile.get_experiment(0).get_msidata(0)
        print "Time for data write:"+str(allRes[tr]['write'] )
        
        #Select x/y
        resXY = np.zeros( repeats , dtype='f' )
        for te in xrange( 0 , repeats ) :
            
            valX = random.randint(0, xdim-sliceWidthX-1 )
            valY = random.randint(0, ydim-sliceWidthY-1 )
            start = time.time()
            d=np.sum( data[valX:(valX+sliceWidthX),valY:(valY+sliceWidthY),:] )
            resXY[te] = (time.time() - start)
        
        np.savetxt( "resXY_"+str(xchunk)+"_"+str(ychunk)+"_"+str(zchunk)+".txt" , resXY )
        allRes[tr]['xy_min'] = np.min( resXY )
        allRes[tr]['xy_avg'] = np.average( resXY )
        allRes[tr]['xy_median'] = np.median( resXY )
        allRes[tr]['xy_max'] = np.max( resXY )
        print "XY-Slicing: "+str(np.max(resXY)) +":"+str(np.average( resXY ))+":"+str(np.median( resXY ))+":"+str(np.min(resXY))
        
        omsiFile.close_file()
        os.remove( omsiOutFile  )
    
    
        
    for tr in xrange( 0 , len(xtest) ) :
         allRes[tr]['write'] =  allRes[tr]['write']/float(numTests)
         
    f = open( resultsFile , 'w' )
    f.write( str(allRes.dtype.names) )
    f.write("\n")
    np.savetxt( f , allRes )
    f.close()
    
    

def generateTestFile( omsiOutFile , xdim , ydim, zdim, xchunk, ychunk , zchunk ) :
    
    #Create the output HDF5 file
    try:
        omsiFile = omsi_file( omsiOutFile )
    except:
        print "Unexpected error creating the output file:", sys.exc_info()[0]
        exit(0)
        
    exp = omsiFile.create_experiment( exp_identifier = "test" )
    #Create an empty method descrition
    sample = exp.create_method_info()
    #Create an empty instrument description
    mzdata = np.ones( zdim )
    instrument = exp.create_instrument_info(instrumentname="undefined" , mzdata=mzdata )
    #Allocate space in the HDF5 file for the img data
    data = exp.create_msidata(data_shape=( xdim , ydim , zdim  ) , data_type = 'uint16' , chunks=(xchunk,ychunk,zchunk))
    itertest=0 
    numChunksX = int( math.ceil( float(xdim)/float(xchunk) ) )
    numChunksY = int( math.ceil( float(ydim)/float(ychunk) ) )
    numChunksZ = int( math.ceil( float(zdim)/float(zchunk) ) )
    print "NumChunks : "+str(numChunksX)+" "+str(numChunksY)+" "+str(numChunksZ)
    numChunks =  numChunksX*numChunksY*numChunksZ
    #Write data one spectrum at a time
    for xi in xrange( 0 , xdim ) :
        sys.stdout.write("[" +str( int( 100.* float(xi)/float(xdim) )) +"%]"+ "\r")
        sys.stdout.flush()
        for yi in xrange( 0 , ydim ) :
            #Save the spectrum to the hdf5 file
            data[xi,yi,:] = (xi*ydim + yi)
    
    #Write data into all the chunks
    """for xt in xrange(0, numChunksX ) :
        xstart = xt*xchunk
        xend = min(  xstart+xchunk , xdim)
        for yt in xrange(0, numChunksY ) :
            ystart = yt*ychunk
            yend = min( ystart+ychunk , ydim )
            for zt in xrange(0, numChunksZ ) :
                zstart = zt*zchunk
                zend = min( zstart+zchunk , zdim )
                #print "Write : "+str(xstart)+" "+str(xend)+" "+str(ystart)+" "+str(yend)+" "+str(zstart)+" "+str(zend)
                data[xstart:xend , ystart:yend, zstart:zend ] = itertest
                itertest+=1
                sys.stdout.write("Generating Data: [" +str( int( 100.* float(itertest)/float(numChunks) )) +"%]"+ "\r")
                sys.stdout.flush()"""
    omsiFile .close_file()
    
def printHelp():
    """Print the help explaining the usage of testHDF5Optimiation"""
    
    print "USAGE: Call \"testHDF5Optimiation zval HDF5File resultsFile\" "


if __name__ == "__main__":
    main() 
