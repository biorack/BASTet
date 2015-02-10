import ctypes as c
import numpy as np
import multiprocessing as mp

from omsi.dataformat.omsi_file import *
import numpy as np
import os
import time
import random
import sys
import math
import sys
from sys import argv,exit

selection = 1 #0=slice , 1= spectra


def main(argv=None):
    """Then main function"""

    global selection
    
    if argv is None:
        argv = sys.argv
        
    #Check for correct usage
    if len(argv) !=4 :
        printHelp()
        exit(0)
        
    omsiFile  = argv[1]
    resultFile = argv[2]
    numWorkers  = int( argv[3] )
    donorFile = "/work2/bowen/11042008_NIMS.h5" #"/project/projectdirs/openmsi/omsi_data/old/TEST.h5"
    repeats = 50
    xdim = 100
    ydim = 100
    zdim = 100000 

    if selection == 0 :
        mp_arr = mp.Array(c.c_uint16, (xdim*ydim) ) # shared, can be used from multiple processes
    else :
        mp_arr = mp.Array(c.c_uint16, (zdim) ) # shared, can be used from multiple processes
   
    allRes = np.zeros( numWorkers , dtype=[ ('x','int32') , ('y','int32') , ('z','int32') , ('numProc' , 'int32') ,  ('select_min','f') , ('select_avg','f') , ('select_median','f') , ('select_max','f') , ('select_std','f') , ('select_var','f')  ]  )

    for numP in range(1, (numWorkers+1) ) :

        #generateBaseTestFile( omsiFile , xdim , ydim, zdim )
        generateChunkedTestFile( omsiFile , xdim , ydim, zdim, 4, 4 , 2048  , True , donorFile )      
        allRes[numP-1]['x'] = 4
        allRes[numP-1]['y'] = 4
        allRes[numP-1]['z'] = 2048
        allRes[numP-1]['numProc'] = numP

        timings = np.zeros( repeats )
        for i in range(0,repeats ) :

            start = time.time()
            procs = map(create_process, map(lambda lst : (omsiFile, numP, mp_arr , xdim, ydim, zdim, lst) ,  range(0,numP) ))  
            map(lambda p : p.start(), procs)  
            map(lambda p : p.join(), procs)  
            timings[i]  = time.time() - start
            print timings[i]
            for p in procs :
                p.terminate()
       
        np.savetxt( "res_shaft_5x5SpectrumMean_"+str(numP)+".txt" , timings )
 
        allRes[numP-1]['select_min'] = np.min( timings) 
        allRes[numP-1]['select_avg'] = np.average( timings)
        allRes[numP-1]['select_median'] = np.median( timings)
        allRes[numP-1]['select_max'] = np.max( timings)
        allRes[numP-1]['select_std'] = np.std( timings)
        allRes[numP-1]['select_var'] = np.var( timings)

        print "----------"+str(numP)+"-----------"
        print "Max:" +str( np.max( timings ) )
        print "Median:"+ str( np.median( timings) )
        print "Mean:" + str( np.mean( timings) )
        print "Std: " + str( np.std( timings ) )
        print "Min: " + str( np.min( timings ) )

        os.remove( omsiFile  )
    
         
    f = open( resultFile , 'w' )
    f.write( str(allRes.dtype.names) )
    f.write("\n")
    np.savetxt( f , allRes )
    f.close()




def create_process( argv ) :  
   #print argv
   global selection
   if selection == 0 :
       p = mp.Process(target=sliceSelect, args=(argv,))  
   else :
       p = mp.Process(target=spectraSelect, args=(argv,))  
   #print p
   return p


def sliceSelect( args ):

    filename, numWorkers, mp_arr, xdim, ydim, zdim,  id = args

    arr = np.frombuffer(mp_arr.get_obj() , dtype = 'uint16'  )  # mp_arr and arr share the same memory
    b = arr.reshape((xdim,ydim)) # b and arr share the same memory

    zrange = 25 # 20000
    zmin = random.randint(0, zdim-zrange-1 )
    zmax = zmin+zrange+1

    #id =   int( mp.current_process()._identity[0] )
    omsiFile = omsi_file( filename , 'r' )
    d = omsiFile.get_experiment(0).get_msidata(0)
    xstep = xdim / numWorkers
    xstart = xstep*(id)
    xend  = xstep*(id+1)
    if id == (numWorkers-1) :
        xend = xdim 
    #print str(xstart) + "  " + str(xend)+"     :"+str(zmin)+" "+str(zmax)
    b[xstart:xend , :  ] = np.var( d[xstart:xend , 0:ydim , zmin:zmax ] , axis=2 )

    omsiFile.close_file()
    sys.stdout.flush()


def spectraSelect( args ):

    filename, numWorkers, mp_arr, xdim, ydim, zdim,  id = args

    arr = np.frombuffer(mp_arr.get_obj() , dtype = 'uint16'  )  # mp_arr and arr share the same memory
    b = arr.reshape( (zdim) ) # b and arr share the same memory


    xmin = random.randint(0, xdim-6 )
    xmax = xmin+5

    ymin = random.randint(0, ydim-6 )
    ymax = ymin+5

    #id =   int( mp.current_process()._identity[0] )
    omsiFile = omsi_file( filename , 'r' )
    d = omsiFile.get_experiment(0).get_msidata(0)
    zstep = zdim / numWorkers
    zmin = zstep*(id)
    zmax  = zstep*(id+1)
    if id == (numWorkers-1) :
        zmax = zdim 
    #print str(id)+": "+str(xmin) + "  " + str(xmax)+" | "+str(ymin) + "  " + str(ymax)+" | "+str(zmin)+" "+str(zmax)
    b[zmin:zmax] = np.mean( np.mean( d[xmin:xmax , ymin:ymax , zmin:zmax ] , axis=0 ) , axis =0 )
    
    omsiFile.close_file()
    sys.stdout.flush()



def generateBaseTestFile( omsiOutFile , xdim , ydim, zdim ) :
    
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
    start = time.time()
    #Allocate space in the HDF5 file for the img data
    data = exp.create_msidata(data_shape=( xdim , ydim , zdim  ) , data_type = 'uint16' , chunks=None )
    #Write data one spectrum at a time
    for xi in xrange( 0 , xdim ) :
        sys.stdout.write("[" +str( int( 100.* float(xi)/float(xdim) )) +"%]"+ "\r")
        sys.stdout.flush()
        for yi in xrange( 0 , ydim ) :
            #Save the spectrum to the hdf5 file
            data[xi,yi,:] = (xi*ydim + yi)

    omsiFile .close_file()
    return ( time.time() - start )





def generateChunkedTestFile( omsiOutFile , xdim , ydim, zdim, xchunk, ychunk , zchunk  , compress=False , donorFile = "/project/projectdirs/openmsi/omsi_data/old/TEST.h5" ) :
    
    writeFullSpectra=False
    useChunking = (xchunk>0 and ychunk>0 and zchunk>0 ) 
    print useChunking
    useDonorFile = True
    if useDonorFile :
        inFile = omsi_file( donorFile )
        inData = inFile.get_experiment(0).get_msidata(0)[:]
        inFile.close_file()


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
    start = time.time()
    if compress : 
        #Use compresion
        data = exp.create_msidata(data_shape=( xdim , ydim , zdim  ) , data_type = 'uint16' ,  chunks=(xchunk,ychunk,zchunk) ,  compression='gzip' , compression_opts=4 )
    elif useChunking :
        #Use chunking
        data = exp.create_msidata(data_shape=( xdim , ydim , zdim  ) , data_type = 'uint16' ,  chunks=(xchunk,ychunk,zchunk) ) #,  compression='gzip' , compression_opts=4 )
    else :
        #Don't use chunking and compression
        data = exp.create_msidata(data_shape=( xdim , ydim , zdim  ) , data_type = 'uint16' )
 
    itertest=0 
    if useChunking :

        numChunksX = int( math.ceil( float(xdim)/float(xchunk) ) )
        numChunksY = int( math.ceil( float(ydim)/float(ychunk) ) )
        numChunksZ = int( math.ceil( float(zdim)/float(zchunk) ) )
        print "NumChunks : "+str(numChunksX)+" "+str(numChunksY)+" "+str(numChunksZ)
        numChunks =  numChunksX*numChunksY*numChunksZ

    else :
       #Write on spectrum at a time if we want a contiguous data layout 
       numChunksX = int( math.ceil( float(xdim)/4.0 ) )
       numChunksY = int( math.ceil( float(ydim)/4.0 ) )
       numChunksZ = 2

    if not useDonorFile :
        if writeFullSpectra :
            print "Writing mxm spectra at a time (artifical data)"
            #Write data one x/y chunk at a time (i.e., multiple z-chunks are writtent at once)
            for xi in xrange( 0 , numChunksX ) :
                sys.stdout.write("[" +str( int( 100.* float(xi)/float(numChunksX) )) +"%]"+ "\r")
                sys.stdout.flush()
                xstart = xi*xchunk
                xend = min(  xstart+xchunk , xdim)
                for yi in xrange( 0 , numChunksY ) :
                    ystart = yi*ychunk
                    yend = min( ystart+ychunk , ydim )
                    #Save the spectrum to the hdf5 file
                    data[xstart:xend , ystart:yend, : ] = (xi*ydim + yi)
    
        else :
            #Write data one x/y/z chunk at a time
            print "Writing one x/y/z chunk at a time (artifical data)"
            for xt in xrange(0 , numChunksX) :
                sys.stdout.write("[" +str( int( 100.* float(xt)/float(numChunksX) )) +"%]"+ "\r")
                sys.stdout.flush()
                xstart = xt*xchunk
                xend = min(  xstart+xchunk , xdim)
                for yt in xrange(0, numChunksY ) :
                    ystart = yt*ychunk
                    yend = min( ystart+ychunk , ydim )
                    for zt in xrange(0, numChunksZ ) :
                        zstart = zt*zchunk
                        zend = min( zstart+zchunk , zdim )
                        data[xstart:xend , ystart:yend, zstart:zend ] =  (xt*ydim*zdim + yt*zdim  +zt)

    else :
        print "Writing one x/y/z chunk at a time (donor dataset)"
        #Write data into all the chunks one x,y,z chunk at a time using the donor file
        for xt in xrange(0, numChunksX ) :
            sys.stdout.write("[" +str( int( 100.* float(xt)/float(numChunksX) )) +"%]"+ "\r")
            sys.stdout.flush()
            xstart = xt*xchunk
            xend = min(  xstart+xchunk , xdim)
            for yt in xrange(0, numChunksY ) :
                ystart = yt*ychunk
                yend = min( ystart+ychunk , ydim )
                for zt in xrange(0, numChunksZ ) :
                    zstart = zt*zchunk
                    zend = min( zstart+zchunk , zdim )
                    #print "Write : "+str(xstart)+" "+str(xend)+" "+str(ystart)+" "+str(yend)+" "+str(zstart)+" "+str(zend)
                    diff = inData.shape[2] - zend
                    myend = zend
                    mystart = zstart
                    if inData.shape[2] < zend :
                        myend = inData.shape[2]
                        mystart = inData.shape[2]-(zend-zstart)
                    a = inData[xstart:xend , ystart:yend, mystart:myend ]
                    #print str(zt)+" : "+str(a.shape)+" : "+str(mystart)+" "+str(myend)
                    data[xstart:xend , ystart:yend, zstart:zend ] = a
        
                    itertest+=1
                    #sys.stdout.write("Generating Data: [" +str( int( 100.* float(itertest)/float(numChunks) )) +"%]"+ "\r")
                    #sys.stdout.flush()

    omsiFile .close_file()
    
    return (time.time() - start)





if __name__ == "__main__":
    main() 

