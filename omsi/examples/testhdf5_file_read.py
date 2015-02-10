"""Simple test script used to test the performance of different HDF5 
   optimizations (using chunking) to improve the performance of 
   hyperslap selections"""

from omsi.dataformat.omsi_file import *
import time
import json

def main(argv=None):
    """Then main function"""

    import sys
    from sys import argv,exit
    
    if argv is None:
        argv = sys.argv
        
    #Check for correct usage
    if len(argv) < 2 :
        printHelp()
        exit(0)
        
    if len(argv) == 8 :
        
        infile = argv[1]
        xmin = int(argv[2])
        xmax = int(argv[3])
        ymin = int(argv[4])
        ymax = int(argv[5])
        zmin = int(argv[6])
        zmax = int(argv[7])
        start = time.time()
        d = omsi_file(infile, 'r').get_experiment(0).get_msidata(0)
        loaddata = d[xmin:xmax , ymin:ymax, zmin:zmax]
        #content = json.dumps( loaddata.tolist() )
        stop = (time.time() - start) 
        print stop
        exit(0)
    
    import numpy as np
    import os
    import random
    import subprocess

    repeats = 50
    outfolder  = argv[1]
    if not outfolder.endswith("/") :
        outfolder = outfolder+"/"
    
    #Baseline filelist
    filelist =  [ "/project/projectdirs/openmsi/manuscript_data/baseline/11042008_NIMS.h5" , "/project/projectdirs/openmsi/manuscript_data/baseline/20121012_lipid_extracts.h5", "/project/projectdirs/openmsi/manuscript_data/baseline/20110929_Tumor624.h5" ,  "/project/projectdirs/openmsi/manuscript_data/baseline/Brain.h5" , "/project/projectdirs/openmsi/manuscript_data/baseline/20111012_Tumor458_50micronSS_D.h5" , "/project/projectdirs/openmsi/manuscript_data/baseline/Microbial_Coculture.h5", "/project/projectdirs/openmsi/manuscript_data/baseline/20111208_KBL_Roots_SmallChip_BigRoot.h5" , "/project/projectdirs/openmsi/manuscript_data/baseline/nimzyme.h5", "/project/projectdirs/openmsi/manuscript_data/baseline/20120801_metabolite_standards.h5", "/project/projectdirs/openmsi/manuscript_data/baseline/20111207_KBL_Roots_BigChip_SmallRoots.h5" ]
    
    #Compressed 4x4x2048 filelist
    #filelist =  [ "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/11042008_NIMS.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20121012_lipid_extracts.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20110929_Tumor624.h5" ,  "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/Brain.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20111012_Tumor458_50micronSS_D.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/Microbial_Coculture.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20111208_KBL_Roots_SmallChip_BigRoot.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/nimzyme.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20120801_metabolite_standards.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_4_4_2048/20111207_KBL_Roots_BigChip_SmallRoots.h5" ]
    
    #Uncompressed 4x4x2048 filelist
    #filelist =  [ "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/11042008_NIMS.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20121012_lipid_extracts.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20110929_Tumor624.h5" ,  "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/Brain.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20111012_Tumor458_50micronSS_D.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/Microbial_Coculture.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20111208_KBL_Roots_SmallChip_BigRoot.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/nimzyme.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20120801_metabolite_standards.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_4_4_2048/20111207_KBL_Roots_BigChip_SmallRoots.h5" ]
    
    #Uncompressed autochunking filelist
    #filelist =  [ "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/11042008_NIMS.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20121012_lipid_extracts.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20110929_Tumor624.h5" ,  "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/Brain.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20111012_Tumor458_50micronSS_D.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/Microbial_Coculture.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20111208_KBL_Roots_SmallChip_BigRoot.h5" , "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/nimzyme.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20120801_metabolite_standards.h5", "/project/projectdirs/openmsi/manuscript_data/uncompressed_autochunking/20111207_KBL_Roots_BigChip_SmallRoots.h5" ]
    
    #Compressed autochunking filelist
    #filelist =  [ "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/11042008_NIMS.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20121012_lipid_extracts.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20110929_Tumor624.h5" ,  "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/Brain.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20111012_Tumor458_50micronSS_D.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/Microbial_Coculture.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20111208_KBL_Roots_SmallChip_BigRoot.h5" , "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/nimzyme.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20120801_metabolite_standards.h5", "/project/projectdirs/openmsi/manuscript_data/compressed_autochunking/20111207_KBL_Roots_BigChip_SmallRoots.h5" ]
    
    
    #filelist = ["/work2/bowen/TEST_DEP3.h5"]  
    
    data_shapes = {}
    results = {}
    for filename in filelist :

        #Initialze data shape
        f = omsi_file( filename , 'r' )
        d = f.get_experiment(0).get_msidata(0)
        data_shapes[filename] = d.shape
        f.close_file()
        #Initialze output storage
        results[filename] = np.zeros( repeats ,  dtype=[ ('mz-slice','f') , ('spectrum','f') , ('xyz-cube','f'), ('mz-slice-all','f') , ('spectrum-all','f') , ('xyz-cube-all','f')  , ('filesize' , 'f') ] )
        results[filename]['filesize'] = os.stat( filename ).st_size

    #Note: We compute each test seperately so that we have touched enough data from 
    #other files to avoid biases due to already cached data already at the beginning
    #of the tests.
    #Note: Depending on the file system, a significant amount of data (and in some 
    #cases complete files) may be cached by the file system itself. This can results
    #in large variations in the times for data acceses. This is particularly the 
    #case when two consecutive accesses happen to by chance access a similar portion
    #of the data. This behavior is expected and is what one expects to happen in real
    #life as well. It is, therefore, often informative to look at the general 
    #variability of results. E.g.,if all data is stored in a single block, then we
    #may see very slow access times at the beginning and then, once, all data has
    #been cached access times drop significantly. For well-chunked data, this
    #variability between access should be much lower.

    #Compute the slice query test
    for filename in filelist :
        print filename+" 25 mz-slices"
        #mz-slice selection 250,000 elements
        sliceWidthZ = 25 #xdim=100 , ydim=100
        for ri in xrange( 0 , repeats ) :

            xmin = 0
            xmax = data_shapes[filename][0]
            ymin = 0
            ymax = data_shapes[filename][1]
            zmin = random.randint(0, data_shapes[filename][2]-sliceWidthZ-1 )
            zmax = zmin + sliceWidthZ
            callCommand = ["python", "testhdf5_file_read.py" , filename , str(xmin) , str(xmax), str(ymin), str(ymax), str(zmin) , str(zmax) ]
            start = time.time()
            p2 = subprocess.Popen(callCommand , stdout = subprocess.PIPE)
            readTime =  float(p2.stdout.read())
            stop = (time.time() - start)
            results[filename]['mz-slice'][ri] = readTime
            results[filename]['mz-slice-all'][ri] = stop
            print str(results[filename]['mz-slice'][ri]) + "   " +str( results[filename]['mz-slice-all'][ri] )+ " " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " " + str(zmin)  + " " + str(zmax)

    #Compute the spectra test
    for filename in filelist :
        print filename+" 3 x 3  spectra"
        #mz-slice selection 250,000 elements
        sliceWidthX = 3 
        sliceWidthY = 3
        for ri in xrange( 0 , repeats ) :

            xmin = random.randint(0, data_shapes[filename][0]-sliceWidthX-1 )
            xmax = xmin + sliceWidthX
            ymin = random.randint(0, data_shapes[filename][1]-sliceWidthY-1 )
            ymax = ymin + sliceWidthY
            zmin = 0
            zmax = data_shapes[filename][2]
            callCommand = ["python", "testhdf5_file_read.py" , filename , str(xmin) , str(xmax), str(ymin), str(ymax), str(zmin) , str(zmax) ]
            start = time.time()
            p2 = subprocess.Popen(callCommand , stdout = subprocess.PIPE)
            readTime =  float(p2.stdout.read())
            stop = (time.time() - start)
            results[filename]['spectrum'][ri] = readTime
            results[filename]['spectrum-all'][ri] = stop
            print str(results[filename]['spectrum'][ri]) + "   " +str( results[filename]['spectrum-all'][ri] )+ " " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " " + str(zmin)  + " " + str(zmax)

    #Compte the cube test 
    for filename in filelist :
        print filename+" 20 x 20 x 1000  cube"
        #mz-slice selection 250,000 elements
        sliceWidthX = 20
        sliceWidthY = 20
        sliceWidthZ = 1000
        for ri in xrange( 0 , repeats ) :

            xmin = random.randint(0, data_shapes[filename][0]-sliceWidthX-1 )
            xmax = xmin + sliceWidthX
            ymin = random.randint(0, data_shapes[filename][1]-sliceWidthY-1 )
            ymax = ymin + sliceWidthY
            zmin = random.randint(0, data_shapes[filename][2]-sliceWidthZ-1 )
            zmax = zmin + sliceWidthZ
            callCommand = ["python", "testhdf5_file_read.py" , filename , str(xmin) , str(xmax), str(ymin), str(ymax), str(zmin) , str(zmax) ]
            start = time.time()
            p2 = subprocess.Popen(callCommand , stdout = subprocess.PIPE)
            readTime =  float(p2.stdout.read())
            stop = (time.time() - start)
            results[filename]['xyz-cube'][ri] = readTime
            results[filename]['xyz-cube-all'][ri] = stop
            print str(results[filename]['xyz-cube'][ri]) + "   " +str( results[filename]['xyz-cube-all'][ri] )+ " " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax) + " " + str(zmin)  + " " + str(zmax)
    
    for filename in filelist :
        
        infilename = os.path.split( filename )[1]
        outfile = outfolder+infilename+"_timings.txt"
        
        f = open( outfile , 'w' )
        for colName in results[filename].dtype.names :
            f.write( colName+" " )
        f.write("\n")
        np.savetxt( f , results[filename] )
        f.close()

    exit(0)


def printHelp():
    """Print the help explaining the usage of testHDF5Optimiation"""
    
    print "USAGE: Call \"testhdf5_file_read resultsFile\" "
    print "Execute query: testhdf5_file_read filename xmin xmax ymin ymax zmin zmax"


if __name__ == "__main__":
    main() 
