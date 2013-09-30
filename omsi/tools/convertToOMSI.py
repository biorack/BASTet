"""Tool used to convert img files to OpenMSI HDF5 files. 

  For usage information execute: python convertToOMSI --help
"""

from omsi.dataformat.bruckerflex_file import bruckerflex_file
from omsi.dataformat.img_file import img_file
from omsi.dataformat.omsi_file import *
from omsi.analysis.multivariate_stats.omsi_nmf import omsi_nmf
from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
from omsi.analysis.findpeaks.omsi_findpeaks_local import omsi_findpeaks_local
from omsi.analysis.omsi_analysis_data import *
import time
import numpy as np
import math
import sys
import os
from sys import argv,exit
try :
    from PIL import Image
except :
    pass


####################################################################
#  Variables used for internal storage                             #
####################################################################
""" :param dataList: List of python dictionaries decribing specific conversion \
             settings for each conversion task. Each dictionary contains the follwing keys:
    
             * 'basename' : Name of the file to be converted 
             * 'format' : File format to be used (see availableFormats) 
             * 'exp' : Indicate the experiment the dataset should be stored with. Valid values are \
                       
                          * 'new' : Generate a new experiment for the dataset
                          * 'previous' : Use the same experiment as used for the previous dataset 
                          * 1,2,3...   : Integer value indicating the index of the experiment to be used. 
            * 'region' : Optional key, indicating the region index to be converted or None if all regions should be merged. 
"""
datasetList = []
omsiFile = None #The openMSI output data file to be used.

####################################################################
#  Define available options for different parameters              ##
####################################################################
availableFormats = ["img", "brukerflex", "auto"] #List of available data formats
availableRegionOptions = ["split" , "merge" , "split+merge"]    #List defining the different options available for handling regions
availableioOptions = ["chunk" , "spectrum", "all" ] #Available options for the data write. One chunk at a time'chunk', one spectrum at a time ('spectrum') or all at one once ('all')

####################################################################
#  Define the input parameters used during the conversion         ##
####################################################################
ioOption = "chunk" #Define how the data should be written to file, one cunk at a time ('chunk'), one spectrum at a time ('spectrum') or all at one once ('all')
formatOption = "auto" #Define which file format reader should be used
regionOption = "split+merge" #Define the region option to be used
auto_chunk = True #Automatically decide which chunking should be used
chunks=(4,4,2048)  #Define the main chunking for the data
user_additional_chunks = [] #User-defined optional additional chunked copies of the data used to optimize data seletion
compression='gzip' #We use gzip by default because szip is not as generally available and LZF is usually only available with h5py.
compression_opts=4 #For gzip this a value between 0-9 (inclusive) indicating how agressive the compression should be. 0=fastest speed, 9=best compression ratio.
suggest_file_chunkings = False #Should we only suggest file chunkings and quit without converting any data

####################################################################
#  Define analysis options and parameter settings                 ##
####################################################################
executeNMF = True  #Define whether NMF should be performed
executeFPG = True  #Define whether global peak finding should be executed
executeFPL = False #Define whether local peak finding should be execuated
generate_thumbnail = True  #Should we generate thumbnail

#Default NMF parameter settings
nmf_numComponents=20 #Number of components for the NMF
nmf_timeout=600 #Timeout for the NMF
nmf_numIter=2000 #Maximum number of iterations for the NMF
nmf_tolerance=0.0001  #Tolerance level for the nmf
nmf_useRawData = False  #Should the NMF be computed from the raw data or from the global peak finding data


def main(argv=None):
    '''The main function defining the control flow for the conversion'''
    #Get the global variables
    global omsiFile
    global datasetList
    global formatOption
    global regionOption

    
    ####################################################################
    #   Determine the settings based on the user input                 #
    ####################################################################
    #Get the user input options
    if argv is None:
        argv = sys.argv
    #Parse the input arguments 
    inputError, inputWarning, omsiOutFile, inputFilenames = parseInputArgs( argv )
    #  Terminate in case an error or warning has occured while processing the user input parameters.
    if inputError :
        print "Terminated. One or more error occured while parsing the command line inputs."
    if inputWarning:
        print "Terminated. Conflicting input parameters found. See WARNINGS above for details."
    if inputWarning or inputError :
        exit()

    ####################################################################
    #   Create the output HDF5 file if needed                          #
    ####################################################################
    try:
        if omsiOutFile is not None : 
            omsiFile = omsi_file( omsiOutFile )
    except:
        print "Unexpected error creating the output file:", sys.exc_info()[0]
        exit(0)
    
    ####################################################################
    # Generate the list of datasets to be converted                    #
    ####################################################################
    datasetList = create_datasetList( inputFilenames=inputFilenames, formatOption=formatOption, regionOption=regionOption)
    print "Number of conversion: "+str(len( datasetList ))

    ####################################################################
    #  Suggest only chunking for the files if requested                #
    ####################################################################
    if suggest_file_chunkings : 
        suggest_chunkings_for_files( datasetList )
        exit()

    ####################################################################
    # Convert all files                                                #
    ####################################################################
    try :
        convert_files()
    except:
        omsiFile.close_file()
        print "ERROR: An error occured during the file conversion. Closing the output file and terminating."
        raise

    ####################################################################
    #  Close the HDF5 file and exit                                    #
    ####################################################################
    omsiFile.close_file()
    exit(0)


def convert_files() :
    """Convert all files in the given list of files with their approbriate conversion options"""

    #Get the global variables
    global availableFormats
    global availableRegionOptions
    global ioOption
    global formatOption
    global regionOption
    global auto_chunk
    global chunks
    global user_additional_chunks
    global compression
    global compression_opts
    global suggest_file_chunkings 
    global executeNMF
    global executeFPG
    global executeFPL
    global generate_thumbnail
    global nmf_numComponents
    global nmf_timeout
    global nmf_numIter
    global nmf_tolerance
    global nmf_useRawData
    global datasetList
    global omsiFile


    ####################################################################
    #  Convert the MSI files and compute the requested analyses       ##
    ####################################################################
    #Iterate over all img files
    for i in datasetList :  # xrange(startIndex,len(argv)-1) :

        ####################################################################
        #  Convert the raw data to HDF5                                   ##
        ####################################################################
        basefile = i['basename'] #argv[i].replace('"' , "")
        print "Converting: " + basefile
        try :
            print "Selected Region: "+str(i['region'])
        except :
            pass
        print "HDF5 chunking: "+str(chunks)
        print "HDF5 compression: "+str(compression)+", "+str(compression_opts)

        #Open the img file
        try:
            currFormat = i['format']
            if currFormat is "img" :
                inputFile = img_file( hdrFile=basefile+".hdr" , t2mFile=basefile+".t2m" , imgFile=basefile+".img" )
            elif currFormat is "bruckerflex" :
                inputFile = bruckerflex_file( spotlist_filename=basefile )
                inputFile.set_region_selection( region_index=i['region'])
            else :
                print "ERROR: The following file will not be converted because the file type could not be determined: "+basefile
                print "INFO: If you know the correct file format then try converting the file using the approbirate --format <formatname> option."
                
            print "In data shape: "+str(inputFile.shape) 
        except:
            print "ERROR: Unexpected error opening the input file:", sys.exc_info()[0]
            print basefile+" failure during converion. Skipping the file and continue conversion of the other files.\n" 
            continue 

        if auto_chunk : 
        
            spec_chunks , img_chunks , tempChunks = suggest_chunking( xsize=inputFile.shape[0], ysize=inputFile.shape[1], mzsize=inputFile.shape[2], dtype=inputFile.data_type, print_results=True)
            chunks = spec_chunks
            additional_chunks = user_additional_chunks + [img_chunks]
            print "Converting data using the following chunking options:"
            print "     - Spectrum chunking: "+str(chunks)
            print "     - Image chunking:    "+str(additional_chunks[0])

	else :
	    additional_chunks = user_additional_chunks

        #Get the mz data
        mzdata = inputFile.mz

        #Define the layout for the img data in the HDF5
        #Create a new experiment for the img file using the basefile as identifer string
        if i['exp'] == 'new': #Create a new experiment for this dataset
            exp = omsiFile.create_exp( exp_identifier = basefile )
            #Create an empty sample descrition
            sample = exp.create_sample_info()
            #Create an empty instrument description
            instrument = exp.create_instrument_info(instrument_name="undefined" , mzdata=mzdata )
        elif i['exp'] == 'previous' : #Store this dataset in combination with the previous one
            exp = omsiFile.get_exp( omsiFile.get_num_exp() -1  )
        else :  #Store this dataset in combination with the given experiment index
            exp = omsiFile.get_exp( i['exp'] )
        
        #Allocate space in the HDF5 file for the img data
        data_dataset, mz_dataset, data_group = exp.create_msidata_full_cube( data_shape=inputFile.shape , data_type = inputFile.data_type , chunks=chunks, compression=compression, compression_opts=compression_opts)
        mz_dataset[:] = mzdata
        data = omsi_file_msidata( data_group=data_group , preload_mz=False, preload_xy_index=False) 
        write_data( inputFile, data , ioOption=ioOption, chunks=chunks)
        omsiFile.flush()
        
        #Generate any additional data copies if requested
        for c in additional_chunks :
            print "Generating optimized data copy: "+str(c)
            tempData = data.create_optimized_chunking(chunks=c, compression=compression, compression_opts=compression_opts , copy_data=False, print_status=True)
            write_data( inputFile, tempData, ioOption=ioOption, chunks=c)
            omsiFile.flush()

        #Close the img file and free up any allocated memory
        inputFile.close_file()
        
        ####################################################################
        #  Execute the requested analyses                                 ##
        ####################################################################
        #Compute the local peak finding and add the peak-cube data to the file
        if executeFPL:
            print "Executing local peak finding"
            #Execute the peak finding
            fpl = omsi_findpeaks_local( nameKey = "omsi_findpeaks_local_"+str(time.ctime()) )
            fpl.execute_peakfinding(data , mzdata, printStatus=True, msidata_dependency=    data)
            #Save the peak-finding to file
            ana , ai = exp.create_analysis( fpl )
        #Compute the global peak finding and add the peak-cube data to the file
        if executeFPG:
            print "Executing global peak finding"
            #Execute the peak finding
            fpg = omsi_findpeaks_global( nameKey = "omsi_findpeaks_global_"+str(time.ctime()) )
            fpg.execute_peakfinding(data , mzdata, msidata_dependency=data)
            #Save the peak-finding to file
            ana , ai = exp.create_analysis( fpg )
            fpgPath = ana.get_h5py_analysisgroup().name
        #Compute the nmf analysis if requested
        if executeNMF :
            print "Executing nmf"
            #Exectue the nmf analysis on the peak-cube or if peak finding was not performed then 
            #execute nmf on the raw data data
            nmf = omsi_nmf( nameKey = "omsi_nmf_"+str(time.ctime()) )
            if executeFPG and not nmf_useRawData :
                print "   Using peak-cube data for NMF"
                fpgData = fpg[ 'peak_cube' ]
                dataPath = fpgPath+"/peak_cube"
                dataMzPath =  fpgPath+"/peak_mz"
                nmf.execute_nmf( fpgData , nmf_numComponents, nmf_timeout, nmf_numIter, nmf_tolerance,  msidata_dependency=ana)
            else :
                print "   Using raw data for NMF"
                nmf.execute_nmf( data , nmf_numComponents, nmf_timeout, nmf_numIter, nmf_tolerance, msidata_dependency=data)
            #Save the nmf results to file
            ana, ai = exp.create_analysis( nmf )
        
        
        ####################################################################
        #  Generate the thumbnail image for the current file              ##
        ####################################################################
        if generate_thumbnail : 
            print "Generating the thumbnail image"
            if executeNMF :
                print "   Generating thumbnail from NMF data"
                #Get the NMF data 
                ho = nmf['ho']
                Nx = ho.shape[0]
                Ny = ho.shape[1]
                #Generate images for the first three NMF components
                d1 = np.log( ho[:,:,0].reshape( (Nx,Ny) )+1)
                d1 = d1 / np.max( d1 )
                im1 = Image.fromarray( d1.astype('float')*255 ).convert('L')
                d2 = np.log( ho[:,:,1].reshape( (Nx,Ny) )+1)
                d2 = d2 / np.max( d2 )
                im2 = Image.fromarray( d2.astype('float')*255 ).convert('L')
                d3 = np.log( ho[:,:,2].reshape( (Nx,Ny) )+1)
                d3 = d3 / np.max( d3 )
                im3 = Image.fromarray( d3.astype('float')*255 ).convert('L')
                #Generate thumbnail by merging the three gray-scale images as an RGB image
                thumbnail = Image.merge( 'RGB' , (im1,im2,im3) )
                expName = str( exp.get_h5py_experimentgroup().name )
                expIndex = expName[7:len(expName)]
                thumbnailFilename = omsiOutFile+"_"+expIndex+".png"
                thumbnail.save( thumbnailFilename , 'PNG' )
            elif executeFPG :
                print "    Generating thumbnail from FPG data"
                #Get the global peak finding data and compute the maximum peak values for each peak
                fpgData = fpg[ 'peak_cube' ][:]
                maxPeakValues = fpgData.max(axis=0).max(axis=0)
                Nx = fpgData.shape[0]
                Ny = fpgData.shape[1]
                #Generate images for the three most intense peaks
                s = np.argsort( maxPeakValues )
                d1 = np.log( fpgData[:,:,s[-1]].reshape( (Nx,Ny) )+1)
                d1 = d1 / np.max( d1 )
                im1 = Image.fromarray( d1.astype('float')*255 ).convert('L')
                d2 = np.log( fpgData[:,:,s[-2]].reshape( (Nx,Ny) )+1)
                d2 = d2 / np.max( d2 )
                im2 = Image.fromarray( d2.astype('float')*255 ).convert('L')
                d3 = np.log( fpgData[:,:,s[-3]].reshape( (Nx,Ny) )+1)
                d3 = d3 / np.max( d3 )
                im3 = Image.fromarray( d3.astype('float')*255 ).convert('L')
                #Generate thumbnail by merging the three gray-scale images as an RGB image
                thumbnail = Image.merge( 'RGB' , (im1,im2,im3) )
                expName = str( exp.get_h5py_experimentgroup().name )
                expIndex = expName[7:len(expName)]
                thumbnailFilename = omsiOutFile+"_"+expIndex+".png"
                thumbnail.save( thumbnailFilename , 'PNG' )
            else :
                print "Generation of thumbnail from raw data is not yet supported. No thumbnail has been generated."
                print "Enable --nmf or --fpg in order to generate a thumbnail image."
                #print "    Generating thumbnail from raw data"
                ##Find three most intense peaks that are at least 1% of the m/z range appart
                #Nx = data.shape[0]
                #Ny = data.shape[1]
                #Nz = data.shape[2]
                #minMzStep = Nz / 100
                #mzImageRange = Nz/4000
                #maxPeakValues = np.zeros( Nz )
                #stepping = np.arange( 0 , Nz-5000 , 5000 )
                #for i in stepping :
                    #maxPeakValues[i:(i+5000)] = data[:,:,i:(i+5000)].std(axis=0).std(axis=0).reshape(5000)
                #maxPeakValues[stepping[-1]:Nz] = data[:,:,stepping[-1]:Nz].std(axis=0).std(axis=0).reshape( Nz - stepping[-1]  )
                #print maxPeakValues.shape
                #s = np.argsort( maxPeakValues )
                #i1 = s[-1]
                #i2 = s[-2]
                #for i in reversed( range(  0 , len(s) ) ) :
                    #if abs(s[i] - i1) > minMzStep :
                        #i2 = s[i]
                        #break
                #i3 = s[-3]
                #for i in reversed( range(  0 , len(s) ) ) :
                    #if abs(s[i] - i1) > minMzStep and abs(s[i]-i2)>minMzStep :
                        #i2 = s[i]
                        #break
                #print minMzStep
                #print str( i1)+" "+str(i2)+" "+str(i3)
                #print str( mzdata[i1])+" "+str(mzdata[i2])+" "+str(mzdata[i3])
                #low = max(0 , s[i1]-mzImageRange)
                #hi = min( Nz , s[i1]+mzImageRange)
                #print str(low) + " " + str(hi)
                #d1 = data[:,:,low:hi].max(axis=2).reshape( (Nx,Ny) )
                #d1 = d1 / np.max( d1 )
                #im1 = Image.fromarray( d1.astype('float')*255 ).convert('L')
                #low = max(0 , s[i2]-mzImageRange)
                #hi = min( Nz , s[i2]+mzImageRange)
                #print str(low) + " " + str(hi)
                #d2 =  data[:,:,low:hi].max(axis=2).reshape( (Nx,Ny) )
                #d2 = d2 / np.max( d2 )
                #im2 = Image.fromarray( d2.astype('float')*255 ).convert('L')
                #low = max(0 , s[i3]-mzImageRange)
                #hi = min( Nz , s[i3]+mzImageRange)
                #print str(low) + " " + str(hi)
                #d3 = data[:,:,low:hi].max(axis=2).reshape( (Nx,Ny) )
                #d3 = d3 / np.max( d3 )
                #im3 = Image.fromarray( d3.astype('float')*255 ).convert('L')
                ##Generate thumbnail by merging the three gray-scale images as an RGB image
                #thumbnail = Image.merge( 'RGB' , (im1,im2,im3) )
                #expName = str( exp.get_h5py_experimentgroup().name )
                #expIndex = expName[7:len(expName)]
                #thumbnailFilename = omsiOutFile+"_"+expIndex+".png"
                #thumbnail.save( thumbnailFilename , 'PNG' )

    
    ####################################################################
    #  Generate the XDMF header file for the HDF5 file                ##
    ####################################################################
    omsiFile.write_XDMF_header( omsiFile.get_filename() +".xdmf" )
    


def check_format( name , formatOption ) :
    """Helper function used to determine the file format that should be used
    
       :param name: Name of the folder/file that we should read
       :param formatOption: String indicating the format-option given by the user. \
                    If the format is not determined (i.e., "auto") then this function \
                    tries to determine the approbriate foramt. Otherwise this option \
                    is returned as is, as the user explicilty said which format should \
                    be used. 
       :returns: String indicating the approbriate format. Returns None in case no valid option was found.
    """
    global availableFormats
    #Option 1: The user told us the format we should use
    if formatOption is not "auto":
        if formatOption in availableFormats :
            return formatOption
        else :
            return None
    #Option 2: We need to determine the format ourselves
    if name.endswith( "Spot List.txt") :
        return "bruckerflex"
    elif os.path.exists( name+".hdr") and os.path.exists( name+".t2m") and os.path.exists( name+".img") :
        return "img"
    else :
        return None


def create_datasetList( inputFilenames, formatOption='auto', regionOption="split+merge" ) :
    """Based on the list of inputFilenames, generate the datasetList, which contains a dictionary describing
       each conversion job
       
       :param inputFilenames: List of names of files to be converted.
       :param formatOption: Define which file-format should be used. Default value is 'auto' indicating the \
                   function should determine for each file the format to be used. See also availableFormats parameter.
       :param regionOption: Define how different regions defined for a file should be handled. E.g., one may want \
                    to split all regions into indiviudal datasets ('split'), merge all regions into a single \
                    dataset ('merge'), or do both ('split+merge'). See also the availableRegionOptions parameter \
                    for details. By default the function will do 'split+merge'.
       
       :returns: List of dictionaries describing the various conversion jobs
    """
    re_datasetList = []
    for i in inputFilenames :
        currDS = {}
        currDS['basename'] = i.rstrip('"').rstrip("'").lstrip('"').lstrip("'") #Remove " in case the user has entered file names with spaces using the "name" syntax
        currDS['format']   = check_format( name=currDS['basename'] , formatOption=formatOption )
        currDS['exp'] = 'new'
        if currDS['format'] is 'bruckerflex' :
            if regionOption == "merge" or regionOption =="split+merge" :
                currDS['region'] = None
                re_datasetList.append( currDS )
        
            if regionOption == 'split' or regionOption == 'split+merge':
                try :
                    tempFile = bruckerflex_file( spotlist_filename=currDS['basename'] , readall=False)
                    for i in xrange(0,tempFile.get_number_of_regions() ) :
                        nDS = {}
                        nDS['basename'] = currDS['basename']
                        nDS['format']   = currDS['format']
                        nDS['region'] = i
                        nDS['exp'] = 'previous'
                        re_datasetList.append( nDS )
                except :
                    print "ERROR: Unexpected error opening the input file:", sys.exc_info()[0]
                    print currDS['basename']+" failure during converion. Skipping the file and continue conversion of the other files.\n"
                    continue

            if regionOption not in availableRegionOptions : #This is just to make sure that we have checked everything
                print "WARNING: Undefined region option. Using the default option split+merge"
        else :
            re_datasetList.append( currDS )
    
    #Return the list of conversion jobs
    return re_datasetList



def suggest_chunkings_for_files( datasetList ) :
    """Helper function used to suggest food chunking strategies for a given set of files.
    
       :param datafiles: Python list of dictionaries describing the settings to be used for the file conversion
       
       :returns: This function simply prints results to standard-out but does not return anything.
    """
    for i in datasetList :
        
        basefile = i["basename"]
        try:
            print "Suggested Chunkings: " + basefile
            currFormat = i["format"]
            if currFormat is "img" :
                inputFile = img_file( hdrFile=basefile+".hdr" , t2mFile=basefile+".t2m" , inputFile=basefile+".img" )
            elif currFormat is "bruckerflex" :
                inputFile = bruckerflex_file( spotlist_filename=basefile , readall=False )
                inputFile.set_region_selection( i["region"] )
            else :
                print "WARNING: Type of file could not be determined for: "+basefile
                continue
            print "In data shape: "+str(inputFile.shape)
            suggest_chunking( xsize=inputFile.shape[0], ysize=inputFile.shape[1], mzsize=inputFile.shape[2], dtype=inputFile.data_type, print_results=True)
        except :
             print "Error while trying to generate chunking suggestion for "+basefile


def suggest_chunking( xsize , ysize, mzsize, dtype, print_results=False ) :
    """Helper function used to suggest god chunking strategies for a given data cube 
    
       :param xsize: Size of the dataset in x.
       :param ysize: Size o the dataset in y.
       :param mzsize: Size of the dataset in mz.
       :param print_results: Print the results to the console. 
       
       :returns: Three tupes:
       
        * ``spectrum_chunk`` : The chunking to be used to optimize selection of spectra.
        * ``slice_chunk`` : The chunking to be used to optimize selection of image slices.
        * ``balanced_chunk`` : The chunking that would provide a good balance in performance for \
                            different selection strategies.
    """
    #Make sure that all sizes are treated as 64bit int to avoid errors due to artificial cutoffs.
    xsize = int(xsize)
    ysize = int(ysize)
    mzsize = int(mzsize)
    imageSize = xsize * ysize
    
    #Variable settings
    numBytes = np.dtype( dtype ).itemsize #Number of bytes per data value
    suggestedNumValues = (1024*64) /  numBytes #Nober of values that fit into a 64Kb chunk
    maxFactor = 4 #How much larger than the suggestNumValues should we allow a chunk to become

    #Define the chunking that would be good for spectra
    spectrum_xchunk = 1
    spectrum_ychunk = 1
    factor1 = math.floor(mzsize/float(suggestedNumValues))
    spectrum_mzchunk1 =  int( math.ceil(mzsize/ float(factor1)) )
    factor1 = math.ceil( mzsize / float(spectrum_mzchunk1) )
    overhead1 = ( (factor1*spectrum_mzchunk1) - mzsize ) *  imageSize * numBytes
    #overhead1 = ( (spectrum_mzchunk1 *  math.ceil (float(mzsize) / float(spectrum_mzchunk1))) - \
    #              (spectrum_mzchunk1 *  math.floor(float(mzsize) / float(spectrum_mzchunk1))) ) * \
    #            imageSize * numBytes
    #overhead1 = (spectrum_mzchunk1 - math.ceil(float(mzsize) % float(spectrum_mzchunk1))) * imageSize * numBytes
    factor2 = math.ceil(mzsize/float(suggestedNumValues))
    spectrum_mzchunk2 =  int( math.ceil(mzsize/ float(factor2)) )
    factor2 = math.ceil( mzsize / float(spectrum_mzchunk2) )
    overhead2 = ( (factor2*spectrum_mzchunk2) - mzsize ) *  imageSize * numBytes
    #overhead2 = (spectrum_mzchunk2 - math.ceil(float(mzsize) % float(spectrum_mzchunk2))) * imageSize *numBytes
    #overhead2 = ( (spectrum_mzchunk2 *  math.ceil (float(mzsize) / float(spectrum_mzchunk2))) - \
    #              (spectrum_mzchunk2 *  math.floor(float(mzsize) / float(spectrum_mzchunk2))) ) * \
    #            imageSize * numBytes
    if overhead1 < overhead2 : 
        spectrum_mzchunk = spectrum_mzchunk1
        spectrum_chunk_overhead = overhead1
    else :
        spectrum_mzchunk = spectrum_mzchunk2
        spectrum_chunk_overhead = overhead2
    spectrum_chunk = (spectrum_xchunk, spectrum_ychunk, spectrum_mzchunk )
    if print_results :
        
        print "Spectrum selection chunking: "+str(spectrum_chunk)
        print "     - Ideal for selection of full spectra."
        print "     - Overhead: "+str(spectrum_chunk_overhead)+" Byte ("+str( int(math.ceil(spectrum_chunk_overhead / (1024.*1024.))))+" MB)"
    
    #Define a chunking that would be good for images
    slice_xchunk = xsize
    slice_ychunk = ysize
    slice_mzchunk = 1
    chunkSize = slice_xchunk*slice_ychunk*slice_mzchunk
    if math.ceil(float(chunkSize) / float(suggestedNumValues)) > maxFactor :
        slice_xchunk = int(math.ceil( xsize / 2. ))
        slice_ychunk = int(math.ceil( ysize / 2. ))
    slice_chunk = ( slice_xchunk , slice_ychunk, slice_mzchunk )
    slice_chunk_xoverhead = ( slice_xchunk * math.ceil( float(xsize) / float(slice_xchunk))  ) - xsize #math.ceil( float(xsize) % float(slice_xchunk) )
    slice_chunk_yoverhead = ( slice_ychunk * math.ceil( float(ysize) / float(slice_ychunk))  ) - ysize #math.ceil( float(ysize) % float(slice_ychunk) )
    slice_chunk_overhead = ( (slice_chunk_xoverhead*ysize) + (slice_chunk_yoverhead*xsize) - (slice_chunk_xoverhead*slice_chunk_yoverhead) ) * mzsize  * numBytes
    if print_results : 
    
        print "Slice selection chunking: "+str(slice_chunk) 
        print "     - Ideal for selection of full image slices."
        print "     - Overhead: "+str(slice_chunk_overhead)+" Byte ("+str( int(math.ceil(slice_chunk_overhead / (1024.*1024.))))+" MB)"
    
    #Define a balanced chunkgin 
    balanced_xchunk  = 4
    balanced_ychunk  = 4 
    balanced_mzchunk = 2048
    balanced_chunk = (balanced_xchunk , balanced_ychunk , balanced_mzchunk )
    if print_results :
        
        print "Balanced chunking: " + str(balanced_chunk)
        print "     - This chunking tries to compromise between selection of slices and spectra."
        
    return spectrum_chunk, slice_chunk, balanced_chunk
    
    
def write_data( inputFile, data , ioOption="spectrum", chunks=None) :
    """Helper function used to implement different data write options.

        :param inputFile: The input img data file
        :param data: The output dataset (either an h5py dataset or omsi_file_msidata object.
        :param ioOption: String indicating the data write method to be used. One of:

            * ``spectrum``: Write the data one spectrum at a time 
            * ``all`` : Write the complete dataset at once.
            * ``chunk`` : Write the data one chunk at a time.

        :param chunks: The chunking used by the data. Needed to decide how the data should be written when a chunk-aligned write is requested.

    """
    if ioOption == "spectrum" or (ioOption == "chunk" and (chunks is None)) :
        for xi in xrange( 0 , inputFile.shape[0] ) :
            sys.stdout.write("[" +str( int( 100.* float(xi)/float(inputFile.shape[0]) )) +"%]"+ "\r")
            sys.stdout.flush()
            for yi in xrange( 0 , inputFile.shape[1] ) :
                #Save the spectrum to the hdf5 file
                data[xi,yi,:] = inputFile[xi , yi, :]
    elif ioOption == "all" :
        data[:] = inputFile[:]
    elif ioOption == "chunk" :
        xdim = inputFile.shape[0]
        ydim = inputFile.shape[1]
        zdim = inputFile.shape[2]
        numChunksX = int( math.ceil( float(xdim)/float(chunks[0]) ) )
        numChunksY = int( math.ceil( float(ydim)/float(chunks[1]) ) )
        numChunksZ = int( math.ceil( float(zdim)/float(chunks[2]) ) )
        numChunks = numChunksX*numChunksY*numChunksZ
        itertest=0 
        for xt in xrange(0, numChunksX ) :
            xstart = xt*chunks[0]
            xend = min(  xstart+chunks[0] , xdim)
            for yt in xrange(0, numChunksY ) :
                ystart = yt*chunks[1]
                yend = min( ystart+chunks[1] , ydim )
                for zt in xrange(0, numChunksZ ) :
                    zstart = zt*chunks[2]
                    zend = min( zstart+chunks[2] , zdim )
                    #print "Write : "+str(xstart)+" "+str(xend)+" "+str(ystart)+" "+str(yend)+" "+str(zstart)+" "+str(zend)
                    data[xstart:xend , ystart:yend, zstart:zend ] = inputFile[xstart:xend , ystart:yend, zstart:zend ]
                    itertest+=1
                    sys.stdout.write("[" +str( int( 100.* float(itertest)/float(numChunks) )) +"%]"+ "\r")
                    sys.stdout.flush()


def parseInputArgs( argv ) :
    """Process input parameters and define the script settings.
    
       :param argv: The list of input arguments 
       
       :returns: This function returns the following four values:
       
            * 'inputError' : Boolean indicating whether an error has occured during the processing of the inputs
            * 'inputWarning' : Boolean indicating whether a warning has been raised during the processing of the inputs
            * 'outputFilename' : Name for the output HDF5 file
            * 'inputFilenames' : List of strings indicating the list of input filenames
    """
    #Get the global variables
    global availableFormats
    global availableRegionOptions
    global ioOption
    global formatOption
    global regionOption
    global auto_chunk
    global chunks
    global user_additional_chunks
    global compression
    global compression_opts
    global suggest_file_chunkings 
    global executeNMF
    global executeFPG
    global executeFPL
    global generate_thumbnail
    global nmf_numComponents
    global nmf_timeout
    global nmf_numIter
    global nmf_tolerance
    global nmf_useRawData
    
    #Initalize the output values to be returned 
    inputError = False          #Parameter used to track whether an error has occured while parsing the input parameters
    inputWarning = False        #Parameter used to indicate conflicting input parameter options
    outputFilename = None       #The output filename
    inputFilenames = []         #The list of input filenames 
    #Basic sanity check
    if len(argv) <3 :
        ar = argv[-1] 
        if ar=="--help" or ar=="--h" or ar=="-help" or ar=="-h":
	    printHelp()
            exit(0)
        else :
            inputError = True
            return inputError, inputWarning, outputFilename, inputFilenames
      


    #Using the inputError and inputWarning paramter we try to keep parsing the input as long as possible to
    #make sure we can inform the user of as many problems as possible 
    startIndex = 1
    i=startIndex
    while i < (len(argv)-1) :
        i = startIndex
        ar = argv[i]
        if ar == "--no-nmf" :
            startIndex = startIndex+1
            executeNMF = False
            print "Disable NMF"
        elif ar == "--nmf" :
            startIndex = startIndex+1
            executeNMF = True
            print "Enable NMF"
        elif ar == "--nmf-nc" :
            startIndex = startIndex+2
            nmf_numComponents = int(argv[i+1])
            print "Set nmf-nc="+str(nmf_numComponents)
        elif ar == "--nmf-timeout" :
            startIndex = startIndex+2
            nmf_timeout = int(argv[i+1])
            print "Set nmf_timeout="+str(nmf_timeout)
        elif ar == "--nmf-niter" :
            startIndex = startIndex+2
            nmf_numIter = int(argv[i+1])
            if nmf_numIter < 2 :
                inputWarning = True 
                print "WARNING: --nfm-niter must be 2 or larger"
            print "Set nmf-niter="+str(nmf_numIter)
        elif ar == "--nmf-tolerance" :
            startIndex = startIndex+2
            nmf_numIter = float(argv[i+1])
            print "Set nmf-tolerance="+str(nmf-tolerance)
        elif ar == "--nmf-raw" :
            startIndex = startIndex+1
            nmf_useRawData = True
            print "Set nmf-raw="+str(nmf-raw)
        elif ar == "--fpg" :
            startIndex = startIndex+1
            executeFPG = True
            print "Enable find peaks global"
        elif ar == "--no-fpg" :
            startIndex = startIndex+1
            executeFPG = False
            print "Disable find peaks global"
            if "--fpg" in argv : 
                inputWarning = True 
                print "WARNING: --no-fpg and --fpg options are conflicting."
        elif ar == "--fpl" :
            startIndex = startIndex+1
            executeFPL = True
            print "Enable find peaks local"
        elif ar == "--no-fpl" :
            startIndex = startIndex+1
            executeFPL = False
            print "Disable find peaks local"
            if "--fpl" in argv : 
                inputWarning = True 
                print "WARNING: --no-fpl and --fpl options are conflicting."
        elif ar == "--auto-chunking" :
            startIndex = startIndex+1
            auto_chunk = True
            if "--chunking" in argv :
                inputWarning = True
                print "WARNING: --chunking options will be ignored due to the use of --auto-chunking"
            if "--no-chunking" in argv :
                inputWarning = True 
                print "WARNING: --no-chunking and --auto-chunking options are conflicting"
            if "--optimized-chunking"  in argv :
                inputWarning = True
                print "WARNING: --optimized-chunking and --auto-chunking options are conflicting"
        elif ar == "--chunking" :
            startIndex = startIndex+4
            try : 
                chunks = ( int(argv[i+1]) , int(argv[i+2]), int(argv[i+3] )   )
            except:
                print "An error accured while parsing the --chunking command. Something may be wrong with the indicated chunk sizes for x y z."
                inputError = True
            if "--auto-chunking" not in argv :
                auto_chunk = False 
            print "Enable chunking: "+str(chunks)
        elif ar == "--no-chunking" :
            startIndex = startIndex+1
            chunks = None
            auto_chunk = False
            print "Disable chunking"
            if "--auto-chunking" in argv or "--chunking" in argv :
                inputWarning = True
                print "WARNGING: --no-chunking option is conflicting with another chunking option"
        elif ar == "--optimized-chunking" :
            startIndex = startIndex+4
            try :
                user_additional_chunks.append( ( int(argv[i+1]) , int(argv[i+2]), int(argv[i+3] )   ) )
            except : 
                print "An error accured while parsing the --optimized-chunking command. Something may be wrong with the indicated chunk sizes for x y z."
                inputError = True
        elif ar == "--compression" :
            startIndex = startIndex+1
            #This is already the default
            print "Enable compression"
        elif ar == "--no-compression" :
            startIndex = startIndex+1
            compression=None
            compression_opts=None
            print "Disable compression"
            if "--compression" in argv : 
                inputWarning = True 
                print "WARNING: --no-compression and --compression options are conflicting."
        elif ar == "--io" :
            startIndex = startIndex + 2
            try :
                ioOption = str(argv[i+1])
                if ioOption!="spectrum" and ioOption!="all" and ioOption != "chunk" :
                    raise ValueError("Invalid io option")
            except:
                print "An error accured while parsing the --io command. Something may be wrong with the indicated io-type."
                inputError = True
        elif ar == "--thumbnail" :
            startIndex = startIndex+1
            generate_thumbnail = True
            print "Enable thumbnail"
        elif ar == "--no-thumbnail" :
            startIndex = startIndex+1
            generate_thumbnail = False
            print "Disable thumbnail"
            if "--thumbnail" in argv : 
                inputWarning = True 
                print "WARNING: --no-thumbnail and --no-thumbnail options are conflicting."
        elif ar=="--help" or ar=="--h" or ar=="-help" or ar=="-h":
            printHelp()
            exit(0)
        elif ar=="--suggest-chunking" :
            startIndex = startIndex+1
            suggest_file_chunkings = True
        elif ar=="--format" :
            startIndex = startIndex+2
            formatOption = str(argv[i+1])
            if formatOption not in availableFormats :
                print "ERROR: The indicated --format option "+formatOption+" is not supported. Available options are:"
                print "     "+str(availableFormats)
                inputError=True
        elif ar=="--regions" :
            startIndex = startIndex+2
            regionOption = str(argv[i+1])
            if regionOption not in availableRegionOptions :
                print "ERROR: The indicated --regions option "+regionOption+" is not supported. Avaiable options are:"
                print "       "+str(availableRegionOptions)
                inputError=True
        elif ar.startswith("--") :
            startIndex = startIndex+1
            inputError = True
            print "Unrecognized input option: "+ar
        else :
            #print "End of input parameters"
            break
    

    #Determine the list of input filenames
    inputFilenames =  [ name.replace('"' , "") for name in argv[ startIndex:(len(argv)-1)] ] #remove " in case the user has entered file names with spaces using the "name" syntax
    #Determine the output filename
    if "--suggest-chunking" in argv :  
        outputFilename = None  #Do not generate an output file if the user just wants to know which chunking to be used
        inputFilenames.append( argv[-1].replace('"' , "")  ) #Check chunking also for the last filename
    else : 
        outputFilename = argv[-1].replace('"' , "")  #Remove " in case the user has entered file names with spaces using the "name" syntax

    #Check whether all packages needed for generating thumbnails are available on the system
    if generate_thumbnail :
        try : 
            from PIL import Image
        except :
            generate_thumbnail = False
            print "PIL not available. Generation of thumbnail images disabled"
    
    #Enable chunking if compression is requested and chunking has been disabled
    if (chunks == None) and (compression is not None) :
        print "WARNING: HDF5 compression is only available with chunking enabled. Do you want to enable chunking? (Y/N)"
        userInput = raw_input()
        if userInput == "Y" or userInput == "y" or userInput=="Yes" or userInput=="yes" or userInput=="YES": 
            chunks=(4,4,2048)
            print "Chunking enabled with (4,4,2048)"
        elif userInput == "N" or userInput == "n" or userInput=="No" or userInput=="no" or userInput=="NO":
            compression=None
            compression_opts=None
            print "Compression disabled"
        else :
            exit()
    
    print "Execute global peak finding (fpg): "+str(executeFPG)
    print "Execute local peak finding (fpl): "+str(executeFPL)
    print "Execute nmf: "+str(executeNMF)
    print "Number of MSI files: "+str(len( inputFilenames ))
    print "Output OMSI file: "+outputFilename

    #Finish and return
    return inputError, inputWarning, outputFilename, inputFilenames




def printHelp():
    """Function used to print the help for this script"""
    
    #Load gloabl variables
    global availableFormats
    global availableRegionOptions
    
    #Print the help explaining the usage of convertToHDF5
    
    print "USAGE: Call \"convertToOMSI [options] imgBaseFile1 imgBaseFile2 ... imgBaseFileN HDF5File\" "
    print ""
    print "This converter script takes the basename (i.e., path+basefilename) of a single"
    print "or multiple MSI files as input and converts them to HDF5. Each MSI file is"
    print "stored as a separate experiment in the output HDF5 file. If an input file"
    print "defines multiple regions, then those regions can either be stored as separate"
    print "datasets of the same experiment and/or merged to a single MSI dataset."
    print "Using the various paramter settings described below, one can define how the"
    print "conversion should be performed, how the data should be stored in HDF5, and"
    print "indicate which analyses should be exectued."
    print ""
    print "===HELPER OPTIONS==="
    print "--suggest-chunking : Iterate over all given input files and suggest a chunking strategy."
    print "                     No data is converted when this option is given, i.e., no name for the"
    print "                     HDF5File should be given, but only input files should be listed."
    print ""
    print "===INPUT DATA OPTIONS==="
    print ""
    print "Default input data options: --format auto --regions split+merge"
    print "--format <option>: Define which file format is used as input. By default the program tries to"
    print "           automatically determine the input format. This option can be used to indicate"
    print "           the format explicitly to in case the auto option fails. Available options are:"
    print "          "+str(availableFormats)
    print "--regions <option>: Some file formats (e.g., bruker) allow multiple regions to be imaged and stored"
    print "           in a single file. This option allows one to specify how these regions should be"
    print "           treated during file conversion. E.g., one may want to store i) each region as a "
    print "           separate dataset in the output file (--regions split), ii) all regions combined "
    print "           in a single dataset (--regions merge), or both (--regions split+merge)"
    print "           Available options are:"
    print "          "+str(availableRegionOptions)
    print ""
    print "===FILE WRITE OPTIONS==="
    print ""
    print "---FILE WRITE OPTIONS: Chunking---"
    print ""
    print "Default HDF5 Chunking options: Enabled by default using --auto-chunking :"
    print "--auto-chunking : Automatically decide which chunking should be used. This option"
    print "                automatically generates two copies of the data, one with a chunking"
    print "                optimized for selection of spectra and another one optimized for "
    print "                selection of ion image slices. All --chunking, --no-chunking, and"
    print "                --optimized-chunking options are ignored if this paramter is given"
    print "--chunking <x y z> : Use chunking when writing the HDF5 file. (DEFAULT, x=4,y=4,z=2048)"
    print "--no-chunking : Disable chunking when writing the HDF5 file. Use in combination with"
    print "                --no-compression since compression depends on chunking and will enable"
    print "                it if compression is used."
    print "--optimized-chunking <x y z> : Use this option to generate additional copies of the data"
    print "                with different chunked data layouts. Generating multiple copies of the"
    print "                data with different chunked data layouts can be help accelerate selective "
    print "                data read opeations. (DEFAULT OFF). We recommend a spectra-aligned chunking"
    print "                for the raw data, e.g., '--chunking 1 1 32768' and an image-aligned chunked"
    print "                secondary copy of the data, e.g., '--optimzied-chunking 20 20 100'."
    print ""
    print "---FILE WRITE OPTIONS: Compression---"
    print "HDF5 Compression: Default ON using (gzip , 4):"
    print "--compression: Enable compression using (gzip,4). NOTE: Compression requires the use of chunking."
    print "--no-compression: Disable the use of compression."
    print ""
    print "===I/O OPTIONS==="
    print "--io <option>: Available options are: "
    print "             i) all : Read the full data in memory and write it at once" 
    print "             ii) spectrum : Read one spectrum at a time and write it to the file. "
    print "             iii) chunk : Read one chunk at a time and write it to the file."
    print ""
    print "---ANALYSIS OPTIONS---"
    print ""
    print "NMF: Default ON: (nc=20, timout=600, niter=2000, tolerance=0.0001, raw=False)"
    print "--nmf : Compute the nmf for all the input data files and store the results in the"
    print "        HDF5 file. NOTE: If global peak-finding (fpg) is performed, then"
    print "        nmf will be performed on the peak-cube, otherwise on the raw data"
    print "--no-nmf: Disable the execution of nmf"
    print "--nmf-nc <nummber>: Number of components to be computed by the NMF. (default nc=20)"
    print "--nmf-timeout <number>: Maximum time in seconds to be used for computing the NMF. (default timeout=600)"
    print "--nmf-niter <number>: Number of iterations (minimum is 2)(default niter=2000)"
    print "--nmf-tolerance <number>: Tolerance value for a relative stopping condition. (default tolerance=0.0001)"
    print "--nmf-raw <number>: Force execution of the NMF on the raw data. By default the results from"
    print "            the global peak finding (--fpg) are used to compute the NMF."
    print ""
    print "Global Peak Finding: Default ON:"
    print "--fpg : Compute the global peak finding for all input data files and save results"
    print "           in the HDF5 file (DEFAULT)"
    print "--no-fpg: Disable the global peak finding"
    print ""
    print "Global Peak Finding: Default OFF:"
    print "--fpl : Compute the local peak finding for all input data files and save results"
    print "        in the HDF5 file"
    print "--no-fpl: Disable the local peak finding (DEFAULT)"
    print ""
    print "---OTHER OPTIONS---"
    print ""
    print "Generate Thunmbnail image: Default ON:"
    print "--thumbnail: Generate thumbnail image for the file based on, in order of avalability:"
    print "             * The frist three components of the NMF"
    print "             * The three most intense peaks from the global peak finding (fpg)"
    print "             * The three most intense peaks in the raw data that are at least 1 percent"
    print "               of the total m/z range apart."
    print "--no-thumbnail: Do not generate a thumbnail image."

if __name__ == "__main__":
    main() 

