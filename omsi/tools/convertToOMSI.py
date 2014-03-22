"""Tool used to convert img files to OpenMSI HDF5 files.

  For usage information execute: python convertToOMSI --help
"""

from omsi.dataformat.bruckerflex_file import bruckerflex_file
from omsi.dataformat.img_file import img_file
from omsi.dataformat.omsi_file import *
from omsi.analysis.multivariate_stats.omsi_nmf import omsi_nmf
from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
from omsi.analysis.findpeaks.omsi_findpeaks_local import omsi_findpeaks_local
import time
import numpy as np
import math
import sys
import os
import warnings
#from sys import exit
try:
    from PIL import Image
    pil_available = True
except:
    try:
        import Image
        pil_available = True
    except:
        pil_available = False
try:
    import urllib2
    #import urllib
except:
    # This is to ensure the script is usable without urllib2 when register to
    # DB is not requested
    pass


####################################################################
#  Variables used for internal storage                             #
####################################################################
""" :param dataList: List of python dictionaries describing specific conversion \
             settings for each conversion task. Each dictionary contains the following keys:
    
             * 'basename' : Name of the file to be converted 
             * 'format' : File format to be used (see available_formats) 
             * 'exp' : Indicate the experiment the dataset should be stored with. Valid values are \
                       
                          * 'new' : Generate a new experiment for the dataset
                          * 'previous' : Use the same experiment as used for the previous dataset 
                          * 1, 2,3...   : Integer value indicating the index of the experiment to be used. 
            * 'region' : Optional key with index of the region to be converted. None if all regions should be merged.
"""
dataset_list = []
omsi_output_file = None  # The openMSI output data file to be used.

####################################################################
#  Define available options for different parameters              ##
####################################################################
# List of available data formats
available_formats = ["img", "bruckerflex", "auto"]
# List defining the different options available for handling regions
available_region_options = ["split", "merge", "split+merge"]
# Available options for the data write. One chunk at a time ('chunk'), one
# spectrum at a time ('spectrum'), or all at one once ('all')
available_io_options = ["chunk", "spectrum", "all"]
available_error_options = ["terminate-and-cleanup",
                           "terminate-only", "continue-on-error"]

####################################################################
#  Define the input parameters used during the conversion         ##
####################################################################
# Define how the data should be written to file, one chunk at a time
# ('chunk'), one spectrum at a time ('spectrum') or all at one once
# ('all')
io_option = "chunk"
format_option = "auto"  # Define which file format reader should be used
region_option = "split+merge"  # Define the region option to be used
auto_chunk = True  # Automatically decide which chunking should be used
chunks = (4, 4, 2048)  # Define the main chunking for the data
# User-defined optional additional chunked copies of the data used to
# optimize data seletion
user_additional_chunks = []
# We use gzip by default because szip is not as generally available and
# LZF is usually only available with h5py.
compression = 'gzip'
# For gzip this a value between 0-9 (inclusive) indicating how aggressive
# the compression should be. 0=fastest speed, 9=best compression ratio.
compression_opts = 4
# Should we only suggest file chunkings and quit without converting any data
suggest_file_chunkings = False

####################################################################
#  Define file add options                                        ##
####################################################################
add_file_to_db = True  # Add the file to the database
# This option is to ensure that when the script is run elsewhere that we
# do not just call add file without actually saving the file at NERSC
check_add_nersc = True
default_db_server_url = "https://openmsi.nersc.gov/"
# Specify the server where the database is registered
db_server_url = default_db_server_url
file_owner = None  # Specify for which owner the file should be registered
# require_login = False #Require login before conversion

####################################################################
#  Define error-handling options                                  ##
####################################################################
error_handling = "terminate-and-cleanup"
recorded_warnings = []

####################################################################
# Define error email notification options                         ##
####################################################################
email_success_recipients = []
email_error_recipients = []


####################################################################
#  Define analysis options and parameter settings                 ##
####################################################################
execute_nmf = True  # Define whether NMF should be performed
execute_fpg = True  # Define whether global peak finding should be executed
execute_fpl = False  # Define whether local peak finding should be executed
generate_thumbnail = True  # Should we generate thumbnail

# Default NMF parameter settings
nmf_num_components = 20  # Number of components for the NMF
nmf_timeout = 600  # Timeout for the NMF
nmf_num_iter = 2000  # Maximum number of iterations for the NMF
nmf_tolerance = 0.0001  # Tolerance level for the nmf
# Should the NMF be computed from the raw data or from the global peak
# finding data
nmf_use_raw_data = False

#TODO Register all print-outs with the email message
def main(argv=None):
    """The main function defining the control flow for the conversion"""
    # Get the global variables
    global omsi_output_file
    global dataset_list
    global format_option
    global region_option
    global add_file_to_db
    global db_server_url
    global file_owner
    global error_handling
    global available_error_options
    global recorded_warnings
    global email_error_recipients
    global email_success_recipients
    #global require_login

    ####################################################################
    #   Determine the settings based on the user input                 #
    ####################################################################
    # Get the user input options
    if argv is None:
        argv = sys.argv
    # Parse the input arguments
    inputError, inputWarning, omsiOutFile, inputFilenames = parse_input_args(argv)
    # Terminate in case an error or warning has occured while processing the
    # user input parameters.
    if inputError:
        emailmsg = "One or more error occurred while parsing the command line inputs."
        send_email(subject="ERROR: Conversion of file failed: " + omsiOutFile,
                   body=emailmsg,
                   email_type='error')
        raise ValueError(emailmsg)
    if inputWarning:
        emailmsg = "Conflicting input parameters found. See WARNINGS above for details."
        send_email(subject="ERROR: Conversion of file failed: " + omsiOutFile,
                   body=emailmsg,
                   email_type='error')
        raise ValueError(emailmsg)

    ####################################################################
    # Generate the list of datasets to be converted                    #
    ####################################################################
    try:
        dataset_list = create_dataset_list(inputFilenames=inputFilenames,
                                           format_type=format_option, data_region_option=region_option)
        print "Number of conversion: " + str(len(dataset_list))
    except:
        emailmsg =  "ERROR: An error occurred during the generation of the input filelist.\n"
        emailmsg +=  "       -- No HDF5 output file has been generated."
        emailmsg +=  "       -- No file has been added to the database."
        emailmsg +=  "       -- Terminating"
        emailmsg += str(sys.exc_info())
        send_email(subject="ERROR: Conversion of file failed: " + omsiOutFile,
                   body=emailmsg,
                   email_type='error')
        print emailmsg
        raise

    ####################################################################
    #  Suggest only chunking for the files if requested                #
    ####################################################################
    if suggest_file_chunkings:
        suggest_chunkings_for_files(dataset_list)
        exit()

    ####################################################################
    #   Create the output HDF5 file if needed                          #
    ####################################################################
    try:
        if omsiOutFile is not None:
            omsi_output_file = omsi_file(omsiOutFile)
    except:
        emailmsg = "Unexpected error creating the output file:", sys.exc_info()[0]
        send_email(subject="ERROR: Conversion of file failed: " + omsiOutFile,
                   body=emailmsg,
                   email_type='error')
        print emailmsg
        raise

    ####################################################################
    # Convert all files                                                #
    ####################################################################
    try:
        #Convert all files and record warnings for reporting purposes
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            convert_files()
            recorded_warnings += w
    except:
        emailmsg =  "ERROR: An error occured during the file conversion."
        print "ERROR: An error occured during the file conversion."
        # Try to close the output file
        try:
            omsi_output_file.close_file()
            emailmsg += "  - Successfully closed output file"
        except:
            emailmsg += "  - Closing of output HDF5 file failed" + str(sys.exc_info())
            pass
        if error_handling == "terminate-and-cleanup":
            emailmsg += "  -The generated HDF5 will not be added to the database."
            print "--The generated HDF5 will not be added to the database."
            add_file_to_db = False
            emailmsg += "  -Attempting to delete the generated HDF5 file."
            print "--Attempting to delete the generated HDF5 file."
            os.remove(omsiOutFile)
            emailmsg += "  -Successfully deleted the generated HDF5 file: " + str(omsiOutFile)
            print "--Successfully deleted the generated HDF5 file: " + str(omsiOutFile)
        if error_handling == "terminate-only" or error_handling == "continue-on-error":
            emailmsg += "  -The generated HDF5 will not be added to the database."
            print "--The generated HDF5 will not be added to the database."
            add_file_to_db = False
            emailmsg += "  -The output HDF5 file (if generate) remains at: " + str(omsiOutFile)
            emailmsg += "  -Output file found: " + str(os.path.exists(omsiOutFile))
            print "--The output HDF5 file (if generate) remains at: " + str(omsiOutFile)
            print "  Output file found: " + str(os.path.exists(omsiOutFile))

        #Add warnings to the email message
        emailmsg += "\n"
        emailmsg += "---------------------------------------------"
        emailmsg += "\n"
        for warn in recorded_warnings:
            emailmsg += warn.message + "\n"

        #Send email notification if needed
        send_email(subject="ERROR: Conversion of file failed: " + omsiOutFile,
                   body=emailmsg,
                   email_type='error')

        # Pass on which-ever error has occurred
        raise

    ####################################################################
    #  Close the HDF5 file and exit                                    #
    ####################################################################
    omsi_output_file.close_file()

    ####################################################################
    #  Register the file with the database                             #
    ####################################################################
    if add_file_to_db:
        status = register_file_with_db(
            filepath=omsiOutFile, db_server=db_server_url, file_owner_name=file_owner)
        print "Registered file with DB: " + str(status)

    ####################################################################
    #  Report all recorded warnings on the command line                #
    ####################################################################
    warningmsg = ""
    if len(recorded_warnings) > 0:
        warningmsg += "WARNING: The following warnings occurred during the conversion process"
        for warn in recorded_warnings:
            warningmsg += warn.message + "\n"
        print warningmsg

    ####################################################################
    #  Send email notification if requested                            #
    ####################################################################
    if len(recorded_warnings) == 0 :
        send_email( recipients=email_success_recipients,
                    subject='Conversion complete: '+str(omsiOutFile)  ,
                    body='Success',
                    email_type='success')
    elif len(email_error_recipients) > 0:
        send_email( recipients=email_success_recipients,
                    subject='Conversion completed with warnings: '+str(omsiOutFile)  ,
                    body=warningmsg,
                    email_type='warning')

    ####################################################################
    #  Exit                                                            #
    ####################################################################
    #exit(0)


def convert_files():
    """Convert all files in the given list of files with their approbriate conversion options"""

    # The following global variables are used here but not changed
    global available_formats
    global available_region_options
    global io_option
    global format_option
    global region_option
    global auto_chunk
    global chunks
    global user_additional_chunks
    global compression
    global compression_opts
    global suggest_file_chunkings
    global execute_nmf
    global execute_fpg
    global execute_fpl
    global generate_thumbnail
    global nmf_num_components
    global nmf_timeout
    global nmf_num_iter
    global nmf_tolerance
    global nmf_use_raw_data
    global dataset_list
    global omsi_output_file
    global error_handling
    global available_error_options


    ####################################################################
    #  Convert the MSI files and compute the requested analyses       ##
    ####################################################################
    # Iterate over all img files
    for i in dataset_list:  # xrange(startIndex,len(argv)-1):

        ####################################################################
        #  Convert the raw data to HDF5                                   ##
        ####################################################################
        basefile = i['basename']  # argv[i].replace('"', "")
        print "Converting: " + basefile
        try:
            print "Selected Region: " + str(i['region'])
        except:
            pass
        if not auto_chunk:
            print "HDF5 chunking: " + str(chunks)
        print "HDF5 compression: " + str(compression) + ", " + str(compression_opts)

        # Open the input file
        try:
            currFormat = i['format']
            print "Input file format: " + str(currFormat)
            if currFormat == "img":
                # hdrFile=basefile+".hdr", t2mFile=basefile+".t2m",
                # imgFile=basefile+".img")
                inputFile = img_file(basename=basefile)
            elif currFormat == "bruckerflex":
                inputFile = bruckerflex_file(spotlist_filename=basefile)
                inputFile.set_region_selection(region_index=i['region'])
            else:
                print "ERROR: The following file will not be converted. File type could not be determined: " + basefile
                print "INFO: If you know the correct file format then try setting appropriate --format option."
                if error_handling == "continue-on-error":
                    pass
                elif error_handling == "terminate-and-cleanup" or error_handling == "terminate-only":
                    raise ValueError("Unrecognized file format.")

            print "In data shape: " + str(inputFile.shape)
        except:
            print "ERROR: Unexpected error opening the input file:", sys.exc_info()[0]
            if error_handling == "continue-on-error":
                print basefile + " failure during conversion. Skipping the file and continue conversion."
                continue
            elif error_handling == "terminate-and-cleanup" or error_handling == "terminate-only":
                raise

        if auto_chunk:

            specChunks, imgChunks, tempChunks = suggest_chunking(
                xsize=inputFile.shape[0],
                ysize=inputFile.shape[1],
                mzsize=inputFile.shape[2],
                dtype=inputFile.data_type,
                print_results=True)
            chunks = specChunks
            additionalChunks = user_additional_chunks + [imgChunks]
            print "Converting data using the following chunking options:"
            print "     - Spectrum chunking: " + str(chunks)
            print "     - Image chunking:    " + str(additionalChunks[0])

        else:
            additionalChunks = user_additional_chunks

        # Get the mz data
        mzdata = inputFile.mz

        # Define the layout for the img data in the HDF5
        # Create a new experiment for the img file using the basefile as
        # identifer string
        if i['exp'] == 'new':  # Create a new experiment for this dataset
            exp = omsi_output_file.create_exp(exp_identifier=basefile)
            # Create an empty sample descrition
            sample = exp.create_sample_info()
            # Create an empty instrument description
            instrument = exp.create_instrument_info(
                instrument_name="undefined", mzdata=mzdata)
        # Store this dataset in combination with the previous one
        elif i['exp'] == 'previous':
            exp = omsi_output_file.get_exp(omsi_output_file.get_num_exp() - 1)
        # Store this dataset in combination with the given experiment index
        else:
            exp = omsi_output_file.get_exp(i['exp'])

        # Allocate space in the HDF5 file for the img data
        dataDataset, mzDataset, dataGroup = exp.create_msidata_full_cube(
            data_shape=inputFile.shape,
            data_type=inputFile.data_type,
            chunks=chunks,
            compression=compression,
            compression_opts=compression_opts)
        mzDataset[:] = mzdata
        data = omsi_file_msidata(
            data_group=dataGroup, preload_mz=False, preload_xy_index=False)
        write_data(inputFile, data, data_io_option=io_option, chunk_shape=chunks)
        omsi_output_file.flush()

        # Generate any additional data copies if requested
        for chunkSpec in additionalChunks:
            print "Generating optimized data copy: " + str(chunkSpec)
            tempData = data.create_optimized_chunking(chunks=chunkSpec,
                                                      compression=compression,
                                                      compression_opts=compression_opts,
                                                      copy_data=False, print_status=True)
            write_data(inputFile, tempData, data_io_option=io_option,
                       chunk_shape=chunkSpec)
            omsi_output_file.flush()

        # Close the input data file and free up any allocated memory
        inputFile.close_file()

        ####################################################################
        #  Execute the requested analyses                                 ##
        ####################################################################
        # Compute the local peak finding and add the peak-cube data to the file
        if execute_fpl:
            print "Executing local peak finding"
            # Execute the peak finding
            fpl = omsi_findpeaks_local(
                nameKey="omsi_findpeaks_local_" + str(time.ctime()))
            fpl.execute(msidata=data, mzdata=mzdata, printStatus=True)
            # Save the peak-finding to file
            ana, anaIndex = exp.create_analysis(fpl)
        # Compute the global peak finding and add the peak-cube data to the
        # file
        if execute_fpg:
            print "Executing global peak finding"
            # Execute the peak finding
            fpg = omsi_findpeaks_global(
                nameKey="omsi_findpeaks_global_" + str(time.ctime()))
            fpg.execute(msidata=data, mzdata=mzdata)
            # Save the peak-finding to file
            ana, anaIndex = exp.create_analysis(fpg)
            fpgPath = ana.get_h5py_analysisgroup().name
        # Compute the nmf analysis if requested
        if execute_nmf:
            print "Executing nmf"
            # Exectue the nmf analysis on the peak-cube or if peak finding was not performed then
            # execute nmf on the raw data data
            nmf = omsi_nmf(nameKey="omsi_nmf_" + str(time.ctime()))
            if execute_fpg and not nmf_use_raw_data:
                print "   Using peak-cube data for NMF"
                fpgData = ana['peak_cube']
                nmf.execute(msidata=fpgData, numComponents=nmf_num_components,
                            timeOut=nmf_timeout, numIter=nmf_num_iter, tolerance=nmf_tolerance)
            else:
                print "   Using raw data for NMF"
                nmf.execute(msidata=data, numComponents=nmf_num_components,
                            timeOut=nmf_timeout, numIter=nmf_num_iter, tolerance=nmf_tolerance)
            # Save the nmf results to file
            ana, anaIndex = exp.create_analysis(nmf)

        ####################################################################
        #  Generate the thumbnail image for the current file              ##
        ####################################################################
        try:
            if generate_thumbnail:
                print "Generating the thumbnail image"
                if execute_nmf:
                    print "   Generating thumbnail from NMF data"
                    # Get the NMF data
                    ho = nmf['ho']
                    numX = ho.shape[0]
                    numY = ho.shape[1]
                    # Generate images for the first three NMF components
                    d1 = np.log(ho[:, :, 0].reshape((numX, numY)) + 1)
                    d1 = d1 / np.max(d1)
                    im1 = Image.fromarray(d1.astype('float') * 255).convert('L')
                    d2 = np.log(ho[:, :, 1].reshape((numX, numY)) + 1)
                    d2 = d2 / np.max(d2)
                    im2 = Image.fromarray(d2.astype('float') * 255).convert('L')
                    d3 = np.log(ho[:, :, 2].reshape((numX, numY)) + 1)
                    d3 = d3 / np.max(d3)
                    im3 = Image.fromarray(d3.astype('float') * 255).convert('L')
                    # Generate thumbnail by merging the three gray-scale images as
                    # an RGB image
                    thumbnail = Image.merge('RGB', (im1, im2, im3))
                    expName = str(exp.get_h5py_experimentgroup().name)
                    expIndex = expName[7:len(expName)]
                    thumbnailFilename = omsi_output_file.hdf_filename + \
                        "_" + expIndex + ".png"
                    thumbnail.save(thumbnailFilename, 'PNG')
                elif execute_fpg:
                    print "    Generating thumbnail from FPG data"
                    # Get the global peak finding data and compute the maximum peak
                    # values for each peak
                    fpgData = fpg['peak_cube'][:]
                    maxPeakValues = fpgData.max(axis=0).max(axis=0)
                    numX = fpgData.shape[0]
                    numY = fpgData.shape[1]
                    # Generate images for the three most intense peaks
                    s = np.argsort(maxPeakValues)
                    d1 = np.log(fpgData[:, :, s[-1]].reshape((numX, numY)) + 1)
                    d1 = d1 / np.max(d1)
                    im1 = Image.fromarray(d1.astype('float') * 255).convert('L')
                    d2 = np.log(fpgData[:, :, s[-2]].reshape((numX, numY)) + 1)
                    d2 = d2 / np.max(d2)
                    im2 = Image.fromarray(d2.astype('float') * 255).convert('L')
                    d3 = np.log(fpgData[:, :, s[-3]].reshape((numX, numY)) + 1)
                    d3 = d3 / np.max(d3)
                    im3 = Image.fromarray(d3.astype('float') * 255).convert('L')
                    # Generate thumbnail by merging the three gray-scale images as
                    # an RGB image
                    thumbnail = Image.merge('RGB', (im1, im2, im3))
                    expName = str(exp.get_h5py_experimentgroup().name)
                    expIndex = expName[7:len(expName)]
                    thumbnailFilename = omsi_output_file.hdf_filename + \
                        "_" + expIndex + ".png"
                    thumbnail.save(thumbnailFilename, 'PNG')
                else:
                    print "Generation of thumbnail from raw data is not yet supported. No thumbnail has been generated."
                    print "Enable --nmf or --fpg in order to generate a thumbnail image."
                    # print "    Generating thumbnail from raw data"
                    # Find three most intense peaks that are at least 1% of the m/z range appart
                    #numX = data.shape[0]
                    #numY = data.shape[1]
                    #numZ = data.shape[2]
                    #minMzStep = numZ / 100
                    #mzImageRange = numZ/4000
                    #maxPeakValues = np.zeros(numZ)
                    #stepping = np.arange(0, numZ-5000, 5000)
                    # for i in stepping:
                        #maxPeakValues[i:(i+5000)] = data[:, :, i:(i+5000)].std(axis=0).std(axis=0).reshape(5000)
                    #maxPeakValues[stepping[-1]:numZ] = data[:, :, stepping[-1]:numZ].std(axis=0).std(axis=0).reshape(numZ - stepping[-1])
                    # print maxPeakValues.shape
                    #s = np.argsort(maxPeakValues)
                    #i1 = s[-1]
                    #i2 = s[-2]
                    # for i in reversed(range(0, len(s))):
                        # if abs(s[i] - i1) > minMzStep:
                            #i2 = s[i]
                            # break
                    #i3 = s[-3]
                    # for i in reversed(range(0, len(s))):
                        # if abs(s[i] - i1) > minMzStep and abs(s[i]-i2)>minMzStep:
                            #i2 = s[i]
                            # break
                    # print minMzStep
                    # print str(i1)+" "+str(i2)+" "+str(i3)
                    # print str(mzdata[i1])+" "+str(mzdata[i2])+" "+str(mzdata[i3])
                    #low = max(0, s[i1]-mzImageRange)
                    #hi = min(numZ, s[i1]+mzImageRange)
                    # print str(low) + " " + str(hi)
                    #d1 = data[:, :,low:hi].max(axis=2).reshape((numX, numY))
                    #d1 = d1 / np.max(d1)
                    #im1 = Image.fromarray(d1.astype('float')*255).convert('L')
                    #low = max(0, s[i2]-mzImageRange)
                    #hi = min(numZ, s[i2]+mzImageRange)
                    # print str(low) + " " + str(hi)
                    #d2 = data[:, :,low:hi].max(axis=2).reshape((numX, numY))
                    #d2 = d2 / np.max(d2)
                    #im2 = Image.fromarray(d2.astype('float')*255).convert('L')
                    #low = max(0, s[i3]-mzImageRange)
                    #hi = min(numZ, s[i3]+mzImageRange)
                    # print str(low) + " " + str(hi)
                    #d3 = data[:, :,low:hi].max(axis=2).reshape((numX, numY))
                    #d3 = d3 / np.max(d3)
                    #im3 = Image.fromarray(d3.astype('float')*255).convert('L')
                    # Generate thumbnail by merging the three gray-scale images as an RGB image
                    #thumbnail = Image.merge('RGB', (im1, im2, im3))
                    #expName = str(exp.get_h5py_experimentgroup().name)
                    #expIndex = expName[7:len(expName)]
                    #thumbnailFilename = omsi_output_file.hdf_filename+"_"+expIndex+".png"
                    #thumbnail.save(thumbnailFilename, 'PNG')
        except ImportError:
            print "ERROR: Thumbnail generation failed. I/O error({0}): {1}".format(e.errno, e.strerror)
            pass
        except:
            print "ERROR: Thumbnail generation failed. Unexpected error:", sys.exc_info()[0]
            pass

    ####################################################################
    #  Generate the XDMF header file for the HDF5 file                ##
    ####################################################################
    omsi_output_file.write_xdmf_header(
        omsi_output_file.get_filename() + ".xdmf")


def check_format(name, format_type):
    """Helper function used to determine the file format that should be used

       :param name: Name of the folder/file that we should read
       :param format_type: String indicating the format-option given by the user. \
                    If the format is not determined (i.e., "auto") then this function \
                    tries to determine the approbriate foramt. Otherwise this option \
                    is returned as is, as the user explicilty said which format should \
                    be used.
       :returns: String indicating the approbriate format. Returns None in case no valid option was found.
    """
    global available_formats
    # Option 1: The user told us the format we should use
    if format_type is not "auto":
        if format_type in available_formats:
            return format_type
        else:
            return None
    # Option 2: We need to determine the format ourselves
    # if os.path.exists(name+".hdr") and os.path.exists(name+".t2m") and os.path.exists(name+".img"):
    #    return "img"
    if img_file.is_img(name):
        return "img"
    elif bruckerflex_file.is_bruckerflex(name):
        return "bruckerflex"

    # elif name.endswith("Spot List.txt"):
    #    return "bruckerflex"
    # elif os.path.isdir(name):
    #    return "bruckerflex"
    else:
        return None


def create_dataset_list(inputFilenames, format_type='auto', data_region_option="split+merge"):
    """Based on the list of inputFilenames, generate the dataset_list, which contains a dictionary describing
       each conversion job

       :param inputFilenames: List of names of files to be converted.
       :param format_type: Define which file-format should be used. Default value is 'auto' indicating the \
                   function should determine for each file the format to be used. See also available_formats parameter.
       :param data_region_option: Define how different regions defined for a file should be handled. E.g., one may \
                    want to split all regions into indiviudal datasets ('split'), merge all regions into a single \
                    dataset ('merge'), or do both ('split+merge'). See also the available_region_options parameter \
                    for details. By default the function will do 'split+merge'.

       :returns: List of dictionaries describing the various conversion jobs
    """
    global available_error_options
    global error_handling

    re_dataset_list = []
    for i in inputFilenames:
        currDS = {}
        # Remove " in case the user has entered file names with spaces using
        # the "name" syntax
        currDS['basename'] = i.rstrip('"').rstrip("'").lstrip('"').lstrip("'")
        currDS['format'] = check_format(name=currDS['basename'],
                                        format_type=format_type)
        currDS['exp'] = 'new'
        if currDS['format'] == 'bruckerflex':
            if data_region_option == "merge" or data_region_option == "split+merge":
                currDS['region'] = None
                re_dataset_list.append(currDS)

            if data_region_option == 'split' or data_region_option == 'split+merge':
                try:
                    tempFile = bruckerflex_file(
                        spotlist_filename=currDS['basename'], readall=False)
                    for ri in xrange(0, tempFile.get_number_of_regions()):
                        nDS = {'basename': currDS['basename'],
                               'format': currDS['format'],
                               'region': ri,
                               'exp': 'previous'}
                        re_dataset_list.append(nDS)
                except:
                    print "ERROR: Unexpected error opening the input file:", sys.exc_info()[0]
                    print currDS['basename'] + " failure during converion. Skipping the file and continue conversion."
                    if error_handling == "continue-on-error":
                        continue
                    elif error_handling == "terminate-and-cleanup" or error_handling == "terminate-only":
                        raise
                    else:
                        raise

            # This is just to make sure that we have checked everything
            if data_region_option not in available_region_options:
                print "WARNING: Undefined region option. Using the default option split+merge"
        else:
            re_dataset_list.append(currDS)

    # Return the list of conversion jobs
    return re_dataset_list


def suggest_chunkings_for_files(in_dataset_list):
    """Helper function used to suggest food chunking strategies for a given set of files.

       :param in_dataset_list: Python list of dictionaries describing the settings to be used for the file conversion

       :returns: This function simply prints results to standard-out but does not return anything.
    """
    for i in in_dataset_list:

        basefile = i["basename"]
        try:
            print "Suggested Chunkings: " + basefile
            currFormat = i["format"]
            if currFormat is "img":
                img_file()
                inputFile = img_file(
                    hdr_filename=basefile + ".hdr", t2m_filename=basefile + ".t2m", img_filename=basefile + ".img")
            elif currFormat is "bruckerflex":
                inputFile = bruckerflex_file(
                    spotlist_filename=basefile, readall=False)
                inputFile.set_region_selection(i["region"])
            else:
                print "WARNING: Type of file could not be determined for: " + basefile
                continue
            print "In data shape: " + str(inputFile.shape)
            suggest_chunking(xsize=inputFile.shape[0], ysize=inputFile.shape[
                             1], mzsize=inputFile.shape[2], dtype=inputFile.data_type, print_results=True)
        except:
            print "Error while trying to generate chunking suggestion for " + basefile


def suggest_chunking(xsize, ysize, mzsize, dtype, print_results=False):
    """Helper function used to suggest god chunking strategies for a given data cube

       :param xsize: Size of the dataset in x.
       :param ysize: Size o the dataset in y.
       :param mzsize: Size of the dataset in mz.
       :param print_results: Print the results to the console.

       :returns: Three tupes:

        * ``spectrumChunk`` : The chunking to be used to optimize selection of spectra.
        * ``sliceChunk`` : The chunking to be used to optimize selection of image slices.
        * ``balancedChunk`` : The chunking that would provide a good balance in performance for \
                            different selection strategies.
    """
    # Make sure that all sizes are treated as 64bit int to avoid errors due to
    # artificial cutoffs.
    xsize = int(xsize)
    ysize = int(ysize)
    mzsize = int(mzsize)
    imageSize = xsize * ysize

    # Variable settings
    numBytes = np.dtype(dtype).itemsize  # Number of bytes per data value
    # Nober of values that fit into a 64Kb chunk
    suggestedNumValues = (1024 * 64) / numBytes
    # How much larger than the suggestNumValues should we allow a chunk to
    # become
    maxFactor = 4

    # Define the chunking that would be good for spectra
    spectrumXChunk = 1
    spectrumYChunk = 1
    factor1 = math.floor(mzsize / float(suggestedNumValues))
    spectrumMzChunk1 = int(math.ceil(mzsize / float(factor1)))
    factor1 = math.ceil(mzsize / float(spectrumMzChunk1))
    overhead1 = ((factor1 * spectrumMzChunk1) - mzsize) * imageSize * numBytes
    # overhead1 = ((spectrumMzChunk1 *  math.ceil (float(mzsize) / float(spectrumMzChunk1))) - \
    #              (spectrumMzChunk1 *  math.floor(float(mzsize) / float(spectrumMzChunk1)))) * \
    #            imageSize * numBytes
    #overhead1 = (spectrumMzChunk1 - math.ceil(float(mzsize) % float(spectrumMzChunk1))) * imageSize * numBytes
    factor2 = math.ceil(mzsize / float(suggestedNumValues))
    spectrumMzChunk2 = int(math.ceil(mzsize / float(factor2)))
    factor2 = math.ceil(mzsize / float(spectrumMzChunk2))
    overhead2 = ((factor2 * spectrumMzChunk2) - mzsize) * imageSize * numBytes
    #overhead2 = (spectrumMzChunk2 - math.ceil(float(mzsize) % float(spectrumMzChunk2))) * imageSize *numBytes
    # overhead2 = ((spectrumMzChunk2 *  math.ceil (float(mzsize) / float(spectrumMzChunk2))) - \
    #              (spectrumMzChunk2 *  math.floor(float(mzsize) / float(spectrumMzChunk2)))) * \
    #            imageSize * numBytes
    if overhead1 < overhead2:
        spectrumMzChunk = spectrumMzChunk1
        spectrumChunkOverhead = overhead1
    else:
        spectrumMzChunk = spectrumMzChunk2
        spectrumChunkOverhead = overhead2
    spectrumChunk = (spectrumXChunk, spectrumYChunk, spectrumMzChunk)
    if print_results:

        print "Spectrum selection chunking: " + str(spectrumChunk)
        print "     - Ideal for selection of full spectra."
        print "     - Overhead: " + str(spectrumChunkOverhead) + \
              " Byte (" + str(int(math.ceil(spectrumChunkOverhead / (1024. * 1024.)))) + " MB)"

    # Define a chunking that would be good for images
    sliceChunkX = xsize
    sliceChunkY = ysize
    sliceChunkMz = 1
    chunkSize = sliceChunkX * sliceChunkY * sliceChunkMz
    if math.ceil(float(chunkSize) / float(suggestedNumValues)) > maxFactor:
        sliceChunkX = int(math.ceil(xsize / 2.))
        sliceChunkY = int(math.ceil(ysize / 2.))
    sliceChunk = (sliceChunkX, sliceChunkY, sliceChunkMz)
    # math.ceil(float(xsize) % float(sliceChunkX))
    sliceChunkOverheadX = (
        sliceChunkX * math.ceil(float(xsize) / float(sliceChunkX))) - xsize
    # math.ceil(float(ysize) % float(sliceChunkY))
    sliceChunkOverheadY = (
        sliceChunkY * math.ceil(float(ysize) / float(sliceChunkY))) - ysize
    sliceChunkOverhead = ((sliceChunkOverheadX * ysize) + (sliceChunkOverheadY * xsize)
                          - (sliceChunkOverheadX * sliceChunkOverheadY)) * mzsize * numBytes
    if print_results:

        print "Slice selection chunking: " + str(sliceChunk)
        print "     - Ideal for selection of full image slices."
        print "     - Overhead: " + str(sliceChunkOverhead) + \
              " Byte (" + str(int(math.ceil(sliceChunkOverhead / (1024. * 1024.)))) + " MB)"

    # Define a balanced chunkgin
    balancedChunkX = 4
    balancedChunkY = 4
    balancedChunkMz = 2048
    balancedChunk = (balancedChunkX, balancedChunkY, balancedChunkMz)
    if print_results:

        print "Balanced chunking: " + str(balancedChunk)
        print "     - This chunking tries to compromise between selection of slices and spectra."

    return spectrumChunk, sliceChunk, balancedChunk


def write_data(inputFile, data, data_io_option="spectrum", chunk_shape=None):
    """Helper function used to implement different data write options.

        :param inputFile: The input img data file
        :param data: The output dataset (either an h5py dataset or omsi_file_msidata object.
        :param data_io_option: String indicating the data write method to be used. One of:

            * ``spectrum``: Write the data one spectrum at a time
            * ``all`` : Write the complete dataset at once.
            * ``chunk`` : Write the data one chunk at a time.

        :param chunk_shape: The chunking used by the data. Needed to decide how the data should \
                            be written when a chunk-aligned write is requested.

    """
    if data_io_option == "spectrum" or (data_io_option == "chunk" and (chunk_shape is None)):
        for xindex in xrange(0, inputFile.shape[0]):
            sys.stdout.write(
                "[" + str(int(100. * float(xindex) / float(inputFile.shape[0]))) + "%]" + "\r")
            sys.stdout.flush()
            for yindex in xrange(0, inputFile.shape[1]):
                # Save the spectrum to the hdf5 file
                data[xindex, yindex, :] = inputFile[xindex, yindex, :]
    elif data_io_option == "all":
        data[:] = inputFile[:]
    elif data_io_option == "chunk":
        xdim = inputFile.shape[0]
        ydim = inputFile.shape[1]
        zdim = inputFile.shape[2]
        numChunksX = int(math.ceil(float(xdim) / float(chunk_shape[0])))
        numChunksY = int(math.ceil(float(ydim) / float(chunk_shape[1])))
        numChunksZ = int(math.ceil(float(zdim) / float(chunk_shape[2])))
        numChunks = numChunksX * numChunksY * numChunksZ
        itertest = 0
        for xChunkIndex in xrange(0, numChunksX):
            xstart = xChunkIndex * chunk_shape[0]
            xend = min(xstart + chunk_shape[0], xdim)
            for yChunkIndex in xrange(0, numChunksY):
                ystart = yChunkIndex * chunk_shape[1]
                yend = min(ystart + chunk_shape[1], ydim)
                for zChunkIndex in xrange(0, numChunksZ):
                    zstart = zChunkIndex * chunk_shape[2]
                    zend = min(zstart + chunk_shape[2], zdim)
                    # print "Write : "+str(xstart)+" "+str(xend)+"
                    # "+str(ystart)+" "+str(yend)+" "+str(zstart)+" "+str(zend)
                    data[xstart:xend, ystart:yend, zstart:zend] = inputFile[
                        xstart:xend, ystart:yend, zstart:zend]
                    itertest += 1
                    sys.stdout.write(
                        "[" + str(int(100. * float(itertest) / float(numChunks))) + "%]" + "\r")
                    sys.stdout.flush()

"""
def loginUser(requestPassword=False):

    global db_server_url 
    import getpass
    #Setup the user login mechanism in urllib2
    if db_server_url.startswith("http:"):
        db_server_url = "https:" + db_server_url.lstrip("http:")
    password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
    username = raw_input("Login: >> ") #Enter username
    password = getpass.getpass()       #Enter password
    password_mgr.add_password(None, db_server_url, username, password)
    handler = urllib2.HTTPBasicAuthHandler(password_mgr)
    opener = urllib2.build_opener(handler)
    urllib2.install_opener(opener)

    #Login the user to the website
    login_data = urllib.urlencode({
            'username' : username,
            'password' : password,
        })
        
    response = urllib2.urlopen(db_server_url, login_data)
    print response
"""


def register_file_with_db(filepath, db_server, file_owner_name):
    """ Function used to register a given file with the database

        :param filepath: Path of the file to be added to the database
        :param db_server: The database server url
        :param file_owner_name: The owner to be used, or None if the owner should be determined based on the file URL.

        :returns: Boolean indicating whether the operation was successful

    """
    global default_db_server_url
    global check_add_nersc

    # Check if the
    if db_server == default_db_server_url and check_add_nersc:
        if (not "/project/projectdirs/openmsi" in filepath) and (not filepath.startswith("/data/openmsi/")):
            print "WARNING: Attempt to add a file to openmsi.nersc.gov that is not in a default location."
            print "Do you want to add the file? (Y/N):"
            numTrys = 3
            for i in range(numTrys):
                userInput = raw_input()
                if userInput == "Y" or userInput == "y" or userInput == "Yes" or \
                        userInput == "yes" or userInput == "YES":
                    break
                elif userInput == "N" or userInput == "n" or userInput == "No" or \
                        userInput == "no" or userInput == "NO":
                    return False
                else:
                    if i == (numTrys - 1):
                        print "Aborting adding file to the database."
                        return False
                    print "Unrecognized repsonse. Do you want to add the file? (Y/N): "

    # If we are at NERSC then set the NERSC Apache permissions
    if 'nersc.gov' in db_server:
        set_apache_acl(filepath)

    # Determine the owner
    currOwner = file_owner_name
    if not currOwner:
        currOwner = os.path.dirname(filepath).split("/")[-1]
    if not currOwner:
        print "ERROR: File could not be added to DB. Owner could not be determined."
        return False

    # Construct the db add-file url
    addFileURL = os.path.join(db_server, "openmsi/resources/addfile")
    addFileURL = addFileURL + "?file=" + \
        os.path.abspath(filepath) + "&owner=" + currOwner

    # Make the url request
    try:
        print "Registering file with DB: " + addFileURL
        urlResponse = urllib2.urlopen(url=addFileURL)
        if urlResponse.code == 200:
            return True
    except urllib2.HTTPError as requestError:
        print "ERROR: File could not be added to DB:"
        print "      Error-code:" + str(requestError.code)
        print "      Error info:" + str(requestError.read())

    return False


def set_apache_acl(filepath):
    """Helper function used to set acl permissions for apache to make the given file accesible
       to Apache at NERSC. This necessary to make the file readable for adding it to the
       database.
    """
    print "Setting NERSC ACL permssiosn for Apache"
    #command = "setfacl -R -m u:apache:rwx "+filepath
    command = "setfacl -R -m u:48:rwx " + filepath
    os.system(command)


def send_email(subject, body, sender='convert@openmsi.nersc.gov', email_type='success'):
    """Send email notification to users.

       :param subject: Subject line of the email
       :param body: Body text of the email.
       :param sender: The originating email address
       :param email_type: One of 'success, 'error', 'warning'. Error messages are sent
                 tp email_error_recipients, success messages to email_success_recipients and
                 warning messages are sent to both lists.

    """

    global email_success_recipients
    global email_error_recipients

    #Define the list of recipients
    if email_type == 'success':
        recipients = email_success_recipients
    elif email_type == 'error':
        recipients = email_error_recipients
    else:
        recipients = email_error_recipients + email_success_recipients
    #Remove duplicates from the list of recipients
    recipients = list(set(recipients))
    #Check if we have any recipients
    if len(recipients) == 0:
        return

    from smtplib import SMTP
    from email.MIMEText import MIMEText
    from email.Header import Header
    from email.Utils import parseaddr, formataddr

    header_charset = 'ISO-8859-1'
    for body_charset in 'US-ASCII', 'ISO-8859-1', 'UTF-8':
        try:
            body.encode(body_charset)
        except UnicodeError:
            pass
        else:
            break

    # Define the sender and recipients
    sender_name, sender_addr = parseaddr(sender)
    sender_name = str(Header(unicode(sender_name), header_charset))
    sender_addr = sender_addr.encode('ascii')

    tostr = ""
    for ri in range(len(recipients)):
        rec = recipients[ri]
        recname, recaddr = parseaddr(rec)
        recname = str(Header(unicode(recname), header_charset))
        recaddr = recaddr.encode('ascii')
        tostr += formataddr((recname, recaddr))
        if ri < (len(recipients)-1):
            tostr += ", "

    # Construct the message
    msg = MIMEText(body.encode(body_charset), 'plain', body_charset)
    msg['From'] = formataddr((sender_name, sender_addr))
    msg['To'] = tostr
    msg['Subject'] = Header(unicode(subject), header_charset)

    # Send the message using sendmail
    try:
        smtp = SMTP("localhost")
        smtp.sendmail(sender, recipients, msg.as_string())
        smtp.quit()
    except:
        Warnings.warn('Email could not be sent' + str(sys.exc_info()))



def parse_input_args(argv):
    """Process input parameters and define the script settings.

       :param argv: The list of input arguments

       :returns: This function returns the following four values:

            * 'inputError' : Boolean indicating whether an error has occured during the processing of the inputs
            * 'inputWarning' : Boolean indicating whether a warning has been raised during the processing of the inputs
            * 'outputFilename' : Name for the output HDF5 file
            * 'inputFilenames' : List of strings indicating the list of input filenames
    """
    # Get the global variables
    global available_formats  # Used but will not be changed
    global available_region_options  # Used but will not be changed
    global available_error_options  # Used but will not be changed
    global io_option
    global format_option
    global region_option
    global auto_chunk
    global chunks
    global user_additional_chunks
    global compression
    global compression_opts
    global suggest_file_chunkings
    global execute_nmf
    global execute_fpg
    global execute_fpl
    global generate_thumbnail
    global nmf_num_components
    global nmf_timeout
    global nmf_num_iter
    global nmf_tolerance
    global nmf_use_raw_data
    global add_file_to_db
    global db_server_url
    global file_owner
    global error_handling
    global check_add_nersc
    global email_error_recipients
    global email_success_recipients
    #global require_login

    # Initalize the output values to be returned
    # Parameter used to track whether an error has occured while parsing the
    # input parameters
    inputError = False
    # Parameter used to indicate conflicting input parameter options
    inputWarning = False
    outputFilename = None  # The output filename
    inputFilenames = []  # The list of input filenames
    # Basic sanity check
    if len(argv) < 3:
        lastArg = argv[-1]
        if lastArg == "--help" or \
           lastArg == "--h" or \
           lastArg == "-help" or \
           lastArg == "-h":
            print_help()
            exit(0)
        else:
            inputError = True
            return inputError, inputWarning, outputFilename, inputFilenames

    # Using the inputError and inputWarning paramter we try to keep parsing the input as long as possible to
    # make sure we can inform the user of as many problems as possible
    startIndex = 1
    i = startIndex
    while i < (len(argv) - 1):
        i = startIndex
        currentArg = argv[i]
        if currentArg == "--no-nmf":
            startIndex += 1
            execute_nmf = False
            print "Disable NMF"
        elif currentArg == "--nmf":
            startIndex += 1
            execute_nmf = True
            print "Enable NMF"
        elif currentArg == "--nmf-nc":
            startIndex += 2
            nmf_num_components = int(argv[i + 1])
            print "Set nmf-nc=" + str(nmf_num_components)
        elif currentArg == "--nmf-timeout":
            startIndex += 2
            nmf_timeout = int(argv[i + 1])
            print "Set nmf_timeout=" + str(nmf_timeout)
        elif currentArg == "--nmf-niter":
            startIndex += 2
            nmf_num_iter = int(argv[i + 1])
            if nmf_num_iter < 2:
                inputWarning = True
                print "WARNING: --nfm-niter must be 2 or larger"
            print "Set nmf-niter=" + str(nmf_num_iter)
        elif currentArg == "--nmf-tolerance":
            startIndex += 2
            nmf_tolerance = float(argv[i + 1])
            print "Set nmf-tolerance=" + str(nmf_tolerance)
        elif currentArg == "--nmf-raw":
            startIndex += 1
            nmf_use_raw_data = True
            print "Set nmf-raw=" + str(nmf_use_raw_data)
        elif currentArg == "--fpg":
            startIndex += 1
            execute_fpg = True
            print "Enable find peaks global"
        elif currentArg == "--no-fpg":
            startIndex += 1
            execute_fpg = False
            print "Disable find peaks global"
            if "--fpg" in argv:
                inputWarning = True
                print "WARNING: --no-fpg and --fpg options are conflicting."
        elif currentArg == "--fpl":
            startIndex += 1
            execute_fpl = True
            print "Enable find peaks local"
        elif currentArg == "--no-fpl":
            startIndex += 1
            execute_fpl = False
            print "Disable find peaks local"
            if "--fpl" in argv:
                inputWarning = True
                print "WARNING: --no-fpl and --fpl options are conflicting."
        elif currentArg == "--auto-chunking":
            startIndex += 1
            auto_chunk = True
            if "--chunking" in argv:
                inputWarning = True
                print "WARNING: --chunking options will be ignored due to the use of --auto-chunking"
            if "--no-chunking" in argv:
                inputWarning = True
                print "WARNING: --no-chunking and --auto-chunking options are conflicting"
            if "--optimized-chunking" in argv:
                inputWarning = True
                print "WARNING: --optimized-chunking and --auto-chunking options are conflicting"
        elif currentArg == "--chunking":
            startIndex += 4
            try:
                chunks = (int(argv[i + 1]), int(argv[i + 2]), int(argv[i + 3]))
            except:
                print "An error accured while parsing the --chunking command. " + \
                      "Something may be wrong with the indicated chunk sizes for x y z."
                inputError = True
            if "--auto-chunking" not in argv:
                auto_chunk = False
            print "Enable chunking: " + str(chunks)
        elif currentArg == "--no-chunking":
            startIndex += 1
            chunks = None
            auto_chunk = False
            print "Disable chunking"
            if "--auto-chunking" in argv or "--chunking" in argv:
                inputWarning = True
                print "WARNGING: --no-chunking option is conflicting with another chunking option"
        elif currentArg == "--optimized-chunking":
            startIndex += 4
            try:
                user_additional_chunks.append(
                    (int(argv[i + 1]), int(argv[i + 2]), int(argv[i + 3])))
            except:
                print "An error accured while parsing the --optimized-chunking command. " + \
                      "Something may be wrong with the indicated chunk sizes for x y z."
                inputError = True
        elif currentArg == "--compression":
            startIndex += 1
            # This is already the default
            print "Enable compression"
        elif currentArg == "--no-compression":
            startIndex += 1
            compression = None
            compression_opts = None
            print "Disable compression"
            if "--compression" in argv:
                inputWarning = True
                print "WARNING: --no-compression and --compression options are conflicting."
        elif currentArg == "--io":
            startIndex += 2
            try:
                io_option = str(argv[i + 1])
                if io_option not in available_io_options:
                    raise ValueError("Invalid io option")
            except:
                print "An error accured while parsing the --io command. Something may be wrong with the io-type."
                inputError = True
        elif currentArg == "--thumbnail":
            startIndex += 1
            generate_thumbnail = True
            print "Enable thumbnail"
        elif currentArg == "--no-thumbnail":
            startIndex += 1
            generate_thumbnail = False
            print "Disable thumbnail"
            if "--thumbnail" in argv:
                inputWarning = True
                print "WARNING: --no-thumbnail and --no-thumbnail options are conflicting."
        elif currentArg == "--help" or currentArg == "--h" or currentArg == "-help" or currentArg == "-h":
            print_help()
            exit(0)
        elif currentArg == "--suggest-chunking":
            startIndex += 1
            suggest_file_chunkings = True
        elif currentArg == "--format":
            startIndex += 2
            format_option = str(argv[i + 1])
            if format_option not in available_formats:
                print "ERROR: Indicated --format option " + format_option + " not supported. Available options are:"
                print "     " + str(available_formats)
                inputError = True
        elif currentArg == "--regions":
            startIndex += 2
            region_option = str(argv[i + 1])
            if region_option not in available_region_options:
                print "ERROR: Indicated --regions option " + region_option + " not supported. Available options are:"
                print "       " + str(available_region_options)
                inputError = True
        elif currentArg == "--add-to-db":
            startIndex += 1
            add_file_to_db = True
            check_add_nersc = False
        elif currentArg == "--no-add-to-db":
            startIndex += 1
            add_file_to_db = False
        elif currentArg == "--db-server":
            startIndex += 2
            db_server_url = str(argv[i + 1])
            check_add_nersc = False
        elif currentArg == "--owner":
            startIndex += 2
            file_owner = str(argv[i + 1])
        elif currentArg == "--email":
            #Consume all email addresses that follow
            for ni in range((i+1),len(argv)):
                if not argv[ni].startswith("--"):
                    email_success_recipients.append(str(argv[ni]))
                    email_error_recipients.append(str(argv[ni]))
                    startIndex = ni+1
        elif currentArg == "--email-success":
            #Consume all email addresses that follow
            for ni in range((i+1),len(argv)):
                if not argv[ni].startswith("--"):
                    email_success_recipients.append(str(argv[ni]))
                    startIndex = ni+1
        elif currentArg == "--email-error":
            #Consume all email addresses that follow
            for ni in range((i+1),len(argv)):
                if not argv[ni].startswith("--"):
                    email_error_recipients.append(str(argv[ni]))
                    startIndex = ni+1
        elif currentArg == "--error-handling":
            startIndex += 2
            errorOption = str(argv[i + 1])
            if errorOption in available_error_options:
                error_handling = errorOption
            else:
                print "ERROR: Indicated --error-handling option " + errorOption + " not supported. Available options:"
                print "     " + str(available_error_options)
                inputError = True
        elif currentArg.startswith("--"):
            startIndex += 1
            inputError = True
            print "Unrecognized input option: " + currentArg
        else:
            # print "End of input parameters"
            break

    # Determine the list of input filenames
    # remove " in case the user has entered file names with spaces using the
    # "name" syntax
    inputFilenames = [name.replace('"', "")
                      for name in argv[startIndex:(len(argv) - 1)]]
    # Determine the output filename
    if "--suggest-chunking" in argv:
        # Do not generate an output file if the user just wants to know which
        # chunking to be used
        outputFilename = None
        # Check chunking also for the last filename
        inputFilenames.append(argv[-1].replace('"', ""))
    else:
        # Remove " in case the user has entered file names with spaces using
        # the "name" syntax
        outputFilename = argv[-1].replace('"', "")

    # Check whether all packages needed for generating thumbnails are
    # available on the system
    if generate_thumbnail:
        if not pil_available:
            generate_thumbnail = False
            print "PIL not available. Generation of thumbnail images disabled"

    # Enable chunking if compression is requested and chunking has been
    # disabled
    if (chunks is None) and (compression is not None):
        print "WARNING: HDF5 compression is only available with chunking enabled. Do you want to enable chunking? (Y/N)"
        userInput = raw_input()
        if userInput == "Y" or userInput == "y" or userInput == "Yes" or userInput == "yes" or userInput == "YES":
            chunks = (4, 4, 2048)
            print "Chunking enabled with (4,4, 2048)"
        elif userInput == "N" or userInput == "n" or userInput == "No" or userInput == "no" or userInput == "NO":
            compression = None
            compression_opts = None
            print "Compression disabled"
        else:
            exit()

    print "Execute global peak finding (fpg): " + str(execute_fpg)
    print "Execute local peak finding (fpl): " + str(execute_fpl)
    print "Execute nmf: " + str(execute_nmf)
    print "Number of MSI files: " + str(len(inputFilenames))
    print "Output OMSI file: " + outputFilename

    # Finish and return
    return inputError, inputWarning, outputFilename, inputFilenames


def print_help():
    """Function used to print the help for this script"""

    # Load gloabl variables
    #global available_formats
    #global available_region_options

    # Print the help explaining the usage of convertToHDF5

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
    print "===HELPER OPTIONS=== "
    print "--suggest-chunking : Iterate over all given input files and suggest a chunking strategy."
    print "                     No data is converted when this option is given, i.e., no name for the"
    print "                     HDF5File should be given, but only input files should be listed."
    print ""
    print "===ERROR HANDLING OPTIONS=== "
    print "--error-handling <options>: Define how errors should be handled. Options are:"
    print "                   i)   terminate-and-cleanup (default) : Terminate the conversion, delete the"
    print "                           the HDF5 file and do not add the file to the database."
    print "                   ii)  terminate-only, : Leave the generated HDF5 output file in place but  do not"
    print "                            add the file to the database."
    print "                   iii) continue-on-error: Ignore errors if possible and continue, even if this"
    print "                            means that some data may be missing from the output."
    print "--email <email1 email2 ...>: Send notification in case of both error or success to the given email address."
    print "--email-success <email1 email2 ...>>: Send notification in case of success to the given email address."
    print "--email-error <email1 email2 ...>>: Send notification in case of error to the given email address."
    print ""
    print "===INPUT DATA OPTIONS=== "
    print ""
    print "Default input data options: --format auto --regions split+merge"
    print "--format <option>: Define which file format is used as input. By default the program tries to"
    print "           automatically determine the input format. This option can be used to indicate"
    print "           the format explicitly to in case the auto option fails. Available options are:"
    print "          " + str(available_formats)
    print "--regions <option>: Some file formats (e.g., brucker) allow multiple regions to be imaged and stored"
    print "           in a single file. This option allows one to specify how these regions should be"
    print "           treated during file conversion. E.g., one may want to store i) each region as a "
    print "           separate dataset in the output file (--regions split), ii) all regions combined "
    print "           in a single dataset (--regions merge), or both (--regions split+merge)"
    print "           Available options are:"
    print "          " + str(available_region_options)
    print ""
    print "===FILE WRITE OPTIONS=== "
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
    print "HDF5 Compression: Default ON using (gzip, 4):"
    print "--compression: Enable compression using (gzip,4). NOTE: Compression requires the use of chunking."
    print "--no-compression: Disable the use of compression."
    print ""
    print "===I/O OPTIONS=== "
    print "--io <option>: Available options are: " + str(available_io_options)
    print "             i) all : Read the full data in memory and write it at once"
    print "             ii) spectrum : Read one spectrum at a time and write it to the file. "
    print "             iii) chunk : Read one chunk at a time and write it to the file."
    print ""
    print "===DATABSE OPTIONS=== "
    print ""
    print "These options control whether the generated output file should be added to a server database"
    print "to manage web file access permissions"
    print "Default options are: --add-to-db --db-server http://openmsi.nersc.gov"
    print "--add-to-db : Add the output HDF5 file to the database."
    print "--no-add-to-db : Disable adding the file to the database."
    print "--db-server : Specify the online server where the file should be registers. Default is"
    print "              http://openmsi.nersc.gov "
    print "--owner : Name of the user that should be assigned as owner. By default the owner is"
    print "          determined automatically based on the file path."
    # print "--login : If set, then the script will ask for login information
    # at the beginning of the script."
    print ""
    print "===ANALYSIS OPTIONS=== "
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
