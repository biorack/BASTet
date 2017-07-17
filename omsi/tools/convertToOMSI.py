"""Tool used to convert img files to OpenMSI HDF5 files.

  For usage information execute: python convertToOMSI --help

"""

# TODO Add write method to optimize the generation of optimized data chunkings by reading many spectra
# TODO Switch job submission over to use Edison
# TODO Test file converter and memory usage
# TODO Make NMF/FPG/FPL etc. seperate workflows
# TODO --io option is not used right now after the changes to convert things iteratively

# TODO ConvertSettings.file_user Get the user from the file-path rather then from the command line
# TODO Register all print-outs with the email message (use logging instead of print)
# TODO add ability have --format option for each dataset to be converted

import sys
import getpass

from omsi.dataformat import *
from omsi.dataformat import omsi_file
from omsi.dataformat import file_reader_base
from omsi.analysis.multivariate_stats.omsi_nmf import omsi_nmf
from omsi.analysis.findpeaks.omsi_findpeaks_global import omsi_findpeaks_global
from omsi.analysis.findpeaks.omsi_findpeaks_local import omsi_findpeaks_local
from omsi.analysis.msi_filtering.omsi_tic_norm import omsi_tic_norm
from omsi.shared.omsi_web_helper import WebHelper
from omsi.shared.omsi_web_helper import UserInput
from omsi.datastructures.dependency_data import dependency_dict
from omsi.datastructures.metadata.metadata_data import metadata_value
from omsi.shared.log import log_helper
import warnings
import os
import numpy as np
import math
import h5py
import time

# Imports for thumbnail image rendering
try:
    from PIL import Image
    pil_available = True
except ImportError:
    try:
        import Image
        pil_available = True
    except ImportError:
        pil_available = False
import urllib2


####################################################################
####################################################################
#  The main function of the command-line tool                      #
####################################################################
####################################################################
def main(argv=None):
    """The main function defining the control flow for the conversion"""

    ####################################################################
    #   Determine the settings based on the user input                 #
    ####################################################################
    # Get the user input options
    if argv is None:
        argv = sys.argv
    # Parse the input arguments
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        input_error, omsi_outfile, input_filenames = ConvertSettings.parse_input_args(argv)
        input_warning = len(w) > 0  # Check if warning occurred
        ConvertSettings.recorded_warnings += w

    # Terminate in case an error or warning has occurred while processing the
    # user input parameters.
    if input_error:
        emailmsg = "One or more errors occurred while parsing the command line inputs."
        for warn in ConvertSettings.recorded_warnings:
            emailmsg += unicode(warn.message) + u"\n"
        WebHelper.send_email(subject="ERROR: Conversion of file failed: " + omsi_outfile,
                             body=emailmsg,
                             email_type='error',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
        if ConvertSettings.job_id is not None:
            try:
                WebHelper.update_job_status(filepath=omsi_outfile,
                                            db_server=ConvertSettings.db_server_url,
                                            jobid=ConvertSettings.job_id,
                                            status='error')
            except ValueError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
        raise ValueError(emailmsg)
    if input_warning:
        emailmsg = "WARNINGS occurred while parsing command line inputs (e.g., conflicting options). \n"
        emailmsg += "See WARNING messages below for details:\n"
        for warn in ConvertSettings.recorded_warnings:
            emailmsg += unicode(warn.message) + u"\n"
        WebHelper.send_email(subject="ERROR: Conversion of file failed: " + omsi_outfile,
                             body=emailmsg,
                             email_type='error',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
        log_helper.warning(__name__, emailmsg)
        if ConvertSettings.job_id is not None:
            try:
                WebHelper.update_job_status(filepath=omsi_outfile,
                                            db_server=ConvertSettings.db_server_url,
                                            jobid=ConvertSettings.job_id,
                                            status='error')
            except ValueError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
        raise ValueError(emailmsg)

    # Notify the system that the job is running
    if ConvertSettings.job_id is not None:
        try:
            WebHelper.update_job_status(filepath=omsi_outfile,
                                        db_server=ConvertSettings.db_server_url,
                                        jobid=ConvertSettings.job_id,
                                        status='running')
        except urllib2.URLError as e:
            log_helper.info(__name__, "Update of job-status to running failed " + unicode(e))
            ConvertSettings.recorded_warnings.append("Update of job-status to running failed " + unicode(e))
        except ValueError as e:
            log_helper.warning(__name__, "Update of job-status to running failed " + unicode(e))
            ConvertSettings.recorded_warnings.append("Update of job-status to running failed " + unicode(e))

    ####################################################################
    # Generate the list of datasets to be converted                    #
    ####################################################################
    try:
        ConvertSettings.dataset_list = ConvertFiles.create_dataset_list(
            input_filenames=input_filenames,
            format_type=ConvertSettings.format_option,
            data_region_option=ConvertSettings.region_option)
        log_helper.info(__name__, "Number of conversion: " + str(len(ConvertSettings.dataset_list)))
    except:
        emailmsg = "ERROR: An error occurred during the generation of the input filelist. \n"
        emailmsg += "       -- No HDF5 output file has been generated. \n"
        emailmsg += "       -- No file has been added to the database. \n"
        emailmsg += "       -- Terminating \n"
        emailmsg += unicode(sys.exc_info())
        WebHelper.send_email(subject="ERROR: Conversion of file failed: " + omsi_outfile,
                             body=emailmsg,
                             email_type='error',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
        log_helper.warning(__name__, emailmsg)
        if ConvertSettings.job_id is not None:
            try:
                WebHelper.update_job_status(filepath=omsi_outfile,
                                            db_server=ConvertSettings.db_server_url,
                                            jobid=ConvertSettings.job_id,
                                            status='error')
            except urllib2.URLError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
                ConvertSettings.recorded_warnings.append("Update of job-status to error failed " + unicode(e))
            except ValueError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
                ConvertSettings.recorded_warnings.append("Update of job-status to error failed " + unicode(e))
        raise

    ####################################################################
    #  Suggest only chunking for the files if requested                #
    ####################################################################
    if ConvertSettings.suggest_file_chunkings:
        ConvertFiles.suggest_chunkings_for_files(ConvertSettings.dataset_list)
        exit()

    ####################################################################
    #   Create the output HDF5 file if needed                          #
    ####################################################################
    try:
        if omsi_outfile is not None:
            ConvertSettings.omsi_output_file = omsi_file.omsi_file(omsi_outfile)
    except:
        emailmsg = "Unexpected error creating the output file: " + unicode(sys.exc_info()[0])
        WebHelper.send_email(subject="ERROR: Conversion of file failed: " + omsi_outfile,
                             body=emailmsg,
                             email_type='error',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
        log_helper.warning(__name__, emailmsg)
        if ConvertSettings.job_id is not None:
            try:
                WebHelper.update_job_status(filepath=omsi_outfile,
                                            db_server=ConvertSettings.db_server_url,
                                            jobid=ConvertSettings.job_id,
                                            status='error')
            except ValueError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
        raise

    ####################################################################
    # Convert all files                                                #
    ####################################################################
    try:
        # Convert all files and record warnings for reporting purposes
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ConvertFiles.convert_files()
            ConvertSettings.recorded_warnings += w
    except:
        emailmsg = "ERROR: An error occurred during the file conversion. \n"
        log_helper.warning(__name__, "ERROR: An error occurred during the file conversion.")
        # Try to close the output file
        try:
            ConvertSettings.omsi_output_file.close_file()
            emailmsg += "  - Successfully closed output file" + " \n"
        except:
            emailmsg += "  - Closing of output HDF5 file failed" + unicode(sys.exc_info()) + " \n"
            pass
        if ConvertSettings.error_handling == "terminate-and-cleanup":
            emailmsg += "  -The generated HDF5 will not be added to the database." + " \n"
            log_helper.warning(__name__, "--The generated HDF5 will not be added to the database.")
            ConvertSettings.add_file_to_db = False
            emailmsg += "  -Attempting to delete the generated HDF5 file." + " \n"
            log_helper.warning(__name__, "--Attempting to delete the generated HDF5 file.")
            os.remove(omsi_outfile)
            emailmsg += "  -Successfully deleted the generated HDF5 file: " + unicode(omsi_outfile) + " \n"
            log_helper.warning(__name__, "--Successfully deleted the generated HDF5 file: " + unicode(omsi_outfile))
        if ConvertSettings.error_handling == "terminate-only" or ConvertSettings.error_handling == "continue-on-error":
            emailmsg += "  -The generated HDF5 will not be added to the database." + " \n"
            log_helper.warning(__name__, "--The generated HDF5 will not be added to the database.")
            ConvertSettings.add_file_to_db = False
            emailmsg += "  -The output HDF5 file (if generate) remains at: " + unicode(omsi_outfile) + " \n"
            emailmsg += "  -Output file found: " + unicode(os.path.exists(omsi_outfile)) + " \n"
            log_helper.warning(__name__, "--The output HDF5 file (if generate) remains at: " + str(omsi_outfile))
            log_helper.warning(__name__, "  Output file found: " + unicode(os.path.exists(omsi_outfile)))
        emailmsg += "\n"
        emailmsg += unicode(sys.exc_info())

        # Add warnings to the email message
        emailmsg += "\n"
        emailmsg += "---------------------------------------------" + " \n"
        emailmsg += "\n"
        for warn in ConvertSettings.recorded_warnings:
            emailmsg += warn.message + "\n"

        # Send email notification if needed
        WebHelper.send_email(subject="ERROR: Conversion of file failed: " + omsi_outfile,
                             body=emailmsg,
                             email_type='error',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
        if ConvertSettings.job_id is not None:
            try:
                WebHelper.update_job_status(filepath=omsi_outfile,
                                            db_server=ConvertSettings.db_server_url,
                                            jobid=ConvertSettings.job_id,
                                            status='error')
            except ValueError as e:
                log_helper.warning(__name__, "Update of job-status to error failed " + unicode(e))
        # Pass on which-ever error has occurred
        raise

    ####################################################################
    #  Close the HDF5 file and exit                                    #
    ####################################################################
    ConvertSettings.omsi_output_file.close_file()

    ####################################################################
    # Update ACL for the raw data files if needed                      #
    ####################################################################
    if ConvertSettings.job_id is not None or ConvertSettings.add_file_to_db:
        for curr_dataset in ConvertSettings.dataset_list:
            basedir = str(curr_dataset['basename'])
            if not os.path.isdir(basedir):
                basedir = os.path.dirname(basedir)
            try:
                WebHelper.set_apache_acl(basedir)
            except:
                ConvertSettings.recorded_warnings += ["Update of raw file ACL failed for" +
                                                      basedir +
                                                      str(sys.exc_info())]

    ####################################################################
    #  Register the file with the database                             #
    ####################################################################
    if ConvertSettings.add_file_to_db and ConvertSettings.job_id is None:
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                status = WebHelper.register_file_with_db(filepath=omsi_outfile,
                                                         db_server=ConvertSettings.db_server_url,
                                                         file_user_name=ConvertSettings.file_user,
                                                         jobid=ConvertSettings.job_id,
                                                         check_add_nersc=ConvertSettings.check_add_nersc)
                ConvertSettings.recorded_warnings += [w]
        except ValueError:
            ConvertSettings.recorded_warnings += [sys.exc_info()]
            status = False
        print "Registered file with DB: " + str(status)

    ####################################################################
    #  Report all recorded warnings on the command line                #
    ####################################################################
    warningmsg = ""
    if len(ConvertSettings.recorded_warnings) > 0:
        warningmsg += "WARNING: The following warnings occurred during the conversion process.\n"
        warningindex = 0
        for warn in ConvertSettings.recorded_warnings:
            warningmsg += u"=========== " + unicode(warningindex) + u" =========== \n"
            try:
                warningmsg += unicode(warn.message) + u"\n"
            except AttributeError:
                warningmsg += unicode(warn) + u"\n"
            warningindex += 1
        log_helper.warning(__name__, warningmsg)

    #####################################################################
    #  Update the status of the job (including registration of the file #
    #####################################################################
    if ConvertSettings.job_id is not None:
        try:
            WebHelper.update_job_status(filepath=omsi_outfile,
                                        db_server=ConvertSettings.db_server_url,
                                        jobid=ConvertSettings.job_id,
                                        status='complete')
        except ValueError as e:
            log_helper.warning(__name__, "ERROR: Update of job status failed: " + e.message)
            ConvertSettings.recorded_warnings += ["Update of job status failed: "+e.message]

    ####################################################################
    #  Send email notification if requested                            #
    ####################################################################
    if len(ConvertSettings.recorded_warnings) == 0:
        emailbody = 'Success'
        if ConvertSettings.job_id is not None:
            emailbody += '\n jobid = %s' % ConvertSettings.job_id

        WebHelper.send_email(subject='Conversion complete: '+str(omsi_outfile),
                             body=emailbody,
                             email_type='success',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)
    elif len(ConvertSettings.email_error_recipients) > 0:
        emailsubject = 'Conversion completed with warnings: '+str(omsi_outfile)
        if ConvertSettings.job_id is not None:
            emailsubject += ' jobid = %s' % ConvertSettings.job_id

        WebHelper.send_email(subject=emailsubject,
                             body=warningmsg,
                             email_type='warning',
                             email_success_recipients=ConvertSettings.email_success_recipients,
                             email_error_recipients=ConvertSettings.email_error_recipients)


####################################################################
####################################################################
#  Define setting for the data conversion                          #
####################################################################
####################################################################
class ConvertSettings(object):
    """This class is used specify the settings for the data conversion"""
    def __init__(self):
        pass

    dataset_list = []
    """ :param dataList: List of python dictionaries describing specific conversion \
                 settings for each conversion task. Each dictionary contains the following keys:

                 * 'basename' : Name of the file to be converted
                 * 'format' : File format to be used (see ConvertSettings.available_formats)
                 * 'exp' : Indicate the experiment the dataset should be stored with. Valid values are \

                              * 'new' : Generate a new experiment for the dataset
                              * 'previous' : Use the same experiment as used for the previous dataset
                              * 1, 2,3...   : Integer value indicating the index of the experiment to be used.
                * 'region' : Optional key with index of the region to be converted. None to merge all regions.
                * 'dataset' : Optional key with index of the dataset to be converted.
                * 'omsi_object': Optional key used to save a pointer to the omsi data object with the converted data
                * 'dependencies': Additional dependencies that should be added for the dataset
    """
    omsi_output_file = None  # The openMSI output data file to be used.

    ####################################################################
    #  Define available options for different parameters              ##
    ####################################################################
    # List of available data formats
    available_formats = file_reader_base.file_reader_base.available_formats()

    # List defining the different options available for handling regions
    available_region_options = ["split",
                                "merge",
                                "split+merge"]
    # Available options for the data write. One chunk at a time ('chunk'), one
    # spectrum at a time ('spectrum'), or all at one once ('all')
    available_io_options = ["chunk", "spectrum", "all", "spectrum_to_image"]
    available_error_options = ["terminate-and-cleanup",
                               "terminate-only",
                               "continue-on-error"]

    ####################################################################
    #  Define the input parameters used during the conversion         ##
    ####################################################################
    # Define how the data should be written to file, one chunk at a time
    # ('chunk'), one spectrum at a time ('spectrum') or all at one once
    # ('all'); read set of spectra and make into image ('spectrum_to_image')
    io_option = "spectrum_to_image"
    # When using the spectrum_to_image io option, what is the maximum block of images we should load in Byte
    io_block_size_limit = 1024 * 1024 * 500  # 500 MB limit
    format_option = None  # Define which file format reader should be used. None=determine automatically
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
    # Specify the server where the database is registered
    db_server_url = WebHelper.default_db_server_url
    file_user = getpass.getuser()  # Specify for which user the file should be registered
    # require_login = False #Require login before conversion
    # Specify the job id
    job_id = None

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
    execute_ticnorm = False  # Define whether tic normalization should be executed
    generate_thumbnail = False  # Should we generate thumbnail
    generate_xdmf = False  # Should we generate an xdmf header file for the file

    # Default NMF parameter settings
    nmf_num_component = 20  # Number of components for the NMF
    nmf_timeout = 600  # Timeout for the NMF
    nmf_num_iter = 2000  # Maximum number of iterations for the NMF
    nmf_tolerance = 0.0001  # Tolerance level for the nmf
    # Should the NMF be computed from the raw data or from the global peak
    # finding data
    nmf_use_raw_data = False

    ####################################################################
    #  Define metadata options                                        ##
    ####################################################################
    metadata = {}   # Dictionary with keys --method, --instrument, --notes for experiment-level metadata and
                    # --method#, instrument#, notes# for dataset-level metadata

    @classmethod
    def parse_input_args(cls, argv):
        """Process input parameters and define the script settings.

           :param argv: The list of input arguments

           :returns: This function returns the following four values:

                * 'input_error' : Boolean indicating whether an error has occurred during the processing of the inputs
                * 'inputWarning' : Boolean indicating whether a warning occurred during the processing of the inputs
                * 'output_filename' : Name for the output HDF5 file
                * 'input_filenames' : List of strings indicating the list of input filenames
        """
        input_error = False    # Parameter used to track whether an error has occurred while parsing the
        output_filename = ""  # The output filename
        input_filenames = []  # The list of input filenames
        helpargs = ["--help", "--h", "-help", "-h"]
        # allargs = ["--no-nmf", "--nmf", "--nmf-nc", "--nmf-timeout",  "--nmf-niter", "--nmf-tolerance",
        #            "--nmf-raw", "--fpg", "--no-fpg", "--fpl", "--no-fpl", "--auto-chunking", "--chunking",
        #            "--no-chunking", "--optimized-chunking", "--compression", "--no-compression",
        #            "--io", "--thumbnail", "--no-thumbnail", "--xdmf", "--no-xdmf",  "--suggest-chunking",
        #            "--format", "--regions", "--add-to-db", "--no-add-to-db", "--db-server", "--user",
        #            "--email", "--email-success", "--email-error", "--error-handling", "--methods",
        #            "--instrument", "--notes", "--jobid", "--ticnorm"] + helpargs
        # Basic sanity check
        if len(argv) < 3:
            last_arg = argv[-1]
            if last_arg in helpargs:
                cls.print_help()
                exit(0)
            else:
                warnings.warn("No command-line options provided. Printing help.")
                cls.print_help()
                return input_error, output_filename, input_filenames

        # Using the input_error and inputWarning parameter we try to keep parsing the input as long as possible to
        # make sure we can inform the user of as many problems as possible
        start_index = 1
        i = start_index
        while i < (len(argv) - 1):
            i = start_index
            current_arg = argv[i]
            # Check if the argument is listed
            # if current_arg.startswith("--") and current_arg not in allargs:
            #     if not (current_arg.startswith("--method") or
            #             current_arg.startswith("--instrument") or
            #             current_arg.startswith("--notes")):
            #         warnings.warn("Argument "+current_arg+" missing in list of all arguments.")
            # Process the argument
            if current_arg == "--no-nmf":
                start_index += 1
                ConvertSettings.execute_nmf = False
                log_helper.info(__name__, "Disable NMF")
            elif current_arg == "--nmf":
                start_index += 1
                ConvertSettings.execute_nmf = True
                log_helper.info(__name__, "Enable NMF")
            elif current_arg == "--nmf-nc":
                start_index += 2
                ConvertSettings.nmf_num_component = int(argv[i + 1])
                log_helper.info(__name__, "Set nmf-nc=" + str(ConvertSettings.nmf_num_component))
            elif current_arg == "--nmf-timeout":
                start_index += 2
                ConvertSettings.nmf_timeout = int(argv[i + 1])
                log_helper.info(__name__, "Set ConvertSettings.nmf_timeout=" + str(ConvertSettings.nmf_timeout))
            elif current_arg == "--nmf-niter":
                start_index += 2
                ConvertSettings.nmf_num_iter = int(argv[i + 1])
                if ConvertSettings.nmf_num_iter < 2:
                    warnings.warn("WARNING: --nfm-niter must be 2 or larger")
                log_helper.info(__name__, "Set nmf-niter=" + str(ConvertSettings.nmf_num_iter))
            elif current_arg == "--nmf-tolerance":
                start_index += 2
                ConvertSettings.nmf_tolerance = float(argv[i + 1])
                log_helper.info(__name__, "Set nmf-tolerance=" + str(ConvertSettings.nmf_tolerance))
            elif current_arg == "--nmf-raw":
                start_index += 1
                ConvertSettings.nmf_use_raw_data = True
                log_helper.info(__name__, "Set nmf-raw=" + str(ConvertSettings.nmf_use_raw_data))
            elif current_arg == "--fpg":
                start_index += 1
                ConvertSettings.execute_fpg = True
                log_helper.info(__name__, "Enable find peaks global")
            elif current_arg == "--no-fpg":
                start_index += 1
                ConvertSettings.execute_fpg = False
                log_helper.info(__name__, "Disable find peaks global")
                if "--fpg" in argv:
                    warnings.warn("WARNING: --no-fpg and --fpg options are conflicting.")
            elif current_arg == "--fpl":
                start_index += 1
                ConvertSettings.execute_fpl = True
                log_helper.info(__name__, "Enable find peaks local")
            elif current_arg == "--no-fpl":
                start_index += 1
                ConvertSettings.execute_fpl = False
                log_helper.info(__name__, "Disable find peaks local")
                if "--fpl" in argv:
                    warnings.warn("WARNING: --no-fpl and --fpl options are conflicting.")
            elif current_arg == "--ticnorm":
                start_index += 1
                ConvertSettings.execute_ticnorm = True
                log_helper.info(__name__, "Enable tic normalization")
            elif current_arg == "--no-ticnorm":
                start_index += 1
                ConvertSettings.execute_ticnorm = False
                log_helper.info(__name__, "Disable tic normalization")
                if "--ticnorm" in argv:
                    warnings.warn("WARNING: --no-ticnorm and --ticnorm options are conflicting")
            elif current_arg == "--auto-chunking":
                start_index += 1
                ConvertSettings.auto_chunk = True
                if "--chunking" in argv:
                    warnings.warn("WARNING: --chunking options will be ignored due to the use of --auto-chunking")
                if "--no-chunking" in argv:
                    warnings.warn("WARNING: --no-chunking and --auto-chunking options are conflicting")
                if "--optimized-chunking" in argv:
                    warnings.warn("WARNING: --optimized-chunking and --auto-chunking options are conflicting")
            elif current_arg == "--chunking":
                start_index += 4
                try:
                    ConvertSettings.chunks = (int(argv[i + 1]), int(argv[i + 2]), int(argv[i + 3]))
                except:
                    log_helper.warning(__name__, "An error occured while parsing the --chunking command. " +
                                                 "Something may be wrong with the indicated chunk sizes for x y z.")
                    input_error = True
                if "--auto-chunking" not in argv:
                    ConvertSettings.auto_chunk = False
                log_helper.info(__name__, "Enable chunking: " + str(ConvertSettings.chunks))
            elif current_arg == "--no-chunking":
                start_index += 1
                ConvertSettings.chunks = None
                ConvertSettings.auto_chunk = False
                log_helper.info(__name__, "Disable chunking")
                if "--auto-chunking" in argv or "--chunking" in argv:
                    warnings.warn("WARNGING: --no-chunking option is conflicting with another chunking option")
            elif current_arg == "--optimized-chunking":
                start_index += 4
                try:
                    ConvertSettings.user_additional_chunks.append(
                        (int(argv[i + 1]), int(argv[i + 2]), int(argv[i + 3])))
                except:
                    log_helper.warning(__name__, "An error accured while parsing the --optimized-chunking command. " +
                                                 "Something may be wrong with the indicated chunk sizes for x y z.")
                    input_error = True
            elif current_arg == "--compression":
                start_index += 1
                # This is already the default
                log_helper.info(__name__, "Enable compression")
            elif current_arg == "--no-compression":
                start_index += 1
                ConvertSettings.compression = None
                ConvertSettings.compression_opts = None
                log_helper.info(__name__, "Disable compression")
                if "--compression" in argv:
                    warnings.warn("WARNING: --no-compression and --compression options are conflicting.")
            elif current_arg == "--io":
                start_index += 2
                try:
                    ConvertSettings.io_option = str(argv[i + 1])
                    if ConvertSettings.io_option not in ConvertSettings.available_io_options:
                        raise ValueError("Invalid io option")
                except:
                    log_helper.warning(__name__, "An error accured while parsing the --io command. " +
                                       "Something may be wrong with the io-type.")
                    input_error = True
            elif current_arg == "--io-block-limit":
                start_index += 2
                try:
                    ConvertSettings.io_block_size_limit = int(argv[i + 1]) * (1024 * 1024)  # MB to Bytes
                    log_helper.info(__name__, "Set io block limit to: " +
                                    str(ConvertSettings.io_block_size_limit) + "Byte")
                except:
                    log_helper.warning(__name__, "An error accured while parsing the --io-block-limit command.")
                    input_error = True
            elif current_arg == "--thumbnail":
                start_index += 1
                ConvertSettings.generate_thumbnail = True
                log_helper.info(__name__, "Enable thumbnail")
            elif current_arg == "--no-thumbnail":
                start_index += 1
                ConvertSettings.generate_thumbnail = False
                log_helper.info(__name__, "Disable thumbnail")
                if "--thumbnail" in argv:
                    warnings.warn("WARNING: --no-thumbnail and --thumbnail options are conflicting.")
            elif current_arg == "--xdmf":
                start_index += 1
                ConvertSettings.generate_xdmf = True
                log_helper.info(__name__, "Enable xdmf")
            elif current_arg == "--no-xdmf":
                start_index += 1
                ConvertSettings.generate_xdmf = False
                log_helper.info(__name__, "Disable xdmf")
                if "--xdmf" in argv:
                    warnings.warn("WARNING: --no-xdmf and --xdmf options are conflicting.")
            elif current_arg in helpargs:
                cls.print_help()
                exit(0)
            elif current_arg == "--suggest-chunking":
                start_index += 1
                ConvertSettings.suggest_file_chunkings = True
            elif current_arg == "--format":
                start_index += 2
                ConvertSettings.format_option = str(argv[i + 1])
                if ConvertSettings.format_option not in ConvertSettings.available_formats.keys():
                    log_helper.warning(__name__,
                                       "ERROR: Indicated --format option " +
                                       ConvertSettings.format_option +
                                       " not supported. Available options are:")
                    log_helper.info(__name__, "     " + str(ConvertSettings.available_formats.keys()))
                    input_error = True
            elif current_arg == "--regions":
                start_index += 2
                ConvertSettings.region_option = str(argv[i + 1])
                if ConvertSettings.region_option not in ConvertSettings.available_region_options:
                    log_helper.warning(__name__,
                                       "ERROR: Indicated --regions option " +
                                       ConvertSettings.region_option +
                                       " not supported. Available options are:")
                    log_helper.warning(__name__, "       " + str(ConvertSettings.available_region_options))
                    input_error = True
            elif current_arg == "--add-to-db":
                start_index += 1
                ConvertSettings.add_file_to_db = True
                ConvertSettings.check_add_nersc = False
            elif current_arg == "--no-add-to-db":
                start_index += 1
                ConvertSettings.add_file_to_db = False
                log_helper.info(__name__, "Disable adding to database")
            elif current_arg == "--db-server":
                start_index += 2
                ConvertSettings.db_server_url = str(argv[i + 1])
                ConvertSettings.check_add_nersc = False
                log_helper.info(__name__, "Set db server to: " + ConvertSettings.db_server_url)
            elif current_arg == "--user":
                start_index += 2
                ConvertSettings.file_user = str(argv[i + 1])
                log_helper.info(__name__, "Set user to: " + ConvertSettings.file_user)
            elif current_arg == "--jobid":
                start_index += 2
                ConvertSettings.job_id = str(argv[i + 1])
                # Discover the job id from the environment variable if requested
                if ConvertSettings.job_id == 'auto':
                    ConvertSettings.job_id = os.environ.get('PBS_JOBID')
                    if ConvertSettings.job_id is None:
                        ConvertSettings.job_id = os.environ.get('SLURM_JOB_ID')
                log_helper.info(__name__, "Set jobid: " + ConvertSettings.job_id)
            elif current_arg == "--email":
                # Consume all email addresses that follow
                for ni in range((i+1), len(argv)):
                    if not argv[ni].startswith("--") and "@" in argv[ni]:
                        ConvertSettings.email_success_recipients.append(str(argv[ni]))
                        ConvertSettings.email_error_recipients.append(str(argv[ni]))
                        start_index = ni+1
                    else:
                        break
            elif current_arg == "--email-success":
                # Consume all email addresses that follow
                for ni in range((i+1), len(argv)):
                    if not argv[ni].startswith("--") and "@" in argv[ni]:
                        ConvertSettings.email_success_recipients.append(str(argv[ni]))
                        start_index = ni+1
                    else:
                        break
            elif current_arg == "--email-error":
                # Consume all email addresses that follow
                for ni in range((i+1), len(argv)):
                    if not argv[ni].startswith("--") and "@" in argv[ni]:
                        ConvertSettings.email_error_recipients.append(str(argv[ni]))
                        start_index = ni+1
                    else:
                        break
            elif current_arg == "--error-handling":
                start_index += 2
                error_option = str(argv[i + 1])
                if error_option in ConvertSettings.available_error_options:
                    ConvertSettings.error_handling = error_option
                else:
                    log_helper.warning(__name__,
                                       "ERROR: Indicated --error-handling option " +
                                       error_option + " invalid. Available options:")
                    log_helper.warning(__name__, "     " + str(ConvertSettings.available_error_options))
                    input_error = True
            elif current_arg.startswith("--methods") or \
                    current_arg.startswith("--instrument") or \
                    current_arg.startswith("--notes"):
                start_index += 2
                ConvertSettings.metadata[current_arg] = unicode(argv[i+1])
            elif current_arg == "--loglevel" or current_arg == "--log":
                start_index += 2
                if str(argv[i+1]) in log_helper.log_levels:
                    log_helper.set_log_level(level=log_helper.log_levels[str(argv[i+1])])
                else:
                    log_helper.warning(__name__, 'Invalid log level ' + str(argv[i+1]))
                    input_error = True
            elif current_arg.startswith("-"):  # Check for single dash to catch the case when a user mistypes the arg
                start_index += 1
                input_error = True
                log_helper.warning(__name__, "Unrecognized input option: " + current_arg)
            else:
                # print "End of input parameters"
                break

        # Determine the list of input filenames
        # remove " in case the user has entered file names with spaces using the
        # "name" syntax
        input_filenames = [name.replace('"', "")
                           for name in argv[start_index:(len(argv) - 1)]]
        # Determine the output filename
        if "--suggest-chunking" in argv:
            # Do not generate an output file if the user just wants to know which
            # chunking to be used
            output_filename = None
            # Check chunking also for the last filename
            input_filenames.append(argv[-1].replace('"', ""))
        else:
            # Remove " in case the user has entered file names with spaces using
            # the "name" syntax
            output_filename = argv[-1].replace('"', "")

        # Check whether all packages needed for generating thumbnails are
        # available on the system
        if ConvertSettings.generate_thumbnail:
            if not pil_available:
                ConvertSettings.generate_thumbnail = False
                log_helper.warning(__name__, "PIL not available. Generation of thumbnail images disabled")

        # Enable chunking if compression is requested and chunking has been
        # disabled
        if (ConvertSettings.chunks is None) and (ConvertSettings.compression is not None):
            timeout = 5*60  # User input timeput after 5 minutes
            num_iter = 3     # Number of tries
            yes_input = ["y", "Y", "Yes", "yes", "YES"]
            no_input = ["n", "N", "No", "no", "NO"]
            for i in range(num_iter):
                print "WARNING: HDF5 compression is only available with chunking enabled. " +\
                      "Do you want to enable chunking? (Y/N)"
                # user_input = raw_input()
                user_input = UserInput.userinput_with_timeout(timeout=timeout, default=None)
                if user_input is None:
                    input_error = True
                    warnings.warn("WARNING: HDF5 compression is only available with chunking enabled." +
                                  " User did not respond to resolve the conflict. Aborting the conversion.")
                elif user_input == yes_input:
                    ConvertSettings.chunks = (4, 4, 2048)
                    log_helper.info(__name__, "Chunking enabled with (4,4, 2048)")
                    break
                elif user_input in no_input:
                    ConvertSettings.compression = None
                    ConvertSettings.compression_opts = None
                    log_helper.info(__name__, "Compression disabled")
                    break
                elif i == (num_iter-1):
                    warnings.warn("User did not respond to resolve the conflict. Aborting the conversion")
                    input_error = True

        if not input_error and not ConvertSettings.suggest_file_chunkings:
            log_helper.info(__name__, "Execute global peak finding (fpg): " + str(ConvertSettings.execute_fpg))
            log_helper.info(__name__, "Execute local peak finding (fpl): " + str(ConvertSettings.execute_fpl))
            log_helper.info(__name__, "Execute nmf: " + str(ConvertSettings.execute_nmf))
            log_helper.info(__name__, "Number of MSI files: " + str(len(input_filenames)))
            log_helper.info(__name__, "Output OMSI file: " + output_filename)

        # Finish and return
        return input_error, output_filename, input_filenames

    @classmethod
    def print_help(cls):
        """Function used to print the help for this script"""
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
        print "indicate which analyses should be executed."
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
        print "--email <email1 email2 ...>: Send notification in case of both error or success to the email address."
        print "--email-success <email1 email2 ...>>: Send notification in case of success to the given email address."
        print "--email-error <email1 email2 ...>>: Send notification in case of error to the given email address."
        print ""
        print "===INPUT DATA OPTIONS=== "
        print ""
        print "Default input data options: --format auto --regions split+merge"
        print "--format <option>: Define which file format is used as input. By default the program tries to"
        print "           automatically determine the input format. This option can be used to indicate"
        print "           the format explicitly to in case the auto option fails. Available options are:"
        print "          " + str(ConvertSettings.available_formats)
        print "--regions <option>: Some file formats (e.g., brucker) allow multiple regions to be imaged and stored"
        print "           in a single file. This option allows one to specify how these regions should be"
        print "           treated during file conversion. E.g., one may want to store i) each region as a "
        print "           separate dataset in the output file (--regions split), ii) all regions combined "
        print "           in a single dataset (--regions merge), or both (--regions split+merge)"
        print "           Available options are:"
        print "          " + str(ConvertSettings.available_region_options)
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
        print "--io <option>: Available options are: " + str(ConvertSettings.available_io_options)
        print "             i) all : Read the full data in memory and write it at once"
        print "             ii) spectrum : Read one spectrum at a time and write it to the file. "
        print "             iii) chunk : Read one chunk at a time and write it to the file."
        print "             The io option applies only for the generation of subsequent chunkings"
        print "             and not the initial iteration over the file to generate the first convert."
        print "             iv) spectrum-to-image: Default option when creating image chunk version from"
        print "             a spectrum-chunk MSI dataset. Read a block of spectra at a time to"
        print "             complete a set of images and then write the block of images at once."
        print "--io-block-limit <MB>: When using spectrum-to-image io (default when using auto-chunking),"
        print "             what should the maximum block in MB that we load into memory. (Default=2000MB)"
        print ""
        print "===DATABSE OPTIONS=== "
        print ""
        print "These options control whether the generated output file should be added to a server database"
        print "to manage web file access permissions"
        print "Default options are: --add-to-db --db-server http://openmsi.nersc.gov"
        print "--add-to-db : Explicitly add the output HDF5 file to the database. This option has no effect"
        print "              if --jobid is set as the file is added through the update of the job status"
        print "              in this case."
        print "--no-add-to-db : Disable adding the file to the database."
        print "--db-server : Specify the online server where the file should be registers. Default is"
        print "              http://openmsi.nersc.gov "
        print "--user : Name of the user that should be assigned as user. By default the user is"
        print "          determined automatically based on the file path."
        print "--jobid : ID of the job. If set to 'auto' then the environment variable PBS_JOBID is used. "
        print "          NOTE: If job ID is set then we assume that the job has been scheduled via the"
        print "          the automated system and that the job is managed. As such the file will be added,"
        print "          to the database by updating the job status and NOT by explicitly adding the file."
        # print "--login : If set, then the script will ask for login information
        # at the beginning of the script."
        print ""
        print "===ANALYSIS OPTIONS=== "
        print ""
        print "NMF: Default ON: (nc=20, timeout=600, niter=2000, tolerance=0.0001, raw=False)"
        print "--nmf : Compute the nmf for all the input data files and store the results in the"
        print "        HDF5 file. NOTE: If global peak-finding (fpg) is performed, then"
        print "        nmf will be performed on the peak-cube, otherwise on the raw data"
        print "--no-nmf: Disable the execution of nmf"
        print "--nmf-nc <number>: Number of components to be computed by the NMF. (default nc=20)"
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
        print "Local Peak Finding: Default OFF:"
        print "--fpl : Compute the local peak finding for all input data files and save results"
        print "        in the HDF5 file"
        print "--no-fpl: Disable the local peak finding (DEFAULT)"
        print ""
        print ""
        print "TIC normalization:"
        print "--ticnorm : Compute tic normalization"
        print "--no-ticnorm : Disable computation of tic normaltization (DEFAULT)"
        print ""
        print "---OTHER OPTIONS---"
        print ""
        print "Generate Thumbnail image: Default OFF:"
        print "--thumbnail: Generate thumbnail image for the file based on, in order of availability:"
        print "             * The first three components of the NMF"
        print "             * The three most intense peaks from the global peak finding (fpg)"
        print "             * The three most intense peaks in the raw data that are at least 1 percent"
        print "               of the total m/z range apart."
        print "--no-thumbnail: Do not generate a thumbnail image."
        print ""
        print "Generate XDMF header file for output file: Default OFF:"
        print "--xdmf: Write XDMF XML-based header-file for the output HDF5 file."
        print "--no-xdmf: Do not generate a XDMF XML-based header for the HDF5 file."
        print ""
        print "===Metadata Options==="
        print ""
        print "NOTE: Input datasets are numbers starting from 0 based on there order on the command line."
        print ""
        print "--methods : JSON describing the experimental methods"
        print "--methods# : JSON describing the experimental methods for input file number #"
        print "--instrument : JSON dictionary describing the instrument"
        print "--instrument# : JSON dictionary describing the instrument for input file number #"
        print "--notes : JSON dictionary with additional user notes about the data"
        print "--notes# : JSON dictionary with additional notes for input file number #"


####################################################################
####################################################################
#    Data conversion                                               #
####################################################################
####################################################################
class ConvertFiles(object):
    """Class providing a number of functions for converting various file types to OMSI,
       including a number of helper functions related to the data conversion.
    """
    def __init__(self):
        pass

    @staticmethod
    def convert_files():
        """Convert all files in the given list of files with their approbriate conversion options"""
        ####################################################################
        #  Convert the MSI files and compute the requested analyses       ##
        ####################################################################
        # Iterate over all input datasets and convert them
        for curr_index, curr_dataset in enumerate(ConvertSettings.dataset_list):  # xrange(startIndex,len(argv)-1):
            ####################################################################
            #  Convert the raw data to HDF5                                   ##
            ####################################################################
            basefile = curr_dataset['basename']  # argv[i].replace('"', "")
            log_helper.info(__name__, "Converting: " + basefile)
            try:
                log_helper.info(__name__, "Selected Region: " + str(curr_dataset['region']))
            except:
                pass
            if not ConvertSettings.auto_chunk:
                log_helper.info(__name__, "HDF5 chunking: " + str(ConvertSettings.chunks))
            log_helper.info(__name__,
                            "HDF5 compression: " + str(ConvertSettings.compression) +
                            ", " + str(ConvertSettings.compression_opts))

            # Open the input file
            try:
                curr_format = curr_dataset['format']
                log_helper.info(__name__, "Input file format: " + str(curr_format))
                if curr_format is not None:
                    input_file = ConvertSettings.available_formats[curr_format](basename=basefile,
                                                                                requires_slicing=False)
                    if input_file.supports_regions():
                        input_file.set_region_selection(region_index=curr_dataset['region'])
                    if input_file.supports_multidata():
                        input_file.set_dataset_selection(dataset_index=curr_dataset['dataset'])
                else:
                    log_helper.warning(__name__,
                                       "ERROR: The following file will not be converted. File type unknown for: " +
                                       basefile)
                    log_helper.info(__name__,
                                    "INFO: If you know the correct file format then try setting " +
                                    "appropriate --format option.")
                    if ConvertSettings.error_handling == "continue-on-error":
                        pass
                    elif ConvertSettings.error_handling == "terminate-and-cleanup" or\
                            ConvertSettings.error_handling == "terminate-only":
                        raise ValueError("Unrecognized file format.")

                log_helper.info(__name__, "In data shape: " + str(input_file.shape))
            except:
                log_helper.warning(__name__, "ERROR: Unexpected error opening the input file:", sys.exc_info()[0])
                if ConvertSettings.error_handling == "continue-on-error":
                    log_helper.info(__name__,
                                    basefile + " failure during conversion. Skipping the file and continue conversion.")
                    continue
                elif ConvertSettings.error_handling == "terminate-and-cleanup" or \
                        ConvertSettings.error_handling == "terminate-only":
                    raise

            if ConvertSettings.auto_chunk:
                spec_chunks, img_chunks, temp_chunks = ConvertFiles.suggest_chunking(
                    xsize=input_file.shape[0],
                    ysize=input_file.shape[1],
                    mzsize=input_file.shape[2],
                    dtype=input_file.data_type,
                    print_results=True)
                ConvertSettings.chunks = spec_chunks
                additional_chunks = ConvertSettings.user_additional_chunks + [img_chunks]
                log_helper.info(__name__, "Converting data using the following chunking options:")
                log_helper.info(__name__, "     - Spectrum chunking: " + str(ConvertSettings.chunks))
                log_helper.info(__name__, "     - Image chunking:    " + str(additional_chunks[0]))

            else:
                additional_chunks = ConvertSettings.user_additional_chunks

            # Get the mz data
            mzdata = input_file.mz

            # Define the layout for the data in the HDF5
            # Create a new experiment for the img file using the basefile as
            # identifer string
            # Create a new experiment for this dataset
            if curr_dataset['exp'] == 'new' or ConvertSettings.omsi_output_file.get_num_experiments() == 0:
                exp = ConvertSettings.omsi_output_file.create_experiment(exp_identifier=basefile)
                # Create an empty method description
                sample = exp.create_method_info()
                # Create an empty instrument description
                instrument = exp.create_instrument_info(instrument_name="undefined", mzdata=mzdata)
            # Store this dataset in combination with the previous one
            elif curr_dataset['exp'] == 'previous':
                exp = ConvertSettings.omsi_output_file.get_experiment(ConvertSettings.omsi_output_file.get_num_experiments() - 1)
            # Store this dataset in combination with the given experiment index
            else:
                exp = ConvertSettings.omsi_output_file.get_experiment(curr_dataset['exp'])

            # Allocate space in the HDF5 file for the img data
            data_dataset, mz_dataset, data_group = exp.create_msidata_full_cube(
                data_shape=input_file.shape,
                data_type=input_file.data_type,
                chunks=ConvertSettings.chunks,
                compression=ConvertSettings.compression,
                compression_opts=ConvertSettings.compression_opts)
            mz_dataset[:] = mzdata
            data = omsi_file.omsi_file_msidata(data_group=data_group,
                                               preload_mz=False,
                                               preload_xy_index=False)
            ConvertFiles.write_data(input_file=input_file,
                                    data=data,
                                    data_io_option='spectrum',  # ConvertSettings.io_option,
                                    chunk_shape=ConvertSettings.chunks,
                                    write_progress=(ConvertSettings.job_id is None))
            ConvertSettings.omsi_output_file.flush()

            # Generate any additional data copies if requested
            for chunkSpec in additional_chunks:
                log_helper.info(__name__, "Generating optimized data copy: " + str(chunkSpec))
                tempdata = data.create_optimized_chunking(chunks=chunkSpec,
                                                          compression=ConvertSettings.compression,
                                                          compression_opts=ConvertSettings.compression_opts,
                                                          copy_data=False,
                                                          print_status=True)

                # Determine the number of images in a block when using spectrum_to_image io
                num_images_block = ConvertSettings.io_block_size_limit / float(input_file.shape[0] *
                                                                               input_file.shape[1] *
                                                                               tempdata.dtype.itemsize)
                num_images_block = int(num_images_block) if num_images_block > 1 else 1
                # Determine the chunking parameter based on the io option
                chunk_shape_write = (input_file.shape[0], input_file.shape[1], num_images_block) if \
                                        ConvertSettings.io_option == "spectrum_to_image" \
                                        else chunkSpec
                # Write the data using the given io strategry
                ConvertFiles.write_data(input_file=data_dataset,
                                        data=tempdata,
                                        data_io_option=ConvertSettings.io_option,
                                        chunk_shape=chunk_shape_write,
                                        write_progress=(ConvertSettings.job_id is None))
                ConvertSettings.omsi_output_file.flush()

            # Save the pointer to the omsi object in the file in the dataset list
            curr_dataset['omsi_object'] = data
            if isinstance(input_file, file_reader_base.file_reader_base_multidata) or \
                    isinstance(input_file, file_reader_base.file_reader_base_with_regions):
                curr_dataset['dependencies'] = input_file.get_dataset_dependencies()

            ####################################################################
            #  Add the metadata from the file-reader for the dataset          ##
            ####################################################################
            dataset_metadata_from_reader = input_file.get_dataset_metadata()
            method_metadata_from_reader = input_file.get_method_metadata()
            instrument_metadata_from_reader = input_file.get_instrument_metadata()
            if dataset_metadata_from_reader is not None:
                curr_dataset['omsi_object'].create_metadata_collection(metadata=dataset_metadata_from_reader)

            ####################################################################
            # Close the input data file and free up any allocated memory      ##
            ####################################################################
            input_file.close_file()

            ####################################################################
            #  Add the user-defined metadata for the dataset                  ##
            ####################################################################
            method_meta_key = '--methods'+str(curr_index)
            notes_meta_key = '--notes'+str(curr_index)
            instrument_meta_key = '--instrument'+str(curr_index)
            has_method_metadata = method_meta_key in ConvertSettings.metadata or\
                                  notes_meta_key in ConvertSettings.metadata or\
                                  len(method_metadata_from_reader) > 0
            has_instrument_metadata = instrument_meta_key in ConvertSettings.metadata or\
                                      len(instrument_metadata_from_reader) > 0
            if has_method_metadata:
                methods_info = data.create_method_info()  # Create will return the existing one if one exists
                if method_meta_key in ConvertSettings.metadata:
                    methods_info.add_metadata(metadata=metadata_value(
                        name='methods',
                        value=unicode(ConvertSettings.metadata[method_meta_key]),
                        description='Description of the experimental methods'))
                if notes_meta_key in ConvertSettings.metadata:
                    methods_info.add_metadata(metadata=metadata_value(
                        name='notes',
                        value=unicode(ConvertSettings.metadata[notes_meta_key]),
                        description='Additional notes about the methods'))
                if len(method_metadata_from_reader) > 0:
                    methods_info.add_metadata(method_metadata_from_reader)

            if has_instrument_metadata:

                parent_instrument_info = exp.get_instrument_info()
                if parent_instrument_info and parent_instrument_info.get_instrument_name():
                    inst_name = parent_instrument_info.get_instrument_name()[0]
                else:
                    inst_name = None
                instrument_info = data.create_instrument_info(instrument_name=inst_name)
                # Create will return the existing one if one exists
                if instrument_meta_key in ConvertSettings.metadata:
                    instrument_info.add_metadata(metadata=metadata_value(
                                                 name='description',
                                                 value=unicode(ConvertSettings.metadata[instrument_meta_key]),
                                                 description='Description of the instrument'))
                if len(instrument_metadata_from_reader) > 0:
                    instrument_info.add_metadata(instrument_metadata_from_reader)

            ####################################################################
            #  Execute the requested analyses                                 ##
            ####################################################################
            # Compute the local peak finding and add the peak-cube data to the file
            try:
                if ConvertSettings.execute_fpl:
                    log_helper.info(__name__, "Executing local peak finding")
                    # Execute the peak finding
                    fpl = omsi_findpeaks_local(name_key="omsi_findpeaks_local_" + str(time.ctime()))
                    fpl.execute(msidata=data, mzdata=mzdata, printStatus=True)
                    # Save the peak-finding to file
                    ana, anaindex = exp.create_analysis(fpl)
            except:
                warnings.warn("Local peak finding failed.")
                e = sys.exc_info()[0]
                log_helper.warning(__name__, "Local peak finding failed. " + str(e))
            # Compute the global peak finding and add the peak-cube data to the
            # file
            try:
                fpg_executed = False
                if ConvertSettings.execute_fpg:
                    log_helper.info(__name__, "Executing global peak finding")
                    # Execute the peak finding
                    fpg = omsi_findpeaks_global(name_key="omsi_findpeaks_global_" + str(time.ctime()))
                    fpg.execute(msidata=data, mzdata=mzdata)
                    # Save the peak-finding to file
                    ana, anaindex = exp.create_analysis(fpg)
                    fpgpath = ana.get_managed_group().name
                    fpg_executed = True
            except:
                warnings.warn("Global peak finding failed")
                e = sys.exc_info()[0]
                log_helper.warning(__name__, "ERROR: Global peak finding failed: " + str(e))
                fpg_executed = False
            try:
                # Compute the nmf analysis if requested
                if ConvertSettings.execute_nmf:
                    log_helper.info(__name__, "Executing nmf")
                    # Exectue the nmf analysis on the peak-cube or if peak finding was not performed then
                    # execute nmf on the raw data data
                    nmf = omsi_nmf(name_key="omsi_nmf_" + str(time.ctime()))
                    if ConvertSettings.execute_fpg and not ConvertSettings.nmf_use_raw_data:
                        if fpg_executed:
                            log_helper.info(__name__, "   Using peak-cube data for NMF")
                            fpgdata = ana['peak_cube']
                            nmf.execute(msidata=fpgdata,
                                        numComponents=ConvertSettings.nmf_num_component,
                                        timeOut=ConvertSettings.nmf_timeout,
                                        numIter=ConvertSettings.nmf_num_iter,
                                        tolerance=ConvertSettings.nmf_tolerance)
                        else:
                            warnings.warn("Skipping NMF due to failed global peak finder.")
                            log_helper.warning(__name__, "Skipping NMF due to failed global peak finder.")
                    else:
                        log_helper.info(__name__, "   Using raw data for NMF")
                        nmf.execute(msidata=data,
                                    numComponents=ConvertSettings.nmf_num_component,
                                    timeOut=ConvertSettings.nmf_timeout,
                                    numIter=ConvertSettings.nmf_num_iter,
                                    tolerance=ConvertSettings.nmf_tolerance)
                    # Save the nmf results to file
                    ana, anaindex = exp.create_analysis(nmf)
            except:
                warnings.warn("NMF failed.")
                e = sys.exc_info()[0]
                log_helper.warning(__name__, "ERROR: NMF failed: " + str(e))

            try:
                # Compute the tic normalization
                if ConvertSettings.execute_ticnorm:
                    log_helper.info(__name__, "Executing tic normalization")
                    tic = omsi_tic_norm(name_key="tic_norm_" + str(time.ctime()))
                    tic.execute(msidata=data, mzdata=mzdata)
                    ana, anaindex = exp.create_analysis(tic)
            except:
                warnings.warn("TIC normalization failed.")
                e = sys.exc_info()[0]
                log_helper.warning(__name__, "ERROR: TIC normalization failed: " + str(e))

            ####################################################################
            #  Generate the thumbnail image for the current file              ##
            ####################################################################
            try:
                if ConvertSettings.generate_thumbnail:
                    log_helper.info(__name__, "Generating the thumbnail image")
                    if ConvertSettings.execute_nmf:
                        log_helper.info(__name__, "   Generating thumbnail from NMF data")
                        # Get the NMF data
                        ho = nmf['ho']
                        numx = ho.shape[0]
                        numy = ho.shape[1]
                        # Generate images for the first three NMF components
                        d1 = np.log(ho[:, :, 0].reshape((numx, numy)) + 1)
                        d1 = d1 / np.max(d1)
                        im1 = Image.fromarray(d1.astype('float') * 255).convert('L')
                        d2 = np.log(ho[:, :, 1].reshape((numx, numy)) + 1)
                        d2 = d2 / np.max(d2)
                        im2 = Image.fromarray(d2.astype('float') * 255).convert('L')
                        d3 = np.log(ho[:, :, 2].reshape((numx, numy)) + 1)
                        d3 = d3 / np.max(d3)
                        im3 = Image.fromarray(d3.astype('float') * 255).convert('L')
                        # Generate thumbnail by merging the three gray-scale images as
                        # an RGB image
                        thumbnail = Image.merge('RGB', (im1, im2, im3))
                        expname = str(exp.get_managed_group().name)
                        expindex = expname[7:len(expname)]
                        thumbnail_filename = ConvertSettings.omsi_output_file.hdf_filename + \
                            "_" + expindex + ".png"
                        thumbnail.save(thumbnail_filename, 'PNG')
                    elif ConvertSettings.execute_fpg:
                        log_helper.info(__name__, "    Generating thumbnail from FPG data")
                        # Get the global peak finding data and compute the maximum peak
                        # values for each peak
                        fpgdata = fpg['peak_cube'][:]
                        max_peak_values = fpgdata.max(axis=0).max(axis=0)
                        numx = fpgdata.shape[0]
                        numy = fpgdata.shape[1]
                        # Generate images for the three most intense peaks
                        s = np.argsort(max_peak_values)
                        d1 = np.log(fpgdata[:, :, s[-1]].reshape((numx, numy)) + 1)
                        d1 = d1 / np.max(d1)
                        im1 = Image.fromarray(d1.astype('float') * 255).convert('L')
                        d2 = np.log(fpgdata[:, :, s[-2]].reshape((numx, numy)) + 1)
                        d2 = d2 / np.max(d2)
                        im2 = Image.fromarray(d2.astype('float') * 255).convert('L')
                        d3 = np.log(fpgdata[:, :, s[-3]].reshape((numx, numy)) + 1)
                        d3 = d3 / np.max(d3)
                        im3 = Image.fromarray(d3.astype('float') * 255).convert('L')
                        # Generate thumbnail by merging the three gray-scale images as
                        # an RGB image
                        thumbnail = Image.merge('RGB', (im1, im2, im3))
                        expname = str(exp.get_managed_group().name)
                        expindex = expname[7:len(expname)]
                        thumbnail_filename = ConvertSettings.omsi_output_file.hdf_filename + \
                            "_" + expindex + ".png"
                        thumbnail.save(thumbnail_filename, 'PNG')
                    else:
                        log_helper.info(__name__,
                                        "Generation of thumbnail from raw data is not yet supported. " +
                                        "Thumbnail not generated.")
                        log_helper.info(__name__, "Enable --nmf or --fpg in order to generate a thumbnail image.")
                        # print "    Generating thumbnail from raw data"
                        # Find three most intense peaks that are at least 1% of the m/z range appart
                        # numx = data.shape[0]
                        # numy = data.shape[1]
                        # numZ = data.shape[2]
                        # minMzStep = numZ / 100
                        # mzImageRange = numZ/4000
                        # max_peak_values = np.zeros(numZ)
                        # stepping = np.arange(0, numZ-5000, 5000)
                        # for i in stepping:
                        #    max_peak_values[i:(i+5000)] = data[:, :, i:(i+5000)].std(axis=0).std(axis=0).reshape(5000)
                        # max_peak_values[stepping[-1]:numZ] = \
                        #     data[:, :, stepping[-1]:numZ].std(axis=0).std(axis=0).reshape(numZ - stepping[-1])
                        # print max_peak_values.shape
                        # s = np.argsort(max_peak_values)
                        # i1 = s[-1]
                        # i2 = s[-2]
                        # for i in reversed(range(0, len(s))):
                        #   if abs(s[i] - i1) > minMzStep:
                        #       i2 = s[i]
                        #       break
                        # i3 = s[-3]
                        # for i in reversed(range(0, len(s))):
                        #   if abs(s[i] - i1) > minMzStep and abs(s[i]-i2)>minMzStep:
                        #       i2 = s[i]
                        #       break
                        # print minMzStep
                        # print str(i1)+" "+str(i2)+" "+str(i3)
                        # print str(mzdata[i1])+" "+str(mzdata[i2])+" "+str(mzdata[i3])
                        # low = max(0, s[i1]-mzImageRange)
                        # hi = min(numZ, s[i1]+mzImageRange)
                        # print str(low) + " " + str(hi)
                        # d1 = data[:, :,low:hi].max(axis=2).reshape((numx, numy))
                        # d1 = d1 / np.max(d1)
                        # im1 = Image.fromarray(d1.astype('float')*255).convert('L')
                        # low = max(0, s[i2]-mzImageRange)
                        # hi = min(numZ, s[i2]+mzImageRange)
                        # print str(low) + " " + str(hi)
                        # d2 = data[:, :,low:hi].max(axis=2).reshape((numx, numy))
                        # d2 = d2 / np.max(d2)
                        # im2 = Image.fromarray(d2.astype('float')*255).convert('L')
                        # low = max(0, s[i3]-mzImageRange)
                        # hi = min(numZ, s[i3]+mzImageRange)
                        # print str(low) + " " + str(hi)
                        # d3 = data[:, :,low:hi].max(axis=2).reshape((numx, numy))
                        # d3 = d3 / np.max(d3)
                        # im3 = Image.fromarray(d3.astype('float')*255).convert('L')
                        # Generate thumbnail by merging the three gray-scale images as an RGB image
                        # thumbnail = Image.merge('RGB', (im1, im2, im3))
                        # expname = str(exp.get_managed_group().name)
                        # expindex = expname[7:len(expname)]
                        # thumbnail_filename = ConvertSettings.omsi_output_file.hdf_filename+"_"+expindex+".png"
                        # thumbnail.save(thumbnail_filename, 'PNG')
            except ImportError:
                log_helper.warning(__name__, "ERROR: Thumbnail generation failed. I/O error.", sys.exc_info()[0])
                warnings.warn("ERROR: Thumbnail generation failed. I/O error. " + str(sys.exc_info()))
                pass
            except:
                log_helper.warning(__name__, "ERROR: Thumbnail generation failed. Unexpected error:", sys.exc_info()[0])
                warnings.warn("ERROR: Thumbnail generation failed. Unexpected error. " + str(sys.exc_info()))
                pass

        ####################################################################
        # Generate any dependencies the datasets may have                  #
        ####################################################################
        for curr_index, curr_dataset in enumerate(ConvertSettings.dataset_list):
            if 'dependencies' in curr_dataset:
                if len(curr_dataset['dependencies']) > 0:
                    for dependency in curr_dataset['dependencies']:
                        # Find the corresponding dataset that we converted if necessary
                        if 'omsi_object' not in dependency or dependency['omsi_object'] is None:
                            for dep_dataset in ConvertSettings.dataset_list:
                                match_found = dep_dataset['basename'] == dependency['basename']
                                if 'dataset' in dependency and dependency['dataset'] is not None:
                                    if 'dataset' in dep_dataset:
                                        match_found &= dep_dataset['dataset'] == dependency['dataset']
                                    else:
                                        continue
                                if 'region' in dependency and dependency['region'] is not None:
                                    if 'region' in dep_dataset:
                                        match_found &= dep_dataset['region'] == dependency['region']
                                    else:
                                        continue
                                if match_found:
                                    dependency['omsi_object'] = dep_dataset['omsi_object']
                                    break

                        # We may not have converted all datasets so we need to check whether we found
                        # the corresponding omsi_object for the current dependency
                        if dependency['omsi_object'] is not None:
                            dependency.pop('basename')
                            dependency.pop('region')
                            dependency.pop('dataset')
                            new_omsi_dependency = dependency_dict()
                            new_omsi_dependency.update(dependency)
                            # Define and save the dependency
                            # new_omsi_dependency = dependency_dict(link_name=dependency['link_name'],
                            #                                      omsi_object=dependency['omsi_object'],
                            #                                      help=dependency['help'],
                            #                                      dependency_type=dependency['dependency_type'])
                            curr_dataset['omsi_object'].add_dependency(new_omsi_dependency)

        ####################################################################
        #  Generate the XDMF header file for the HDF5 file                ##
        ####################################################################
        if ConvertSettings.generate_xdmf:
            ConvertSettings.omsi_output_file.write_xdmf_header(
                ConvertSettings.omsi_output_file.get_filename() + ".xdmf")

    @staticmethod
    def check_format(name, format_type):
        """Helper function used to determine the file format that should be used

           :param name: Name of the folder/file that we should read
           :param format_type: String indicating the format-option given by the user.
                        If the format is not determined (i.e., "auto") then this function
                        tries to determine the appropriate format. Otherwise this option
                        is returned as is, as the user explicitly said which format should
                        be used.
           :returns: String indicating the appropriate format. Returns None in case no valid option was found.
        """
        # Option 1: The user told us the format we should use
        if format_type is not None:
            if format_type in ConvertSettings.available_formats.keys():
                return format_type
            else:
                return None
        # Option 2: We need to determine the format ourselves
        # if os.path.exists(name+".hdr") and os.path.exists(name+".t2m") and os.path.exists(name+".img"):
        #    return "img"
        for formatname, formatclass in ConvertSettings.available_formats.items():
            try:
                if formatclass.is_valid_dataset(name):
                    return formatname
            except:
                pass
        return None

    @staticmethod
    def create_dataset_list(input_filenames, format_type=None, data_region_option="split+merge"):
        """Based on the list of input_filenames, generate the ConvertSettings.dataset_list, which contains a
           dictionary describing each conversion job

           :param input_filenames: List of names of files to be converted.
           :param format_type: Define which file-format should be used. Default value is 'auto' indicating the
                       function should determine for each file the format to be used.
                       See also ConvertSettings.available_formats parameter.
           :param data_region_option: Define how different regions defined for a file should be handled. E.g., one may
                        want to split all regions into indiviudal datasets ('split'), merge all regions into a single
                        dataset ('merge'), or do both ('split+merge'). See also the
                        ConvertSettings.available_region_options parameter for details. By default the function will
                        do 'split+merge'.

           :returns: List of dictionaries describing the various conversion jobs. Each job is described by a dict
                with the following keys:

                * 'basename' : The base name of the file or directory with the data
                * 'format' : The data format to be used
                * 'dataset' : The index of the dataset to be converted if the input stores multiple data cubes
                * 'region' : The index of the region to be converted if the input defines multiple regions
                * 'exp' : One of 'previous' or 'new', defining whether a new experiment should be created or whether \
                          the experiment from the previous conversion(s) should be reused.

        """
        re_dataset_list = []
        for i in input_filenames:
            # Remove " in case the user has entered file names with spaces using
            # the "name" syntax
            currds = {'basename': i.rstrip('"').rstrip("'").lstrip('"').lstrip("'"),
                      'exp': 'new'}
            currds['format'] = ConvertFiles.check_format(name=currds['basename'],
                                                         format_type=format_type)
            supports_regions = False
            supports_multidata = False
            if currds['format'] is not None:
                supports_regions = ConvertSettings.available_formats[currds['format']].supports_regions()
                supports_multidata = ConvertSettings.available_formats[currds['format']].supports_multidata()
            if supports_multidata:
                try:
                    tempfile = ConvertSettings.available_formats[currds['format']](
                        basename=currds['basename'],
                        requires_slicing=False)
                    for dataset_index in xrange(0, tempfile.get_number_of_datasets()):
                        nds = {'basename': currds['basename'],
                               'format': currds['format'],
                               'dataset': dataset_index,
                               'exp': 'previous'}
                        re_dataset_list.append(nds)
                except:
                    log_helper.info(__name__, "ERROR: Unexpected error opening the input file:", sys.exc_info()[0])
                    if ConvertSettings.error_handling == "continue-on-error":
                        log_helper.info(__name__, currds['basename'] + " failure during conversion. Skipping the file and continue.")
                        continue
                    elif ConvertSettings.error_handling == "terminate-and-cleanup" or \
                            ConvertSettings.error_handling == "terminate-only":
                        raise
                    else:
                        raise
            elif supports_regions:
                if data_region_option == "merge" or data_region_option == "split+merge":
                    currds['region'] = None
                    re_dataset_list.append(currds)

                if data_region_option == 'split' or data_region_option == 'split+merge':
                    try:
                        tempfile = ConvertSettings.available_formats[currds['format']](
                            basename=currds['basename'],
                            requires_slicing=False)
                        for ri in xrange(0, tempfile.get_number_of_regions()):
                            nds = {'basename': currds['basename'],
                                   'format': currds['format'],
                                   'region': ri,
                                   'exp': 'previous'}
                            re_dataset_list.append(nds)
                    except:
                        log_helper.warning(__name__,
                                           "ERROR: Unexpected error opening the input file:" + str(sys.exc_info()[0]))
                        if ConvertSettings.error_handling == "continue-on-error":
                            log_helper.info(__name__,
                                            currds['basename'] +
                                            " failure during conversion. Skipping the file and continue.")
                            continue
                        elif ConvertSettings.error_handling == "terminate-and-cleanup" or \
                                ConvertSettings.error_handling == "terminate-only":
                            raise
                        else:
                            raise

                # This is just to make sure that we have checked everything
                if data_region_option not in ConvertSettings.available_region_options:
                    warnings.warn("WARNING: Undefined region option. Using the default option split+merge")
            else:
                re_dataset_list.append(currds)

        # Return the list of conversion jobs
        return re_dataset_list

    @staticmethod
    def suggest_chunkings_for_files(in_dataset_list):
        """Helper function used to suggest good chunking strategies for a given set of files.

           :param in_dataset_list: Python list of dictionaries describing the settings to be used
                                   for the file conversion

           :returns: This function simply prints results to standard-out but does not return anything.
        """
        spectrum_chunk = [()] * len(in_dataset_list)
        slice_chunk = [()] * len(in_dataset_list)
        balanced_chunk = [()] * len(in_dataset_list)
        for di in range(len(in_dataset_list)):

            currdataset = in_dataset_list[di]
            basefile = currdataset["basename"]
            try:
                print "Suggested Chunkings: " + basefile
                currformat = currdataset["format"]
                inputfile = ConvertSettings.available_formats[currformat](basename=basefile,
                                                                          requires_slicing=False)
                if inputfile.supports_regions():
                    inputfile.set_region_selection(currdataset["region"])
                elif inputfile.supports_multidata():
                    inputfile.set_dataset_selection(currdataset['dataset'])
                else:
                    warnings.warn("WARNING: Type of file could not be determined for: " + basefile)
                    continue
                print "In data shape: " + str(inputfile.shape)
                spectrum_chunk[di], slice_chunk[di], balanced_chunk[di] = ConvertFiles.suggest_chunking(
                    xsize=inputfile.shape[0],
                    ysize=inputfile.shape[1],
                    mzsize=inputfile.shape[2],
                    dtype=inputfile.data_type,
                    print_results=True)
            except:
                warnings.warn("Error while trying to generate chunking suggestion for " +
                              basefile + " " + str(sys.exc_info()))
        return spectrum_chunk, slice_chunk, balanced_chunk

    @staticmethod
    def suggest_chunking(xsize, ysize, mzsize, dtype, print_results=False):
        """Helper function used to suggest god chunking strategies for a given data cube

           :param xsize: Size of the dataset in x.
           :param ysize: Size o the dataset in y.
           :param mzsize: Size of the dataset in mz.
           :param print_results: Print the results to the console.

           :returns: Three tuples:

            * ``spectrum_chunk`` : The chunking to be used to optimize selection of spectra.
            * ``slice_chunk`` : The chunking to be used to optimize selection of image slices.
            * ``balanced_chunk`` : The chunking that would provide a good balance in performance for \
                                different selection strategies.
        """
        # Make sure that all sizes are treated as 64bit int to avoid errors due to
        # artificial cutoffs.
        xsize = int(xsize)
        ysize = int(ysize)
        mzsize = int(mzsize)
        imagesize = xsize * ysize

        # Variable settings
        numbytes = np.dtype(dtype).itemsize  # Number of bytes per data value
        # Nober of values that fit into a 64Kb chunk
        suggested_num_values = (1024 * 64) / numbytes
        # How much larger than the suggestNumValues should we allow a chunk to
        # become
        max_factor = 4

        # Define the chunking that would be good for spectra
        spectrum_x_chunk = 1
        spectrum_y_chunk = 1
        factor1 = math.floor(mzsize / float(suggested_num_values))
        if factor1 == 0:
            factor1 = 1
        spectrum_mz_chunk1 = int(math.ceil(mzsize / float(factor1)))
        factor1 = math.ceil(mzsize / float(spectrum_mz_chunk1))
        overhead1 = ((factor1 * spectrum_mz_chunk1) - mzsize) * imagesize * numbytes
        # overhead1 = ((spectrum_mz_chunk1 *  math.ceil (float(mzsize) / float(spectrum_mz_chunk1))) - \
        #              (spectrum_mz_chunk1 *  math.floor(float(mzsize) / float(spectrum_mz_chunk1)))) * \
        #            imagesize * numbytes
        # overhead1 = (spectrum_mz_chunk1 - math.ceil(float(mzsize) % float(spectrum_mz_chunk1))) * imagesize * numbytes
        factor2 = math.ceil(mzsize / float(suggested_num_values))
        if factor2 == 0:
            factor2 = 1
        spectrum_mz_chunk2 = int(math.ceil(mzsize / float(factor2)))
        factor2 = math.ceil(mzsize / float(spectrum_mz_chunk2))
        overhead2 = ((factor2 * spectrum_mz_chunk2) - mzsize) * imagesize * numbytes
        # overhead2 = (spectrum_mz_chunk2 - math.ceil(float(mzsize) % float(spectrum_mz_chunk2))) * imagesize *numbytes
        # overhead2 = ((spectrum_mz_chunk2 *  math.ceil (float(mzsize) / float(spectrum_mz_chunk2))) - \
        #              (spectrum_mz_chunk2 *  math.floor(float(mzsize) / float(spectrum_mz_chunk2)))) * \
        #            imagesize * numbytes
        if overhead1 < overhead2:
            spectrum_mz_chunk = spectrum_mz_chunk1
            spectrum_chunk_overhead = overhead1
        else:
            spectrum_mz_chunk = spectrum_mz_chunk2
            spectrum_chunk_overhead = overhead2
        spectrum_chunk = (spectrum_x_chunk, spectrum_y_chunk, spectrum_mz_chunk)
        if print_results:

            print "Spectrum selection chunking: " + str(spectrum_chunk)
            print "     - Ideal for selection of full spectra."
            print "     - Overhead: " + str(spectrum_chunk_overhead) + \
                  " Byte (" + str(int(math.ceil(spectrum_chunk_overhead / (1024. * 1024.)))) + " MB)"

        # Define a chunking that would be good for images
        slice_chunk_x = xsize
        slice_chunk_y = ysize
        slice_chunk_mz = 1
        chunk_size = slice_chunk_x * slice_chunk_y * slice_chunk_mz
        if math.ceil(float(chunk_size) / float(suggested_num_values)) > max_factor:
            slice_chunk_x = int(math.ceil(xsize / 2.))
            slice_chunk_y = int(math.ceil(ysize / 2.))
        slice_chunk = (slice_chunk_x, slice_chunk_y, slice_chunk_mz)
        # math.ceil(float(xsize) % float(slice_chunk_x))
        slice_chunk_overhead_x = (slice_chunk_x * math.ceil(float(xsize) / float(slice_chunk_x))) - xsize
        # math.ceil(float(ysize) % float(slice_chunk_y))
        slice_chunk_overhead_y = (slice_chunk_y * math.ceil(float(ysize) / float(slice_chunk_y))) - ysize
        slice_chunk_overhead = ((slice_chunk_overhead_x * ysize) + (slice_chunk_overhead_y * xsize)
                                - (slice_chunk_overhead_x * slice_chunk_overhead_y)) * mzsize * numbytes
        if print_results:

            print "Slice selection chunking: " + str(slice_chunk)
            print "     - Ideal for selection of full image slices."
            print "     - Overhead: " + str(slice_chunk_overhead) + \
                  " Byte (" + str(int(math.ceil(slice_chunk_overhead / (1024. * 1024.)))) + " MB)"

        # Define a balanced chunkgin
        balanced_chunk_x = 4
        balanced_chunk_y = 4
        balanced_chunk_mz = 2048
        balanced_chunk = (balanced_chunk_x, balanced_chunk_y, balanced_chunk_mz)
        if print_results:

            print "Balanced chunking: " + str(balanced_chunk)
            print "     - This chunking tries to compromise between selection of slices and spectra."

        return spectrum_chunk, slice_chunk, balanced_chunk

    @staticmethod
    def write_data(input_file, data, data_io_option="spectrum", chunk_shape=None, write_progress=True):
        """Helper function used to implement different data write options.

            :param input_file: The input data file
            :param data: The output dataset (either an h5py dataset or omsi_file_msidata object.
            :param data_io_option: String indicating the data write method to be used. One of:

                * ``spectrum``: Write the data one spectrum at a time
                * ``all`` : Write the complete dataset at once.
                * ``chunk`` : Write the data one chunk at a time.

            :param chunk_shape: The chunking used by the data. Needed to decide how the data should \
                                be written when a chunk-aligned write is requested.
            :param write_progress: Write progress in % to standard out while data is being written.
            :type write_progress: bool

        """
        if data_io_option == "spectrum" or (data_io_option == "chunk" and (chunk_shape is None)):
            num_spectra = float(input_file.shape[0] * input_file.shape[1])
            if isinstance(input_file, file_reader_base.file_reader_base):
                spectrum_index = 0
                for spectrum in input_file.spectrum_iter():
                    xindex = spectrum[0][0]
                    yindex = spectrum[0][1]
                    vals = spectrum[1]
                    data[xindex, yindex, :] = vals
                    spectrum_index += 1
                    if write_progress:
                        try:
                            sys.stdout.write("[" + str(int(100. * spectrum_index / num_spectra)) + "%]" + "\r")
                            sys.stdout.flush()
                        except ValueError:
                            write_progress = False
            else:
                for xindex in xrange(0, input_file.shape[0]):
                    if write_progress:
                        sys.stdout.write("[" + str(int(100. * float(xindex) / float(input_file.shape[0]))) + "%]" + "\r")
                        sys.stdout.flush()
                    for yindex in xrange(0, input_file.shape[1]):
                        # Save the spectrum to the hdf5 file
                        data[xindex, yindex, :] = input_file[xindex, yindex, :]
        elif data_io_option == "all":
            data[:] = input_file[:]
        elif data_io_option == "spectrum_to_image":
            # Determine the I/O settings
            log_helper.info(__name__, "Spectrum-to-image I/O. Write %s blocks at a time" % (chunk_shape, ))
            xdim = input_file.shape[0]
            ydim = input_file.shape[1]
            zdim = input_file.shape[2]
            num_chunks_x = int(math.ceil(float(xdim) / float(chunk_shape[0])))
            num_chunks_y = int(math.ceil(float(ydim) / float(chunk_shape[1])))
            num_chunks_z = int(math.ceil(float(zdim) / float(chunk_shape[2])))
            dtype = input_file.data_type if isinstance(input_file, file_reader_base.file_reader_base) else input_file.dtype

            # Convert the data in image block in z-chunks of size chunk_shape[2]
            for zChunkIndex in xrange(0, num_chunks_z):
                # Determine the start and stop of the current z block
                zstart = zChunkIndex * chunk_shape[2]
                zend = min(zstart + chunk_shape[2], zdim)
                # Reserve space in memory for the current block of full images
                temp_data = np.zeros(shape=(xdim, ydim, zend-zstart), dtype=dtype)
                # Load the images for the block into memory in block of x/y chunks
                for yChunkIndex in xrange(0, num_chunks_y):
                    ystart = yChunkIndex * chunk_shape[1]
                    yend = min(ystart + chunk_shape[1], ydim)
                    for xChunkIndex in xrange(0, num_chunks_x):
                        xstart = xChunkIndex * chunk_shape[0]
                        xend = min(zstart + chunk_shape[0], xdim)
                        temp_data[xstart:xend, ystart:yend, :] = input_file[xstart:xend, ystart:yend, zstart:zend]
                # Write the data for the current z block the file
                data[:, :, zstart:zend] = temp_data[:]
                # Update the progress indicator
                if write_progress:
                    sys.stdout.write("[" + str(int(100. * float(zChunkIndex) / float(num_chunks_z))) + "%]" + "\r")
                    sys.stdout.flush()
                # Clean-up
                if isinstance(data, h5py.Dataset):
                    data.file.flush()
                del temp_data

        elif data_io_option == "chunk":
            xdim = input_file.shape[0]
            ydim = input_file.shape[1]
            zdim = input_file.shape[2]
            num_chunks_x = int(math.ceil(float(xdim) / float(chunk_shape[0])))
            num_chunks_y = int(math.ceil(float(ydim) / float(chunk_shape[1])))
            num_chunks_z = int(math.ceil(float(zdim) / float(chunk_shape[2])))
            num_chunks = num_chunks_x * num_chunks_y * num_chunks_z
            itertest = 0
            for xChunkIndex in xrange(0, num_chunks_x):
                xstart = xChunkIndex * chunk_shape[0]
                xend = min(xstart + chunk_shape[0], xdim)
                for yChunkIndex in xrange(0, num_chunks_y):
                    ystart = yChunkIndex * chunk_shape[1]
                    yend = min(ystart + chunk_shape[1], ydim)
                    for zChunkIndex in xrange(0, num_chunks_z):
                        zstart = zChunkIndex * chunk_shape[2]
                        zend = min(zstart + chunk_shape[2], zdim)
                        # print "Write : "+str(xstart)+" "+str(xend)+"
                        # "+str(ystart)+" "+str(yend)+" "+str(zstart)+" "+str(zend)
                        data[xstart:xend, ystart:yend, zstart:zend] = input_file[
                            xstart:xend, ystart:yend, zstart:zend]
                        itertest += 1
                        if write_progress:
                            sys.stdout.write("[" + str(int(100. * float(itertest) / float(num_chunks))) + "%]" + "\r")
                            sys.stdout.flush()


if __name__ == "__main__":
    main()
