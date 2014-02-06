Converting and Accessing Files
==============================

.. code-block:: none
    
    ssh edison.nersc.gov
    cd /project/projectdirs/openmsi/devel/convert
    source setupEnvironment.csh
    python convertToOMSI.py <infile1> <output HDF5File>

Note if you use the bash shell then use ``setupEnvironment.bash`` instead. 


Making a converted file accesilbe to OpenMSI (Private)
------------------------------------------------------

By default, the conversion script will automatically register new files with the OpenMSI site and
assign them private to a single user if the output file placed in the OpenMSI private data location:

/project/projectdirs/openmsi/omsi_data/private/<username>/<filename>

The username in the path will determine the user the file is assigned to. The filename will also be
the name used in the listing on the site. E.g.:

.. code-block:: none

    python convertToOMSI.py <infile1> /project/projectdirs/openmsi/omsi_data/private/<username>/<filename>
    
Changing file permissions:
--------------------------

Using the OpenMSI website the owner of the file can assign permissions to files online:

https://openmsi.nersc.gov/openmsi/resources/filemanager?file=<filename>


    


``convertToOMSI``: Usage and Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: In order to view a current, complete list of conversions options use:

.. code-block:: bash

    python convertToOMSI.py --help

.. code-block:: none

    USAGE: Call "convertToOMSI [options] imgBaseFile1 imgBaseFile2 ... imgBaseFileN HDF5File" 

    This converter script takes the basename (i.e., path+basefilename) of a single
    or multiple MSI files as input and converts them to HDF5. Each MSI file is
    stored as a separate experiment in the output HDF5 file. If an input file
    defines multiple regions, then those regions can either be stored as separate
    datasets of the same experiment and/or merged to a single MSI dataset.
    Using the various paramter settings described below, one can define how the
    conversion should be performed, how the data should be stored in HDF5, and
    indicate which analyses should be exectued.

    ===HELPER OPTIONS===
    --suggest-chunking : Iterate over all given input files and suggest a chunking strategy.
                         No data is converted when this option is given, i.e., no name for the
                         HDF5File should be given, but only input files should be listed.

    ===ERROR HANDLING OPTIONS===
    --error-handling <options>: Define how errors should be handled. Options are:
                       i)   terminate-and-cleanup (default) : Terminate the conversion, delete the
                               the HDF5 file and do not add the file to the database.
                       ii)  terminate-only,: Leave the generated HDF5 output file in place but  do not
                                add the file to the database.
                       iii) continue-on-error: Ignore errors if possible and continue, even if this
                                means that some data may be missing from the output.

    ===INPUT DATA OPTIONS===

    Default input data options: --format auto --regions split+merge
    --format <option>: Define which file format is used as input. By default the program tries to
               automatically determine the input format. This option can be used to indicate
               the format explicitly to in case the auto option fails. Available options are:
              ['img', 'bruckerflex', 'auto']
    --regions <option>: Some file formats (e.g., brucker) allow multiple regions to be imaged and stored
               in a single file. This option allows one to specify how these regions should be
               treated during file conversion. E.g., one may want to store i) each region as a 
               separate dataset in the output file (--regions split), ii) all regions combined 
               in a single dataset (--regions merge), or both (--regions split+merge)
               Available options are:
              ['split', 'merge', 'split+merge']

    ===FILE WRITE OPTIONS===

    ---FILE WRITE OPTIONS: Chunking---

    Default HDF5 Chunking options: Enabled by default using --auto-chunking :
    --auto-chunking : Automatically decide which chunking should be used. This option
                    automatically generates two copies of the data, one with a chunking
                    optimized for selection of spectra and another one optimized for 
                    selection of ion image slices. All --chunking, --no-chunking, and
                    --optimized-chunking options are ignored if this paramter is given
    --chunking <x y z> : Use chunking when writing the HDF5 file. (DEFAULT, x=4,y=4,z=2048)
    --no-chunking : Disable chunking when writing the HDF5 file. Use in combination with
                    --no-compression since compression depends on chunking and will enable
                    it if compression is used.
    --optimized-chunking <x y z> : Use this option to generate additional copies of the data
                    with different chunked data layouts. Generating multiple copies of the
                    data with different chunked data layouts can be help accelerate selective 
                    data read opeations. (DEFAULT OFF). We recommend a spectra-aligned chunking
                    for the raw data, e.g., '--chunking 1 1 32768' and an image-aligned chunked
                    secondary copy of the data, e.g., '--optimzied-chunking 20 20 100'.

    ---FILE WRITE OPTIONS: Compression---
    HDF5 Compression: Default ON using (gzip , 4):
    --compression: Enable compression using (gzip,4). NOTE: Compression requires the use of chunking.
    --no-compression: Disable the use of compression.

    ===I/O OPTIONS===
    --io <option>: Available options are: 
                 i) all : Read the full data in memory and write it at once
                 ii) spectrum : Read one spectrum at a time and write it to the file. 
                 iii) chunk : Read one chunk at a time and write it to the file.

    ===DATABSE OPTIONS===

    These options control whether the generated output file should be added to a server database
    to manage web file access permissions
    Default options are: --add-to-db --db-server http://openmsi.nersc.gov
    --add-to-db : Add the output HDF5 file to the database.
    --no-add-to-db : Disable adding the file to the database.
    --db-server : Specify the online server where the file should be registers. Default is
                  http://openmsi.nersc.gov 
    --owner : Name of the user that should be assigned as owner. By default the owner is
              determined automatically based on the file path.

    ===ANALYSIS OPTIONS===

    NMF: Default ON: (nc=20, timout=600, niter=2000, tolerance=0.0001, raw=False)
    --nmf : Compute the nmf for all the input data files and store the results in the
            HDF5 file. NOTE: If global peak-finding (fpg) is performed, then
            nmf will be performed on the peak-cube, otherwise on the raw data
    --no-nmf: Disable the execution of nmf
    --nmf-nc <nummber>: Number of components to be computed by the NMF. (default nc=20)
    --nmf-timeout <number>: Maximum time in seconds to be used for computing the NMF. (default timeout=600)
    --nmf-niter <number>: Number of iterations (minimum is 2)(default niter=2000)
    --nmf-tolerance <number>: Tolerance value for a relative stopping condition. (default tolerance=0.0001)
    --nmf-raw <number>: Force execution of the NMF on the raw data. By default the results from
                the global peak finding (--fpg) are used to compute the NMF.

    Global Peak Finding: Default ON:
    --fpg : Compute the global peak finding for all input data files and save results
               in the HDF5 file (DEFAULT)
    --no-fpg: Disable the global peak finding

    Global Peak Finding: Default OFF:
    --fpl : Compute the local peak finding for all input data files and save results
            in the HDF5 file
    --no-fpl: Disable the local peak finding (DEFAULT)

    ---OTHER OPTIONS---

    Generate Thunmbnail image: Default ON:
    --thumbnail: Generate thumbnail image for the file based on, in order of avalability:
                 * The frist three components of the NMF
                 * The three most intense peaks from the global peak finding (fpg)
                 * The three most intense peaks in the raw data that are at least 1 percent
                   of the total m/z range apart.
    --no-thumbnail: Do not generate a thumbnail image.
