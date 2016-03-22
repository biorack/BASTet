.. _converting-files:

Converting Files at NERSC and Making them Accessible
====================================================

To convert a mass spectrometry imaging file, e.g., in img or bruckerflex format, to HDF5 do the following:

.. code-block:: none

    ssh cori.nersc.gov
    cd /project/projectdirs/openmsi/devel/convert
    source setupEnvironment.csh
    python convertToOMSI.py <infile1> <output HDF5File>

Note if you use the bash shell then use ``setupEnvironment.bash`` instead.

Note, if you want to use the output file from openmsi.nersc.gov then the output file path should be:

.. code-block:: none

    /project/projectdirs/openmsi/omsi_data_private/<username>/<filename>

where username is the name of the primary user that owns the file.

Making a converted file accessible to OpenMSI (Private)
-------------------------------------------------------

The conversion script will by default automatically try to register new files with the OpenMSI site and
assign them private to a single user if the output file is placed in the OpenMSI private data location:


.. code-block:: none

    python convertToOMSI.py <infile1> /project/projectdirs/openmsi/omsi_data_private/<username>/<filename>

The username in the path will determine the user the file is assigned to. E.g: The filename will also be
the name used in the listing on the site. In order to generate the HDF5 file without adding to the
database, use the ``--no-add-to-db`` command line option, e.g,:


Changing file permissions:
--------------------------

Using the OpenMSI website the owner of the file can assign permissions to files online:

https://openmsi.nersc.gov/openmsi/resources/filemanager?file=<filename>



``convertToOMSI``: Usage and Options
------------------------------------

NOTE: In order to view a current, complete list of conversions options use:

.. code-block:: bash

    python convertToOMSI.py --help

.. code-block:: none

    python convertToOMSI.py --help
    USAGE: Call "convertToOMSI [options] imgBaseFile1 imgBaseFile2 ... imgBaseFileN HDF5File"

    This converter script takes the basename (i.e., path+basefilename) of a single
    or multiple MSI files as input and converts them to HDF5. Each MSI file is
    stored as a separate experiment in the output HDF5 file. If an input file
    defines multiple regions, then those regions can either be stored as separate
    datasets of the same experiment and/or merged to a single MSI dataset.
    Using the various paramter settings described below, one can define how the
    conversion should be performed, how the data should be stored in HDF5, and
    indicate which analyses should be executed.

    ===HELPER OPTIONS===
    --suggest-chunking : Iterate over all given input files and suggest a chunking strategy.
                         No data is converted when this option is given, i.e., no name for the
                         HDF5File should be given, but only input files should be listed.

    ===ERROR HANDLING OPTIONS===
    --error-handling <options>: Define how errors should be handled. Options are:
                       i)   terminate-and-cleanup (default) : Terminate the conversion, delete the
                               the HDF5 file and do not add the file to the database.
                       ii)  terminate-only, : Leave the generated HDF5 output file in place but  do not
                                add the file to the database.
                       iii) continue-on-error: Ignore errors if possible and continue, even if this
                                means that some data may be missing from the output.
    --email <email1 email2 ...>: Send notification in case of both error or success to the email address.
    --email-success <email1 email2 ...>>: Send notification in case of success to the given email address.
    --email-error <email1 email2 ...>>: Send notification in case of error to the given email address.

    ===INPUT DATA OPTIONS===

    Default input data options: --format auto --regions split+merge
    --format <option>: Define which file format is used as input. By default the program tries to
               automatically determine the input format. This option can be used to indicate
               the format explicitly to in case the auto option fails. Available options are:
              {'imzml_file': <class 'omsi.dataformat.imzml_file.imzml_file'>, 'bruckerflex_file': <class 'omsi.dataformat.bruckerflex_file.bruckerflex_file'>, 'img_file': <class 'omsi.dataformat.img_file.img_file'>, 'mzml_file': <class 'omsi.dataformat.mzml_file.mzml_file'>}
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
    HDF5 Compression: Default ON using (gzip, 4):
    --compression: Enable compression using (gzip,4). NOTE: Compression requires the use of chunking.
    --no-compression: Disable the use of compression.

    ===I/O OPTIONS===
    --io <option>: Available options are: ['chunk', 'spectrum', 'all']
                 i) all : Read the full data in memory and write it at once
                 ii) spectrum : Read one spectrum at a time and write it to the file.
                 iii) chunk : Read one chunk at a time and write it to the file.
                 The io option applies only for the generation of subsequent chunkings
                 and not the initial iteration over the file to generate the first convert.
                 iv) spectrum-to-image: Default option when creating image chunk version from
                 a spectrum-chunk MSI dataset. Read a block of spectra at a time to
                 complete a set of images and then write the block of images at once.
    --io-block-limit <MB>: When using spectrum-to-image io (default when using auto-chunking),
                 what should the maximum block in MB that we load into memory. (Default=2000MB)

    ===DATABSE OPTIONS===

    These options control whether the generated output file should be added to a server database
    to manage web file access permissions
    Default options are: --add-to-db --db-server http://openmsi.nersc.gov
    --add-to-db : Explicitly add the output HDF5 file to the database. This option has no effect
                  if --jobid is set as the file is added through the update of the job status
                  in this case.
    --no-add-to-db : Disable adding the file to the database.
    --db-server : Specify the online server where the file should be registers. Default is
                  http://openmsi.nersc.gov
    --user : Name of the user that should be assigned as user. By default the user is
              determined automatically based on the file path.
    --jobid : ID of the job. If set to 'auto' then the environment variable PBS_JOBID is used.
              NOTE: If job ID is set then we assume that the job has been scheduled via the
              the automated system and that the job is managed. As such the file will be added,
              to the database by updating the job status and NOT by explicitly adding the file.

    ===ANALYSIS OPTIONS===

    NMF: Default ON: (nc=20, timeout=600, niter=2000, tolerance=0.0001, raw=False)
    --nmf : Compute the nmf for all the input data files and store the results in the
            HDF5 file. NOTE: If global peak-finding (fpg) is performed, then
            nmf will be performed on the peak-cube, otherwise on the raw data
    --no-nmf: Disable the execution of nmf
    --nmf-nc <number>: Number of components to be computed by the NMF. (default nc=20)
    --nmf-timeout <number>: Maximum time in seconds to be used for computing the NMF. (default timeout=600)
    --nmf-niter <number>: Number of iterations (minimum is 2)(default niter=2000)
    --nmf-tolerance <number>: Tolerance value for a relative stopping condition. (default tolerance=0.0001)
    --nmf-raw <number>: Force execution of the NMF on the raw data. By default the results from
                the global peak finding (--fpg) are used to compute the NMF.

    Global Peak Finding: Default ON:
    --fpg : Compute the global peak finding for all input data files and save results
               in the HDF5 file (DEFAULT)
    --no-fpg: Disable the global peak finding

    Local Peak Finding: Default OFF:
    --fpl : Compute the local peak finding for all input data files and save results
            in the HDF5 file
    --no-fpl: Disable the local peak finding (DEFAULT)


    TIC normalization:
    --ticnorm : Compute tic normalization
    --no-ticnorm : Disable computation of tic normaltization (DEFAULT)

    ---OTHER OPTIONS---

    Generate Thumbnail image: Default OFF:
    --thumbnail: Generate thumbnail image for the file based on, in order of availability:
                 * The first three components of the NMF
                 * The three most intense peaks from the global peak finding (fpg)
                 * The three most intense peaks in the raw data that are at least 1 percent
                   of the total m/z range apart.
    --no-thumbnail: Do not generate a thumbnail image.

    Generate XDMF header file for output file: Default OFF:
    --xdmf: Write XDMF XML-based header-file for the output HDF5 file.
    --no-xdmf: Do not generate a XDMF XML-based header for the HDF5 file.

    ===Metadata Options===

    NOTE: Input datasets are numbers starting from 0 based on there order on the command line.

    --methods : JSON describing the experimental methods
    --methods# : JSON describing the experimental methods for input file number #
    --instrument : JSON dictionary describing the instrument
    --instrument# : JSON dictionary describing the instrument for input file number #
    --notes : JSON dictionary with additional user notes about the data
    --notes# : JSON dictionary with additional notes for input file number #


