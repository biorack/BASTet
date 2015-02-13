Developing a new Analysis for OpenMSI
=====================================

Overview:
---------

The OpenMSI Toolkit includes a basic template for implementing new analyses as part of OpenMSI. The template is located in ``omsi.templates.omsi_analysis_template.py``. The template provides step-by-step instructions on how to implement a new analysis. Simply search top-to-bottom for EDIT_ME markers to find locations that need to be edited and what changes need to be made.

The implementation of a new analysis is divided into three main steps. We here provide a brief overview of these steps. A detailed walk-through the required implementation is provided in the following sections.

**Step 1) Basic integration of your analysis with OpenMSI (Required)**

 The basic integration is simple and should requires only minimal additional effort. The basic integration with the OpenMSI provides full integration of the analysis with the OpenMSI file format and API and provide full support for OpenMSI's data provenance capabilities. It also provides basic integration of the analysis with the website, in that a user will be able to browse the analysis in the online file browser. The basic integration automatically provides full support for the ``qmetadata`` and ``qcube`` URL data access patterns. The basic integration provides limited support for the ``qslice``, ``qspectrum``, and ``qmz`` patterns, in that it automatically exposes all dependencies of the analysis that support these patterns but it does not expose the data of the analysis itself. This is part of step 2.

**Step 2)  Integrating your analysis with the OpenMSI web-based viewer (Recommended)**

Once the basic integration is complete, you may want integrate your analysis fully with the OpenMSI online viewer, in order to make your analysis easily accesible to the OpenMSI user community. This step requires the implementation of the ``qslice``, ``qspectrum``, and ``qmz`` URL patterns for the analysis. This step completes the integration with the OpenMSI framework itself.

**Step 3) Making your analysis self-sufficient (Recommended)**

This step makes your analysis "self-sufficient" in that it allows you to execute your analysis from the command-line.


Some important features of ``omsi_analysis_base``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``omsi.analysis.omsi_analysis_base`` is the base class for all omsi analysis functionality. The class provides a large set of functionality designed to facilitate i) storage of analysis data in the omsi HDF5 file format and ii) intergration of new analysis capabilities with the OpenMSI web API and the OpenMSI web-based viewer (see Viewer functions below for details), and iii) support for data provenance.

Slicing
"""""""
``omsi_analysis_base`` implements basic slicing to access data stored in the main member variables. By default the data is retrieved from __data_list by the __getitem__(key) function, which implements the [..] operator, i.e., the functions returns __data_list[key]['data']. The key is a string indicating the name of the paramter to be retrieved. If the key is not found in the __data_list then the function will try to retrieve the data from __parameter_list instead. By adding "parameter/key" or "dependency/key" one may also explicitly retrieve values from the __parameter_list and __dependency_list.

Important Member Variables
""""""""""""""""""""""""""

* ``analysis_identifier`` defines the name for the analysis used as key in search operations.
* ``__data_list`` defines a dictonary of ``omsi_analysis_data`` to be written to the HDF5 file. Derived classes need to add all data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See ``omsi.analysis.omsi_analysis_data`` for details.
* ``__parameter_list``  defines a dictonary of ``omsi_analysis_data`` to be written to the HDF5 file. Derived classes need to add all parameter data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See omsi.analysis.omsi_analysis_data for details.
* ``__dependency_list`` defines a dictonary of ``omsi_dependency`` to be written to the HDF5 file. Derived classes need to add all dependencies data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See ``omsi.analysis.omsi_analysis_data`` for details.

I/O functions
"""""""""""""

These functions can be optionally overwritten to control how the analysis data should be written/read from the omsi HDF5 file. Default implementations are provided here, which should be sufficient for most cases.

* ``write_to_omsi_file``: The default implementation is empty as the default data write is  managed by the omsi_file_experiment.create_analysis() function.  Overwrite this function, in case that the analysis needs to write data to the HDF5 omsi file beyond what the defualt omsi data API does.

* ``read_from_omsi_file``: The default implementation tries to reconstruct the original data as far  as possible, however, in particular in case that a custom write_to_omsi_file            funtion has been implemented, the default implementation may not be sufficien. The default implementation reconstructs: i) analysis_identifier and reads all custom data into ii)__data_list. Note, an error will be raised in case that the analysis type specified in the HDF5 file does not match the analysis type specified by get_analysis_type(). This function can be optionally overwritten to implement a custom data read.

Web API Functions
"""""""""""""""""

Several convenient functions are used to allow the OpenMSI online viewer to interact with the analysis and to visualize it. The default implementations provided here simply indicate that the analysis does not support the data access operations required by the online viewer. Overwrite these functions in the derived analysis classes in order to interface them with the viewer. All viewer-related functions start with ``v\_...`` .

NOTE: the default implementation of the viewer functions defined in ``omsi_analysis_base`` are designed to take care of the common requirement for providing viewer access to data from all depencies of an analysis. In many cases, the default implementation is often sill called at the end of custom viewer functions.

NOTE: The viewer functions typically support a viewerOption parameter. viewerOption=0 is expected to refer to the analysis itself.

* ``v_qslice``: Retrieve/compute data slices as requested via qslice URL requests. The corrsponding view of the DJANGO data access server already translates all input parameters and takes care of generating images/plots if needed. This function is only responsible for retrieving the data.
* ``v_qspectrum``: Retrieve/compute spectra as requested via qspectrum URL requests. The corrsponding view of the DJANGO data access server already translates all input parameters and takes care of generating images/plots if needed. This function is only responsible for retrieving the data.
* ``v_qmz``: Define the m/z axes for image slices and spectra as requested by qspectrum URL requests.
* ``v_qspectrum_viewer_options``: Define a list of strings, describing the different viewer options available for the analysis for qspectrum requests (i.e., ``v_qspectrum``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.
* ``v_qslice_viewer_options``: Define a list of strings, describing the different viewer options available for the analysis for qslice requests (i.e., ``v_qslice``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.


Executing, saving, and restoring an analysis object
"""""""""""""""""""""""""""""""""""""""""""""""""""

Any analysis based on the infrastructure provided by ``omsi_analysis_base`` is fully integrated with OpenMSI file API provided by``omsi.dataformat.omsi_file``. This means the analysis can be directly  API saved to an OMSI HDF5 file and  the saved analysis can be restored from file. In OMSI files, analyses are generally associated with experiments, so that we use the ``omsi.dataformat.omsi_file.omsi_file_experiment`` API here.

.. code-block:: python
    :linenos:
    :emphasize-lines: 12, 21,22,23,24

    # Open the MSI file and get the desired experiment
    from omsi.dataformat.omsi_file import *
    f = omsi_file( filename, 'a' )
    e = f.get_experiment(0)

    # Execute the analysis
    d = e.get_msidata(0)
    a = omsi_myanalysis()
    a.execute(msidata=d, integration_width=10, msidata_dependency=d)

    # Save the analysis object.
    analysis_object , analysis_index = exp.create_analysis( a )

    # This single line is sufficient to store the complete analysis to the omsi file.
    # By default the call will block until the write is complete. Setting the
    # parameter flushIO=False enables buffered write, so that the call will
    # return once all data write operations have been scheduled. Here we get
    # an omsi.dataformat.omsi_file.omsi_file_analysis
    # object for management of the data stored in HDF5 and the integer index of the analysis.

    # To restore an analysis from file, i.e., read all the analysis data from file
    # and store it in a corresponding analysis object we can do the following.
    # Similar to the read_from_omsi_file(...) function of omsi_analysis_base
    # mentioned below, we can decide via parameter settings of the function,
    # which portions of the analysis should be loaded into memory
    a2 = analysis_object.restore_analysis()

    # If we want to now re-execute the same analysis we can simply call
    a2.execute()

    # If we want to rerun the analysis but change one or more parameter settings,
    # then we can simply change those parameters when calling the execute function
    d2 = e.get_msidata(1)   # Get another MSI dataset
    a2.execute(msidata=d2)  # Execute the analysis on the new MSI dataset

    # The omsi_file_analysis class also provides a convenient function that allows us
    # to recreate, i.e., restore and run the analysis, in a single function call
    a3 = analysis_object.recreate_analysis()

    # The recreate_analysis(...) function supports additional keyword arguments
    # which will be passed to the execute(...) call of the analysis, so that we
    # can change parameter settings for the analysis also when using the
    # recreate analysis call.

    # If we know the type of analysis object (which we can also get from file), then we
    # naturally also restore the analysis from file ourselves via
    a4 = omsi_myanalysis().read_from_omsi_file(analysisGroup=analysis_object, \
                                               load_data=True, \
                                               load_parameters=True,\
                                               load_runtime_data=True, \
                                               dependencies_omsi_format=True )
    # By setting load_data and/or load_parameters to False, we create h5py instead of
    # numpy objects, avoiding the actual load of the data. CAUTION: To avoid the accidental
    # overwrite of data we recommend to use load_data and load_parameters as False only
    # when the file has been opened in read-only mode 'r'.

    # Rerunning the same analysis again
    a4.execute()





Integrating a new Analysis using the OpenMSI Analysis Template
--------------------------------------------------------------

Step 1) Basic integration
^^^^^^^^^^^^^^^^^^^^^^^^^

The simple steps outlined below provide you now with full integration of your analysis with the OpenMSI file format and API and full support for OpenMSI's data provenance capabilities. It also provides basic integration of your analysis with the OpenMSI website, in that a user will be able to browse your analysis in the online file browser. The basic integration also automatically provides full support for the ``qmetadata`` and ``qcube`` URL data access patterns, so that you can start to program against your analysis remotely. The basic integration provides limited support for the ``qslice``, ``qspectrum``, and ``qmz`` patterns, in that it automatically exposes all dependencies of the analysis that support these patterns but it does not expose the data of your analysis itself. This is part of step 2. Once you have completed the basic integation yout final analysis code should look something like this:

.. code-block:: python
    :linenos:
    :emphasize-lines: 6,7,27,28,29

    class omsi_mypeakfinder(omsi_analysis_base) :

        def __init__(self, name_key="undefined"):
            """Initalize the basic data members"""

            super(omsi_mypeakfinder,self).__init__()
            self.parameter_names = [ 'msidata' , 'mzdata', 'integration_width', 'peakheight' ]
            self.data_names = [ 'peak_cube' , 'peak_mz' ]
            self.analysis_identifier = name_key

        def execute_analysis(self) :
            """..."""

            #Set default parameter values for optional parameters
            if not self['integration_width'] :
                self['integration_width']=10

            #Copy parameters to local variables for convenience
            msidata = self['msidata']
            mzdata = self['mzdata']
            integration_width = self['integration_width'][0]
            peakheight = self['peakheight'][0]

            #Implementation of my peakfinding algorithm

            ...

            #Save output data
            self['peak_cube'] = peakCube
            self['peak_mz']   = peakMZ

        ...


1.1 Creating a new analysis skeleton
""""""""""""""""""""""""""""""""""""

- Copy the analysis template to the appropriate location where your analysis should live. Any new anlysis should be located in a submodule of the ``omsi.analysis.`` module. E.g., if you implement a new peak finding algorithm, it should be placed in omsi/analysis/findpeaks. For example:

.. code-block:: none

    cp omsi/templates/omsi_analysis_template.py openmsi-tk/omsi/analysis/findpeaks/omsi_mypeakfinder.py


- Replace all occurrences of ``omsi_analysis_template`` in the file with the name of your analysis class, e.g, omsi_mypeakfinder. You can do this easily using "Replace All" feature of most text editors.  or on most Unix systems  (e.g, Linux or MacOS) on the commandline via:

.. code-block:: none

    cd openmsi-tk/omsi/analysis/findpeaks
    sed -i.bak 's/omsi_analysis_template/omsi_mypeakfinder/' omsi_mypeakfinder.py
    rm omsi_mypeakfinder.py.bak

- Add your analysis to the ``__init__.py`` file of the python module where your analysis lives. In the ``__init__.py`` file you need to add the name of your analysis class to the ``all__`` list and add a an import of your class, e.g,  ``from omsi_mypeakfinder import *`` . For example:

.. code-block:: python
    :linenos:
    :emphasize-lines: 1,7

    all__ = [ "omsi_mypeakfinder",  "omsi_findpeaks_global" , ...]
    from omsi_findpeaks_global import *
    from omsi_findpeaks_local import *
    from omsi_lpf import *
    from omsi_npg import *
    from omsi_peakcube import *
    from omsi_mypeakfinder import *

1.2 Specifying analysis inputs and outputs
""""""""""""""""""""""""""""""""""""""""""


In the ``__init__`` function specify the names of the input parameters of your analysis as well as the names of the ouput data generated by your analysis. E.g.,

.. code-block:: python
    :linenos:
    :emphasize-lines: 5,6

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""

        super(omsi_mypeakfinder,self).__init__()
        self.parameter_names = [ 'msidata' , 'mzdata', 'integration_width', 'peakheight' ]
        self.data_names = [ 'peak_cube' , 'peak_mz' ]
        self.analysis_identifier = name_key

1.3: Implementing the ``execute_analysis`` function
"""""""""""""""""""""""""""""""""""""""""""""""""""

**1.3.1** Document your execute_analysis function. OpenMSI typically uses Sphynx notation in the doc-string.

.. code-block:: python
    :linenos:
    :emphasize-lines: 2-11

    def execute_analysis(self) :
        """This analysis computes global peaks in MSI data...

           :param msidata: The input MSI data
           :param mzdata: The mz axis information for the MSI dataset
           :param integration_width: The integration width to be used
           :param peakheight: Minimum peak height threshold.

           :returns: The funtion generated a 'peak_cube' dataset of all
                     global peaks and 'peak_mz' with the m/z values.

        """

**1.3.2** Define default values for any optional input parameters and implement your analysis. s

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,5

    def execute_analysis(self) :
        """..."""

        if not self['integration_width'] :
            self['integration_width']=10



**1.3.3**  Implement your analysis. For conveniece it is often useful to assign the your parameters to local variables, although, this is by no means required. Note, all values are stored as 1D+ numpy arrays. I.e., for scalar parameters we need to access the [0] value. E.g.:

.. code-block:: python

    integration_width = self['integration_width'][0]

**1.3.4** Save the results of your analysis to the respective output variables. This will allow users to conveniently access your results and it enables the OpenMSI file API to save your results to file. We here automatically convert single scalars to 1D numpy arrays to ensure consistency. Although, the data write function can handle a large rangeof python built_in types by automaticaly converting them to numpy for storage in HDF5, we generally recommend to convert use numpy directly here to save your data.

.. code-block:: python

    self['peak_cube'] = peakCube
    self['peak_mz']   = peakMZ

With this you have now completed the basic integration of your analysis with the OpenMSI framework.



Optional: Custom data save
^^^^^^^^^^^^^^^^^^^^^^^^^^

In most cases the default data save and restore functions should be sufficient. However, the ``omsi_analysis_base`` API also supports implementation of custom HDF5 write. To extend the existing data write code, simple implement the following function provided by ``omsi_analysis_base`` .

.. code-block:: python
    :linenos:
    :emphasize-lines: 1

    def write_to_omsi_file(self , analysisGroup) :
        """This function can be optionally overwritten to implement a custom data write
           function for the analysis to be used by the omsi_file API.

           Note, this function should be used only to add additional data to the analysis
           group. The data that is written by default is typically still written by the
           omsi_file_experiment.create_analysis() function, i.e., the following data is
           wirtten by default: i) analysis_identifier ,ii) get_analysis_type,
           iii)__data_list, iv) __parameter_list , v) __dependency_list. Since the
           omsi_file.experiment.create_analysis() functions takes care of setting up the
           basic structure of the analysis storage (included the subgroubs for storing
           parameters and data dependencies) this setup can generally be assumed to exist
           before this function is called. This function is called automatically at the
           end omsi_file.experiment.create_analysis() (i.e, actually
           omsi_file_analysis.__populate_analysis__(..)) so that this function does not need
           to be called explicitly.

           Keyword Arguments:

           :param analysisGroup: The omsi_file_analysis object of the group for the
                                 analysis that can be used for writing.

           """
        pass

Optional: Custom analysis restore
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similarly in order implement custom data restore behavior we can overwrite the default implementation of ``omsi_analysis_base.read_from_omsi_file)`` . In this case one will usually call the default implementation via ``super(omsi_myanalysis,self).read_from_omsi_file(...)`` first and then add any additional behavior.

Step 2) Integrating the Analysis with the OpenMSI Web API:
----------------------------------------------------------

Once the analysis is stored in the OMSI file format, integration with ``qmetadata`` and ``qcube`` calls of the web API is automatic. The ``qmetadata`` and ``qcube`` functions provide general purpose access to the data so that we can immediatly start to program against our analysis.

Some applications---such as the OpenMSI web-based viewer---utilize the simplified, special data access patterns ``qslice``, ``qspectrum``, and ``qmz`` in order to interact with the data. The default implementation of these function available in ``omsi.analysis.omsi_analysis_base`` exposes the data from all depencdencies of the analysis that support these patterns. For full integration with the web API, however, we need to implement this functionality in our analysis class. The ``qmz`` pattern in particular is relevant to both the ``qslice`` and ``qspectrum`` pattern and should be always implemented as soon as one of the two patterns is defined.

2.1 Implementing the ``qslice`` pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,5,18,20,21,22,23,30,31,32,33,38,39,42,44,45

    class omsi_myanalysis(omsi_analysis_base) :
        ...

        @classmethod
        def v_qslice(cls , anaObj , z , viewer_option=0) :
            """Implement support for qslice URL requests for the viewer

               anaObj: The omsi_file_analysis object for which slicing should be performed.
               z: Selection string indicting which z values should be selected.
               viewer_option: An analysis can provide different default viewer behaviors
                             for how slice operation should be performed on the data.
                             This is a simple integer indicating which option is used.

               :returns: numpy array with the data to be displayed in the image slice
                         viewer. Slicing will be performed typically like [:,:,zmin:zmax].

            """
            from omsi.shared.omsi_data_selection import *
            #Implement custom analysis viewer options
            if viewer_option == 0 :
                dataset =  anaObj[ 'labels' ] #We assume labels was a 3D image cube of labels
                zselect = selection_string_to_object(z)
                return dataset[ : , :, zselect ]

            #Expose recursively the slice options for any data dependencies. This is useful
            #to allow one to trace back data and generate comlex visualizations involving
            #multiple different data sources that have some from of dependency in that they
            #led to the generation of this anlaysis. This behavior is already provided by
            #the default implementation of this function ins omsi_analysis_base.
            elif viewer_option >= 0 :
                #Note, the base class does not know out out viewer_options so we need to adjust
                #the vieweOption accordingly by substracting the number of our custom options.
                return super(omsi_myanalysis,cls).v_qslice( anaObj , z, viewer_option-1)
            #Invalid viewer_option
            else :
                return None

         @classmethod
            def v_qslice_viewer_options(cls , anaObj ) :
                """Define which viewer_options are supported for qspectrum URL's"""
                #Get the options for all data dependencies
                dependent_options = super(omsi_findpeaks_global,cls).v_qslice_viewer_options(anaObj)
                #Define our custom viewer options
                re = ["Labels"] + dependent_options
                return re



NOTE: We here convert the selection string to a python selection (i.e., a list, slice, or integer) object using the ``omsi.shared.omsi_data_selection.check_selection_string(...)`` . This has the advantage that we can use the given selection directly in our code and avoids the use of a potentially dangerous ``eval`` , e.g., ``return eval("dataset[:,:, %s]" %(z,))`` . While we can also check the validity of the  seletion string  using  ``omsi.shared.omsi_data_selection.check_selection_string(...)`` , it is recommened to convert the string to a valid python selection to avoid possible attacks.


2.2 Implementing the ``qspectrum`` pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:
    :emphasize-lines: 3,4,38,39,40,41,47,48,49,54,55,56,58,64,67,72,75,77,78

    class omsi_myanalysis(omsi_analysis_base) :
        ...
        @classmethod
        def v_qspectrum( cls, anaObj , x, y , viewer_option=0) :
            """Implement support for qspectrum URL requests for the viewer.

               anaObj: The omsi_file_analysis object for which slicing should be performed
               x: x selection string
               y: y selection string
               viewer_option: If multiple default viewer behaviors are available for a given
                             analysis then this option is used to switch between them.

               :returns: The following two elemnts are expected to be returned by this function :

                    1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The spectrum axis,
                    e.g., mass (m/z), must be the last axis. For index selection x=1,y=1 a 1D array
                    is usually expected. For indexList selections x=[0]&y=[1] usually a 2D array
                    is expected. For range selections x=0:1&y=1:2 we one usually expect a 3D array.
                    This behavior is consistent with numpy and h5py.

                    2) None in case that the spectra axis returned by v_qmz are valid for the
                    returned spectrum. Otherwise, return a 1D numpy array with the m/z values
                    for the spectrum (i.e., if custom m/z values are needed for interpretation
                    of the returned spectrum).This may be needed, e.g., in cases where a
                    per-spectrum peak analysis is performed and the peaks for each spectrum
                    appear at different m/z values.

                Developer Note: h5py currently supports only a single index list. If the user provides
                an index-list for both x and y, then we need to construct the proper merged list and
                load the data manually, or, if the data is small enough, one can load the full data
                into a numpy array which supports mulitple lists in the selection. This, however, is
                only recommended for small datasets.

            """

            data = None
            customMZ = None
            if viewer_option == 0 :
                from omsi.shared.omsi_data_selection import *
                dataset =  anaObj[ 'labels' ]
                if (check_selection_string(x) == selection_type['indexlist']) and  \
                   (check_selection_string(y) == selection_type['indexlist']) :
                    #Assuming that the data is small enough, we can handle the multiple list
                    #selection case here just by loading the full data use numpy to do the
                    #subselection. Note, this version would work for all selection types but
                    #we would like to avoid loading the full data if we don't have to.
                    xselect = selection_string_to_object(x)
                    yselect = selection_string_to_object(y)
                    data = dataset[:][xselect,yselect,:]
                    #Since we alredy confirmed that both selection strings are index lists we could
                    #also just do an eval as follows.
                    #data = eval("dataset[:][%s,%s, :]" %(x,y))
                else :
                    xselect = selection_string_to_object(x)
                    yselect = selection_string_to_object(y)
                    data = dataset[xselect,yselect,:]
                #Return the spectra and indicate that no customMZ data values (i.e. None) are needed
                return data, None
            #Expose recursively the slice options for any data dependencies. This is useful
            #to allow one to trace back data and generate comlex visualizations involving
            #multiple different data sources that have some from of dependency in that they
            #led to the generation of this anlaysis. This behavior is already provided by
            #the default implementation of this function ins omsi_analysis_base.
            elif viewer_option >= 0 :
                #Note, the base class does not know out out viewer_options so we need to adjust
                #the vieweOption accordingly by substracting the number of our custom options.
                return super(omsi_findpeaks_global,cls).v_qspectrum( anaObj , x , y, viewer_option-1)

            return data , customMZ

        @classmethod
        def v_qspectrum_viewer_options(cls , anaObj ) :
            """Define which viewer_options are supported for qspectrum URL's"""
            #Get the options for all data dependencies
            dependent_options = super(omsi_findpeaks_global,cls).v_qspectrum_viewer_options(anaObj)
            #Define our custom viewer options
            re = ["Labels"] + dependent_options
            return re

2.3 Implementing the ``qmz`` pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:
    :emphasize-lines: 5,6,41,42,43,45,46,47,48,50,51,52,53,54,55,57,58,59,60,61,62,64

    class omsi_myanalysis(omsi_analysis_base) :
        ...

        @classmethod
            def v_qmz(cls, anaObj, qslice_viewer_option=0, qspectrum_viewer_option=0) :
                """Implement support for qmz URL requests for the viewer.

                    Get the mz axes for the analysis

                    anaObj: The omsi_file_analysis object for which slicing should be performed.
                    qslice_viewer_option: If multiple default viewer behaviors are available for
                                a given analysis then this option is used to switch between them
                                for the qslice URL pattern.
                    qspectrum_viewer_option: If multiple default viewer behaviors are available
                                for a given analysis then this option is used to switch between
                                them for the qspectrum URL pattern.

                    :returns: The following four arrays are returned by the analysis:

                      - mzSpectra : 1D numpy array with the static mz values for the spectra.
                      - labelSpectra : String with lable for the spectral mz axis
                      - mzSlice : 1D numpy array of the static mz values for the slices or
                                  None if identical to the mzSpectra array.
                      - labelSlice : String with label for the slice mz axis or None if
                                     identical to labelSpectra.

                     Developer Note: Here we need to handle the different possible combinations
                     for the differnent viewer_option patterns. It is in general safe to populate
                     mzSlice and lableSlice also if they are identical with the spectrum settings,
                     however, this potentially has a significant overhead when the data is transfered
                     via a slow network connection, this is why we allow those values to be None
                     in case that they are identical.

                """
                #The four values to be returned
                mzSpectra =  None
                labelSpectra = None
                mzSlice = None
                labelSlice = None
                #Both qslice and qspectrum here point to our custom analysis
                if qspectrum_viewer_option == 0 and qslice_viewer_option==0: #Loadings
                    mzSpectra =  anaObj[ 'labels' ][:]
                    labelSpectra = "Labels"
                #Both viewer_options point to a data dependency
                elif qspectrum_viewer_option > 0 and qslice_viewer_option>0 :
                    mzSpectra, labelSpectra, mzSlice, labelSlice = \
                           super(omsi_findpeaks_global,cls).v_qmz( anaObj, \
                                 qslice_viewer_option-1 , qspectrum_viewer_option-1)
                #Only the a qlice options point to a data dependency
                elif qspectrum_viewer_option == 0 and qslice_viewer_option>0 :
                    mzSpectra =  anaObj[ 'peak_mz' ][:]
                    labelSpectra = "m/z"
                    tempA, tempB, mzSlice, labelSlice = \
                            super(omsi_findpeaks_global,cls).v_qmz( anaObj, \
                                  qslice_viewer_option-1 , 0)
                #Onlye the qlice option points to a data dependency
                elif qspectrum_viewer_option > 0 and qslice_viewer_option==0 :
                    mzSlice =  anaObj[ 'peak_mz' ][:]
                    labelSlice = "m/z"
                    mzSpectra, labelSpectra, tempA, tempB = \
                            super(omsi_findpeaks_global,cls).v_qmz( anaObj, \
                                  0 , qspectrum_viewer_option-1)

                return mzSpectra, labelSpectra, mzSlice, labelSlice
