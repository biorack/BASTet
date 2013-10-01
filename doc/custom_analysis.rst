Developing a new Analysis for OpenMSI
=====================================

Basic Structure
---------------

Any new anlysis that should be part of the OpenMSI Toolkit should be located in a submodule of the omsi.analysis. The anlysis class itself should inherit from ``omsi.analysis.omsi_analysis_base``


.. code-block:: python
    :linenos:
    :emphasize-lines: 1,2,3,5,10

    from omsi.analysis.omsi_analysis_base import omsi_analysis_base
    from omsi.analysis.omsi_analysis_data import omsi_analysis_data
    from omsi.shared.omsi_dependency import *

    class omsi_myanalysis(omsi_analysis_base) :
    
        def __init__(self, nameKey="undefined"):
            """Initalize the basic data members and define any settings of omsi_analysi"""
            super(omsi_myanalysis,self).__init__()
            self.analysis_identifier = nameKey
            
        def execute( self ) :
            #Code to execute the analysis 
            
        def main(argv=None) :
            """To make the code immediately useful to others it is highly recommended
               to provide a main function that allows execution of the analysis on a 
               fiven omsi HDF5 data file"""
            import sys
            if argv is None:
                argv = sys.argv
            omsiInFile = argv[1]
            omsiFile =  omsi_file( omsiInFile , 'r' )
            ...

    if __name__ == "__main__":
        main()
   
Some important features of ``omsi_analysis_base``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``omsi.analysis.omsi_analysis_base`` is the base class for all omsi analysis functionality. The class provides a large set of functionality designed to facilitate i) storage of analysis data in the omsi HDF5 file format and ii) intergration of new analysis capabilities with the OpenMSI web API and the OpenMSI web-based viewer (see Viewer functions below for details).
    
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
* ``v_qspectrum_viewerOptions``: Define a list of strings, describing the different viewer options available for the analysis for qspectrum requests (i.e., ``v_qspectrum``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.
* ``v_qslice_viewerOptions``: Define a list of strings, describing the different viewer options available for the analysis for qslice requests (i.e., ``v_qslice``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.


        
Integrating the Analysis with the OMSI File Format
--------------------------------------------------

Making an analysis object saveable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the follwing we are extending the execute function of our analysis to record parameters, dependencies and outputs of the analysis into to make them accesible for automatic addition to an OMSI HDF5 file. At this point the data is added to a set of structured dictionaries that are defined in ``omsi_analysis_base`` , i.e, no data is actually saved to an HDF5 file yet, but we are preparing the object so that it can be easily saved.


.. code-block:: python
    :linenos:
    :emphasize-lines: 12,13,14,19,24,25,33

    def execute(self, msidata, integration_width=10, msidata_dependency=None) :
        
        #1) msidata is a raw data input for which we only want to store the dependencies. 
        #   We expect that the user specfies this dependency in the input parameter 
        #   msidata_dependency. One can also try to infer the dependency from msidata
        #   directly but that may not always be possible.
        #2) integration_width and number_of_clusters are other paramters which we want to store
        #3) We assume the code generates two numpy arrays as outputs: labels, centroids
        
        #Clear any previously stored analysis data. This is needed in case the analysis
        #has been excuted before.
        self.clear_analysis_data()
        self.clear_parameter_data()
        self.clear_dependency_data()
        
        #Save the analysis parameters to the __parameter_list so that the data can be 
        #saved automatically vie the omsi file API
        p1 = np.asarray( [integration_width] ) #Convert the object to numpy
        self.add_parameter_data( name='integration_width' , data=iw , dtype=str(iw.dtype) ) 
        #Do the same for any other parameters that need to be recored
        
        #Save the analysis output data to the __data_list so that the data can be saved 
        #automatically vie the omsi file API
        self.add_analysis_data( name='labels' , data=labels , dtype=str(labels.dtype) ) 
        self.add_analysis_data( name='centroids' , data=centroids , dtype=str(centroids.dtype) ) 
        
        #Save the analysis dependencies to the __dependency_list so that the data can be saved 
        #automatically by the omsi HDF5 file API
        if msidata_dependency is not None :
            #If the user has provided a full definition of the dependency
            if isinstance( msidata_dependency , omsi_dependency ) :
                #Add the dependency as given by the user
                self.add_dependency_data( msidata_dependency )
            #The user only gave us the object that we depend on
            #so we need to construct the dependency object
            else :
                #Assuming that we have  an h5py object or an instance of an
                #omsi_file object management class given we can construct the dependency
                #directly. This will result in a ValueError in case an incompatible object
                #is given as input.
                self.add_dependency_data( omsi_dependency( param_name = 'msidata', \
                                                           link_name='msidata', \
                                                           omsi_object=msidata_dependency, \ 
                                                           selection=None ) )
                
Executing, saving, and restoring an analysis object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that the analysis specifies what data needs to be stored we can directly use the infrastructure provided by ``omsi_analysis_base`` and the ``omsi.dataformat.omsi_file`` API to save the data from the analysis to an OMSI HDF5 file. In OMSI files, analyses are generally associated with experiments, so that we use the ``omsi.dataformat.omsi_file.omsi_file_experiment`` API here.

.. code-block:: python
    :linenos:
    :emphasize-lines: 12, 21,22,23,24

    #Open the MSI file and get the desired experiment
    from omsi.dataformat.omsi_file import *
    f = omsi_file( filename, 'a' )
    e = f.get_exp(0)
    
    #Execute the analysis
    d = e.get_msidata(0)
    a = omsi_myanalysis()
    a.execute(msidata=d, integration_width=10, msidata_dependency=d) 
    
    #Save the analysis object. 
    analysis_object , analysis_index = exp.create_analysis( a )
    #This single line is sufficient to store the complete analysis to the omsi file.
    #By default the call will block until the write is complete. Setting the
    #parameter flushIO=False enables buffered write, so that the call will 
    #return once all data write operations have been scheduled. Here we get
    #an omsi.dataformat.omsi_file.omsi_file_analysis
    #object for management of the data stored in HDF5 and the integer index of the analysis.
    
    #Restoring the analysis from file. Here we can decide which data should be loaded.
    a2 = omsi_myanalysis().read_from_omsi_file(analysisGroup=analysis_object, \
                                               load_data=True, \
                                               load_parameters=True,\
                                               dependencies_omsi_format=True ) 
    #By setting load_data and/or load_parameters to False, we create h5py instead of
    #numpy objects, avoiding the actual load of the data. CAUTION: To avoid the accidental
    #overwrite of data we recommend to use load_data and load_parameters as False only
    #when the file has been opend in read-only mode 'r'.
    
Custom data save
^^^^^^^^^^^^^^^^

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
           omsi_file_analysis.__create_analysis__(..)) so that this function does not need 
           to be called explicitly.
           
           Keyword Arguments:

           :param analysisGroup: The omsi_file_analysis object of the group for the
                                 analysis that can be used for writing.

           """
        pass
    
Custom analysis restore
^^^^^^^^^^^^^^^^^^^^^^^
    
Similarly in order implement custom data restore behavior we can overwrite the default implementation of ``omsi_analysis_base.read_from_omsi_file)`` . In this case one will usually call the default implementation via ``super(omsi_myanalysis,self).read_from_omsi_file(...)`` first and then add any additional behavior.
    
Integrating the Analysis with the OpenMSI Web API:
--------------------------------------------------

Once the analysis is stroed in the OMSI file format, integration with ``qmetadata`` and ``qcube`` calls of the web API is automatic. These ``qmetadata`` and ``qcube`` functions provide general purpose access to the data so that we can immediatly start to program against our analysis.

Some applications---such as the OpenMSI web-based viewer---utilize the simplified, special data access patterns ``qslice``, ``qspectrum``, and ``qmz`` in order to interact with the data. The default implementation of these function available in ``omsi.analysis.omsi_analysis_base`` exposes the data from all depencdencies of the analysis that support these patterns. For full integration with the web API, however, we need to implement this functionality in our analysis class. The ``qmz`` pattern in particular is relevant to both the ``qslice`` and ``qspectrum`` pattern and should be always implemented as soon as one of the two patterns is defined. 

Implementing the ``qslice`` pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:
    :emphasize-lines: 4,5,18,20,21,22,23,30,31,32,33,38,39,42,44,45

    class omsi_myanalysis(omsi_analysis_base) :
        ...

        @classmethod
        def v_qslice(cls , anaObj , z , viewerOption=0) :
            """Implement support for qslice URL requests for the viewer
            
               anaObj: The omsi_file_analysis object for which slicing should be performed. 
               z: Selection string indicting which z values should be selected.
               viewerOption: An analysis can provide different default viewer behaviors 
                             for how slice operation should be performed on the data. 
                             This is a simple integer indicating which option is used. 
           
               :returns: numpy array with the data to be displayed in the image slice 
                         viewer. Slicing will be performed typically like [:,:,zmin:zmax].
            
            """ 
            from omsi.shared.omsi_data_selection import *
            #Implement custom analysis viewer options 
            if viewerOption == 0 :
                dataset =  anaObj[ 'labels' ] #We assume labels was a 3D image cube of labels
                zselect = selection_string_to_object(z) 
                return dataset[ : , :, zselect ]
            
            #Expose recursively the slice options for any data dependencies. This is useful
            #to allow one to trace back data and generate comlex visualizations involving 
            #multiple different data sources that have some from of dependency in that they 
            #led to the generation of this anlaysis. This behavior is already provided by 
            #the default implementation of this function ins omsi_analysis_base.
            elif viewerOption >= 0 :
                #Note, the base class does not know out out viewerOptions so we need to adjust 
                #the vieweOption accordingly by substracting the number of our custom options.
                return super(omsi_myanalysis,cls).v_qslice( anaObj , z, viewerOption-1)
            #Invalid viewerOption
            else :
                return None
             
         @classmethod
            def v_qslice_viewerOptions(cls , anaObj ) :
                """Define which viewerOptions are supported for qspectrum URL's"""
                #Get the options for all data dependencies
                dependent_options = super(omsi_findpeaks_global,cls).v_qslice_viewerOptions(anaObj)
                #Define our custom viewer options
                re = ["Labels"] + dependent_options
                return re

                
                
NOTE: We here convert the selection string to a python selection (i.e., a list, slice, or integer) object using the ``omsi.shared.omsi_data_selection.check_selection_string(...)`` . This has the advantage that we can use the given selection directly in our code and avoids the use of a potentially dangerous ``eval`` , e.g., ``return eval("dataset[:,:, %s]" %(z,))`` . While we can also check the validity of the  seletion string  using  ``omsi.shared.omsi_data_selection.check_selection_string(...)`` , it is recommened to convert the string to a valid python selection to avoid possible attacks. 


Implementing the ``qspectrum`` pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:
    :emphasize-lines: 3,4,38,39,40,41,47,48,49,54,55,56,58,64,67,72,75,77,78

    class omsi_myanalysis(omsi_analysis_base) :
        ...
        @classmethod
        def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
            """Implement support for qspectrum URL requests for the viewer.
           
               anaObj: The omsi_file_analysis object for which slicing should be performed 
               x: x selection string
               y: y selection string
               viewerOption: If multiple default viewer behaviors are available for a given 
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
            if viewerOption == 0 :
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
            elif viewerOption >= 0 :
                #Note, the base class does not know out out viewerOptions so we need to adjust 
                #the vieweOption accordingly by substracting the number of our custom options.
                return super(omsi_findpeaks_global,cls).v_qspectrum( anaObj , x , y, viewerOption-1)
                
            return data , customMZ
        
        @classmethod
        def v_qspectrum_viewerOptions(cls , anaObj ) :
            """Define which viewerOptions are supported for qspectrum URL's"""
            #Get the options for all data dependencies
            dependent_options = super(omsi_findpeaks_global,cls).v_qspectrum_viewerOptions(anaObj)
            #Define our custom viewer options
            re = ["Labels"] + dependent_options
            return re
            
Implementing the ``qmz`` pattern
--------------------------------

.. code-block:: python
    :linenos:
    :emphasize-lines: 5,6,41,42,43,45,46,47,48,50,51,52,53,54,55,57,58,59,60,61,62,64

    class omsi_myanalysis(omsi_analysis_base) :
        ...
        
        @classmethod
            def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
                """Implement support for qmz URL requests for the viewer.
                
                    Get the mz axes for the analysis

                    anaObj: The omsi_file_analysis object for which slicing should be performed.
                    qslice_viewerOption: If multiple default viewer behaviors are available for
                                a given analysis then this option is used to switch between them 
                                for the qslice URL pattern.
                    qspectrum_viewerOption: If multiple default viewer behaviors are available 
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
                     for the differnent viewerOption patterns. It is in general safe to populate 
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
                if qspectrum_viewerOption == 0 and qslice_viewerOption==0: #Loadings
                    mzSpectra =  anaObj[ 'labels' ][:]
                    labelSpectra = "Labels"
                #Both viewerOptions point to a data dependency
                elif qspectrum_viewerOption > 0 and qslice_viewerOption>0 :
                    mzSpectra, labelSpectra, mzSlice, labelSlice = \ 
                           super(omsi_findpeaks_global,cls).v_qmz( anaObj, \ 
                                 qslice_viewerOption-1 , qspectrum_viewerOption-1)
                #Only the a qlice options point to a data dependency
                elif qspectrum_viewerOption == 0 and qslice_viewerOption>0 :
                    mzSpectra =  anaObj[ 'peak_mz' ][:]
                    labelSpectra = "m/z"
                    tempA, tempB, mzSlice, labelSlice = \ 
                            super(omsi_findpeaks_global,cls).v_qmz( anaObj, \ 
                                  qslice_viewerOption-1 , 0)
                #Onlye the qlice option points to a data dependency
                elif qspectrum_viewerOption > 0 and qslice_viewerOption==0 :
                    mzSlice =  anaObj[ 'peak_mz' ][:]
                    labelSlice = "m/z"
                    mzSpectra, labelSpectra, tempA, tempB = \
                            super(omsi_findpeaks_global,cls).v_qmz( anaObj, \
                                  0 , qspectrum_viewerOption-1)
                
                return mzSpectra, labelSpectra, mzSlice, labelSlice
