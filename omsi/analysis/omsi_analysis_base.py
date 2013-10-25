from omsi.dataformat.omsi_file import omsi_file_analysis, omsi_file_msidata
from omsi.analysis.omsi_analysis_data import *
from omsi.shared.omsi_dependency import *

class omsi_analysis_base(object) :
    """Base class for omsi analysis functionality. The class provides a large set of functionality designed
       to facilitate storage of analysis data in the omsi HDF5 file format. The class also provides a set
       of functions to enable easy intergration of new analysis with the OpenMSI web-based viewer (see
       Viewer functions below for details).
    
       Slicing :
       
            This class supports basic slicing to access data stored in the main member variables. By default the data is retrieved from __data_list and the __getitem__(key) function. which implemtent the [..] operator, returns __data_list[key]['data']. The key is a string indicating the name of the paramter to be retrieved. If the key is not found in the __data_list then the function will try to retrieve the data from __parameter_list instead. By adding "parameter/key" or "dependency/key" one may also explicitly retrieve values from the __parameter_list and __dependency_list.
    
       Member Variables:
           
           :param analysis_identifier: Define the name for the analysis used as key in search operations
           :param __data_list: Dictonary of omsi_analysis_data to be written to the HDF5 file. Derived classes need to add all data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See omsi.analysis.omsi_analysis_data for details.
           :param __parameter_list: Dictonary of omsi_analysis_data to be written to the HDF5 file. Derived classes need to add all parameter data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See omsi.analysis.omsi_analysis_data for details.
           :param __dependency_list: Dictonary of omsi_dependency to be written to the HDF5 file. Derived classes need to add all dependencies data that should be saved for the analysis in the omsi HDF5 file to this dictionary. See omsi.analysis.omsi_analysis_data for details.
    
       I/O functions:
       
            These functions can be optionally overwritten to control how the analysis data should be written/read from the omsi HDF5 file. Default implementations are provided here, which should be sufficient for most cases. 
            
            * ``write_to_omsi_file``: The default implementation is empty as the default data write is  managed by the omsi_file_experiment.create_analysis() function.  Overwrite this function, in case that the analysis needs to write data to the HDF5 omsi file beyond what the defualt omsi data API does.
            
            * ``read_from_omsi_file``: The default implementation tries to reconstruct the original data as far  as possible, however, in particular in case that a custom write_to_omsi_file            funtion has been implemented, the default implementation may not be sufficien. The default implementation reconstructs: i) analysis_identifier and reads all custom data into ii)__data_list. Note, an error will be raised in case that the analysis type specified in the HDF5 file does not match the analysis type specified by get_analysis_type(). This function can be optionally overwritten to implement a custom data read.
            
        Viewer functions:
        
            Several convenient functions are used to allow the OpenMSI online viewer to interact with the analysis and to visualize it. The default implementations provided here simply indicate that the analysis does not support the data access operations required by the online viewer. Overwrite these functions in the derived analysis classes in order to interface them with the viewer. All viewer-related functions start with ``v\_...`` .
            
            NOTE: the default implementation of the viewer functions defined in ``omsi_analysis_base`` are designed to take care of the common requirement for providing viewer access to data from all depencies of an analysis. In many cases, the default implementation is often sill called at the end of custom viewer functions.
            
            NOTE: The viewer functions typically support a viewerOption parameter. viewerOption=0 is expected to refer to the analysis itself.
            
            * ``v_qslice``: Retrieve/compute data slices as requested via qslice URL requests. The corrsponding view of the DJANGO data access server already translates all input parameters and takes care of generating images/plots if needed. This function is only responsible for retrieving the data.
            * ``v_qspectrum``: Retrieve/compute spectra as requested via qspectrum URL requests. The corrsponding view of the DJANGO data access server already translates all input parameters and takes care of generating images/plots if needed. This function is only responsible for retrieving the data.
            * ``v_qmz``: Define the m/z axes for image slices and spectra as requested by qspectrum URL requests.
            * ``v_qspectrum_viewerOptions``: Define a list of strings, describing the different viewer options available for the analysis for qspectrum requests (i.e., ``v_qspectrum``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.
            * ``v_qslice_viewerOptions``: Define a list of strings, describing the different viewer options available for the analysis for qslice requests (i.e., ``v_qslice``). This feature allows the analysis developer to define multiple different visualization modes for the analysis. For example, when performing a data rediction (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are most interested in.
    """
    def __init__(self):
        """Initalize the basic data members"""
        self.analysis_identifier = "undefined"
        self.__data_list = []
        self.__parameter_list = []
        self.__dependency_list = []
        
    def __getitem__(self, key) :
        """This class supports basic slicing to access data stored in the main member variables. 
           By default the data is retrieved from __data_list and the __getitem__(key) function.
           which implemtent the [..] operator, returns __data_list[key]['data']. The key is 
           a string indicating the name of the paramter to be retrieved. If the key is not 
           found in the __data_list then the function will try to retrieve the data from 
           __parameter_list instead. By adding "parameter/key" or "dependency/key" 
           one may also explicitly retrieve values from the __parameter_list and __dependency_list.
        """
        if isinstance( key , str ) :
            
            try :
                for i in self.__data_list :
                    if i['name'] == key : 
                        return i['data']
            except :
                pass
            try:
                for i in self.__parameter_list :
                    if i['name'] == key : 
                        return i['data']
            except :
                pass
            if key.startswith("parameter/" ) :
                key = key.lstrip('parameter/')
                try : 
                    for i in self.__parameter_list :
                        if i['name'] == key : 
                            return i['data']
                except:
                    pass
            if key.startswith("dependency/" ) :
                key = key.lstrip('dependency/')
                try : 
                    for i in self.__dependency_list :
                        if i['name'] == key : 
                            return i['omsi_object']
                except :
                    pass
            return None
        else :
            return None
                

    @classmethod
    def v_qslice(cls , anaObj , z , viewerOption=0) :
        """Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer
        
           :param anaObj: The omsi_file_analysis object for which slicing should be performed 
           :param z: Selection string indicting which z values should be selected.
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them.
           
           :returns: numpy array with the data to be displayed in the image slice viewer. Slicing will be performed typically like [:,:,zmin:zmax].
           
           :raises: NotImplementedError in case that v_qslice is not supported by the analsis.
        """
        from omsi.analysis.omsi_viewer_helper import omsi_viewer_helper
        from omsi.shared.omsi_data_selection import check_selection_string, selection_type, selection_to_indexlist
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex = cls.__construct_dependent_viewer_options__(anaObj)
        #Check whether the given selection is valid
        zType = check_selection_string(z)
        if zType == selection_type['invalid'] :
            return None 
        if isinstance( re_slicedata[viewerOption] , omsi_file_msidata ) :
            
            try :
                data = eval("re_slicedata[viewerOption][:,:, %s]" %(z,))
                return data 
            except :
                return None 
        
        elif isinstance( re_slicedata[viewerOption] , omsi_file_analysis ) : 
        
            analysisType = str(re_slicedata[viewerOption].get_analysis_type()[0])
            return omsi_viewer_helper.__string_to_class__( analysisType ).v_qslice(anaObj=re_slicedata[viewerOption], z=z, viewerOption=re_slice_optionIndex[viewerOption])
            
        else:
            return None
        
        
        #raise NotImplementedError( "Implement me: The v_qslice function should return a numpy array with the requested image slice data" )
    
   
    @classmethod
    def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
        """Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer
        
           Developer Note: h5py currently supports only a single index list. If the user provides an index-list for both
                           x and y, then we need to construct the proper merged list and load the data manually, or if
                           the data is small enough, one can load the full data into a numpy array which supports 
                           mulitple lists in the selection. 
        
           :param anaObj: The omsi_file_analysis object for which slicing should be performed 
           :param x: x selection string
           :param y: y selection string
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them.
           
           :returns: The following two elemnts are expected to be returned by this function :
           
                1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The mass (m/z) axis must be the last axis. For index selection x=1,y=1 a 1D array is usually expected. For indexList selections x=[0]&y=[1] usually a 2D array is expected. For ragne selections x=0:1&y=1:2 we one usually expects a 3D arrya.
                2) None in case that the spectra axis returned by v_qmz are valid for the returned spectrum. Otherwise, return a 1D numpy array with the m/z values for the spectrum (i.e., if custom m/z values are needed for interpretation of the returned spectrum).This may be needed, e.g., in cases where a per-spectrum peak analysis is performed and the peaks for each spectrum appear at different m/z values. 
        """
        import numpy as np
        from omsi.analysis.omsi_viewer_helper import omsi_viewer_helper
        from omsi.shared.omsi_data_selection import check_selection_string, selection_type, selection_to_indexlist,  selection_string_to_object
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex = cls.__construct_dependent_viewer_options__(anaObj)
        #Check whether the given selection is valid
        xType = check_selection_string(x)
        yType = check_selection_string(y)
        xSel = selection_string_to_object( y )
        ySel = selection_string_to_object( x )
        
        if (xType == selection_type['invalid']) or (yType == selection_type['invalid']) :
            return None , None
        if isinstance( re_spectrumdata[viewerOption] , omsi_file_msidata ) :
            if xType == selection_type['invalid'] or yType == selection_type['invalid'] :
                return None , None
            if xType == selection_type['indexlist'] and yType == selection_type['indexlist'] :
                #We now need to match up the index lists and load all the individial values
                #ToDo: This is can be very inefficient for large lists
                xSize = len(xSel)
                ySize = len(ySel)
                if xSize != ySize :
                    raise KeyError( "Selection lists don't match" )
                dataSet = re_spectrumdata[viewerOption]
                zSize = dataSet.shape[2]
                #Allocate the required memory
                data = np.zeros( (xSize, zSize) , dtype=dataSet.dtype )
                for i in xrange(0,xSize) :
                    data[i,:] = dataSet[ xSel[i] , ySel[i] , : ]
            else :
                data = re_spectrumdata[viewerOption][ xSel , ySel, : ]
            
            return data ,None
            
            #try :
            #    data = eval("re_spectrumdata[viewerOption][%s,%s, :]" %(x,y))
            #    return data , None
            #except :
            #    import sys
            #    print str(sys.exc_info())
            #    return None , None
        
        elif isinstance( re_spectrumdata[viewerOption] , omsi_file_analysis ) : 
        
            analysisType = str(re_spectrumdata[viewerOption].get_analysis_type()[0])
            return omsi_viewer_helper.__string_to_class__( analysisType ).v_qspectrum(anaObj=re_spectrumdata[viewerOption], x=x, y=y, viewerOption=re_spectrum_optionIndex[viewerOption])
        
        else:
            return None , None
        
        #raise NotImplementedError( "Implement me: The v_spectrum function should return: i) a numpy array with the requested spectrum data and ii) information about the m/z axis for the spectrum")  
    
    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
        """ Get the mz axes for the analysis
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed
            :param qslice_viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them for the qslice URL pattern.
            :param qspectrum_viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them for the qspectrum URL pattern.
        
            :returns: The following four arrays are returned by the analysis:
            
                - mzSpectra : Array with the static mz values for the spectra.
                - labelSpectra : Lable for the spectral mz axis 
                - mzSlice : Array of the static mz values for the slices or None if identical to the mzSpectra.
                - labelSlice : Lable for the slice mz axis or None if identical to labelSpectra.
        """
        from omsi.analysis.omsi_viewer_helper import omsi_viewer_helper
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex = cls.__construct_dependent_viewer_options__(anaObj)
        mzSpectra=None
        labelSpectra=None
        mzSlice=None 
        labelSlice = None
        #Determine the spectra mz axis
        if len(re_spectrum) > 0 :
            if isinstance( re_spectrumdata[qspectrum_viewerOption] , omsi_file_msidata ) :
                mzSpectra = re_spectrumdata[qspectrum_viewerOption].mz[:]
                labelSpectra = "m/z"
            elif isinstance( re_spectrumdata[qspectrum_viewerOption]  , omsi_file_analysis ) :
                mzSpectra, labelSpectra, tempA, tempB = omsi_viewer_helper.get_axes( re_spectrumdata[qspectrum_viewerOption], qslice_viewerOption=re_slice_optionIndex[qslice_viewerOption], qspectrum_viewerOption=re_spectrum_optionIndex[qspectrum_viewerOption])
            else :
                mzSpectra = None
                labelSpectra = None
        #Determine the slice mz axis
        if len(re_slice) > 0 :
            if re_spectrumdata[qspectrum_viewerOption] != re_slicedata[qslice_viewerOption] :
                if isinstance( re_slicedata[qslice_viewerOption] , omsi_file_msidata ) :
                    mzSlice = re_slicedata[qslice_viewerOption].mz[:]
                    labelSlice= "m/z"
                elif isinstance( re_slicedata[qslice_viewerOption]  , omsi_file_analysis ) :
                    tempA, tempB, mzSlice, labelSlice = omsi_viewer_helper.get_axes( re_slicedata[qslice_viewerOption], qslice_viewerOption=re_slice_optionIndex[qslice_viewerOption], qspectrum_viewerOption=re_spectrum_optionIndex[qspectrum_viewerOption] )
                else :
                    mzSlice=None 
                    labelSlice = None
                
        return mzSpectra, labelSpectra, mzSlice, labelSlice
        #raise NotImplementedError( "Implement me: The v_qmz function should return i) array with the statix mz values for the spectra, ii) the lable to be used for the spectra and iii) array with the static mz values for the slices (or None if the same as the spectra), and iv) the lable to be used to the slices (or None if the same as for the spectra)" )  
        #return None, None, None, None
    
    @classmethod
    def v_qspectrum_viewerOptions(cls , anaObj ) :
        """Get a list of strings describing the different default viewer options for the analysis for qspectrum. 
           The default implementation tries to take care of handling the spectra retrieval for all the depencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed.  For most cases this is not needed here as the support for slice operations is usually a static decission based on the class type, however, in some cases additional checks may be needed (e.g., ensure that the required data is available).
        
            :returns: List of strings indicating the different available viewer options. The list should be empty if the analysis does not support qspectrum requests (i.e., v_qspectrum(...) is not available).
        """
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex = cls.__construct_dependent_viewer_options__(anaObj)
        return re_spectrum
        #raise NotImplementedError( "Implement me: The v_qspectrum_viewerOptions function should return i) an empty list of v_qspectrum is not available, ii) a list of strings indicating the different types of data for which spectra can be retrieved for this analysis" )  
        
    @classmethod
    def v_qslice_viewerOptions(cls , anaObj ) :
        """Get a list of strings describing the different default viewer options for the analysis for qslice.
           The default implementation tries to take care of handling the spectra retrieval for all the depencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed.  For most cases this is not needed here as the support for slice operations is usually a static decission based on the class type, however, in some cases additional checks may be needed (e.g., ensure that the required data is available).
        
            :returns: List of strings indicating the different available viewer options. The list should be empty if the analysis does not support qslice requests (i.e., v_qslice(...) is not available).
        """
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex = cls.__construct_dependent_viewer_options__(anaObj)
        return re_slice
        
        #raise NotImplementedError( "Implement me: The v_qslice_viewerOptions function should return i) an empty list of v_qslice is not available, ii) a list of strings indicating the different types of data for which image slices can be retrieved for this analysis" )  
       
    @classmethod
    def __construct_dependent_viewer_options__(cls, anaObj) :
        """Internal helper function to construc the viewer options for analysis depencies.
        
          :returns: The following 4 lists:
          
            * ``re_slice`` : List of names for the slice options
            * ``re_spectrum`` : List of names for the spectrum options
            * ``re_slicedata`` : omsi_file_* API object associated with the slice options. These are either omis_file_analysis or omsi_file_msidata objects.
            * ``re_spectrumdata`` : omsi_file_* API object associated with the slice options. Theser are either omis_file_analysis or omsi_file_msidata objects.
            * ``re_slice_optionIndex``: List of integrers indicating for the given entry the viewerOption to be used with the re_slicedata object for the given option.
            * ``re_spectrum_optionIndex``: List of integrers indicating for the given entry the viewerOption to be used with the re_spectrumdata object for the given option. 
        """
        
        from omsi.analysis.omsi_viewer_helper import omsi_viewer_helper
        re_slice = []
        re_slicedata = [] 
        re_slice_optionIndex = []
        re_spectrum = []
        re_spectrumdata = []
        re_spectrum_optionIndex = []
        
        dependencies = anaObj.get_all_dependency_data_recursive()
        for di in dependencies :
            #Check if we can slice the data
            if isinstance( di['omsi_object'] , omsi_file_msidata ) :
                re_spectrum.append( "Raw Data: "+di['link_name'] )
                re_slice.append( "Raw Data: "+di['link_name'] )
                re_slice_optionIndex.append(0)
                re_spectrumdata.append( di['omsi_object'])
                re_slicedata.append( di['omsi_object'] )
                re_spectrum_optionIndex.append(0)
            elif isinstance( di['omsi_object'] , omsi_file_analysis ) :
                slice_options    = omsi_viewer_helper.get_qslice_viewerOptions( di['omsi_object'] )
                spectrum_options = omsi_viewer_helper.get_qspectrum_viewerOptions( di['omsi_object'] )
                for si in range(0,len(slice_options) ) :
                    re_slice.append( slice_options[si] )
                    re_slicedata.append( di['omsi_object'] )
                    re_slice_optionIndex.append( si )
                for si in range(0,len(spectrum_options) ) :
                    re_spectrum.append( spectrum_options[si] )
                    re_spectrumdata.append( di['omsi_object'] )
                    re_spectrum_optionIndex.append( si )    
                #analysisType = str(anaObj.get_analysis_type()[0])
                #if omsi_viewer_helper.supports_slice( di['omsi_object']) :
                #    re_slice.append( "Analysis: "+str(di['omsi_object'].get_analysis_identifier()[0]) )
                #    re_slicedata.append( di['omsi_object'] )
                #if omsi_viewer_helper.supports_spectra( di['omsi_object'] ) :
                #    re_spectrum.append( "Analysis: "+str(di['omsi_object'].get_analysis_identifier()[0]) )
                #    re_spectrumdata.append( di['omsi_object'] )
            else :
                print "Unknown dependency"
        return re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_optionIndex, re_spectrum_optionIndex
        
        
    def get_analysis_type(self) :
        """Return a string indicating the type of analysis performed"""
        return self.__module__ #self.__class__.__name__
        #raise NotImplementedError( "Implement me: The get_analysis_type function should return a string indicating the analysis type" )    
    
    def get_analysis_data_names(self) :
        """Get a list of all analysis dataset names."""
        
        return self.__data_list.keys()
    
    def get_parameter_data_names(self) :
        """Get a list of all parameter dataset names."""
        
        return self.__parameter_list.keys()
    
    def get_dependency_data_names(self) :
        """Get a list of all parameter dataset names."""
        
        return self.__dependency_list.keys()
    
    def get_analysis_data(self, index):
        """Given the index return the associated dataset to be written to the HDF5 file

          :param index : Retrun the index entry of the private member __data_list.
        """
        return self.__data_list[index]
        
    def get_parameter_data(self, index):
        """Given the index return the associated dataset to be written to the HDF5 file

          :param index : Retrun the index entry of the private member __parameter_list.
        """
        return self.__parameter_list[index]
        
    def get_dependency_data(self, index):
        """Given the index return the associated dataset to be written to the HDF5 file

          :param index : Retrun the index entry of the private member __parameter_list.
        """
        return self.__dependency_list[index]
        
    def get_analysis_data_by_name(self, dataname) :
        """Given the key name of the data return the associated omsi_analysis_data object.

          :param dataname: Name of the analysis data requested from the private __data_list member.
        """
        for i in self.__data_list :
            if i['name'] == dataname : 
                return i
                
    def get_parameter_data_by_name(self, dataname) :
        """Given the key name of the data return the associated omsi_analysis_data object.

          :param dataname: Name of the parameter requested from the private __parameter_list member.
        """
        for i in self.__parameter_list :
            if i['name'] == dataname : 
                return i
                
    def get_dependency_data_by_name(self, dataname) :
        """Given the key name of the data return the associated omsi_analysis_data object.

          :param dataname: Name of the dependency requested from the private __dependency_list member.
        """
        for i in self.__dependency_list :
            if i['name'] == dataname : 
                return i
        
    def get_all_analysis_data(self) :
        """Get the complete list of all analysis datasets to be written to the HDF5 file"""
        return self.__data_list
        
    def get_all_parameter_data(self) :
        """Get the complete list of all parameter datasets to be written to the HDF5 file"""
        return self.__parameter_list
        
    def get_all_dependency_data(self) :
        """Get the complete list of all direct dependencies to be written to the HDF5 file
        
           NOTE: These are only the direct dependencies as sepecified by the analysis itself. \
           Use  get_all_dependency_data_recursive(..) to also get the indirect depencies of \
           the analysis due to dependencies of the depencies themself.
        """
        return self.__dependency_list
        
    def get_num_analysis_data(self) :
        """Retrun the number of analysis datasets to be wirtten to the HDF5 file"""
        return len(self.__data_list)
        
    def get_num_parameter_data(self) :
        """Retrun the number of parameter datasets to be wirtten to the HDF5 file"""
        return len(self.__parameter_list)
        
    def get_num_dependency_data(self) :
        """Retrun the number of dependencies to be wirtten to the HDF5 file"""
        return len(self.__dependency_list)
        
    def clear_analysis_data(self):
        """Clear the list of analysis data"""
        self.__data_list = []
        
    def clear_parameter_data(self) : 
        """Clear the list of parameter data"""
        self.__parameter_list = []
        
    def clear_dependency_data(self) : 
        """Clear the list of parameter data"""
        self.__dependency_list = []
        
    def add_analysis_data(self , name, data, dtype ) : 
        """Add a new dataset to the list of data to be written to the HDF5 file

          The input parameters will be transformed into a omsi_analysis_data dictionary.

          :param name: The name for the dataset in the HDF5 format
          :param data: The numpy array to be written to HDF5. The data write function 
                omsi_file_experiment.create_analysis used for writing of the data to file can
                in principal also handel other primitive data types by explicitly converting them
                to numpy. However, in this case the dtype is determined based on the numpy conversion
                and correct behavior is not guaranteed. I.e., even single scalars should be stored as
                a 1D numpy array here.
          :param dtype: The data type to be used during writing. For standard numpy data types this is just 
                 the dtype  of the dataset, i.e., ['data'].dtype. Other allowed datatypes are: 
                 
                 * For string:  omsi_format.str_type (omsi_format is located in omsi.dataformat.omsi_file )
                 * To generate data links: ana_hdf5link   (omsi_analysis_data)
        """
        self.__data_list.append( omsi_analysis_data(name=name, data=data, dtype=dtype) )
        
    def add_parameter_data(self , name, data, dtype ) :
        """Add a new analysis parameter dataset to the list of data to be written to the HDF5 file

           The input parameters will be transformed into a omsi_analysis_data dictionary.

          :param name: The name for the dataset in the HDF5 format
          :param data: The numpy array to be written to HDF5. The data write function 
                omsi_file_experiment.create_analysis used for writing of the data to file can
                in principal also handel other primitive data types by explicitly converting them
                to numpy. However, in this case the dtype is determined based on the numpy conversion
                and correct behavior is not guaranteed. I.e., even single scalars should be stored as
                a 1D numpy array here.
          :param dtype: The data type to be used during writing. For standard numpy data types this is just 
                 the dtype  of the dataset, i.e., ['data'].dtype. Other allowed datatypes are: 
                 
                 * For string:  omsi_format.str_type (omsi_format is located in omsi.dataformat.omsi_file )
                 * To generate data links: ana_hdf5link   (omsi_analysis_data)
          
        """
        self.__parameter_list.append( omsi_analysis_data(name=name, data=data, dtype=dtype) )
        
        
    def add_dependency_data(self , dependency ) :
        """Add a new analysis parameter dataset to the list of data to be written to the HDF5 file

          :param dependency: The data that should be added to __dependency_list member to be written to HDF5.
          :type data: omsi_dependency 

          :raises: ValueError in case data is not of omsi_dependency object. 
        """
        if isinstance( dependency , omsi_dependency ) :
            self.__dependency_list.append( dependency )
        else :
             raise ValueError( "Invalid input for add_data_dependency function" )
            
        #from omsi.dataformat.omsi_file import omsi_file, omsi_file_analysis, omsi_file_experiment, omsi_file_instrument, omsi_file_sample, omsi_file_msidata
        #if isinstance( link
        #
        #if isinstance( link_path , str ) :
        #    self.__dependency_list.append( omsi_analysis_data( name=name , data=link_path  , dtype=omsi_analysis_data.ana_hdf5link  ) )
        #elif isinstance( link_path , omsi_file_analysis ) or isinstance( link_path, omsi_file_msidata ) \
        #     or isinstance( link_path, omsi_file_experiment ) or isinstance( link_path, omsi_file_sample ) \
        #     or isinstance( link_path, omsi_file_instrument ) or isinstance( link_path, omsi_file ) :
        #
        #    self.__dependency_list.append( omsi_analysis_data( name=name , data=link_path.name , dtype=omsi_analysis_data.ana_hdf5link  ) )
        #else :
        #    errorMessage = "Could not add dependency data to the analysis as no valid HDF5 link could be identified."
        #    raise ValueError( errorMessage )
        

    def write_to_omsi_file(self , analysisGroup) :
        """This function can be optionally overwritten to implement a custom data write 
           function for the analysis to be used by the omsi_file API.
           
           Note, this function should be used only to add additional data to the analysis
           group. The data that is written by default is typically still written by the omsi_file_experiment.create_analysis() 
           function, i.e., the following data is wirtten by default: i) analysis_identifier ,ii) get_analysis_type, iii)__data_list,
           iv) __parameter_list , v) __dependency_list. Since the omsi_file.experiment.create_analysis() functions takes care
           of setting up the basic strucutre of the analysis storage (included the subgroubs for storing parameters and data
           dependencies) this setup can generally be assumed to exist before this function is called. This function is called
           automatically at the end omsi_file.experiment.create_analysis() (i.e, actually omsi_file_analysis.__create_analysis__(..)
           so that this function typically does not need to be called explicitly.
           
           Keyword Arguments:

           :param analysisGroup: The omsi_file_analysis object of the group for the analysis that can be used for writing.

           """
        pass
    
    def read_from_omsi_file(self, analysisGroup, load_data=True, load_parameters=True, dependencies_omsi_format=True ) :
        """This function can be optionally overwritten to implement a custom data read.
           
           The default implementation tries to reconstruct the original data as far
           as possible, however, in particular in case that a custom write_to_omsi_file
           funtion has been implemented, the default implementation may not be sufficient.
           The default implementation reconstructs: i) analysis_identifier and reads all
           custom data into iii)__data_list. Note, an error will be raised in case that
           the analysis type specified in the HDF5 file does not match the analysis type
           specified by get_analysis_type()
           
           Keyword Parameters:

           :param analysisGroup: The omsi_file_analysis object associated with the hdf5 data group with the anlaysis data_list
           :param load_data: Should the analysis data be loaded from file (default) or just stored as h5py data objects
           :param load_parameters: Should parameters be loaded from file (default) or just stored as h5py data objects. 
           :param dependencies_omsi_format: Should dependencies be loaded as omsi_file_ API objects (default) or just as h5py objects. 
           
           Returns

           :returns bool: Boolean indicating whether the data was read succesfully
           
           Raises:

           :raise: TypeError : A type error will be raised in case that the analysis type specified \
                       by the file does not match the analysis type provided by self.get_analysis_type()
           
           
        """
        if str(analysisGroup.get_analysis_type()[0]) != str(self. get_analysis_type()) :
            errorMessage = "The type of the analysis specified in the omsi data file does not match the analysis type of the object "
            errorMessage = errorMessage + str(analysisGroup.get_analysis_type()[0]) + " != " + str(self.get_analysis_type())
            raise TypeError( errorMessage )
        
        identifier = analysisGroup.get_analysis_identifier()
        if identifier is not None :
            self.analysis_identifier = identifier[0]
        else :
            print "The analysis identifier could not be read from the omsi file"
        
        self.__data_list = analysisGroup.get_all_analysis_data(load_data=load_data)
        self.__parameter_list = analysisGroup.get_all_parameter_data(load_data= load_parameters)
        self.__dependency_list = analysisGroup.get_all_dependency_data(omsi_format= dependencies_omsi_format)
        
        return True
        
       
        
    def get_analysis_identifier(self) :
        """Return the name of the analysis used as key when searching for a particular analysis"""
        return self.analysis_identifier
        
    def set_analysis_identifier(self , newName) :
        """Set the name of the analysis to newName

           :param newName: The new analysis identifier string to be used (should be unique)
           :type newName: str
        """
        self.analysis_identifier = newName
        


