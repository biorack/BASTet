"""Helper functions and classes for interfacing different anlaysis algorithms with the webbased viewer"""
import numpy as np
from omsi.dataformat.omsi_file import *
from omsi.shared.omsi_data_selection import *
from omsi.analysis import *
import types


class omsi_viewer_helper(object) :
    """Helper class for interfacing different anlaysis algorithms with the webbased viewer"""
         
    def __init__(self):
        """Initalize the basic data members"""
        pass
    
    
    @classmethod
    def get_slice(cls, anaObj , z ,  normalize=-1, reduction=None, axis=0 , reduction2=None, axis2=0, viewerOption=0) :
        """Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer
        
           :param anaObj: The omsi_file_analysis object for which slicing should be performed 
           :param z: Selection string indicting which z values should be selected. 
           :param normalize: Should the returned data be normalized by deviding by the max value
           :param reduction: Reduction opteration to be performed for the slice
           :param axis: Axis along which the reduction should be performed
           :param reduction2: Secondary reduction to be applied to the data 
           :param axis2: Axis along which the secondary reduction should be performed. (Note, remember that the dimensionality of the data has already been reduced by 1 by the first reduction operation)
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them.
           
           :returns: numpy array with the data to be displayed in the image slice viewer. Slicing will be performed typically like [:,:,zmin:zmax].
        """
        try :
            analysisType = anaObj.get_analysis_type()[0].lstrip('omsi.analysis')
            data = cls.__string_to_class__( analysisType ).v_qslice(anaObj, z, viewerOption)
            if data is None :
                return None
        except :
            return None
        
        #We expect a 3D dataset (x,y,m/z) so if only a single 2D slice is returned then we reshape the data first
        if len(data.shape) == 2 :
            data = data.reshape( (data.shape[0], data.shape[1] , 1 ) )
        
        try:
            #Perform the requested redution operation
            if reduction is not None:
                data = perform_reduction( data , reduction , axis )
                
            if (reduction2 is not None) and ( not isinstance(data, HttpResponse ) ) :
                data = perform_reduction( data , reduction2 , axis2 )
            
            if data is None :
                return None
        except :
            return None
        
        try:
            #Normalize the data if requested
            if int(normalize)>0 :
                dm = float (np.max(data) )
                if dm != 0 :
                    data = data/dm
        except:
            return None

        return data
        
        
    @classmethod
    def get_spectra( cls, anaObj , x, y , normalize=-1, reduction=None, axis=0 , reduction2=None, axis2=0, viewerOption=0) :
        """Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer
        
           Developer Note: h5py currently supports only a single index list. If the user provides an index-list for both
                           x and y, then we need to construct the proper merged list and load the data manually, or if
                           the data is small enough, one can load the full data into a numpy array which supports 
                           mulitple lists in the selection. 
        
           :param anaObj: The omsi_file_analysis object for which slicing should be performed 
           :param x: x selection string
           :param y: y selection string
           :param normalize: Should the returned data be normalized by deviding by the max value
           :param reduction: Reduction opteration to be performed for the slice
           :param axis: Axis along which the reduction should be performed
           :param reduction2: Secondary reduction to be applied to the data 
           :param axis2: Axis along which the secondary reduction should be performed. (Note, remember that the dimensionality of the data has already been reduced by 1 by the first reduction operation)
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis then this option is used to switch between them.
           
           :returns: 2D or 3D numpy array of the requested spectra. The mass (m/z) axis must be the last axis.
        
        """
        try :
            analysisType = str(anaObj.get_analysis_type()[0]).lstrip('omsi.analysis')
            data , customMZ = cls.__string_to_class__( analysisType ).v_qspectrum(anaObj, x, y, viewerOption)
            if data is None :
                return None, None
        except :
            return None , None
        
        try:
            #Perform the requested redution operation
            if reduction is not None and axis<len(data.shape):
                data = perform_reduction( data , reduction , axis )
            
            if reduction2 is not None and axis<len(data.shape):
                data = perform_reduction( data , reduction2 , axis2 )
            
            if data is None :
                return None , None
        except :
            return None , None
        
        try:
            #Normalize the data if requested
            if int(normalize)>0 :
                dm = float (np.max(data) )
                if dm != 0 :
                    data = data/dm
        except:
            return None , None

        return data, customMZ 
        
        
    @classmethod
    def get_axes(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
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
        mzSpectra =  None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        try :
            analysisType = str(anaObj.get_analysis_type()[0]).lstrip('omsi.analysis')
            mzSpectra, labelSpectra, mzSlice, labelSlice = cls.__string_to_class__( analysisType ).v_qmz(anaObj, qslice_viewerOption, qspectrum_viewerOption)
        except :
            pass
        return mzSpectra, labelSpectra, mzSlice, labelSlice
      

    @classmethod
    def supports_slice(cls, anaObj ) :
        """ Get whether a default slice selection behavior is defined for the analysis.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed
        
            :return: Boolean indicating whether get_slice(...) is defined for the analysis object.
        """
        supportsSlice = len( cls.get_qslice_viewerOptions(anaObj) ) > 0
        return supportsSlice
        

    @classmethod
    def supports_spectra(cls, anaObj ) :
        """ Get wheter a default spectra selection behavior is defined for the analysis.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed.
            
            :return: Boolean indicating whether get_spectra(...) is defined for the analysis object.
        """
        supportsSpectra = len( cls.get_qspectrum_viewerOptions(anaObj) ) > 0
        return supportsSpectra
            
            
    @classmethod
    def get_qspectrum_viewerOptions(cls , anaObj ) :
        """Get a list of strings describing the different default viewer options for qspectrum.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed.
        
            :returns: Array of strings indicating the different available viewer options. The array may be empty if now viewerOptions are available, i.e., get_slice and get_spectrum are undefined for the given analysis.
        """
        viewerOptions = []
        try :
            analysisType = str(anaObj.get_analysis_type()[0]).lstrip('omsi.analysis')
            viewerOptions = cls.__string_to_class__( analysisType ).v_qspectrum_viewerOptions(anaObj)
        except :
            message = "An error occured while trying to check which qspectrum viewerOptions are available for the analysis: "+str(sys.exc_info())
            print message
        return viewerOptions
        
    @classmethod
    def get_qslice_viewerOptions(cls , anaObj ) :
        """Get a list of strings describing the different default viewer options for qslice.
        
            :param anaObj: The omsi_file_analysis object for which slicing should be performed.
        
            :returns: Array of strings indicating the different available viewer options. The array may be empty if now viewerOptions are available, i.e., get_slice and get_spectrum are undefined for the given analysis.
        """
        viewerOptions = []
        try :
            analysisType = str(anaObj.get_analysis_type()[0]).lstrip('omsi.analysis')
            viewerOptions = cls.__string_to_class__( analysisType ).v_qslice_viewerOptions(anaObj)
        except :
            message = "An error occured while trying to check which qslice viewerOptions are available for the analysis: "+str(sys.exc_info())
            print message
        return viewerOptions


    @classmethod
    def __string_to_class__(cls , className ) :
        """Convert the given string indicating the class to a python class"""
        try:
            classObj = getattr(sys.modules[__name__], className)
            if isinstance(classObj , (types.ClassType, types.TypeType)):
                return classObj
        except AttributeError:
            raise NameError(className+" doesn't exist or is not a class.")



