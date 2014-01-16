from omsi.analysis.omsi_analysis_base import *

class omsi_analysis_generic(omsi_analysis_base) :
    """This analysis class is used if the specific anlaysis type is unknown, e.g., when loading 
    custom user-defined analysis data that may have not be available in the standard
    omsi package used.
    """
    def __init__(self , nameKey="undefined"):
        """Initalize the basic data members

           Keyword arguments:

           :param nameKey: The name for the analysis
        """
        super(omsi_analysis_generic,self).__init__()
        self.analysis_identifier = nameKey
        self.real_analysis_type = None #This is the analysis type indicated in the HDF5 file
    
    @classmethod
    def v_qslice(cls , anaObj , z , viewerOption=0) :
        """Implement support for qslice URL requests for the viewer""" 
        return None
    
    @classmethod
    def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
        """Implement support for qspectrum URL requests for the viewer"""
        super(omsi_analysis_generic,cls).v_qspectrum( anaObj , x , y, viewerOption )
        
    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
        """Implement support for qmz URL requests for the viewer"""
        super(omsi_analysis_generic,cls).v_qmz( anaObj , qslice_viewerOption ,qspectrum_viewerOption )
    
    @classmethod
    def v_qspectrum_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        return super(omsi_analysis_generic,cls).v_qspectrum_viewerOptions( anaObj )
            

    @classmethod
    def v_qslice_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        return super(omsi_analysis_generic,cls).v_qslice_viewerOptions( anaObj )
    
    
    def read_from_omsi_file(self, analysisGroup , load_data=True, load_parameters=True, dependencies_omsi_format=True ) :
        """Overwrite the default implementation to remove the type check as this class should
        be able to store any unknown types of analysis in a generic fashion.

           Keyword Parameters:

           :param analysisGroup: The omsi_file_analysis object associated with the hdf5 data group with the anlaysis data_list
           :param load_data: Should the analysis data be loaded from file (default) or just stored as h5py data objects
           :param load_parameters: Should parameters be loaded from file (default) or just stored as h5py data objects. 
           :param dependencies_omsi_format: Should dependencies be loaded as omsi_file API objects (default) or just as h5py objects. 
           
           Returns

           :returns bool: Boolean indicating whether the data was read succesfully
           


        """
        self.real_analysis_type = analysisGroup.get_analysis_type() 
        identifier = analysisGroup.get_analysis_identifier()
        if identifier is not None :
            self.analysis_identifier = identifier[0]
        else :
            print "The analysis identifier could not be read from the omsi file"
        
        self.__data_list = analysisGroup.get_all_analysis_data(load_data=load_data)
        self.__parameter_list = analysisGroup.get_all_parameter_data(load_data= load_parameters)
        self.__dependency_list = analysisGroup.get_all_dependency_data(omsi_format= dependencies_omsi_format)
        
        return True
        
    def get_real_analysis_type() :
        """This class is designed to handle generic (including unkown) types of analysis.
        In cases, e.g., were this class is used to store analysis data from an HDF5
        file we may have an actual analysis type available even if we do not have
        a special analysis class may not be available in the current intallation"""
        return self.real_analysis_type

    @classmethod
    def get_analysis_type(self) :
        """Return a string indicating the type of analysis performed"""
        return "generic"
