class omsi_dependency(dict) :
    """Define a dependency between an analysis and another omsi data object
    
        **Required Keyword Arguments**:

           :param param_name: The name of the parameter that has the depency
           :param link_name: The name of for the link to be created in the HDF5 file.
           :param omsi_object: The object to which a link should be established to. This must be either an h5py.Dataset or the omsi_file_analysis or omsi_file_msidata or any of the other omsi_file API  interface ojects.
           :param selection: Optional string type parameter indicating a python selection for the dependency

    """
    def __init__(self, param_name=None , link_name=None, omsi_object=None, selection=None ) :
        """Initatlize the allowed set of keys"""
        super(omsi_dependency,self).__init__()
        dict.__setitem__( self, 'param_name'  , param_name )
        dict.__setitem__( self, 'link_name'  , link_name )
        dict.__setitem__( self, 'omsi_object' , None )
        if omsi_object is not None :
            dict.__setitem__( self, 'omsi_object' , omsi_object )
        dict.__setitem__( self ,'selection' , selection )
        
    def __setitem__(self, key, value ) :
        """Overwrite the __setitem__ function inheritied from dict to ensure that only elements with a specific
           set of keys can be modified"""
        if self.has_key(key) :
            if key == "omsi_object" :
                from omsi.dataformat.omsi_file import omsi_file, omsi_file_analysis, omsi_file_experiment, omsi_file_instrument, omsi_file_sample, omsi_file_msidata
                import h5py
                if isinstance( value , omsi_file ) or isinstance( value, omsi_file_experiment ) or \
                   isinstance( value , omsi_file_sample) or isinstance( value , omsi_file_instrument ) or \
                   isinstance( value , omsi_file_analysis ) or isinstance( value , omsi_file_msidata) or \
                   isinstance( value , h5py.Dataset ) :
                   
                     dict.__setitem__(self, key , value )
            elif key == 'selection' :
                from omsi.shared.omsi_data_selection import check_selection_string
                if check_selection_string( str(value) ) :
                    
                    dict.__setitem__(self, key , str(value) )
                else :
                    raise ValueError( str(value) + " invalid selection string given to omsi_analysis_dependcy.")
            else :
                dict.__setitem__(self, key , value )
        else :
            raise KeyError("\'"+str(key)+'\' key not in default key set of omsi_analysis_dependcy' )
        
        
        