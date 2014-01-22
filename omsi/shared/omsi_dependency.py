from omsi.dataformat.omsi_file import *
import h5py  

class omsi_dependency(dict) :
    """Define a dependency between an analysis and another omsi data object
    
        **Required Keyword Arguments**:

           :param param_name: The name of the parameter that has the depency
           :param link_name: The name of for the link to be created in the HDF5 file.
           :param omsi_object: The object to which a link should be established to. This must be either an h5py.Dataset or the omsi_file_analysis or omsi_file_msidata or any of the other omsi_file API  interface ojects.
           :param selection: Optional string type parameter indicating a python selection for the dependency
           :param dataname: String indicating the dataset within the omsi_object. If the omsi_object is an h5py object within a managed Group, then the omsi_object is automatically split up into the parent object and dataname. 
           :param _data: Private key used to store the data associated with the dependency object. 

    """
    def __init__(self, param_name=None , link_name=None, omsi_object=None, selection=None , dataname=None) :
        """Initatlize the allowed set of keys"""
        super(omsi_dependency,self).__init__()
        #Add all the keys first
        dict.__setitem__( self, 'param_name'  , param_name )
        dict.__setitem__( self, 'link_name'  , link_name )
        #if 'omsi_object' not in self.keys() :
        dict.__setitem__( self, 'omsi_object' , None )
        #if 'dataname' not in self.keys() :
        dict.__setitem__(self, 'dataname' , None )
        dict.__setitem__( self ,'selection' , selection )
        dict.__setitem__( self, '_data' , None)
        
        #Set all the dataname and omsi_object keys to their approbirate values
        if omsi_object is not None :
            self.__setitem__( 'omsi_object' , omsi_object )
        if dataname is not None :
            self.__setitem__(self, 'dataname' , dataname)
        
    def __setitem__(self, key, value ) :
        """Overwrite the __setitem__ function inheritied from dict to ensure that only elements with a specific
           set of keys can be modified"""
        
        if self.has_key(key) :
            if key == "omsi_object" :
                if omsi_file.is_managed( value ) :
                    dict.__setitem__(self, key , omsi_file.get_omsi_object( value ) )
                elif isinstance( value , h5py.Dataset ) or isinstance( value, h5py.Group) : 
                    parent = value.parent
                    if omsi_file.is_managed( parent ) :
                        dict.__setitem__( self, 'omsi_object' , omsi_file.get_omsi_object( parent ) )
                        dict.__setitem__( self, 'dataname'  , unicode(value.name.split('/')[-1])  )
                        #print super(omsi_dependency,self).__str__()
                    else :
                        print "WARNING: The generated dependency does not point to a managed object."
                        dict.__setitem__( self, omsi_object , omsi_file.get_omsi_object( parent ) ) 
                    dict.__setitem__( self, '_data' , None) #Any previously loaded date may be invalide (delete)
                else :
                    raise ValueError( str(value) + " invalid omsi_object parameter for omsi_dependcy without valid data dependency.")
            elif key == 'selection' :
                from omsi.shared.omsi_data_selection import check_selection_string
                if check_selection_string( str(value) ) :
                    dict.__setitem__(self, key , unicode(value) )
                    dict.__setitem__( self, '_data' , None) #Any previously loaded date may be invalide (delete)
                else :
                    raise ValueError( str(value) + " invalid selection string given to omsi_analysis_dependcy.")
            elif key == 'dataname' :
                dict.__setitem__(self, 'dataname' , unicode(value))
                dict.__setitem__( self, '_data' , None) #Any previously loaded date may be invalide (delete)
            else :
                dict.__setitem__(self, key , value )
            #print super(omsi_dependency,self).__str__()
        else :
            raise KeyError("\'"+str(key)+'\' key not in default key set of omsi_analysis_dependcy' )

    def __getitem__( self, key ) :
        """Custom slicing. Return the value associated with the given key if it is one of our predefined keys.
           Otherwise, assume that the user wants to slice into the data assoicated with the dependency and 
           return get_data()[key] instead.
           
           :param key: The key to be used for slicing 
           
           :retuns: Value.
        """
        if key in self.keys() :
            return dict.__getitem__(self, key)
        else :
            return self.get_data()[key]

    def get_data(self) :
        """Get the data associated with the dependency.
        
           :returns: For dependencies that support data load (e.g., h5py.Dataset, omsi_file_msidata) 
                     the associate data is loaded and the numpy data array is returned.
                     Otherwise the ['omsi_object'] is returned. 
        """
        if self['_data'] :
            return self['_data']
        else :
            if self['dataname'] :
                dataObj = self['omsi_object'][ self['dataname'] ]
            else :
                dataObj = self['omsi_object']
            try :
                if not self['selection'] :
                    data = dataObj[:]
                    self['_data'] = data
                else :
                    from omsi.shared.omsi_data_selection import check_selection_string
                    selection_str = self['selection']
                    if check_selection_string( str(selection_str) ) :
                        data = eval( "dataObj[ "+selection_str+" ]" )
                        self['_data'] = data 
                    else :
                        self['data'] = None
                        raise ValueError('Invalid selection string')
                return self['_data']
            except :
                import sys
                print "ERROR: "+str(sys.exc_info())
                return dataObj
                

        
        
        