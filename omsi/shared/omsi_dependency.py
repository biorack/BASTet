"""Define a dependency to another omsi object"""

# from omsi.dataformat.omsi_file import *
from omsi.dataformat.omsi_file.common import omsi_file_common
import h5py
import warnings


class omsi_dependency(dict):
    """
    Define a dependency to another omsi file-based data object

    **Required Keyword Arguments**:

    :ivar param_name: The name of the parameter that has the depency
    :ivar link_name: The name of for the link to be created in the HDF5 file.
    :ivar omsi_object: The object to which a link should be established to. This
        must be either an h5py.Dataset or the omsi_file_analysis or omsi_file_msidata
        or any of the other omsi_file API  interface ojects.
    :ivar selection: Optional string type parameter indicating a python selection for the dependency
    :ivar dataname: String indicating the dataset within the omsi_object. If the omsi_object
        is an h5py object within a managed Group, then the omsi_object is automatically split
        up into the parent object and dataname.
    :ivar _data: Private key used to store the data associated with the dependency object.

    """
    def __init__(self,
                 param_name=None,
                 link_name=None,
                 omsi_object=None,
                 selection=None,
                 dataname=None):
        """
        Initialize the allowed set of keys.

        :param param_name: The name of the parameter that has the dependency
        :param link_name: The name of for the link to be created in the HDF5 file.
        :param omsi_object: The object to which a link should be established to. This
            must be either an h5py.Dataset or the omsi_file_analysis or omsi_file_msidata
            or any of the other omsi_file API  interface ojects.
        :param selection: Optional string type parameter indicating a python selection for the dependency
        :param dataname: String indicating the dataset within the omsi_object. If the omsi_object
            is an h5py object within a managed Group, then the omsi_object is automatically split
            up into the parent object and dataname.

        """
        super(omsi_dependency, self).__init__()
        # Add all the keys first
        dict.__setitem__(self, 'param_name', param_name)
        dict.__setitem__(self, 'link_name', link_name)
        # if 'omsi_object' not in self.keys():
        dict.__setitem__(self, 'omsi_object', None)
        # if 'dataname' not in self.keys():
        dict.__setitem__(self, 'dataname', None)
        dict.__setitem__(self, 'selection', None)
        dict.__setitem__(self, '_data', None)

        # Set all the dataname and omsi_object keys to their appropriate values
        if omsi_object is not None:
            self.__setitem__('omsi_object', omsi_object)
        if dataname is not None:
            self.__setitem__('dataname', dataname)
        if selection is not None and len(selection) > 0:  # Ignore None and empty string
            self.__setitem__('selection', selection)

    def __setitem__(self,
                    key,
                    value):
        """Overwrite the __setitem__ function inherited from dict to ensure that only elements with a specific
           set of keys can be modified"""

        if key in self:
            if key == "omsi_object":
                if omsi_file_common.is_managed(value):
                    dict.__setitem__(self, key, omsi_file_common.get_omsi_object(value))
                elif isinstance(value, h5py.Dataset) or isinstance(value, h5py.Group):
                    parent = value.parent
                    if omsi_file_common.is_managed(parent):
                        dict.__setitem__(self, 'omsi_object', omsi_file_common.get_omsi_object(parent))
                        dict.__setitem__(self, 'dataname', unicode(value.name.split('/')[-1]))
                        # print super(omsi_dependency,self).__str__()
                    else:
                        print "WARNING: The generated dependency does not point to a managed object."
                        dict.__setitem__(self, 'omsi_object', omsi_file_common.get_omsi_object(parent))
                    dict.__setitem__(self, '_data', None)  # Any previously loaded date may be invalid (delete)
                else:
                    raise ValueError(str(value) +
                                     " invalid omsi_object parameter for "
                                     + "omsi_dependency without valid data dependency.")
            elif key == 'selection':
                if value is None or (isinstance(value, basestring) and len(value) == 0):
                    new_value = None
                else:
                    from omsi.shared.omsi_data_selection import selection_to_string
                    new_value = unicode(selection_to_string(selection=value))
                dict.__setitem__(self, key, new_value)
                dict.__setitem__(self, '_data', None)  # Any previously loaded date may be invalid (delete)
            elif key == 'dataname':
                if not isinstance(value, basestring):
                    raise ValueError('Dataname must be a string')
                dict.__setitem__(self, 'dataname', unicode(value))
                dict.__setitem__(self, '_data', None)  # Any previously loaded date may be invalid (delete)
            elif key == 'param_name':
                if not isinstance(value, basestring):
                    raise ValueError('param_name must be a string')
                dict.__setitem__(self, 'param_name', unicode(value))
            elif key == 'link_name':
                if not isinstance(value, basestring):
                    raise ValueError('link_name must be a string')
                dict.__setitem__(self, 'link_name', unicode(value))
            elif key == '_data':
                raise KeyError('_data key is managed by omsi_dependency. Explicit definition of _data not permitted.')
            else:
                dict.__setitem__(self, key, value)
            # print super(omsi_dependency,self).__str__()
        else:
            raise KeyError("\'"+str(key)+'\' key not in default key set of omsi_dependency')

    def __getitem__(self,
                    key):
        """Custom slicing. Return the value associated with the given key if it is one of our predefined keys.
           Otherwise, assume that the user wants to slice into the data associated with the dependency and
           return get_data()[key] instead.

           :param key: The key to be used for slicing

           :returns: Value.
        """
        if key in self.keys():
            return dict.__getitem__(self, key)
        else:
            return self.get_data()[key]

    def get_data(self):
        """Get the data associated with the dependency.

           :returns: If a selection is applied and the dependency object supports
                     array data load (e.g., h5py.Dataset, omsi_file_msidata), then
                     the selected data will be loaded and returned as numpy array.
                     Otherwise the ['omsi_object'] is returned.
        """
        if self['_data'] is not None:
            return self['_data']
        else:
            if self['dataname']:
                data_object = self['omsi_object'][self['dataname']]
            else:
                data_object = self['omsi_object']
            try:
                if self['selection'] is None:
                    # data = data_object[:]
                    # self['_data'] = data_object
                    return data_object
                else:
                    from omsi.shared.omsi_data_selection import selection_string_to_object
                    current_selection = selection_string_to_object(self['selection'])
                    if current_selection is not None:
                        self['_data'] = data_object[current_selection]
                    else:
                        raise ValueError('Invalid selection string')
                return self['_data']
            except:
                import sys
                warnings.warn("ERROR: "+str(sys.exc_info()))
                return data_object
