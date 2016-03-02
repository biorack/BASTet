"""Define a dependency to another omsi object"""

# from omsi.dataformat.omsi_file import *
from omsi.dataformat.omsi_file.common import omsi_file_common
import h5py
import warnings
from omsi.shared.log import log_helper


class dependency_dict(dict):
    """
    Define a dependency to another omsi file-based data object or in-memory analysis_base object

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

    **Optional Keyword arguments**:

    :ivar dependency_type: The type of the dependency being modeled. If not defined then
        the default value of 'parameter' is assumed.

    """
    dependency_types = {'parameter': 'parameter',      # Default value, defining that this dependency is a parameter
                        'link': 'link',                # A not further defined link to a related dataset
                        'co_modality': 'co_modality',  # Data acquired from a related data modality
                        'subset': 'subset',            # The source is a subset of the object we link to
                        'contains': 'contains',        # The dataset we link to is contained in the source
                        'undefined': None              # Undefined dependency type
                        }

    def __init__(self,
                 param_name=None,
                 link_name=None,
                 omsi_object=None,
                 selection=None,
                 dataname=None,
                 help=None,
                 dependency_type=None):
        """
        Initialize the allowed set of keys.

        :param param_name: The name of the parameter that has the dependency
        :param link_name: The name of for the link to be created in the HDF5 file.
        :param omsi_object: The object to which a link should be established to. This
            must be either an h5py.Dataset or the omsi_file_analysis or omsi_file_msidata
            or any of the other omsi_file API  interface objects.
        :param selection: Optional string type parameter indicating a python selection for the dependency
        :param dataname: String indicating the dataset within the omsi_object. If the omsi_object
            is an h5py object within a managed Group, then the omsi_object is automatically split
            up into the parent object and dataname.
        :param help: Optional string describing the object

        """
        super(dependency_dict, self).__init__()
        # Add all the keys first
        dict.__setitem__(self, 'param_name', param_name)
        dict.__setitem__(self, 'link_name', link_name)
        # if 'omsi_object' not in self.keys():
        dict.__setitem__(self, 'omsi_object', None)
        # if 'dataname' not in self.keys():
        dict.__setitem__(self, 'dataname', None)
        dict.__setitem__(self, 'selection', None)
        dict.__setitem__(self, '_data', None)
        dict.__setitem__(self, 'help', '')
        dict.__setitem__(self, 'dependency_type', dependency_type)

        # Set all the dataname and omsi_object keys to their appropriate values
        if omsi_object is not None:
            self.__setitem__('omsi_object', omsi_object)
        if dataname is not None:
            self.__setitem__('dataname', dataname)
        if selection is not None and len(str(selection)) > 0:  # Ignore None and empty string
            self.__setitem__('selection', selection)
        if help is not None:
            self.__setitem__('help', help)

    def copy(self):
        """
        Return a new dependency_dict object with the same data as stored in the current object

        :return: dependency_dict object
        """
        new_dependency = dependency_dict()
        for key, value in self.iteritems():
            if key == '_data':
                new_dependency._force_set_data(value)
            else:
                # Force setting all keys when we copy as we may need to copy a partially defined link
                super(dependency_dict, new_dependency).__setitem__(key, value)
                # new_dependency[key] = value
        return new_dependency

    def __setitem__(self,
                    key,
                    value):
        """Overwrite the __setitem__ function inherited from dict to ensure that only elements with a specific
           set of keys can be modified"""
        from omsi.analysis.base import analysis_base
        from omsi.dataformat.file_reader_base import file_reader_base
        if key in self:
            if key == "omsi_object":
                if omsi_file_common.is_managed(value):
                    dict.__setitem__(self, key, omsi_file_common.get_omsi_object(value))
                elif isinstance(value, h5py.Dataset) or isinstance(value, h5py.Group):
                    parent = value.parent
                    if omsi_file_common.is_managed(parent):
                        dict.__setitem__(self, 'omsi_object', omsi_file_common.get_omsi_object(parent))
                        dict.__setitem__(self, 'dataname', unicode(value.name.split('/')[-1]))
                        # print super(dependency_dict,self).__str__()
                    else:
                        warnings.warn("The generated dependency does not point to a managed object.")
                        dict.__setitem__(self, 'omsi_object', omsi_file_common.get_omsi_object(value))
                    dict.__setitem__(self, '_data', None)  # Any previously loaded date may be invalid (delete)
                elif isinstance(value, analysis_base):
                    dict.__setitem__(self, 'omsi_object', value)
                else:
                    raise ValueError(str(value) +
                                     " invalid omsi_object parameter for " +
                                     "dependency_dict without valid data dependency.")
            elif key == 'selection':
                if value is None or (isinstance(value, basestring) and len(value) == 0):
                    new_value = None
                else:
                    from omsi.shared.data_selection import selection_to_string
                    new_value = unicode(selection_to_string(selection=value))
                dict.__setitem__(self, key, new_value)
                dict.__setitem__(self, '_data', None)  # Any previously loaded data may be invalid (delete)
            elif key == 'dataname':
                if not isinstance(value, basestring):
                    raise ValueError('Dataname must be a string')
                dict.__setitem__(self, 'dataname', unicode(value))
                dict.__setitem__(self, '_data', None)  # Any previously loaded data may be invalid (delete)
            elif key == 'param_name':
                if not isinstance(value, basestring):
                    raise ValueError('param_name must be a string')
                dict.__setitem__(self, 'param_name', unicode(value))
            elif key == 'link_name':
                if not isinstance(value, basestring):
                    raise ValueError('link_name must be a string')
                dict.__setitem__(self, 'link_name', unicode(value))
            elif key == '_data':
                raise KeyError('_data key is managed by dependency_dict. Explicit definition of _data not permitted.')
            elif key == 'help':
                if isinstance(value, basestring):
                    dict.__setitem__(self, 'help', unicode(value))
            elif key == 'dependency_type':
                if value in self.dependency_types.values():
                    dict.__setitem__(self, 'dependency_type', value)
                else:
                    raise ValueError('Unknown dependency type specified. Valid types are: ' +
                                     str(self.dependency_types))
            else:
                dict.__setitem__(self, key, value)
            # print super(dependency_dict,self).__str__()
        else:
            raise KeyError("\'"+str(key)+'\' key not in default key set of dependency_dict')

    def __getitem__(self,
                    key):
        """Custom slicing. Return the value associated with the given key if it is one of our predefined keys.
           Otherwise, assume that the user wants to slice into the data associated with the dependency and
           return get_data()[key] instead. If get_data() is not ready (i.e., creates a dependency), then a
           new dependency_dict is returned that adds the selection.

           :param key: The key to be used for slicing

           :returns: Value. May return a dependency_dict if the selection refers to the data object and the
                dependency cannot be resolved yet.
        """
        if key in self.keys():
            return dict.__getitem__(self, key)
        else:
            data_ref = self.get_data()
            if isinstance(data_ref, dependency_dict):
                if self['selection'] is not None:
                    log_helper.error(__name__, "The current dependency already has a selection. Refinement of " +
                                     "existing selections is not yet supported. A new dependency with the full " +
                                     "current selection will be used instead.")
                copy_ref = self.copy()
                copy_ref['selection'] = key
                return copy_ref
            else:
                return self.get_data()[key]

    def _force_set_data(self, data):
        """
        Force setting of the _data key to the given data value. Use with great caution. This function
        is used internally when a dependency is created for which we have a pointer to the data ready.

        :param data: The data object to be used as the data value for the dependency.
        """
        dict.__setitem__(self, '_data', data)

    def get_data(self):
        """Get the data associated with the dependency.

           :returns: If a selection is applied and the dependency object supports
                     array data load (e.g., h5py.Dataset, omsi_file_msidata), then
                     the selected data will be loaded and returned as numpy array.
                     Otherwise the ['omsi_object'] is returned.
        """
        # Return preloaded data if available
        if self['_data'] is not None:
            return self['_data']
        # Check if we can access the data object
        else:
            # Retrieve the data object
            if self['dataname']:
                data_object = self['omsi_object'][self['dataname']]
                # Ensure that the dependency can actually be resolved. E.g, if the data the dependency points to
                # is not ready yet then we may get a dependency back that points to the same object as we do,
                # which in turn could result in an endless recursion
                if isinstance(data_object, dependency_dict):
                    if data_object['omsi_object'] is self['omsi_object']:
                        return self
            else:
                data_object = self['omsi_object']
            # Resolve any data seletions
            try:
                if self['selection'] is None:
                    # data = data_object[:]
                    # self['_data'] = data_object
                    return data_object
                else:
                    from omsi.shared.data_selection import selection_string_to_object
                    current_selection = selection_string_to_object(self['selection'])
                    if current_selection is not None:
                        dict.__setitem__(self, '_data', data_object[current_selection])
                    else:
                        raise ValueError('Invalid selection string')
                return self['_data']
            except:
                raise
                import sys
                log_helper.error(__name__, "Application of data selection failed. " + str(sys.exc_info()))
                return data_object
