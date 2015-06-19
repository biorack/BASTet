"""
Helper module with data structures for managing analysis-related data.
"""

import numpy as np
from omsi.shared.dependency_data import dependency_dict

import warnings
import sys
import ast
import h5py


class analysis_dtypes(dict):
    """
    Class specifying basic function for specifying common
    data types used as part of an analysis.
    """
    @staticmethod
    def get_dtypes():
        """
        Get a list of available data type specifications
        """
        dtypes = {'int': int,
                  'float': float,
                  'long': long,
                  'complex': complex,
                  'bool': analysis_dtypes.bool_type,
                  'str': str,
                  'unicode': unicode,
                  'ndarray': analysis_dtypes.ndarray}
        return dtypes

    @staticmethod
    def bool_type(argument):
        """
        Implement conversion of boolean input parameters since
        arparse (or bool, depending on the point of view), do not
        handle bool as a type in an intuitive fashion.

        :param argument: The argument to be parsed to a boolean
        :return: The converted value
        """
        try:
            bool(int(argument))
        except ValueError:
            if argument in ('TRUE', 'true', 'True', 't', 'T'):
               return True
            elif argument in ('FALSE', 'false', 'False', 'f', 'F'):
               return False
            else:
               raise ValueError('Parameter could not be converted to type bool')

    @staticmethod
    def ndarray(argument):
        """
        This dtype may be used to indicate numpy ndarrays as
        well as h5py arrays or omsi_dependencies

        :param argument: The argument to be parsed to ndarray

        :return: The converted ndarray
        """
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        from omsi.dataformat.omsi_file.msidata import omsi_file_msidata
        from omsi.dataformat.omsi_file.common import omsi_file_common
        if isinstance(argument, basestring):
            try:
                return np.asarray(ast.literal_eval(argument))
            except (ValueError, SyntaxError):
                omsi_out_object = omsi_file_common.get_omsi_object(h5py_object=argument)
                if omsi_out_object is not None:
                    return omsi_out_object
                else:
                    raise ValueError('String could not be converted to valid ndarray. This may be ' +
                                     'due to, e.g., a syntax error or the file may not exists')
        elif isinstance(argument, dependency_dict) or \
                isinstance(argument, h5py.Dataset) or \
                isinstance(argument, omsi_file_analysis) or \
                isinstance(argument, omsi_file_msidata):
            return argument
        elif argument is None:
            return None
        return np.asarray(argument)

    # @staticmethod
    # def numpy_ndarray(argument):
    #     """
    #     Use this dtypen to indicate that only numpy ndarrays are allowed, i.e.,
    #     h5py datasets, omsi_dependencies etc. are forbidden but only user-defined
    #     numpy arrays are allowed.
    #     """
    #     if isinstance(argument, basestring):
    #         try:
    #             outdata = ast.literal_eval(argument)
    #             return np.asarray(outdata)
    #         except ValueError:
    #             pass
    #         except SyntaxError:
    #             pass
    #     if isinstance(argument, dependency_dict):
    #         return argument.get_data()
    #     if isinstance(argument, omsi_file_analysis) or isinstance(argument, omsi_file_msidata):
    #         return argument[:]
    #     return np.asarray(argument)


class analysis_data(dict):
    """
    Define an output dataset for the analysis that should be written to the omsi HDF5 file
    """
    ana_hdf5link = -1
    """
    Value used to indicate that a hard link to another dataset should be created when saving an analysis object
    """

    def __init__(self, name="undefined", data=None, dtype='float32'):
        """
        The class can be used like a dictionary but restricts the set of keys that can be used
        to the following required keys which should be provided during initalization.

        **Required Keyword Arguments**:

        :param name: The name for the dataset in the HDF5 format
        :param data: The numpy array to be written to HDF5. The data write function
            omsi_file_experiment.create_analysis used for writing of the data to file can
            in principal also handel other primitive data types by explicitly converting them
            to numpy. However, in this case the dtype is determined based on the numpy conversion
            and correct behavior is not guaranteed. I.e., even single scalars should be stored as
            a 1D numpy array here. Default value is None which is mapped to np.empty( shape=(0) , dtype=dtype)
            in __init__
        :param dtype: The data type to be used during writing. For standard numpy data types this is just
             the dtype  of the dataset, i.e., ['data'].dtype. Other allowed datatypes are:

             - For string:  omsi_format.str_type (omsi_format is located in omsi.dataformat.omsi_file )
             - To generate data links: ana_hdf5link   (analysis_data)

        """
        super(analysis_data, self).__init__()
        if data is None:
            data = np.empty(shape=(0,), dtype=dtype)
        dict.__setitem__(self, 'name', name)
        dict.__setitem__(self, 'data', data)
        dict.__setitem__(self, 'dtype', dtype)

    def __setitem__(self, key, value):
        """
        Overwrite the __setitem__ function inherited from dict to ensure that only elements with a specific
        set of keys can be modified
        """
        if key in self:
            dict.__setitem__(self, key, value)
        else:
            raise KeyError("\'"+str(key)+'\' key not in default key set of analysis_data')


class parameter_data(dict):
    """
    Define a single input parameter for an analysis.

    :ivar default_keys: List of allowed dictionary keys:

    Required keys:

        * `name` : The name of the parameter
        * `help` : Help string describing the parameter
        * `type` : Optional type. Default is None, indicating a dynamically typed dataset that the analysis will convert
        * `required` : Boolean indicating whether the parameter is required (True) or optional (False). Default False
        * `default` : Optional default value for the parameter. Default None.
        * `choices` : Optional list of choices with allowed data values. Default None, indicating no choices set.
        * `data` : The data assigned to the parameter. None by default.
        * 'group' : Optional group string used to organize parameters. This may also be a dict of \
                   {'name':<group>, 'description':<description>}

    In the context of the argparse package the default keys have the following mapping:

        * `argparse.name` = `name`
        * `argparse.action` --> The action is constant and set to save value
        * `argparse.nargs`  --> Left as default
        * `argparse.const   --> Not used as action is always save value
        * `argparse.type` = `type`
        * `argparse.choices = `choices`
        * `argparse.required = `required`
        * `argparse.help = `help`
        * `argparse.metavar  --> Not used. Positional arguments are not allowed for analyses
        * `argparse.destination --> Automatically determined by the `name` of the parameter
        * `argparse.add_argument_group(...) --> Automatically determined based on the required parameter and \
           the `group` parameter if set.

    """
    default_keys = ['name',
                    'default',
                    'dtype',
                    'choices',
                    'required',
                    'help',
                    'data',
                    'group']
    """
    List of allowed keys for the parameter dict.
    """

    def __init__(self,
                 name,
                 help='',
                 dtype=None,
                 required=False,
                 default=None,
                 choices=None,
                 data=None,
                 group=None):
        """
        Initialize a new parameter description.

        :param name: Required name for the parameter
        :param help: Required help string for the parameter
        :param dtype: Type argument. Default unicode.
        :param required: Boolean indicating whether the parameter is required (default=True)
        :param default: Optional default value for the parameter. Default None.
        :param choices: Optional list of choices with allowed data values. Default None, indicating no choices set.
        :param data: The data assigned to the parameter. None by default.
        :param group: The parameter group to be used. None by default.

        """
        # Initialize the dict
        super(parameter_data, self).__init__()

        # Assign required keyword arguments
        self['name'] = name
        self['help'] = help
        self['dtype'] = dtype
        self['required'] = required
        self['default'] = default
        self['choices'] = choices
        self['data'] = data
        self['group'] = group

    def __setitem__(self, key, value):
        """
        Overwrite the default dict assignment to ensure that only valid keys are set.
        :param key:
        :param value:

        :raises: KeyError is raised in case that an invalid key is used.
        :raises: ValueError is raised in case that an invalid value is provided
        """
        if key in self.default_keys:
            if key == 'group' and value is not None:
                if not isinstance(value, basestring) and not isinstance(value, dict):
                    raise ValueError('Invalid group description for omsi_analysis_parameter')
                if isinstance(value, dict) and 'name' not in value:
                    raise ValueError('Invalid group description for omsi_analysis_parameter')
            dict.__setitem__(self, key, value)
        else:
            raise KeyError(unicode(key) +
                           " not in the list of allowed keys. Allowed keys are: " +
                           unicode(self.default_keys))

    def copy(self):
        """
        Return a new parameter_data object with the same data as stored in the current object

        :return: dependency_dict object
        """
        new_parameter = parameter_data('')
        new_parameter.update(self)
        return new_parameter

    def data_set(self):
        """
        Check if a data has been assigned for the parameter.
        """
        return self['data'] is not None

    def data_ready(self):
        """
        This function check if the data points to a dependency and if so, then check if the dependency can be
        resolved or not
        """
        if self.is_dependency():
            if isinstance(self['data'].get_data(), dependency_dict):
                return False
        return True

    def is_dependency(self):
        """
        Check whether the parameter defines a dependency.

        :return: Boolean indicating whether the parameter defines a dependency.
        """
        return isinstance(self['data'], dependency_dict)

    def get_group_name(self):
        """
        Get the name of the group to be used.

        :return: String with the name of the group of None if not set
        """
        if isinstance(self['group'], dict) and 'name' in self['group']:
            return self['group']['name']
        else:
            return self['group']

    def get_group_description(self):
        """
        Get the description for the group if available.

        :return: String with the group description or None.
        """
        if isinstance(self['group'], dict) and 'description' in self['group']:
            return self['group']['description']
        else:
            return None

    def get_data_or_default(self):
        """
        Get the data of the parameter if set, otherwise get the default value if available.

        :return: The data to be used for the parameter.

        :raises: KeyError is raised in case that neither 'default' nor 'data' are available.
            This should never be the case if the object was created properly.
        """
        get_data = True  # Should we get the 'data' (True) or 'default' key (False)
        if 'data' in self.keys():
            if 'default' in self.keys() and self['default'] is not None and self['data'] is None:
                get_data = False
        elif 'default' in self.keys():
            get_data = False
        else:
            raise KeyError('No data or default setting available for the parameter ' + self['name'])

        if get_data:
            if isinstance(self['data'], dependency_dict):
                outdata = self['data'].get_data()
            else:
                outdata = self['data']
        else:
            outdata = self['default']

        if not self['required'] and outdata is None:
            return outdata

        # Convert the data to the proper type if needed
        # Convert the data to a scalar type if needed
        scalar_types = [int, float, long, complex, bool, str, unicode]
        if self['dtype'] in scalar_types:
            curr_dtype = scalar_types[scalar_types.index(self['dtype'])]
            try:
                outdata = curr_dtype(outdata)
            except:
                # raise
                warnings.warn('Conversion of parameter data to the expected dtype failed. ' +
                              self['name'] + "  " + str(outdata) + "  " + unicode(sys.exc_info()))

        # Convert to ndarray if needed
        try:
            is_array_type = self['dtype'] == analysis_dtypes.ndarray
        except:
            try:
                is_array_type = self['dtype'] == 'ndarray'
            except:
                is_array_type = False
        if is_array_type:
            outdata = self['dtype'](outdata)

        # Retrieve scalar from numpy scalar array
        if isinstance(outdata, np.ndarray) and outdata.shape == ():
            outdata = outdata[()]

        return outdata

    def clear_data(self):
        """
        Remove the currently assigned data.
        """
        self['data'] = None
