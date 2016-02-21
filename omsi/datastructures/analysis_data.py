"""
Helper module with data structures for managing analysis-related data.
"""

import warnings
import sys
import ast

import numpy as np
import h5py

from omsi.datastructures.dependency_data import dependency_dict
from omsi.shared.log import log_helper


#########################################################
#             data_types                                #
#########################################################
class data_dtypes(dict):
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
                  'bool': data_dtypes.bool_type,
                  'str': str,
                  'unicode': unicode,
                  'ndarray': data_dtypes.ndarray}
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
                isinstance(argument, h5py.Dataset) or isinstance(argument, h5py.Group) or \
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


#########################################################
#             analysis_data                             #
#########################################################
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
        if key in self or key in ['name', 'data', 'dtype']:  # The second check is to allow for pickle
            dict.__setitem__(self, key, value)
        else:
            raise KeyError("\'"+str(key)+'\' key not in default key set of analysis_data')


#########################################################
#             parameter_data                            #
#########################################################
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
                try:
                    warnings.warn('Conversion of parameter data to the expected dtype failed. ' +
                                  self['name'] + "  " + unicode(outdata) + "  " +
                                  unicode(self['dtype']) + "  " + unicode(sys.exc_info()))
                except (UnicodeDecodeError, UnicodeEncodeError):
                    warnings.warn('Conversion of parameter data to the expected dtype failed. ' + self['name']
                                  + "  " + unicode(self['dtype']) + "  " + unicode(sys.exc_info()))

        # Convert to ndarray if needed
        try:
            is_array_type = self['dtype'] == data_dtypes.ndarray
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


#########################################################
#             parameter_manager                         #
#########################################################
class parameter_manager(object):
    """
    Base class for objects that manage their own parameters.

    Parameters are set and their values retrieved by name using dict-like slicing. Derived classes
    may overwrite __getitem__ and __setitem__ to implement their own behavior but
    we exepct that the functionality of the interface is preserved, i.e., others should
    still be able set parameter value and retrieve values via dict slicing.
    """
    def __init__(self):
        """

        """
        super(parameter_manager, self).__init__()
        self.parameters = []

    def __len__(self):
        return self.get_num_parameter_data()

    def __getitem__(self, item):
        """
        Convenience function used to access analysis objects directly.
        Same as self.analysis_tasks.__getitem__
        :param item:
        :return: Output of self.analysis_tasks.__getitem__ implemented by omsi.workflow.common
        """
        if isinstance(item, basestring):
            for param in self.parameters:
                if param['name'] == item:
                    return param.get_data_or_default()

        raise KeyError('Invalid parameter key')

    def __setitem__(self, key, value):
        """
        Set parameter options directly via slicing

        Overwrite this function in child classes to implement custom setting behavior, e.g., error
        checking for valid values before setting a non-standard parameter.

        :param key: name of the parameters
        :param value: new value

        :raise: ValueError if an invalid value is given
        :raise: KeyError if an invalid key is given
        """
        # Check if we have a valid key
        param_set = False
        if isinstance(key, basestring):
            for param in self.parameters:
                if param['name'] == key:
                    log_helper.debug(__name__, "Setting parameter " + key)
                    param['data'] = value
                    param_set = True
        if not param_set:
            raise KeyError('Invalid parameter key')

    def keys(self):
        """
        Get a list of all valid keys, i.e., a list of all parameter names.

        :return: List of strings with all input parameter and output names.
        """
        return self.get_parameter_names()

    def get_all_parameter_data(self,
                               exclude_dependencies=False):
        """
        Get the complete list of all parameter datasets to be written to the HDF5 file

        :param exclude_dependencies: Boolean indicating whether we should exclude parameters
            that define dependencies from the list
        """
        if exclude_dependencies:
            return [param for param in self.parameters if not param.is_dependency()]
        else:
            return self.parameters

    def get_parameter_data(self,
                           index):
        """
        Given the index return the associated dataset to be written to the HDF5 file

        :param index : Return the index entry of the private member parameters. If a
            string is given, then get_parameter_data_by_name(...) will be used instead.

        :raises: IndexError is raised when the index is out of bounds
        """
        if isinstance(index, basestring):
            return self.get_parameter_data_by_name(index)
        else:
            return self.get_all_parameter_data()[index]

    def get_parameter_data_by_name(self,
                                   dataname):
        """
        Given the key name of the data return the associated parameter_data object.

        :param dataname: Name of the parameter requested from the parameters member.

        :returns: The parameter_data object or None if not found
        """
        for i in self.parameters:
            if i['name'] == dataname:
                return i
        return None

    def get_parameter_names(self):
        """
        Get a list of all parameter dataset names (including those that may define
        dependencies.
        """
        return [param['name'] for param in self.parameters]

    def get_num_parameter_data(self):
        """Return the number of parameter datasets to be wirtten to the HDF5 file"""
        return len(self.parameters)

    def get_num_dependency_data(self):
        """Return the number of dependencies defined as part of the parameters"""
        return len(self.get_all_dependency_data())

    def get_all_dependency_data(self):
        """
        Get the complete list of all direct dependencies to be written to the HDF5 file

        NOTE: These are only the direct dependencies as specified by the analysis itself.
        Use  get_all_dependency_data_recursive(..) to also get the indirect dependencies of
        the analysis due to dependencies of the dependencies themselves.

        :returns: List of parameter_data objects that define dependencies.

        """
        dependency_list = []
        for param in self.parameters:
            if param.is_dependency():
                dependency_list.append(param)
        return dependency_list

    def define_missing_parameters(self):
        """
        Set any required parameters that have not been defined to their respective default values.

        This function may be overwritten in child classes to customize
        the definition of default parameter values and to apply any
        modifications (or checks) of parameters before the analysis is executed.
        Any changes applied here will be recorded in the parameter of the analysis.
        """
        log_helper.debug(__name__, "Define missing parameters to default")
        for param in self.parameters:
            if param['required'] and not param.data_set():
                param['data'] = param['default']

    def add_parameter(self,
                      name,
                      help,
                      dtype=unicode,
                      required=False,
                      default=None,
                      choices=None,
                      data=None,
                      group=None):
        """
        Add a new parameter for the analysis. This function is typically used in the constructor
        of a derived analysis to specify the parameters of the analysis.

        :param name: The name of the parameter
        :param help: Help string describing the parameter
        :param dtype: Optional type. Default is string.
        :param required: Boolean indicating whether the parameter is required (True) or optional (False). Default False.
        :param default: Optional default value for the parameter. Default None.
        :param choices: Optional list of choices with allowed data values. Default None, indicating no choices set.
        :param data: The data assigned to the parameter. None by default.
        :param group: Optional group string used to organize parameters. Default None, indicating that
            parameters are automatically organized by driver class (e.g. in required and optional parameters)

        :raises: ValueError is raised if the parameter with the given name already exists.
        """
        log_helper.debug(__name__, "Add parameter " + str(name))
        if self.get_parameter_data_by_name(name) is not None:
            raise ValueError('A parameter with the name ' + unicode(name) + " already exists.")
        self.parameters.append(parameter_data(name=name,
                                              help=help,
                                              dtype=dtype,
                                              required=required,
                                              default=default,
                                              choices=choices,
                                              data=data,
                                              group=group))

    def clear_parameter_data(self):
        """Clear the list of parameter data"""
        log_helper.debug(__name__, "Clearing parameter data")
        for param in self.parameters:
            param.clear_data()

    def set_parameter_default_value(self, name, value):
        """
        Set the default value of the parameter with the given name

        :param name: Name of the parameter
        :param value: New value

        :raises: KeyError if parameter not found
        """
        log_helper.debug(__name__, "Setting default value of " +str(name) + " to " + str(value))
        param = self.get_parameter_data_by_name(dataname=name)
        if isinstance(param, parameter_data):
            param['default'] = value
        else:
            raise KeyError('Unknown parameter ' + str(name))



