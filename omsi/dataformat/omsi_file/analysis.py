"""
Module for managing custom analysis data in OMSI HDF5 files.
"""
import numpy as np
import h5py

from omsi.dataformat.omsi_file.format import \
    omsi_format_common, \
    omsi_format_analysis, \
    omsi_format_dependencies
from omsi.dataformat.omsi_file.dependencies import omsi_dependencies_manager
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager
import warnings


class omsi_analysis_manager(omsi_file_object_manager):
    """
    Analysis manager helper class used to define common functionality needed for analysis-related data.
    Usually, a class that defines a format that contains an omsi_file_analysis object
    will inherit from this class (in addition to omsi_file_common) to acquire the common
    features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar analysis_parent: The h5py.Group parent object containing the instrument object to be managed.

    """
    def __init__(self, analysis_parent):
        super(omsi_analysis_manager, self).__init__(analysis_parent)
        self.analysis_parent = analysis_parent

    def create_analysis(self, analysis, flush_io=True):
        """
        Add a new group for storing derived analysis results for the current experiment

        Create the analysis group using omsi_file_analysis.__create___ which in turn uses
        omsi_file_analysis.__populate_analysis__(...) to populate the group with the appropriate data.

        :param analysis: Instance of omsi.analysis.omsi_analysis_base defining the analysis
        :type analysis: omsi.analysis.omsi_analysis_base:
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                       so that all data has been written to file

        :returns: The omsi_file_analysis object for the newly created analysis group and the
                integer index of the analysis
        """
        return omsi_file_analysis.__create__(parent_group=self.analysis_parent,
                                             analysis=analysis,
                                             analysis_index=None,  # Same as self.get_num_analysis()
                                             flush_io=flush_io)

    def get_num_analysis(self):
        """
        Get the number of raw mass spectrometry images stored for a given experiment

        :returns: Integer indicating the number of analyses available for the experiment.
        """
        return omsi_file_common.get_num_items(self.analysis_parent, omsi_format_analysis.analysis_groupname)

    def get_analysis_identifiers(self):
        """
        Get a list of all identifiers for all analysis stored for the experiment

        :returns: List of strings of analysis identifiers.
        """
        output_list = []
        # Iterate through all groups of the root folder
        for item_obj in self.analysis_parent.items():
            if item_obj[0].startswith(omsi_format_analysis.analysis_groupname):
                cur_ana_id = omsi_file_analysis(
                    self.analysis_parent[item_obj[0]]).get_analysis_identifier()
                if cur_ana_id is not None:
                    output_list.append(cur_ana_id[0])
                else:
                    output_list.append("")

        return output_list

    def get_analysis(self, analysis_index):
        """
        Get the omsi_format_analysis analysis object for the experiment with
        the given index.

        :param analysis_index: The index of the analysis
        :type analysis_index: Unsigned integer

        :returns: omsi_file_analysis object for the requested analysis. The function
                      returns None in case the analysis object was not found.
        """
        analysis_name = unicode(omsi_format_analysis.analysis_groupname + str(analysis_index))
        if analysis_name in self.analysis_parent.keys():
            return omsi_file_analysis(self.analysis_parent[analysis_name])
        else:
            return None

    def get_analysis_by_identifier(self, analysis_identifier_string):
        """
        Get the omsi_format_analysis analysis object for the the analysis with the given identifier.

        :param analysis_identifier_string: The string used as identifier for the analysis.
        :type analysis_identifier_string: string

        :returns: h5py obejct of the analysis or None in case the analysis is not found.
        """

        # Iterate through all groups of the root folder
        for item_obj in self.analysis_parent.items():
            if item_obj[0].startswith(omsi_format_analysis.analysis_groupname):
                cur_ana_id = omsi_file_analysis(self.analysis_parent[item_obj[0]]).get_analysis_identifier()
                if cur_ana_id is not None:
                    if cur_ana_id[0] == analysis_identifier_string:
                        return omsi_file_analysis(self.analysis_parent[item_obj[0]])
        return None


class omsi_file_analysis(omsi_dependencies_manager,
                         omsi_file_common):
    """
    Class for managing analysis specific data in omsi hdf5 files
    """
    @classmethod
    def __create__(cls,
                   parent_group,
                   analysis,
                   analysis_index=None,
                   flush_io=True):
        """
        Add a new group for storing derived analysis results for the given parent object

        Create the analysis group and use omsi_file_analysis.__populate_analysis__(...) to populate the
        group with the appropriate data.

        :param parent_group: The h5py.Group where the method object should be created
        :type parent_group: h5py.Group
        :param analysis_index: The integer index of the analysis to be created. Set to None if the
                         function should determine the index automatically.
        :param analysis: Instance of omsi.analysis.omsi_analysis_base defining the analysis
        :type analysis: omsi.analysis.omsi_analysis_base:
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                       so that all data has been written to file
        :type flush_io: True

        :returns: i) the omsi_file_analysis object for the newly created analysis group ii) the integer index
                  for the analysis.

        :raises: IndexError is raised in case that an analysis with the given index already exists.
        """
        # 1 Import required modules
        from omsi.analysis.omsi_analysis_base import omsi_analysis_base
        import time

        # 2 Check if input values are correct
        # 2.1 Confirm that we have a valid analysis object
        if not isinstance(analysis, omsi_analysis_base):

            errormessage = "Could not write the analysis. The input object was of type " + \
                str(type(analysis)) + " not of omsi_analysis_base as expected."
            raise NameError(errormessage)

        # 2.2 Create a name for the group and check if it is valid
        # 2.2.1 Determine the current number of analyses
        if analysis_index is None:
            analysis_index = omsi_file_common.get_num_items(parent_group,
                                                            omsi_format_analysis.analysis_groupname)
        # 2.2.2 Create the name of the group and check if it already exists
        analysis_group_name = omsi_format_analysis.analysis_groupname + str(analysis_index)
        if analysis_group_name in parent_group.keys():
            raise IndexError("An analysis with the requested index already exists.")

        # 3 Create the group and define the required attributes
        analysis_group = parent_group.require_group(analysis_group_name)
        analysis_group.attrs[omsi_format_common.type_attribute] = "omsi_file_analysis"
        analysis_group.attrs[omsi_format_common.version_attribute] = omsi_format_analysis.current_version
        analysis_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())

        # Populate the group with data and return the omsi_file_analysis object
        analysis_object = omsi_file_analysis.__populate_analysis__(analysis_group, analysis)

        # Flush I/O if necessary
        if flush_io:
            parent_group.file.flush()

        return analysis_object, analysis_index

    @classmethod
    def __populate_analysis__(cls,
                              analysis_group,
                              analysis):
        """
        Populate the given h5py group with the analysis data.

        NOTE: This is a private helper function. Use the corresponding create_analysis function
        of omsi_file_experiment to create a completely new analysis.

        :param analysis_group: h5py group in which the analysis data should be stored.
        :param analysis: Instance of omsi.analysis.omsi_analysis_base defining the analysis
        :type analysis: omsi.analysis.omsi_analysis_base:

        :returns: The omsi_file_analysis object for the newly created analysis group. The analysis data is
                  automatically written to file by this function so no addition work is required.

        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies

        # Write the analysis name
        analysis_identifier_data = analysis_group.require_dataset(
            name=unicode(omsi_format_analysis.analysis_identifier),
            shape=(1,),
            dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysis_identifier_data[0] = analysis.get_analysis_identifier()
        else:
            analysis_identifier_data[0] = str(analysis.get_analysis_identifier())

        # Write the analysis type
        analysis_type_data = analysis_group.require_dataset(name=unicode(omsi_format_analysis.analysis_type),
                                                            shape=(1,),
                                                            dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysis_type_data[0] = analysis.get_analysis_type()
        else:
            analysis_type_data[0] = str(analysis.get_analysis_type())

        # Write the analysis data
        for ana_data in analysis.get_all_analysis_data():
            cls.__write_omsi_analysis_data__(analysis_group, ana_data)

        # Write all the parameters
        parameter_group = analysis_group.require_group(omsi_format_analysis.analysis_parameter_group)
        for param_data in analysis.get_all_parameter_data(exclude_dependencies=True):
            if param_data['required'] or \
                    param_data['default'] is not None or\
                    param_data.data_set():
                cls.__write_omsi_analysis_data__(parameter_group, param_data)
                # Try to add the help string attribute
                try:
                    help_attr = omsi_format_analysis.analysis_parameter_help_attr
                    parameter_group[param_data['name']].attrs[help_attr] = param_data['help']
                except KeyError:
                    pass

        # Write all the runtime execution information
        runinfo_group = analysis_group.require_group(omsi_format_analysis.analysis_runinfo_group)
        for run_info_key, run_info_value in analysis.get_all_run_info().items():
            # Generate an omsi_analysis_data object in order to use the
            # __write_omsi_analysis_data function to write the data
            if isinstance(run_info_value, unicode) or isinstance(run_info_value, str):
                anadata = omsi_analysis_data(name=unicode(run_info_key),
                                             data=run_info_value,
                                             dtype=omsi_format_common.str_type)
            else:
                dat = np.asarray(run_info_value)
                if len(dat.shape) == 0:
                    dat = np.asarray([run_info_value])
                anadata = omsi_analysis_data(name=unicode(run_info_key),
                                             data=dat,
                                             dtype=dat.dtype)
            cls.__write_omsi_analysis_data__(runinfo_group, anadata)

        # Write all data dependencies
        dependencies = [dep['data'] for dep in analysis.get_all_dependency_data()]
        omsi_file_dependencies.__create__(parent_group=analysis_group,
                                          dependencies_data_list=dependencies)

        # Execute the custom data write for the analysis
        analysis.write_to_omsi_file(analysis_group)

        return omsi_file_analysis(analysis_group)

    @classmethod
    def __write_omsi_analysis_data__(cls,
                                     data_group,
                                     ana_data):
        """
        Private helper function used to write the data defined by a omsi_analysis_data object to HDF5.

        :param data_group: The h5py data group to which the data should be written to.
        :param ana_data: The omsi_analysis_data object with the description of the data to be written.
        :type ana_data: omsi.analysis.omsi_analysis_data
        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data, omsi_parameter_data
        # Create link in HDF5 to an existing dataset within the file
        if isinstance(ana_data, omsi_analysis_data) and isinstance(ana_data['dtype'], int):
            if ana_data['dtype'] == ana_data.ana_hdf5link:
                linkobject = data_group.file.get(ana_data['data'])
                data_group[ana_data['name']] = linkobject
                omsiobj = omsi_file_common.get_omsi_object(linkobject)
                try:
                    # Check if we already have a type attribute
                    _ = data_group[ana_data['name']].attrs[omsi_format_common.type_attribute]
                except:
                    # Generate the type attribute from scratch
                    if omsiobj is not None:
                        omsiobjtype = omsiobj.__class__.__name__
                    else:
                        omsiobjtype = ""
                    data_group[ana_data['name']].attrs[
                        omsi_format_common.type_attribute] = omsiobjtype
        # Create a new string-type dataset
        elif (ana_data['dtype'] == omsi_format_common.str_type) or (ana_data['dtype'] == h5py.special_dtype(vlen=str)):
            tempdata = data_group.require_dataset(name=unicode(ana_data['name']),
                                                  shape=(1,),
                                                  dtype=omsi_format_common.str_type)
            if len(unicode(ana_data['data'])) > 0:
                if omsi_format_common.str_type_unicode:
                    tempdata[0] = ana_data['data']
                else:
                    tempdata[0] = str(ana_data['data'])
            else:
                warnings.warn("WARNING: " + ana_data['name'] +
                              " dataset generated but not written. The given dataset was empty.")
        # Create a new dataset to store the current numpy-type dataset
        elif 'numpy' in str(type(ana_data['data'])):
            # Decide whether we want to enable chunking for the current
            # analysis dataset
            chunks = None
            if ana_data['data'].size > 1000:
                chunks = True
            # Write the current analysis dataset
            tempdata = data_group.require_dataset(name=ana_data['name'],
                                                  shape=ana_data['data'].shape,
                                                  dtype=ana_data['dtype'],
                                                  chunks=chunks)
            if ana_data['data'].size > 0:
                tempdata[:] = ana_data['data']
            else:
                warnings.warn("WARNING: " + ana_data['name'] +
                              " dataset generated but not written. The given dataset was empty.")
        # Unkown dtype. Attempt to convert the dataset to numpy and write it to
        # file.
        else:
            # Savely convert scalars to numpy but warn in case we see something else
            if ana_data['dtype'] not in [int, float, long, complex, bool, str, unicode,
                                         'int', 'float', 'long', 'complex', 'bool', 'str', 'unicode']:
                warnings.warn("WARNING: " + str(ana_data['name']) + \
                              ": The data specified by the analysis object is not " + \
                              "in numpy format. Attempting to convert the data to numpy")
            try:
                dat = np.asarray(ana_data['data'])
                if len(dat.shape) == 0:
                    dat = np.asarray([ana_data['data']])
                tempdata = data_group.require_dataset(
                    name=ana_data['name'], shape=dat.shape, dtype=str(dat.dtype))
                if dat.size > 0:
                    tempdata[:] = dat
                else:
                    warnings.warn(ana_data['name'] + " dataset generated but not written. The given dataset was empty.")
            except:
                warnings.warn("ERROR: " + str(ana_data['name']) +
                              ": The data specified by the analysis could not be " +
                              "converted to numpy for writing to HDF5")

    def __init__(self,
                 analysis_group):
        """
        Initialize the analysis object given the h5py object of the analysis group.

        :param analysis_group: The h5py object with the analysis group of the omsi hdf5 file.

        """
        import warnings
        super(omsi_file_analysis, self).__init__(analysis_group)
        # The following tasks are performed by the super call
        # self.managed_group = analysis_group
        # self.dependencies = ...
        if self.dependencies is None:
            warnings.warn("No dependencies defined for analysis.")
        self.parameter = self.managed_group[unicode(omsi_format_analysis.analysis_parameter_group)]
        self.runinfo_group = self.managed_group[unicode(omsi_format_analysis.analysis_runinfo_group)]
        self.analysis_omsi_object = None
        self.name = self.managed_group.name

    def __getitem__(self,
                    key):
        """
        Support direct read interaction with the analysis h5py group
        """
        try:
            return self.managed_group[key]
        except:
            try:
                return self.parameter[key]
            except:
                try:
                    return self.dependencies[key]
                except:
                    pass
        # No errors have obscured but no matching object has been found either
        return None

    def __setitem__(self,
                    key,
                    value):
        """
        Support direct write interaction with the analysis h5py group
        """
        try:
            super(omsi_file_analysis, self).__setitem__(key, value)
        except:
            try:
                self.parameter[key] = value
            except:
                try:
                    if key in self.dependencies.items():
                        raise KeyError(
                            "Assignment to dependcies is not permitted via this mechanism")
                    else:
                        raise KeyError(
                            "Requested assignment operation failed")
                except:
                    raise KeyError("Requested assignment operation failed")

    def get_analysis_identifier(self):
        """
        Get the identifier name of the analysis.

        Use get_analysis_identifier()[...] to retrive the identifier string.

        :returns: h5py object for the dataset with the identifier string.
                 Returns None, in case no identifer exisits. This should
                 not be the case for a valid OpenMSI file.
        """
        if self.managed_group is None:
            return None

        try:
            return self.managed_group[unicode(omsi_format_analysis.analysis_identifier)]
        except:
            return None

    def get_analysis_type(self):
        """
        Get the type for the analysis.

        Use get_analysis_type()[...] tor retrieve the type string.

        :returns: h5py object with the dataset of the analysis string. Returns,
                 None in case no analysis type exists. This should not be the
                 case in a valid omsi file.
        """
        if self.managed_group is None:
            return None
        try:
            return self.managed_group[unicode(omsi_format_analysis.analysis_type)]
        except:
            return None

    def get_analysis_data_names(self):
        """
        This function returns all dataset names (and groups) that are custom
        to the analysis, i.e., that are not part of the omsi file standard.

        :returns: List of analysis-specific dataset names.
        """
        output_list = []
        if self.managed_group is not None:
            for item_obj in self.managed_group.items():
                if item_obj[0] != omsi_format_analysis.analysis_identifier and \
                        item_obj[0] != omsi_format_analysis.analysis_type and \
                        item_obj[0] != omsi_format_analysis.analysis_parameter_group and \
                        item_obj[0] != omsi_format_dependencies.dependencies_groupname and \
                        item_obj[0] != omsi_format_analysis.analysis_runinfo_group:
                    output_list.append(item_obj[0])
        return output_list

    def get_analysis_data_shapes_and_types(self):
        """
        This function returns two dictionaries with all dataset names (and groups)
        that are custom to the analysis, i.e., that are not part of the omsi
        file standard, and idenifies the shape of the analysis data objects.

        :returns: Dictonary indicating for each analysis-specific dataset its name (key) and
                 shape (value). And a second dictionariy indicating the name (key) and dtype
                 of the dataset.
        """
        output_shape_dict = {}
        output_type_dict = {}
        if self.managed_group is not None:
            for item_obj in self.managed_group.items():
                if item_obj[0] != omsi_format_analysis.analysis_identifier and \
                        item_obj[0] != omsi_format_analysis.analysis_type and \
                        item_obj[0] != omsi_format_analysis.analysis_parameter_group and \
                        item_obj[0] != omsi_format_dependencies.dependencies_groupname and \
                        item_obj[0] != omsi_format_analysis.analysis_runinfo_group:
                    try:
                        output_shape_dict[item_obj[0]] = self.managed_group[unicode(item_obj[0])].shape
                    except:
                        output_shape_dict[item_obj[0]] = (0,)
                    try:
                        output_type_dict[item_obj[0]] = self.managed_group[unicode(item_obj[0])].dtype
                    except:
                        output_type_dict[item_obj[0]] = "Unkown type"
        return output_shape_dict, output_type_dict

    def get_all_analysis_data(self,
                              load_data=False):
        """
        Get all analysis data associated with the analysis.

        :param load_data: load_data: Should the data be loaded or just the h5py objects be
                         stored in the dictionary.

        :returns: List of omsi_analysis_data objects with the names and h5py or numpy objects.
                 Access using [index]['name'] and [index]['data'].
        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data
        output_list = []
        if self.managed_group is not None:
            for item_obj in self.managed_group.items():
                if item_obj[0] != omsi_format_analysis.analysis_identifier and \
                        item_obj[0] != omsi_format_analysis.analysis_type and \
                        item_obj[0] != omsi_format_analysis.analysis_parameter_group and \
                        item_obj[0] != omsi_format_dependencies.dependencies_groupname and \
                        item_obj[0] != omsi_format_analysis.analysis_runinfo_group:
                    output_list.append(omsi_analysis_data())
                    output_list[-1]['name'] = str(item_obj[0])
                    if load_data:
                        output_list[-1]['data'] = self.managed_group[unicode(item_obj[0])][:]
                    else:
                        output_list[-1]['data'] = self.managed_group[unicode(item_obj[0])]
                    output_list[-1]['dtype'] = str(output_list[-1]['data'].dtype)
        return output_list

    def get_all_parameter_data(self,
                               load_data=False,
                               exclude_dependencies=False):
        """
        Get all parameter data associated with the analysis.

        :param load_data: Should the data be loaded or just the h5py objects be stored in the dictionary.

        :returns: List of omsi_parameter_data objects with names and h5py or numpy object. Access using
                 [index]['name'] and [index]['data'].
        """
        from omsi.analysis.omsi_analysis_data import omsi_parameter_data
        from omsi.shared.omsi_dependency import omsi_dependency
        output_list = []
        if self.parameter is not None:
            for item_obj in self.parameter.items():
                curr_parameter = omsi_parameter_data(name=unicode(item_obj[0]),
                                                     help='')
                curr_parameter_dataset = self.parameter[unicode(item_obj[0])]
                if omsi_format_analysis.analysis_parameter_help_attr in curr_parameter_dataset.attrs:
                    curr_parameter['help'] = unicode(
                        curr_parameter_dataset.attrs[omsi_format_analysis.analysis_parameter_help_attr])
                if load_data:
                    curr_parameter['data'] = curr_parameter_dataset[:]
                else:
                    curr_parameter['data'] = curr_parameter_dataset
                curr_parameter['dtype'] = unicode(curr_parameter['data'].dtype)
                output_list.append(curr_parameter)
        if not exclude_dependencies:
            dependency_data = self.get_all_dependency_data(omsi_dependency_format=False)
            for dep in dependency_data:
                curr_parameter = omsi_parameter_data(name=dep.get_parameter_name(),
                                                     help='')
                curr_parameter['data'] = dep.get_omsi_dependency()
                curr_parameter['dtype'] = omsi_dependency
                output_list.append(curr_parameter)

        return output_list

    def get_all_runinfo_data(self,
                             load_data=False):
        """
        Get a dict of all runtime information stored in the file

        :return: Dict with the runtime information restored.
        """
        output_dict = {}
        if self.runinfo_group is not None:
            for item_obj in self.runinfo_group.items():
                object_name = unicode(item_obj[0])
                if load_data:
                    output_dict[object_name] = self.runinfo_group[object_name][:]
                    if output_dict[object_name].shape == (1,):
                        curr_dtype = self.runinfo_group[object_name].dtype
                        if curr_dtype == h5py.special_dtype(vlen=unicode) or \
                                curr_dtype == h5py.special_dtype(vlen=unicode):
                            output_dict[object_name] = unicode(output_dict[object_name][0])
                        else:
                            output_dict[object_name] = output_dict[object_name][0]
                else:
                    output_dict[object_name] = self.runinfo_group[object_name]
        return output_dict

    def restore_analysis(self,
                         load_data=True,
                         load_parameters=True,
                         load_runtime_data=True,
                         dependencies_omsi_format=True):
        """
        Load an analysis from file and create an instance of the appropriate
        analysis object defined by the analysis type (i.e., a derived class
        of `omsi.analysis. omsi_analysis_base`)

        :param load_data: Should the analysis data be loaded from file (default) or just stored as h5py data objects
        :param load_parameters: Should parameters be loaded from file (default) or just stored as h5py data objects.
        :param load_runtime_data: Should runtime data be loaded from file (default) or just stored as h5py data objects.
        :param dependencies_omsi_format: Should dependencies be loaded as omsi_file API objects (default)
            or just as h5py objects.

        :return: Instance of the specific analysis object (e.g, omsi_nmf) that inherits from
                 omsi.analysis.omsi_analysis_base with the input parameters, output result,
                 and dependencies restored. We can call execute(..) on the returned object
                 to rerun an analysis. May return omsi_analysis_generic in case that the
                 specific analysis is not known.
        """
        from omsi.analysis.omsi_viewer_helper import omsi_viewer_helper
        from omsi.analysis.omsi_analysis_generic import omsi_analysis_generic
        try:
            analysis_class = omsi_viewer_helper.analysis_name_to_class(self.get_analysis_type()[0])
            ignore_type_conflict = False
        except:
            analysis_class = omsi_analysis_generic
            ignore_type_conflict = True
        analysis_instance = analysis_class()
        analysis_instance.read_from_omsi_file(analysis_object=self,
                                              load_data=load_data,
                                              load_parameters=load_parameters,
                                              load_runtime_data=load_runtime_data,
                                              dependencies_omsi_format=dependencies_omsi_format,
                                              ignore_type_conflict=ignore_type_conflict)
        return analysis_instance

    def recreate_analysis(self, **kwargs):
        """
        Load an analysis from file and re-execute it.
        This is equivalent to omsi_analysis.base.restore_analysis().execute()

        :param kwargs: Additional keyword arguments to be passed to the execute function of the analysis

        :return: Instance of the specific analysis object (e.g, omsi_nmf) that inherits from
                 omsi.analysis.omsi_analysis_base with the input parameters and
                 dependencies restored from file. The output, however, is the result from
                 re-executing the analysis. None is returned in case the analysis object
                 cannot be created.
        """
        from omsi.analysis import omsi_analysis_generic
        analysis_instance = self.restore_analysis(**kwargs)
        if isinstance(analysis_instance, omsi_analysis_generic):
            return None
        else:
            analysis_instance.execute()
            return analysis_instance

