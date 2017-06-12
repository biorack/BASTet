"""
Module for managing custom analysis data in OMSI HDF5 files.
"""
import os
import warnings

import numpy as np
import h5py

from omsi.dataformat.omsi_file.format import \
    omsi_format_common, \
    omsi_format_analysis, \
    omsi_format_dependencies
from omsi.dataformat.omsi_file.dependencies import omsi_dependencies_manager
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager
from omsi.datastructures.run_info_data import run_info_dict
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.log import log_helper

#try:
#    import cloudpickle   # Use the version of cloud-pickle installed on the system
#    log_helper.debug(__name__, "Using system cloudpickle module")
#except ImportError:
#    try:
import omsi.shared.third_party.cloudpickle as cloudpickle
log_helper.debug(__name__, "Using fallback cloudpickle version")
#    except ImportError:
#        log_helper.warning(__name__, "cloudpickle could not be imported. Using standard pickle instead. " +
#                           " Some features may not be available.")
#        import pickle as cloudpickle
import pickle

# TODO create_analysis_static(...) and other create functions need to handle the case when a file is opened with the MPI I/O backend. Currently we assume a serial write from root

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
        try:
            # In most cases all members in the inheritance hierarchy are expected to accept the parent group as input
            super(omsi_analysis_manager, self).__init__(analysis_parent)
        except TypeError:
            # We are initialized separate from an omis_file_common hierarchy and our super is just object
            super(omsi_analysis_manager, self).__init__()
        self.analysis_parent = analysis_parent

    @staticmethod
    def create_analysis_static(analysis_parent,
                               analysis,
                               flush_io=True,
                               force_save=False,
                               save_unsaved_dependencies=True,
                               mpi_root=0,
                               mpi_comm=None):
        """
        Same as create_analysis(...) but instead of relying on object-level, this function
        allows additional parameters (specifically the analysis_parent) to be provided as
        input, rather than being determined based on self

        :param analysis_parent: The h5py.Group object or omsi.dataformat.omsi_file.common.omsi_file_common object
            where the analysis should be created
        :param kwargs: Additional keyword arguments for create_analysis(...). See create_analysis(...) for details.

        :return: The output of create_analysis
        """
        if mpi_helper.get_rank(comm=mpi_comm) == mpi_root:
            if isinstance(analysis_parent, h5py.Group):
                parent_group = analysis_parent
            elif isinstance(analysis_parent, omsi_file_common):
                parent_group = analysis_parent.managed_group
            else:
                log_helper.error(__name__, 'Illegal analysis_parent type. Expected h5py.Group or omsi_file_common')
                raise ValueError("Illegal value for analysis parent")
            return omsi_file_analysis.__create__(parent_group=parent_group,
                                                 analysis=analysis,
                                                 analysis_index=None,  # Same as self.get_num_analysis()
                                                 flush_io=flush_io,
                                                 force_save=force_save,
                                                 save_unsaved_dependencies=save_unsaved_dependencies)
        else:
            try:
                analysis.write_analysis_data()
                return None
            except NotImplementedError:
                pass


    def create_analysis(self,
                        analysis,
                        flush_io=True,
                        force_save=False,
                        save_unsaved_dependencies=True,
                        mpi_root=0,
                        mpi_comm=None):
        """
        Add a new group for storing derived analysis results for the current experiment

        Create the analysis group using omsi_file_analysis.__create___ which in turn uses
        omsi_file_analysis.__populate_analysis__(...) to populate the group with the appropriate data.

        NOTE: Dependencies are generally resolved to point to file objects. However, if save_unsaved_dependencies
        is set to False and a given in-memory dependency has not been saved yet, then the value associated with
        that dependency will be saved instead as part of the parameters and, hence, only the value of the
        dependency is persevered in that case and not the full dependency chain.

        NOTE: Dependencies if they only exists in memory are typically saved recursively unless
        save_unsaved_dependencies is set to False. I.e, calling create_analysis may result in the
        creating of multiple other dependent analyses if they have not been saved before.


        :param analysis: Instance of omsi.analysis.analysis_base defining the analysis
        :type analysis: omsi.analysis.analysis_base:

        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                       so that all data has been written to file
        :type flush_io: bool

        :param force_save: Should we save the analysis even if it has been saved in the same location before?
            If force_save is False (default) and the self.omsi_analysis_storage parameter of the analysis object
            contains a matching storage location---i.e., same file and experiment---, then the analysis will
            not be saved again, but the object will only be retrieved from file.
            If force_save is True, then the analysis will be saved either way and the self.omsi_analysis_storage
            parameter will be extended.
        :type force_save: bool

        :param save_unsaved_dependencies: If there are unsaved (in-memory) dependencies, then should those be
            saved to file as well? Default value is True, i.e, by default all in-memory dependencies that have
            not been saved yet, i.e, for which the self.omsi_analysis_storage of the corresponding omsi_analysis_
            base object is empty, are saved as well. If in-memory dependencies have been saved before, then a
            link to those dependencies will be established, rather than re-saving the dependency.
        :type save_unsaved_dependencies: bool

        :param mpi_root: The root MPI process that should perform the writing. This is to allow all analyses
            to call the function and have communication in the analysis.write_analysis_data function be
            handled.

        :param mpi_comm: The MPI communicator to be used. None if default should be used (ie., MPI.COMM_WORLD

        :returns: The omsi_file_analysis object for the newly created analysis group and the
                integer index of the analysis. NOTE: if force_save is False (default), then the group
                returned may not be new but may be simply the first entry in the list of existing
                storage locations for the given analysis. NOTE: If we are in MPI parallel and we are on a
                core that does not write any data, then None is returned instead.
        """
        return self.create_analysis_static(analysis_parent=self.analysis_parent,
                                           analysis=analysis,
                                           flush_io=flush_io,
                                           force_save=force_save,
                                           save_unsaved_dependencies=save_unsaved_dependencies,
                                           mpi_root=mpi_root,
                                           mpi_comm=mpi_comm)

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

        :returns: omsi_file_analysis object for the requested analysis. The function
                      returns None in case the analysis object was not found.
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
                   flush_io=True,
                   force_save=False,
                   save_unsaved_dependencies=True):
        """
        Add a new group for storing derived analysis results for the given parent object

        Create the analysis group and use omsi_file_analysis.__populate_analysis__(...) to populate the
        group with the appropriate data.

        :param parent_group: The h5py.Group where the method object should be created
        :type parent_group: h5py.Group

        :param analysis_index: The integer index of the analysis to be created. Set to None (default) if the
                         function should determine the index automatically.
        :type analysis_index: int or None

        :param analysis: Instance of omsi.analysis.analysis_base defining the analysis
        :type analysis: omsi.analysis.analysis_base:

        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                       so that all data has been written to file
        :type flush_io: True

        :param force_save: Should we save the analysis even if it has been saved in the same location before?
            If force_save is False (default) and the self.omsi_analysis_storage parameter of the analysis object
            contains a matching storage location---i.e., same file and experiment---, then the analysis will
            not be saved again, but the object will only be retrieved from file.
            If force_save is True, then the analysis will be saved either way and the self.omsi_analysis_storage
            parameter will be extended.
        :type force_save: bool

        :param save_unsaved_dependencies: If there are unsaved (in-memory) dependencies, then should those be
            saved to file as well? Default value is True, i.e, by default all in-memory dependencies that have
            not been saved yet, i.e, for which the self.omsi_analysis_storage of the corresponding omsi_analysis_
            base object is empty, are saved as well. If in-memory dependencies have been saved before, then a
            link to those dependencies will be established, rather than re-saving the dependency.
        :type save_unsaved_dependencies: bool

        :returns: i) the omsi_file_analysis object for the newly created analysis group ii) the integer index
                  for the analysis.

        :raises: IndexError is raised in case that an analysis with the given index already exists.
        """
        # 1 Import required modules
        from omsi.analysis.base import analysis_base
        import time

        # 2 Check if input values are correct
        # 2.1 Confirm that we have a valid analysis object
        if not isinstance(analysis, analysis_base):
            errormessage = "Could not write the analysis. The input object was of type " + \
                str(type(analysis)) + " not of analysis_base as expected."
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

        # 3. Check if we need to save the analysis or whether it has already been saved
        if not(force_save or len(analysis.get_omsi_analysis_storage()) == 0):
             # Check if the analysis has been saved before for the same parent_group
            parent_filename = os.path.abspath(parent_group.file.filename)
            for prior_data_store in analysis.get_omsi_analysis_storage():
                prior_group_name = os.path.split(prior_data_store.name)[0]
                prior_filename = os.path.abspath(prior_data_store.managed_group.file.filename)
                if prior_group_name == parent_group.name and prior_filename == parent_filename:
                    return prior_data_store, prior_data_store.get_analysis_index()

        # 4. Create the group and define the required attributes
        analysis_group = parent_group.require_group(analysis_group_name)
        analysis_group.attrs[omsi_format_common.type_attribute] = "omsi_file_analysis"
        analysis_group.attrs[omsi_format_common.version_attribute] = omsi_format_analysis.current_version
        analysis_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())

        # 5. Check if there are any unsaved in-memory dependencies. NOTE: All dependencies are assumed to have
        if save_unsaved_dependencies:
            save_extra_dependencies = []
            for dependency in analysis.get_all_dependency_data():
                if isinstance(dependency['data']['omsi_object'], analysis_base):
                    if not dependency['data']['omsi_object'].has_omsi_analysis_storage():
                        save_extra_dependencies.append(dependency)
            if len(save_extra_dependencies) > 0:  # Same as 'if save_extra_dependencies:'
                temp_analysis_manager = omsi_analysis_manager(analysis_parent=parent_group)
                for dependency in save_extra_dependencies:
                    temp_analysis_manager.create_analysis(analysis=dependency['data']['omsi_object'],
                                                          flush_io=flush_io,
                                                          force_save=False,
                                                          save_unsaved_dependencies=save_unsaved_dependencies)

        # 6. Populate the group with data and return the omsi_file_analysis object
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

        NOTE: At this point we assume that all in-memory dependencies have been resolved. If not,
        then the raw data associated with the given parameter will be saved instead.

        :param analysis_group: h5py group in which the analysis data should be stored.
        :param analysis: Instance of omsi.analysis.analysis_base defining the analysis
        :type analysis: omsi.analysis.analysis_base:

        :returns: The omsi_file_analysis object for the newly created analysis group. The analysis data is
                  automatically written to file by this function so no addition work is required.

        """
        from omsi.datastructures.analysis_data import analysis_data
        from omsi.dataformat.omsi_file.dependencies import omsi_file_dependencies
        from omsi.analysis.base import analysis_base

        # 1. Write the analysis name
        analysis_identifier_data = analysis_group.require_dataset(
            name=unicode(omsi_format_analysis.analysis_identifier),
            shape=(1,),
            dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysis_identifier_data[0] = analysis.get_analysis_identifier()
        else:
            analysis_identifier_data[0] = str(analysis.get_analysis_identifier())

        # 2. Write the analysis type
        analysis_type_data = analysis_group.require_dataset(name=unicode(omsi_format_analysis.analysis_type),
                                                            shape=(1,),
                                                            dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysis_type_data[0] = analysis.get_analysis_type()
        else:
            analysis_type_data[0] = str(analysis.get_analysis_type())

        # 3. Write the analysis data
        try:
            analysis.write_analysis_data(analysis_group=analysis_group)
        except NotImplementedError:
            for ana_data in analysis.get_all_analysis_data():
                cls.__write_omsi_analysis_data__(analysis_group, ana_data)

        # 4. Determine all dependencies and parameters that we need to write
        dependencies = []  # [dep['data'] for dep in analysis.get_all_dependency_data()]
        parameters = []
        # 4.1 Resolve in-memory dependencies if possible
        for dependent_parameter in analysis.get_all_dependency_data():
            # 4.1.1 We have an in-memory dependency
            if isinstance(dependent_parameter['data']['omsi_object'], analysis_base):
                # 4.1.1.1 We can resolve the dependency to an object in an HDF5 file
                if dependent_parameter['data']['omsi_object'].has_omsi_analysis_storage():
                    # Create a new dependency that points to the approbriate file location
                    # NOTE: We do not modify the dependency in the analysis object that we save
                    #       but we only change it for the purpose of storage
                    new_dep = dependent_parameter['data'].copy()
                    new_dep_omsi_object = None
                    # Check if we can find an analysis data store within the same parent (or at least file)
                    parent_filename = os.path.abspath(analysis_group.file.filename)
                    for analysis_store in dependent_parameter['data']['omsi_object'].get_omsi_analysis_storage():
                        analysis_store_filename = os.path.abspath(analysis_store.managed_group.file.filename)
                        if analysis_store.name == analysis_group.parent.name and \
                                analysis_store_filename == parent_filename:
                            new_dep_omsi_object = analysis_store
                            break
                        elif analysis_store_filename == parent_filename:
                            new_dep_omsi_object = analysis_store

                    # We could not find a prior data store within the same file so use one from another file
                    if new_dep_omsi_object is None:
                        dep_object = dependent_parameter['data']['omsi_object']
                        new_dep['omsi_object'] = dep_object.get_omsi_analysis_storage()[0]
                    else:
                        new_dep['omsi_object'] = new_dep_omsi_object
                    # Append it to the list of dependencies
                    dependencies.append(new_dep)
                # 4.1.1.2  We cannot resolve the dependency and need to store it as an parameter instead
                else:
                    # Replace the dependency with the actual data and save it as a parameter instead
                    new_param = dependent_parameter.copy()
                    new_param['data'] = new_param['data'].get_data()
                    parameters.append(new_param)

            # 4.1.2 We have a file-based dependencies so keep it as is and add it to the list of dependencies
            else:
                dependencies.append(dependent_parameter['data'])

        # 4.2 Add all regular parameters to the list of parameters
        parameters += analysis.get_all_parameter_data(exclude_dependencies=True)

        # 5. Write all the parameters
        parameter_group = analysis_group.require_group(omsi_format_analysis.analysis_parameter_group)
        for param_data in parameters:
            if param_data['required'] or param_data.data_set() or param_data['default'] is not None:
                temp_data = param_data.get_data_or_default()
                if temp_data is not None:
                    anadata = analysis_data(name=param_data['name'],
                                            data=param_data.get_data_or_default(),
                                            dtype=param_data['dtype'])
                    cls.__write_omsi_analysis_data__(parameter_group, anadata)
                    # Try to add the help string attribute
                    try:
                        help_attr = omsi_format_analysis.analysis_parameter_help_attr
                        parameter_group[param_data['name']].attrs[help_attr] = param_data['help']
                    except KeyError:
                        pass

        # 6. Write all the runtime execution information
        runinfo_group = analysis_group.require_group(omsi_format_analysis.analysis_runinfo_group)
        for run_info_key, run_info_value in analysis.get_all_run_info().items():
            # Generate an analysis_data object in order to use the
            # __write_omsi_analysis_data function to write the data
            if isinstance(run_info_value, unicode) or isinstance(run_info_value, str):
                anadata = analysis_data(name=unicode(run_info_key),
                                        data=run_info_value,
                                        dtype=omsi_format_common.str_type)
            else:
                dat = np.asarray(run_info_value)
                if len(dat.shape) == 0:
                    dat = np.asarray([run_info_value])
                anadata = analysis_data(name=unicode(run_info_key),
                                        data=dat,
                                        dtype=dat.dtype)
            cls.__write_omsi_analysis_data__(runinfo_group, anadata)

        # 7. Write all dependencies
        omsi_file_dependencies.__create__(parent_group=analysis_group,
                                          dependencies_data_list=dependencies)

        # 8. Execute the custom data write for the analysis
        analysis.add_custom_data_to_omsi_file(analysis_group)

        # 9. Create the output object
        re = omsi_file_analysis(analysis_group)

        # 10. Save the output object in the ist of omsi analysis data stores as part of the analysis object
        analysis.omsi_analysis_storage.append(re)

        # 11. Check if we need to pickle and save the analysis class in case this is a custom class that is not part of BASTet
        try:
            from omsi.analysis.analysis_views import analysis_views
            _ = analysis_views.analysis_name_to_class(analysis.get_analysis_type())
        except NameError:
            class_pickle = cloudpickle.dumps(analysis.__class__)
            # Convert the pickle string to an uint8 array to avoid problems
            # with storing string with NULL characters in HDF5
            class_pickle_arr = np.fromstring(class_pickle,
                                            dtype=omsi_format_analysis.analysis_class_pickle_np_dtype)
            analysis_group[unicode(omsi_format_analysis.analysis_class)] = class_pickle_arr
        except:
            log_helper.warning(__name__, "Could not save the analysis class.")
            pass

        # 12. Retrun the new omsi_file_analysis object
        return re

    @classmethod
    def __write_omsi_analysis_data__(cls,
                                     data_group,
                                     ana_data):
        """
        Private helper function used to write the data defined by a analysis_data object to HDF5.

        :param data_group: The h5py data group to which the data should be written to.
        :param ana_data: The analysis_data object with the description of the data to be written.
        :type ana_data: omsi.analysis.analysis_data
        """
        from omsi.datastructures.analysis_data import analysis_data, data_dtypes
        curr_dtype = ana_data['dtype']
        try:
            if curr_dtype == data_dtypes.get_dtypes()['ndarray']:
                curr_dtype = ana_data['data'].dtype
        except TypeError:
            pass
        try:
            if curr_dtype == data_dtypes.get_dtypes()['bool']:
                curr_dtype = bool
        except TypeError:
            pass
        try:
            if curr_dtype == data_dtypes.get_dtypes()['str']:
                curr_dtype =  omsi_format_common.str_type
        except TypeError:
            pass

        # Create link in HDF5 to an existing dataset within the file
        if isinstance(ana_data, analysis_data) and isinstance(ana_data['dtype'], int):
            if curr_dtype == ana_data.ana_hdf5link:
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
        elif (curr_dtype == omsi_format_common.str_type) or (curr_dtype == h5py.special_dtype(vlen=str)):
            tempdata = data_group.require_dataset(name=unicode(ana_data['name']),
                                                  shape=(1,),
                                                  dtype=omsi_format_common.str_type)
            if len(unicode(ana_data['data'])) > 0:
                if omsi_format_common.str_type_unicode:
                    tempdata[0] = unicode(ana_data['data'])
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
            if ana_data['data'].dtype.type in [np.string_, np.unicode_]:
                tempdata = data_group.require_dataset(name=ana_data['name'],
                                                      shape=ana_data['data'].shape,
                                                      dtype=omsi_format_common.str_type,
                                                      chunks=chunks)
            else:
                tempdata = data_group.require_dataset(name=ana_data['name'],
                                                      shape=ana_data['data'].shape,
                                                      dtype=ana_data['data'].dtype,
                                                      chunks=chunks)
            if ana_data['data'].size > 0:
                try:
                    tempdata[:] = ana_data['data']
                except TypeError:
                    tempdata[()] = ana_data['data']
            else:
                warnings.warn("WARNING: " + ana_data['name'] +
                              " dataset generated but not written. The given dataset was empty.")
        # Unknown dtype. Attempt to convert the dataset to numpy and write it to
        # file.
        else:
            # Safely convert scalars to numpy but warn in case we see something else
            from omsi.datastructures.analysis_data import data_dtypes
            default_dtypes = data_dtypes.get_dtypes()
            if ana_data['dtype'] not in default_dtypes.keys() and ana_data['dtype'] not in default_dtypes.values():
                warnings.warn("WARNING: " + str(ana_data['name']) +
                              ": The data specified by the analysis object is not " +
                              "in numpy format. Attempting to convert the data to numpy")
            try:
                dat = np.asarray(ana_data['data'])
                if len(dat.shape) == 0:
                    dat = dat[np.newaxis]  # np.asarray([ana_data['data']])
                try:
                    tempdata = data_group.require_dataset(name=ana_data['name'],
                                                          shape=dat.shape,
                                                          dtype=str(dat.dtype))
                except TypeError:  # Some Unicode types are not well-understood by h5py
                    if 'U' in str(dat.dtype) or 'S' in str(dat.dtype):
                        tempdata = data_group.require_dataset(name=ana_data['name'],
                                                              shape=dat.shape,
                                                              dtype=omsi_format_common.str_type)
                    else:
                        raise
                if dat.size > 0:
                    try:
                        tempdata[:] = dat
                    except TypeError:
                        tempdata[()] = dat
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
        try:
            self.runinfo_group = self.managed_group[unicode(omsi_format_analysis.analysis_runinfo_group)]
        except:
            self.runinfo_group = None
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

    def get_analysis_index(self):
        """
        Based on the name of the group, get the index of the analysis.

        :return: Integer index of the analysis in the file.
        """
        return int(self.managed_group.name.split('_')[-1])

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

    def get_analysis_class(self):
        """
        Get the analysis class corresponding to the save analysis. This may be
        analysis_generic in case that the analysis type is not known or generic
        :return:
        """
        from omsi.analysis.analysis_views import analysis_views
        from omsi.analysis.generic import analysis_generic
        try:
            # Get the analysis from BASTet itself
            analysis_class = analysis_views.analysis_name_to_class(self.get_analysis_type()[0])
        except NameError:
            # Try to restore the analysis pickle if we saved it to file
            if omsi_format_analysis.analysis_class in self.managed_group.keys():
                class_pickle = self.managed_group[unicode(omsi_format_analysis.analysis_class)][:]
                # Restore the pickle string we stored as uint8 array to avoid problems with NULL
                class_pickle_string = class_pickle.tostring()
                analysis_class = pickle.loads(class_pickle_string)
            else:
                analysis_class = analysis_generic
        except:
            analysis_class = analysis_generic
        return analysis_class

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

        :returns: List of analysis_data objects with the names and h5py or numpy objects.
                 Access using [index]['name'] and [index]['data'].
        """
        from omsi.datastructures.analysis_data import analysis_data
        output_list = []
        if self.managed_group is not None:
            for item_obj in self.managed_group.items():
                if item_obj[0] != omsi_format_analysis.analysis_identifier and \
                        item_obj[0] != omsi_format_analysis.analysis_type and \
                        item_obj[0] != omsi_format_analysis.analysis_parameter_group and \
                        item_obj[0] != omsi_format_dependencies.dependencies_groupname and \
                        item_obj[0] != omsi_format_analysis.analysis_runinfo_group:
                    output_list.append(analysis_data())
                    output_list[-1]['name'] = str(item_obj[0])
                    if load_data:
                        if len(self.managed_group[unicode(item_obj[0])].shape) == 0:  # Scalar dataset
                            output_list[-1]['data'] = self.managed_group[unicode(item_obj[0])][()]
                        else:
                            output_list[-1]['data'] = self.managed_group[unicode(item_obj[0])][:]
                    else:
                        output_list[-1]['data'] = self.managed_group[unicode(item_obj[0])]
                    if not isinstance(output_list[-1]['data'], basestring):
                        output_list[-1]['dtype'] = str(output_list[-1]['data'].dtype)
                    else:
                        output_list[-1]['dtype'] = unicode
        return output_list

    def get_all_parameter_data(self,
                               load_data=False,
                               exclude_dependencies=False):
        """
        Get all parameter data associated with the analysis.

        :param load_data: Should the data be loaded or just the h5py objects be stored in the dictionary.

        :returns: List of parameter_data objects with names and h5py or numpy object. Access using
                 [index]['name'] and [index]['data'].
        """
        from omsi.datastructures.analysis_data import parameter_data
        from omsi.datastructures.dependency_data import dependency_dict
        output_list = []
        if self.parameter is not None:
            for item_obj in self.parameter.items():
                curr_parameter = parameter_data(name=unicode(item_obj[0]),
                                                     help='')
                curr_parameter_dataset = self.parameter[unicode(item_obj[0])]
                if omsi_format_analysis.analysis_parameter_help_attr in curr_parameter_dataset.attrs:
                    curr_parameter['help'] = unicode(
                        curr_parameter_dataset.attrs[omsi_format_analysis.analysis_parameter_help_attr])
                if load_data:
                    curr_parameter['data'] = curr_parameter_dataset[:]
                    # try to restore string parameters as stings
                    try:
                        if curr_parameter['data'].size == 1:
                            if isinstance(curr_parameter['data'][0], basestring):
                                curr_parameter['data'] = curr_parameter['data'][0]
                    except:
                        raise
                else:
                    curr_parameter['data'] = curr_parameter_dataset
                if not isinstance(curr_parameter['data'], basestring):
                    curr_parameter['dtype'] = unicode(curr_parameter['data'].dtype)
                else:
                    curr_parameter['dtype'] = unicode
                output_list.append(curr_parameter)
        if not exclude_dependencies:
            dependency_data = self.get_all_dependency_data(omsi_dependency_format=False)
            for dep in dependency_data:
                curr_parameter = parameter_data(name=dep.get_parameter_name(),
                                                     help='')
                curr_parameter['data'] = dep.get_omsi_dependency()
                curr_parameter['dtype'] = dependency_dict
                output_list.append(curr_parameter)

        return output_list

    def get_all_runinfo_data(self,
                             load_data=False):
        """
        Get a dict of all runtime information stored in the file

        :return: omsi.shared.run_info_data.run_info_dict type python dict with the runtime information restored.
        """
        output_dict = run_info_dict()
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
        of `omsi.analysis. analysis_base`)

        :param load_data: Should the analysis data be loaded from file (default) or just stored as h5py data objects
        :param load_parameters: Should parameters be loaded from file (default) or just stored as h5py data objects.
        :param load_runtime_data: Should runtime data be loaded from file (default) or just stored as h5py data objects.
        :param dependencies_omsi_format: Should dependencies be loaded as omsi_file API objects (default)
            or just as h5py objects.

        :return: Instance of the specific analysis object (e.g, omsi_nmf) that inherits from
                 omsi.analysis.analysis_base with the input parameters, output result,
                 and dependencies restored. We can call execute(..) on the returned object
                 to rerun an analysis. May return analysis_generic in case that the
                 specific analysis is not known.
        """
        from omsi.analysis.generic import analysis_generic
        analysis_class = self.get_analysis_class()
        analysis_instance = analysis_class()
        analysis_instance.read_from_omsi_file(analysis_object=self,
                                              load_data=load_data,
                                              load_parameters=load_parameters,
                                              load_runtime_data=load_runtime_data,
                                              dependencies_omsi_format=dependencies_omsi_format,
                                              ignore_type_conflict=isinstance(analysis_instance, analysis_generic))
        return analysis_instance

    def recreate_analysis(self, **kwargs):
        """
        Load an analysis from file and re-execute it.
        This is equivalent to omsi_analysis.base.restore_analysis().execute()

        :param kwargs: Additional keyword arguments to be passed to the execute function of the analysis

        :return: Instance of the specific analysis object (e.g, omsi_nmf) that inherits from
                 omsi.analysis.analysis_base with the input parameters and
                 dependencies restored from file. The output, however, is the result from
                 re-executing the analysis. None is returned in case the analysis object
                 cannot be created.
        """
        #from omsi.analysis import analysis_generic
        analysis_instance = self.restore_analysis(load_runtime_data=False)  # We don't need the runtime data
        #if isinstance(analysis_instance, analysis_generic):
        #    return None
        #else:
        analysis_instance.execute(**kwargs)
        return analysis_instance

