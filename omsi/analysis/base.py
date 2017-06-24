"""
Module specifying the base analysis API for integrating new analysis with the toolkit and the
OpenMSI science gateway.
"""

# TODO Rename self.driver to self.executor
# TODO Add get_all_dependency_data_recursive back in (see all omsi.datastructures.analysis_data)
# TODO Separate qmz, qslice, qspectrum viewer functionality into a separate base class
# TODO Make list of viewer options unique (if we have a dependency graph rather than a tree) than the same option can occure multiple times
# TODO Add ability to lable inidivdual outputs as ready or not to support per-output-based interactivity

import warnings
import weakref
from collections import OrderedDict

import numpy as np

from omsi.workflow.executor.base import workflow_executor_base
from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
from omsi.dataformat.omsi_file.msidata import omsi_file_msidata
from omsi.datastructures.analysis_data import analysis_data, parameter_data, data_dtypes
from omsi.datastructures.dependency_data import dependency_dict
from omsi.datastructures.analysis_data import parameter_manager
import omsi.shared.mpi_helper as mpi_helper
from omsi.datastructures.run_info_data import run_info_dict
from omsi.shared.log import log_helper


class AnalysisReadyError(Exception):
    """
    Custom exception used to indicate that an analysis is not ready to execute.
    """
    def __init__(self, value, params=None):
        """
        Initialize the AnalysisReadyError

        :param value: Error message string
        :param params: Optional list of dependent parameters that are not ready to be used.
        """
        self.value = value
        self.params = params
        message = unicode(self.value)
        if self.params is not None and len(self.params) > 0:
            message += " The following parameters are not ready: "
            for param in self.params:
                try:
                    message += " " + param['name']
                except KeyError:
                    pass
        super(AnalysisReadyError, self).__init__(repr(message))


class analysis_base(parameter_manager):
    """
    Base class for omsi analysis functionality. The class provides a large set of functionality designed
    to facilitate storage of analysis data in the omsi HDF5 file format. The class also provides a set
    of functions to enable easy intergration of new analysis with the OpenMSI web-based viewer (see
    Viewer functions below for details).

    **Slicing:**

    This class supports basic slicing to access data stored in the main member variables. By
    default the data is retrieved from __data_list and the __getitem__(key) function. which implements
    the [..] operator, returns __data_list[key]['data']. The key is a string indicating the name of
    the parameter to be retrieved. If the key is not found in the __data_list then the function will
    try to retrieve the data from self.parameters list instead. By adding "parameter/key" or "dependency/key"
    one may also explicitly retrieve values from the parameters.

    **Instance Variables:**

    :ivar analysis_identifier: Define the name for the analysis used as key in search operations
    :ivar __data_list: List of analysis_data to be written to the HDF5 file. Derived classes
        need to add all data that should be saved for the analysis in the omsi HDF5 file to this dictionary.
        See omsi.analysis.analysis_data for details.
    :ivar parameters: List of parameter_data objects of all  analysis parameters
         (including those that may have dependencies).
    :ivar data_names: List of strings of all names of analysis output datasets. These are the
         target keys for __data_list.
    :ivar profile_time_and_usage: Boolean indicating whether we should profile the execute_analysis(...) function
        when called as part of the execute(...) function. The default value is false. Use the
        enable_time_and_usage_profiling(..) function to determine which profiling should be performed. The time_and
        _usage profile uses pythons cProfile (or Profile) to monitor how often and for how long particular parts
        of the analysis code executed.
    :ivar profile_memory: Boolean indicating whether we should monitor memory usage (line-by-line) when
        executing the execute_analysis(...) function. The default value is false. Use the
        enable_time_and_usage_profiling(..) function to determine which profiling should be performed.
    :ivar omsi_analysis_storage: List of omsi_file_analysis object where the analysis is stored. The list may be empty.
    :ivar mpi_comm: In case we are running with MPI, this is the MPI communicator used for runnign the analysis.
        Default is MPI.Comm_world/
    :ivar mpi_root: In case we are running with MPI, this is the root rank where data is collected to (e.g., runtime
        data and analysis results)
    :ivar update_analysis: If the value is True, then we should execute the analysis before using the outputs.
        If False, then the analysis has been executed with the current parameter settings.
    :ivar driver: Workflow driver to be used when executing multiple analyses, e.g., via execute_recursive or
        execute_all. Default value is None in which case a new default driver will be used each time we
        execute a workflow.


    **Execution Functions:**

    * ``execute`` : Then main function the user needs to call in order to execute the analysis
    * ``execute_analysis: This function needs to be implemented by child classes of `analysis_base` \
        to implement the specifics of executing the analysis.

    **I/O functions:**

    These functions can be optionally overwritten to control how the analysis data should be written/read
    from the omsi HDF5 file. Default implementations are provided here, which should be sufficient for most cases.

    * ``add_custom_data_to_omsi_file``: The default implementation is empty as the default data write is  managed by \
    the `omsi_file_experiment.create_analysis()` function.  Overwrite this function, in case that the analysis needs \
    to write data to the HDF5 omsi file beyond what the defualt omsi data API does.

    * ``read_from_omsi_file``: The default implementation tries to reconstruct the original data as far as possible, \
    however, in particular in case that a custom add_custom_data_to_omsi_file function has been implemented, the \
    default implementation may not be sufficien. The default implementation reconstructs: i) analysis_identifier \
    and reads all custom data into ii)__data_list. Note, an error will be raised in case that the analysis type \
    specified in the HDF5 file does not match the analysis type specified by get_analysis_type(). This function \
    can be optionally overwritten to implement a custom data read.

    **Viewer functions:**

    Several convenient functions are used to allow the OpenMSI online viewer to interact with the analysis \
    and to visualize it. The default implementations provided here simply indicate that the analysis does not \
    support the data access operations required by the online viewer. Overwrite these functions in the derived \
    analysis classes in order to interface them with the viewer. All viewer-related functions start with ``v\_...`` .

    NOTE: the default implementation of the viewer functions defined in ``analysis_base`` are \
    designed to take care of the common requirement for providing viewer access to data from all dependencies \
    of an analysis. In many cases, the default implementation is often sill called at the end of custom \
    viewer functions.

    NOTE: The viewer functions typically support a viewer_option parameter. viewer_option=0 is expected to  \
    refer to the analysis itself.

    * ``v_qslice``: Retrieve/compute data slices as requested via qslice URL requests. The corresponding view \
    of the DJANGO data access server already translates all input parameters and takes care of generating images/plots \
    if needed. This function is only responsible for retrieving the data.
    * ``v_qspectrum``: Retrieve/compute spectra as requested via qspectrum URL requests. The corresponding view of \
    the DJANGO data access server already translates all input parameters and takes care of generating images/plots \
    if needed. This function is only responsible for retrieving the data.
    * ``v_qmz``: Define the m/z axes for image slices and spectra as requested by qspectrum URL requests.
    * ``v_qspectrum_viewer_options``: Define a list of strings, describing the different viewer options available \
    for the analysis for qspectrum requests (i.e., ``v_qspectrum``). This feature allows the analysis developer \
    to define multiple different visualization modes for the analysis. For example, when performing a data reduction \
    (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum \
    view (v_qspectrum). By providing different viewer options we allow the user to decide which option they are \
    most interested in.
    * ``v_qslice_viewer_options``: Define a list of strings, describing the different viewer options available for \
    the analysis for qslice requests (i.e., ``v_qslice``). This feature allows the analysis developer to define \
    multiple different visualization modes for the analysis. For example, when performing a data reduction \
    (e.g., PCA or NMF) one may want to show the raw spectra or the loadings vector of the projection in the spectrum \
    view (v_qspectrum). By providing different viewer options we allow the user to decide which option they \
    are most interested in.
    """

    _analysis_instances = OrderedDict()
    """
    Class variable used to track all instances of analysis_base. Only the keys of the
    OrderedDict are used to store weak references. The values are set to None in
    all cases. I.e., we are taking advantage of the uniqueness of dict keys and the
    order-preserving feature of the OrderedDict to ensure that analyses are not
    duplicated and that references are retrieved in order of creation.
    """

    def __init__(self):
        """Initialize the basic data members"""
        super(analysis_base, self).__init__()
        self.analysis_identifier = "undefined"
        self.__data_list = []
        self.parameters = []  # Inherited from parent class parameter_data
        self.data_names = []
        self.run_info = run_info_dict()
        self.omsi_analysis_storage = []
        self._analysis_instances[weakref.ref(self)] = None  # Register the object with analysis_base._analysis_instance
        self.mpi_comm = mpi_helper.get_comm_world()
        self.mpi_root = 0
        self.update_analysis = True
        self.driver = None
        self.continue_analysis_when_ready = False  # If we have tasks waiting for us then continue the those tasks when we are done with our execution

        # Add common analysis parameters
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='profile_time_and_usage',
                           help='Enable/disable profiling of time and usage of the analysis function',
                           required=False,
                           default=False,
                           dtype=dtypes['bool'],
                           group=groups['profile'])
        self.add_parameter(name='profile_memory',
                           help='Enable/disable profiling of memory usage of the analysis function',
                           required=False,
                           default=False,
                           dtype=dtypes['bool'],
                           group=groups['profile'])

    def __call__(self, **kwargs):
        """
        Execute the analysis.

        Same a self.execute().

        :param kwargs: Keyword arguments to set parameters.

        :return: The output of the execute_analysis(...) function
        """
        return self.execute(**kwargs)

    def __getitem__(self,
                    key):
        """
        This class supports basic slicing to access data stored in the main member variables.
        By default the data is retrieved from __data_list and the __getitem__(key) function.
        which implemtent the [..] operator, returns __data_list[key]['data']. The key is
        a string indicating the name of the parameter to be retrieved. If the key is not
        found in the __data_list then the function will try to retrieve the data from
        the self.parameters list.
        """
        if isinstance(key, str) or isinstance(key, unicode):
            for i in self.__data_list:
                if i['name'] == key:
                    return i['data']
            for i in self.parameters:
                if i['name'] == key:
                    return i.get_data_or_default()
            if key in self.data_names:
                return dependency_dict(param_name=None,
                                       link_name=None,
                                       dataname=key,
                                       omsi_object=self,
                                       selection=None,
                                       help=None)
        elif isinstance(key, tuple):
            if len(key) == 2:
                if isinstance(key[0], str) or isinstance(key[0], unicode):
                    if key[0] in self.data_names:
                        return dependency_dict(param_name=None,
                                               link_name=None,
                                               dataname=key[0],
                                               omsi_object=self,
                                               selection=key[1],
                                               help=None)

        return None

    def __setitem__(self,
                    key,
                    value):
        """
        Set values in the __data, or parameters lists.
        If the given key is found in the parameters list then it is assigned
        to the paramerters/dependencies, otherwise the key is assumed to be a
        an output that needs to be added to the __data_list
        """
        if key in self.get_parameter_names():
            self.set_parameter_values(**{key: value})
            self.omsi_analysis_storage = []
        elif key in self.data_names:
            if isinstance(value, dependency_dict):
                ana_data = analysis_data(name=key,
                                         data=value)
            elif 'numpy' not in str(type(value)):
                temp_value = np.asarray(value)
                ana_data = analysis_data(name=key,
                                         data=temp_value,
                                         dtype=temp_value.dtype)
            else:
                ana_data = analysis_data(name=key,
                                         data=value,
                                         dtype=value.dtype)
            index = None
            for obj_index, ana_data_obj in enumerate(self.__data_list):
                if ana_data_obj['name'] == key:
                    index = obj_index
                    break
            if index is None:
                self.__data_list.append(ana_data)
            else:
                self.__data_list[index] = ana_data
            self.omsi_analysis_storage = []
        else:
            raise KeyError('Invalid key. The given key was not found as part of the analysis parameters nor output.')

    @classmethod
    def get_analysis_instances(cls):
        """
        Generator function used to iterate through all instances of analysis_base.
        The function creates references for all weak references stored in cls._analysis_instances
        and returns the references if it exists and cleans up the any invalid references after the
        iteration is complete.
        :return: References to analysis_base objects
        """
        # 1) Yield all valid references in order of creation
        invalid_references = set()   # Collect invalid references for clean-up
        for ref in cls._analysis_instances.keys():
            obj = ref()              # Retrieve the reference
            if obj is not None:
                yield obj            # Yield the reference if it is valid
            else:
                invalid_references.add(ref)   # Record invalid reference for removal
        # 2) Remove all invalid references
        for ref in invalid_references:
            cls._analysis_instances.pop(ref)

    @classmethod
    def locate_analysis(cls,
                        data_object,
                        include_parameters=False):
        """
        Given a data_object try to locate the analysis that creates the object as an
        output of its execution (and optionally analyses that have the object as an input).

        :param data_object: The data object of interest.
        :param include_parameters: Boolean indicating whether also input parameters should be considered
            in the search in addition to the outputs of an analysis

        :return: dependency_dict pointing to the relevant object or None in case the
            object was not found.
        """
        for ana_obj in cls.get_analysis_instances():
            ana_params_and_outputs = ana_obj.get_analysis_data_names()
            if include_parameters:
                ana_params_and_outputs += ana_obj.get_parameter_names()
            for dataname in ana_params_and_outputs:
                obj = ana_obj[dataname]
                if obj is data_object:
                    # These types are hashed by Python so we cannot uniquely decide from which analysis they come from
                    # e.g., if two analysis return the string 'test', then the condition of obj is data_object will be
                    # True in both cases.
                    if type(obj) not in (float, int, bool, long, complex, str, unicode):
                        return dependency_dict(param_name=None,
                                               link_name=None,
                                               dataname=dataname,
                                               omsi_object=ana_obj,
                                               selection=None,
                                               help=None)

        return None

    def update_analysis_parameters(self, **kwargs):
        """
        Record the analysis parameters passed to the execute() function.

        The default implementation simply calls the set_parameter_values(...) function.
        This function may be overwritten to customize the behavior of how parameters
        are recorded by the execute function.

        :param kwargs: Dictionary of keyword arguments with the parameters passed to the execute(..) function

        """
        self.set_parameter_values(**kwargs)

    def define_missing_parameters(self):
        """
        Called by the execute function before self.update_analysis_parameters
        to set any required parameters that have not been defined to their respective default values.

        This function may be overwritten in child classes to customize
        the definition of default parameter values and to apply any
        modifications (or checks) of parameters before the analysis is executed.
        Any changes applied here will be recorded in the parameter of the analysis.
        """
        super(analysis_base, self).define_missing_parameters()

    def check_ready_to_execute(self):
        """
        Check if all inputs are ready to determine if the analysis is ready to run.

        :return: List of omsi_analysis_parameter objects that are not ready. If the
            returned list is empty, then the analysis is ready to run.
        """
        pending_inputs = []
        for param in self.parameters:
            if not param.data_ready():
                pending_inputs.append(param)
        return pending_inputs

    def continue_workflow_when_ready(self, executor):
        """
        Continue the workflow defined by the given executor when ready
        :param executor:
        :return:
        """
        self.driver = executor
        self.continue_analysis_when_ready = True

    def outputs_ready(self):
        """
        Call this function to indicate that outputs are ready to use. This function is only used
        in case not all outputs are defined by the execute_analysis function but instead generate
        other dependencies. This may be the case, e.g., when we generate interactive apps where
        outputs depend on user inputs. In this case, this function must be called manually by
        the user to indicate that the outputs are now ready.
        """
         # Continue running tasks that are waiting on us
        if self.continue_analysis_when_ready:
            log_helper.info(__name__, "Outputs of " + str(self) + "are ready. Trying to continue blocked tasks.",
                                     root=self.mpi_root, comm=self.mpi_comm)
            self.driver.execute()

    def record_execute_analysis_outputs(self, analysis_output):
        """
        Function used internally by execute to record the output
        of the custom execute_analysis(...) function to the __data_list.

        This function may be overwritten in child classes in order to
        customize the behavior for recording data outputs. Eg., for some
        analyses one may only want to record a particular set of outputs,
        rather than all outputs generated by the analysis.

        :param analysis_output: The output of the execute_analysis(...) function to be recorded
        """
        # Record the analysis output so that we can save it to file
        if analysis_output is not None:
            if len(self.data_names) == 1:   # We need this case, because analysis_output is not a tuple we can slice
                log_helper.debug(__name__, "Recording output " + str(self.data_names[0]),
                                 root=self.mpi_root, comm=self.mpi_comm)
                self[self.data_names[0]] = analysis_output
            else:
                for data_index, data_name in enumerate(self.data_names):
                    log_helper.debug(__name__, "Recording output " + str(data_name),
                                     root=self.mpi_root, comm=self.mpi_comm)
                    self[data_name] = analysis_output[data_index]

    def execute(self, **kwargs):
        """
        Use this function to run the analysis.

        :param kwargs: Parameters to be used for the analysis. Parameters may also be set using
            the __setitem__ mechanism or as batches using the set_parameter_values function.

        :returns: This function returns the output of the execute analysis function.

        :raises: AnalysisReadyError in case that the analysis is not ready to be executed. This may be
            the case, e.g, when a dependent input parameter is not ready to be used.

        """
        log_helper.debug(__name__, "Execute analysis. " + str(self), root=self.mpi_root, comm=self.mpi_comm)

        log_helper.debug(__name__, "Initializing analysis parameters and environment. " + str(self),
                         root=self.mpi_root, comm=self.mpi_comm)
        # 1) Remove the saved analysis object since we are running the analysis again
        self.omsi_analysis_storage = []

        # 2) Define all parameters and make sure that they are ready
        # 2.1) Set any parameters that are given to the execute function
        self.update_analysis_parameters(**kwargs)

        # 2.2) Set the parameters that are required that have a default value but that have not been initialized
        self.define_missing_parameters()

        # 2.3) Check that all parameters are ready to be used
        pending_params = self.check_ready_to_execute()
        if len(pending_params) > 0:
            raise AnalysisReadyError("The analysis is not ready.", pending_params)

        # 3) Record basic execution provenance while running the analysis
        log_helper.debug(__name__, "Run and profile the analysis: " + str(self),
                         root=self.mpi_root, comm=self.mpi_comm)

        if self.run_info is not None:
            self.run_info.mpi_comm = self.mpi_comm
            self.run_info.mpi_root = self.mpi_root
            self.enable_time_and_usage_profiling(self['profile_time_and_usage'])
            self.enable_memory_profiling(self['profile_memory'])
            analysis_output = self.run_info(self.execute_analysis)()
        else:
            analysis_output = self.execute_analysis()

        log_helper.debug(__name__, "Finished the analysis. " + str(self), root=self.mpi_root, comm=self.mpi_comm)

        # Record the analysis output
        self.record_execute_analysis_outputs(analysis_output=analysis_output)

        # Indicate the analysis is up-to-date
        self.update_analysis = False

        # Return the output of the analysis
        return analysis_output

    def execute_analysis(self):
        """
        Implement this function to implement the execution of the actual analysis.

        This function may not require any input parameters. All input parameters are
        recorded in the parameters and dependencies lists and should be retrieved
        from there, e.g, using basic slicing self[ paramName ]

        Input parameters may be added for internal use ONLY. E.g, we may add parameters that
        are used internally to help with parallelization of the execute_analysis function.
        Such parameters are not recorded and must be strictly optional so that analysis_base.execute(...)
        can call the function.

        :returns: This function may return any developer-defined data. Note, all
                 output that should be recorded must be put into the data list.

        """
        raise NotImplementedError("Implement execute_analysis in order to be able to run the analysis.")

    def execute_recursive(self,
                          **kwargs):
        """
        Recursively execute this analysis and all its dependencies if necessary

        We use a workflow driver to control the execution. To define the workflow driver we can set
        the self.driver variable. If no workflow driver is given (i.e, self.driver==None), then the
        default driver will be created. To change the default driver,
        see `omsi.workflow.base.workflow_executor_base.DEFAULT_EXECUTOR_CLASS`

        :param kwargs: Parameters to be used for the analysis. Parameters may also be set using
            the __setitem__ mechanism or as batches using the set_parameter_values function.

        :return: Same as execute
        """
        log_helper.info(__name__, "Executing the analysis and all its dependencies. " + str(self),
                        root=self.mpi_root, comm=self.mpi_comm)
        self.update_analysis = True   # Force update of the current analysis
        if len(kwargs) > 0:
            self.update_analysis_parameters(**kwargs)
        if self.driver is not None:
            log_helper.debug(__name__, "Using user-defined driver.",  root=self.mpi_root, comm=self.mpi_comm)
            self.driver.clear()
            self.driver.add_analysis(self)
            self.driver.execute()
        else:
            log_helper.debug(__name__, "Creating default driver and workflow to run the analysis.",
                             root=self.mpi_root, comm=self.mpi_comm)
            default_driver = workflow_executor_base.get_default_executor(analysis_objects=self)
            default_driver.execute()
        log_helper.debug(__name__, "Compiling outputs and return", root=self.mpi_root, comm=self.mpi_comm)
        outputs = [self[name] for name in self.data_names]
        return tuple(outputs)

    @classmethod
    def execute_all(cls,
                    force_update=False,
                    executor=None):
        """
        Execute all analysis instances that are currently defined.

        :param force_update: Boolean indicating whether we should force that all analyses are
            executed again, even if they have already been run with the same settings before.
            False by default.
        :param executor: Optional workflow executor to be used for the execution of all analyses.
            The executor will be cleared and then all analyses will be added to executor. Default
            value is None, in which case the function creates a default executor to be used.

        :return: The workflow executor used
        """
        root = 0 if executor is None else executor.mpi_root
        comm = mpi_helper.get_comm_world() if executor is None else executor.mpi_comm
        log_helper.debug(__name__, "Execute all analyses", root=root, comm=comm)
        log_helper.log_var(__name__, force_update=force_update, level='DEBUG', root=root, comm=comm)
        if executor is not None:
            log_helper.debug(__name__, "Using user-defined executor.",  root=root, comm=comm)
            executor.clear()
        else:
            log_helper.debug(__name__, "Creating the default executor and adding all analysis objects",
                             root=root, comm=comm)
            executor = workflow_executor_base.get_default_executor()
        for ana_obj in cls.get_analysis_instances():
            if force_update:
                ana_obj.update_analysis = True
            executor.add_analysis(ana_obj)
        log_helper.debug(__name__, "Execute the workflow using the default executor", root=root, comm=comm)
        executor.execute()
        return executor

    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """
        Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param z: Selection string indicting which z values should be selected.
        :param viewer_option: If multiple default viewer behaviors are available for a given analysis
                            then this option is used to switch between them.

        :returns: numpy array with the data to be displayed in the image slice viewer. Slicing will be
                 performed typically like [:,:,zmin:zmax].

        :raises: NotImplementedError in case that v_qslice is not supported by the analysis.
        """
        from omsi.analysis.analysis_views import analysis_views
        from omsi.shared.data_selection import check_selection_string, \
            selection_type, \
            selection_string_to_object
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index = \
            cls.__construct_dependent_viewer_options__(analysis_object)

        # Check whether the given selection is valid
        z_type = check_selection_string(z)
        if z_type == selection_type['invalid']:
            return None
        if isinstance(re_slicedata[viewer_option], omsi_file_msidata):
            try:
                z_select = selection_string_to_object(selection_string=z)
                data = re_slicedata[viewer_option][:, :, z_select]
                return data
            except:
                return None

        elif isinstance(re_slicedata[viewer_option], omsi_file_analysis):
            current_analysis_type = str(re_slicedata[viewer_option].get_analysis_type()[0])
            return analysis_views.analysis_name_to_class(current_analysis_type).v_qslice(
                analysis_object=re_slicedata[viewer_option],
                z=z,
                viewer_option=re_slice_option_index[viewer_option])
        else:
            return None

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """
        Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer

        **Developer Note:** h5py currently supports only a single index list. If the user provides an index-list
        for both x and y, then we need to construct the proper merged list and load the data manually, or if
        the data is small enough, one can load the full data into a numpy array which supports
        multiple lists in the selection.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param x: x selection string
        :param y: y selection string
        :param viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them.

        :returns: The following two elements are expected to be returned by this function :

                1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The mass (m/z) axis must be \
                the last axis. For index selection x=1,y=1 a 1D array is usually expected. For indexList \
                selections x=[0]&y=[1] usually a 2D array is expected. For ragne selections x=0:1&y=1:2 we \
                one usually expects a 3D array.
                2) None in case that the spectra axis returned by v_qmz are valid for the returned spectrum. \
                Otherwise, return a 1D numpy array with the m/z values for the spectrum (i.e., if custom m/z \
                values are needed for interpretation of the returned spectrum).This may be needed, e.g., in \
                cases where a per-spectrum peak analysis is performed and the peaks for each spectrum appear \
                at different m/z values.
        """
        from omsi.analysis.analysis_views import analysis_views
        from omsi.shared.data_selection import \
            check_selection_string, \
            selection_type, \
            selection_string_to_object   # , selection_to_indexlist
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index = \
            cls.__construct_dependent_viewer_options__(analysis_object)
        # Check whether the given selection is valid
        x_type = check_selection_string(x)
        y_type = check_selection_string(y)
        x_selection = selection_string_to_object(x)
        y_selection = selection_string_to_object(y)

        if (x_type == selection_type['invalid']) or (y_type == selection_type['invalid']):
            return None, None
        if isinstance(re_spectrumdata[viewer_option], omsi_file_msidata):
            if x_type == selection_type['invalid'] or y_type == selection_type['invalid']:
                return None, None
            if x_type == selection_type['indexlist'] and y_type == selection_type['indexlist']:
                # We now need to match up the index lists and load all the individial values
                # ToDo: This is can be very inefficient for large lists
                x_size = len(x_selection)
                y_size = len(y_selection)
                if x_size != y_size:
                    raise KeyError("Selection lists don't match")
                current_dataset = re_spectrumdata[viewer_option]
                z_size = current_dataset.shape[2]
                # Allocate the required memory
                data = np.zeros((x_size, z_size), dtype=current_dataset.dtype)
                for i in xrange(0, x_size):
                    data[i, :] = current_dataset[x_selection[i], y_selection[i], :]
            else:
                data = re_spectrumdata[viewer_option][x_selection, y_selection, :]

            return data, None
        elif isinstance(re_spectrumdata[viewer_option], omsi_file_analysis):
            current_analysis_type = str(re_spectrumdata[viewer_option].get_analysis_type()[0])
            return analysis_views.analysis_name_to_class(current_analysis_type).v_qspectrum(
                analysis_object=re_spectrumdata[viewer_option],
                x=x,
                y=y,
                viewer_option=re_spectrum_option_index[viewer_option])
        else:
            return None, None

    @classmethod
    def v_qmz(cls,
              analysis_object,
              qslice_viewer_option=0,
              qspectrum_viewer_option=0):
        """
        Get the mz axes for the analysis

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param qslice_viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them for the qslice URL pattern.
        :param qspectrum_viewer_option: If multiple default viewer behaviors are available for a
            given analysis then this option is used to switch between them for the qspectrum URL pattern.

        :returns: The following four arrays are returned by the analysis:

            - mzSpectra : Array with the static mz values for the spectra.
            - labelSpectra : Label for the spectral mz axis
            - mzSlice : Array of the static mz values for the slices or None if identical to the mzSpectra.
            - labelSlice : Label for the slice mz axis or None if identical to labelSpectra.
            - values_x: The values for the x axis of the image (or None)
            - label_x: Label for the x axis of the image
            - values_y: The values for the y axis of the image (or None)
            - label_y: Label for the y axis of the image
            - values_z: The values for the z axis of the image (or None)
            - label_z: Label for the z axis of the image

        """
        from omsi.analysis.analysis_views import analysis_views
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index = \
            cls.__construct_dependent_viewer_options__(analysis_object)
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None
        valuesX = None
        labelX = None
        valuesY = None
        labelY = None
        valuesZ = None
        labelZ = None
        # Determine the spectra mz axis
        if len(re_spectrum) > 0:
            if isinstance(re_spectrumdata[qspectrum_viewer_option], omsi_file_msidata):
                mz_spectra = re_spectrumdata[qspectrum_viewer_option].mz[:]
                label_spectra = "m/z"
                valuesX = range(0, re_spectrumdata[qspectrum_viewer_option].shape[0])
                valuesY = range(0, re_spectrumdata[qspectrum_viewer_option].shape[1])
                is3Dimage = len(re_spectrumdata[qspectrum_viewer_option].shape)
                valuesZ = None if is3Dimage <= 3 else range(0, re_spectrumdata[qspectrum_viewer_option].shape[2])
                labelX = 'pixel index X'
                labelY = 'pixel index Y'
                labelZ = None if is3Dimage else 'pixel index Z'
            elif isinstance(re_spectrumdata[qspectrum_viewer_option], omsi_file_analysis):
                mz_spectra, label_spectra, temp_a, temp_b, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                    analysis_views.get_axes(re_spectrumdata[qspectrum_viewer_option],
                                            qslice_viewer_option=re_slice_option_index[qslice_viewer_option],
                                            qspectrum_viewer_option=re_spectrum_option_index[qspectrum_viewer_option])
            else:
                mz_spectra = None
                label_spectra = None
        # Determine the slice mz axis
        if len(re_slice) > 0:
            if isinstance(re_slicedata[qslice_viewer_option], omsi_file_msidata):
                mz_slice = re_slicedata[qslice_viewer_option].mz[:]
                label_slice = "m/z"
            elif isinstance(re_slicedata[qslice_viewer_option], omsi_file_analysis):
                temp_a, temp_b, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                    analysis_views.get_axes(
                        re_slicedata[qslice_viewer_option],
                        qslice_viewer_option=re_slice_option_index[qslice_viewer_option],
                        qspectrum_viewer_option=re_spectrum_option_index[qspectrum_viewer_option])
            else:
                mz_slice = None
                label_slice = None

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """
        Get a list of strings describing the different default viewer options for the analysis for qspectrum.
        The default implementation tries to take care of handling the spectra retrieval for all the dependencies
        but can naturally not decide how the qspectrum should be handled by a derived class. However, this
        implementation is often called at the end of custom implementations to also allow access to data from
        other dependencies.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.
            For most cases this is not needed here as the support for slice operations is usually a static decision
            based on the class type, however, in some cases additional checks may be needed (e.g., ensure that the
            required data is available).

        :returns: List of strings indicating the different available viewer options. The list should be empty if
            the analysis does not support qspectrum requests (i.e., v_qspectrum(...) is not available).

        """
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index = \
            cls.__construct_dependent_viewer_options__(analysis_object)
        return re_spectrum

    @classmethod
    def v_qslice_viewer_options(cls,
                                analysis_object):
        """
        Get a list of strings describing the different default viewer options for the analysis for qslice.
        The default implementation tries to take care of handling the spectra retrieval for all the dependencies
        but can naturally not decide how the qspectrum should be handled by a derived class. However, this
        implementation is often called at the end of custom implementations to also allow access to data from
        other dependencies.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.  For most cases
            this is not needed here as the support for slice operations is usually a static decission based on
            the class type, however, in some cases additional checks may be needed (e.g., ensure that the required
            data is available).

        :returns: List of strings indicating the different available viewer options. The list should be empty
            if the analysis does not support qslice requests (i.e., v_qslice(...) is not available).

        """
        re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index = \
            cls.__construct_dependent_viewer_options__(analysis_object)
        return re_slice

    @classmethod
    def __construct_dependent_viewer_options__(cls,
                                               analysis_object):
        """
        Internal helper function to construct the viewer options for analysis dependencies.

        :returns: The following 4 lists:

        * ``re_slice`` : List of names for the slice options
        * ``re_spectrum`` : List of names for the spectrum options
        * ``re_slicedata`` : omsi_file_* API object associated with the slice options. These \
        are either omis_file_analysis or omsi_file_msidata objects.
        * ``re_spectrumdata`` : omsi_file_* API object associated with the slice options. These \
        are either omis_file_analysis or omsi_file_msidata objects.
        * ``re_slice_optionIndex``: List of integrers indicating for the given entry the viewer_option \
        to be used with the re_slicedata object for the given option.
        * ``re_spectrum_optionIndex``: List of integrers indicating for the given entry the viewer_option \
        to be used with the re_spectrumdata object for the given option.
        """
        from omsi.analysis.analysis_views import analysis_views
        re_slice = []
        re_slicedata = []
        re_slice_option_index = []
        re_spectrum = []
        re_spectrumdata = []
        re_spectrum_option_index = []

        # We don't need to use get_all_dependency_data_recursive here, because when we call
        # analysis_views.get_qslice_ .. (spectrum etc.) the recursion to dependent options
        # occurs automatically.
        all_dependencies = analysis_object.get_all_dependency_data()
        for di in all_dependencies:
            # Check if we can slice the data
            if isinstance(di['omsi_object'], omsi_file_msidata):
                re_spectrum.append("Raw Data: " + di['link_name'])
                re_slice.append("Raw Data: " + di['link_name'])
                re_slice_option_index.append(0)
                re_spectrumdata.append(di['omsi_object'])
                re_slicedata.append(di['omsi_object'])
                re_spectrum_option_index.append(0)
            elif isinstance(di['omsi_object'], omsi_file_analysis):
                slice_options = analysis_views.get_qslice_viewer_options(di['omsi_object'])
                spectrum_options = analysis_views.get_qspectrum_viewer_options(di['omsi_object'])
                for sloption_index in range(0, len(slice_options)):
                    re_slice.append(slice_options[sloption_index])
                    re_slicedata.append(di['omsi_object'])
                    re_slice_option_index.append(sloption_index)
                for sloption_index in range(0, len(spectrum_options)):
                    re_spectrum.append(spectrum_options[sloption_index])
                    re_spectrumdata.append(di['omsi_object'])
                    re_spectrum_option_index.append(sloption_index)
                # analysisType = str(anaObj.get_analysis_type()[0])
                # if analysis_views.supports_slice( di['omsi_object']) :
                #     re_slice.append( "Analysis: "+str(di['omsi_object'].get_analysis_identifier()[0]) )
                #     re_slicedata.append( di['omsi_object'] )
                # if analysis_views.supports_spectra( di['omsi_object'] ) :
                #     re_spectrum.append( "Analysis: "+str(di['omsi_object'].get_analysis_identifier()[0]) )
                #     re_spectrumdata.append( di['omsi_object'] )
            else:
                warnings.warn("Unknown dependency")
        return re_slice, re_spectrum, re_slicedata, re_spectrumdata, re_slice_option_index, re_spectrum_option_index

    @staticmethod
    def get_default_dtypes():
        """
        Get a list of available default dtypes used for analyses.
        Same as `data_dtypes.get_dtypes()`.
        """
        return data_dtypes.get_dtypes()

    @staticmethod
    def get_default_parameter_groups():
        """
        Get a list of commonly used parameter groups and associated descriptions.

        Use of default groups provides consistency and allows other system to
        design custom behavior around the semantic of parameter groups

        :return: Dictionary where the keys are the short names of the groups and the
            values are dicts with following keys:value pairs: 'name' , 'description'.
            Use the 'name' to define the group to be used.
        """
        return {'stop': {'name': 'stop conditions',
                         'description': 'Termination criteria for the analysis'},
                'input': {'name': 'input data',
                          'description': 'Input data to be analyzed'},
                'settings': {'name': 'analysis settings',
                             'description': 'Analysis settings'},
                'parallel': {'name': 'parallel settings',
                             'description': 'Parallel execution settings'},
                'profile': {'name': 'profile analysis',
                            'description': 'Profile the exeution of the analysis'}}

    def enable_time_and_usage_profiling(self, enable=True):
        """
        Enable or disable profiling of time and usage of code parts of execute_analysis.

        :param enable: Enable (True) or disable (False) profiling
        :type enable: bool

        :raises: ImportError is raised if a required package for profiling is not available.
        """
        if enable != self.run_info.get_profile_time_and_usage():
            self.run_info.enable_profile_time_and_usage(enable=enable)
        if self.run_info.get_profile_time_and_usage() != self['profile_time_and_usage']:
            # NOTE when assignign the paramter this function will be called again by the set_parameter_values function
            self['profile_time_and_usage'] = self.run_info.get_profile_time_and_usage()

    def enable_memory_profiling(self, enable=True):
        """
        Enable or disable line-by-line profiling of memory usage of execute_analysis.

        :param enable_memory: Enable (True) or disable (False) line-by-line profiling of memory usage
        :type enable_memory: bool

        :raises: ImportError is raised if a required package for profiling is not available.
        """
        if enable != self.run_info.get_profile_memory():
            self.run_info.enable_profile_memory(enable=enable)
        if self.run_info.get_profile_memory() != self['profile_memory']:
            # NOTE when assignign the paramter this function will be called again by the set_parameter_values function
            self['profile_memory'] = self.run_info.get_profile_memory()

    def results_ready(self):
        """
        Check whether the results of the analysis are ready to be used
        :return: Boolean
        """
        if not self.update_analysis:
            for output in self.data_names:
                if isinstance(self[output], dependency_dict):
                    return False
            return True
        return not self.update_analysis

    def get_memory_profile_info(self):
        """
        Based on the memory profile of the execute_analysis(..) function get
        the string describing the line-by-line memory usage.

        :return: String describing the memory usage profile. None is returned in case that
            no memory profiling data is available.
        """
        if 'profile_mem_stats' in self.run_info:
            return self.run_info['profile_mem_stats']
        else:
            return None

    def get_profile_stats_object(self, consolidate=True, stream=None):
        """
        Based on the execution profile of the execute_analysis(..) function get
        ``pstats.Stats`` object to help with the interpretation of the data.

        :param consolidate: Boolean flag indicating whether multiple stats (e.g., from multiple cores)
            should be consolidated into a single stats object. Default is True.
        :param stream: The optional stream parameter to be used fo the pstats.Stats object.

        :return: A single pstats.Stats object if consolidate is True. Otherwise the function
            returns a list of pstats.Stats objects, one per recorded statistic. None is returned
            in case that the stats objects cannot be created or no profiling data is available.
        """
        return self.run_info.get_profile_stats_object(consolidate=consolidate, stream=stream)

    def get_help_string(self):
        """
        Get a string describing the analysis.

        :return: Help string describing the analysis and its parameters
        """
        from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
        temp_driver = cl_analysis_driver(analysis_class=self.__class__)
        temp_driver.initialize_argument_parser()
        return temp_driver.parser.format_help()

    def get_omsi_analysis_storage(self):
        """
        Get a list of known locations where this analysis has been saved.

        :return: List of `omsi.dataformat.omsi_file.analysis. omsi_file_analysis` objects where the analysis is saved.
        """
        return self.omsi_analysis_storage

    def has_omsi_analysis_storage(self):
        """
        Check whether a storage location is known where the anlaysis has been saved.

        :return: Boolean indicating whether self.omsi_analysis_storage is not empty
        """
        return len(self.omsi_analysis_storage) > 0

    def get_analysis_type(self):
        """
        Return a string indicating the type of analysis performed
        """
        if self.__module__ != '__main__':
            return self.__module__
        else:
            return self.__class__.__name__

    def get_analysis_data_names(self):
        """
        Get a list of all analysis dataset names.
        """
        return self.data_names

    def get_parameter_names(self):
        """
        Get a list of all parameter dataset names (including those that may define
        dependencies.
        """
        return super(analysis_base, self).get_parameter_names()

    def get_analysis_data(self,
                          index):
        """
        Given the index return the associated dataset to be written to the HDF5 file

        :param index : Retrun the index entry of the private member __data_list.
        """
        return self.__data_list[index]

    def get_parameter_data(self,
                           index):
        """
        Given the index return the associated dataset to be written to the HDF5 file

        :param index : Return the index entry of the private member parameters.
        """
        return super(analysis_base, self).get_parameter_data(index=index)

    def get_analysis_data_by_name(self,
                                  dataname):
        """
        Given the key name of the data return the associated analysis_data object.

        :param dataname: Name of the analysis data requested from the private __data_list member.

        :returns: The analysis_data object or None if not found.
        """
        for i in self.__data_list:
            if i['name'] == dataname:
                return i
        return None

    def get_parameter_data_by_name(self,
                                   dataname):
        """
        Given the key name of the data return the associated parameter_data object.

        :param dataname: Name of the parameter requested from the parameters member.

        :returns: The parameter_data object or None if not found
        """
        return super(analysis_base, self).get_parameter_data_by_name(dataname=dataname)

    def get_all_run_info(self):
        """Get the dict with the complete info about the last run of the analysis"""
        return self.run_info

    def get_all_analysis_data(self):
        """Get the complete list of all analysis datasets to be written to the HDF5 file"""
        return self.__data_list

    def get_all_parameter_data(self,
                               exclude_dependencies=False):
        """
        Get the complete list of all parameter datasets to be written to the HDF5 file

        :param exclude_dependencies: Boolean indicating whether we should exclude parameters
            that define dependencies from the list
        """
        return super(analysis_base, self).get_all_parameter_data(exclude_dependencies=exclude_dependencies)

    def get_all_dependency_data(self):
        """
        Get the complete list of all direct dependencies to be written to the HDF5 file

        NOTE: These are only the direct dependencies as specified by the analysis itself.
        Use  get_all_dependency_data_recursive(..) to also get the indirect dependencies of
        the analysis due to dependencies of the dependencies themselves.

        :returns: List of parameter_data objects that define dependencies.

        """
        return super(analysis_base, self).get_all_dependency_data()

    def get_num_analysis_data(self):
        """Retrun the number of analysis datasets to be wirtten to the HDF5 file"""
        return len(self.__data_list)

    def get_num_parameter_data(self):
        """Return the number of parameter datasets to be wirtten to the HDF5 file"""
        return super(analysis_base, self).get_num_parameter_data()

    def get_num_dependency_data(self):
        """Return the number of dependencies to be wirtten to the HDF5 file"""
        return super(analysis_base, self).get_num_dependency_data()

    def keys(self):
        """
        Get a list of all valid keys, i.e., a combination of all input parameter and output names.

        :return: List of strings with all input parameter and output names.
        """
        return self.get_analysis_data_names() + self.get_parameter_names()

    def clear_analysis_data(self):
        """Clear the list of analysis data"""
        log_helper.debug(__name__, "Clearing analysis data. ", root=self.mpi_root, comm=self.mpi_comm)
        self.__data_list = []
        self.update_analysis = True

    def clear_parameter_data(self):
        """Clear the list of parameter data"""
        log_helper.debug(__name__, "Clearing parameter data. ", root=self.mpi_root, comm=self.mpi_comm)
        for param in self.parameters:
            param.clear_data()
        self.update_analysis = True

    def clear_run_info_data(self):
        """Clear the runtime information data"""
        log_helper.debug(__name__, "Clearing runtime information data. ", root=self.mpi_root, comm=self.mpi_comm)

    def clear_analysis(self):
        """Clear all analysis data---i.e., parameter, dependency data, output results, runtime data"""
        log_helper.debug(__name__, "Clearing the analysis. ", root=self.mpi_root, comm=self.mpi_comm)
        self.clear_analysis_data()
        self.clear_parameter_data()
        self.clear_run_info_data()
        self.update_analysis = True

    def clear_and_restore(self, analysis_manager=None, resave=False):
        """
        Clear all analysis data and restore the results from file

        :param analysis_manager: Instance of omsi_analysis_manager (e.g., an omsi_file_experiment) where the
            analysis should be saved.
        :param resave: Boolean indicating whether the analysis should be saved again, even if it has been
            saved before. This parameter only has effect if analysis_manager is given.

        :return: self, i.e., the updated analysis object with all data replaced with HDF5 references
        """
        log_helper.debug(__name__, "Clearing and restoring the analysis", root=self.mpi_root, comm=self.mpi_comm)
        from tempfile import NamedTemporaryFile
        from omsi.dataformat.omsi_file.main_file import omsi_file

        # 1) Write the analysis to file if necessary
        if not self.has_omsi_analysis_storage() or (analysis_manager is not None and resave):
            # 1.1) Create a temporary file for storage if needed
            named_temp_file = None
            if analysis_manager is None:
                named_temp_file = NamedTemporaryFile(suffix=".h5")
                temp_omsi_file = omsi_file(named_temp_file.name, 'a')
                analysis_manager = temp_omsi_file.create_experiment()
                log_helper.debug(__name__, "Created temporary file store " + named_temp_file.name,
                                 root=self.mpi_root, comm=self.mpi_comm)
            # 1.2) Save the analysis to file
            ana_obj, ana_index = analysis_manager.create_analysis(analysis=self,
                                             flush_io=True,
                                             force_save=False,
                                             save_unsaved_dependencies=True,
                                             mpi_root=self.mpi_root,
                                             mpi_comm=self.mpi_comm)
            # 1.3) Make sure that the named temporary file we created does not go out of
            # scope until our data store is deleted for the analysis, by attaching the file to the object
            setattr(ana_obj, named_temp_file.name, named_temp_file)

        # 2) Clear the in-memory data
        self.clear_analysis()
        # 3) Restore the data from file
        print self.get_omsi_analysis_storage()[0].managed_group.items()
        self.read_from_omsi_file(analysis_object=self.get_omsi_analysis_storage()[0],
                                 load_data=True,
                                 load_parameters=True,
                                 load_runtime_data=False,
                                 dependencies_omsi_format=False,
                                 ignore_type_conflict=False)
        # 4) Return the updated object
        return self

    def set_parameter_values(self,
                             **kwargs):
        """
        Set all parameters given as input to the function. The inputs
        are placed in the self.parameters list. If the parameter refers
        to an existing h5py.Dataset, h5py.Group,  managed h5py object,
        or is an instance of an existing omis_analysi_base object, then
        a dependency_dict will be created and stored as value instead.

        :param kwargs: Dictionary of keyword arguments. All keys are
               expected to be strings. All values are expected to be
               either i) numpy arrays, ii) int, float, str or unicode
               variables, iii) h5py.Dataset or  h5py.Group, iv) or any
               the omsi_file API class objects. For iii) and iv) one
               may provide a tuple consisting of the dataobject t[0] and
               an additional selection string t[1].
        """
        log_helper.debug(__name__, "Setting analysis parameters. " + str(kwargs.keys()),
                         root=self.mpi_root, comm=self.mpi_comm)
        import h5py
        from omsi.dataformat.omsi_file.common import omsi_file_common
        for k, v in kwargs.items():
            name = unicode(k)
            value = v
            selection = None
            curr_parameter = self.get_parameter_data_by_name(name)
            dtype = curr_parameter['dtype'] if curr_parameter is not None else unicode  # unicode
            if isinstance(v, tuple):
                value = v[0]
                selection = v[1]
            if isinstance(value, h5py.Dataset) or \
                    isinstance(value, h5py.Group) or \
                    omsi_file_common.is_managed(value) or \
                    isinstance(value, analysis_base):
                value = dependency_dict(param_name=name,
                                        link_name=name,
                                        omsi_object=value,
                                        selection=selection,
                                        help=curr_parameter['help'] if curr_parameter is not None else '',
                                        dependency_type=dependency_dict.dependency_types['parameter'])
                # dtype = data_dtypes.get_dtypes()['ndarray']
            elif isinstance(value, dependency_dict):
                # Set any possibly missing parameters
                value['param_name'] = name
                value['link_name'] = name
                value['dependency_type'] = dependency_dict.dependency_types['parameter']
                value['help'] = curr_parameter['help'] if curr_parameter is not None else ''
                value['dependency_type'] = dependency_dict.dependency_types['parameter']
            else:
                # Try to locate the input parameter to see if it is an output of another analysis
                check_param = self.locate_analysis(value)
                # The given parameter value is the output of another analysis so we create a link to it
                if check_param is not None:
                    check_param['param_name'] = name  # Define the name of the dependent parameter
                    check_param['link_name'] = name   # Define the name of the object when stored to file
                    check_param['dependency_type'] = dependency_dict.dependency_types['parameter']
                    # If the analysis has been saved to file already, then change the dependency to point
                    # to the corresponding file object instead
                    # if len(check_param['omsi_object'].get_omsi_analysis_storage()) > 0:
                    #    check_param['omsi_object'] = check_param['omsi_object'].omsi_analysis_storage[0]

                    # Since we have the data still in memory, add a reference to the data to the
                    # dependency to use the in-memory data rather then forcing a data load from file
                    check_param._force_set_data(value)
                    value = check_param
                # If the input could not be located as an output of another analysis
                else:
                    try:
                        dtype = value.dtype  # if the object specifies a valid numpy dtype
                    except AttributeError:
                        if isinstance(value, float) or \
                                isinstance(value, int) or \
                                isinstance(value, bool) or \
                                isinstance(value, long) or \
                                isinstance(value, complex):
                            value = np.asarray([value])
                            # dtype = value.dtype
                        elif isinstance(value, str) or isinstance(value, unicode):
                            # dtype = omsi_format_common.str_type
                            pass
                        else:
                            value = np.asarray(value)
                            # dtype = value.dtype

            # Parameter set
            if curr_parameter is not None:
                # curr_parameter['dtype'] = dtype
                curr_parameter['data'] = value
                self.update_analysis = True
            # Missing parameter. Add the parameter to ensure that we don't omit any data.
            else:
                # Using the standard implementation path, this should never happen. However, when dynamically
                # wrapping and restoring functions on-the-fly we cannot determine all parameters until later
                # (e.g., when use analysis_generic) so that we need to cover this case here.
                log_helper.warning(__name__, "Parameter " + name +
                                   " not found in analysis_base.set_parameter_values(). Adding a new parameter.",
                                   root=self.mpi_root, comm=self.mpi_comm)
                self.parameters.append(parameter_data(name=name,
                                                      help='',
                                                      dtype=dtype,
                                                      required=False,
                                                      data=value))
                self.update_analysis = True

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
        :param type: Optional type. Default is string.
        :param required: Boolean indicating whether the parameter is required (True) or optional (False). Default False.
        :param default: Optional default value for the parameter. Default None.
        :param choices: Optional list of choices with allowed data values. Default None, indicating no choices set.
        :param data: The data assigned to the parameter. None by default.
        :param group: Optional group string used to organize parameters. Default None, indicating that
            parameters are automatically organized by driver class (e.g. in required and optional parameters)

        :raises: ValueError is raised if the parameter with the given name already exists.
        """
        log_helper.debug(__name__, "Adding parameter. " + str(name), root=self.mpi_root, comm=self.mpi_comm)
        super(analysis_base, self).add_parameter(name=name,
                                                 help=help,
                                                 dtype=dtype,
                                                 required=required,
                                                 default=default,
                                                 choices=choices,
                                                 data=data,
                                                 group=group)
        self.update_analysis = True

#    def add_analysis_data(self , name, data, dtype ) :
#        """Add a new dataset to the list of data to be written to the HDF5 file
#
#          The input parameters will be transformed into a analysis_data dictionary.
#
#          :param name: The name for the dataset in the HDF5 format
#          :param data: The numpy array to be written to HDF5. The data write function
#                omsi_file_experiment.create_analysis used for writing of the data to file can
#                in principal also handel other primitive data types by explicitly converting them
#                to numpy. However, in this case the dtype is determined based on the numpy conversion
#                and correct behavior is not guaranteed. I.e., even single scalars should be stored as
#                a 1D numpy array here.
#          :param dtype: The data type to be used during writing. For standard numpy data types this is just
#                 the dtype  of the dataset, i.e., ['data'].dtype. Other allowed datatypes are:
#
#                 * For string:  omsi_format.str_type (omsi_format is located in omsi.dataformat.omsi_file )
#                 * To generate data links: ana_hdf5link   (analysis_data)
#        """
#        self.__data_list.append( analysis_data(name=name, data=data, dtype=dtype) )

#    def add_parameter_data(self , name, data, dtype ) :
#        """Add a new analysis parameter dataset to the list of data to be written to the HDF5 file
#
#           The input parameters will be transformed into a analysis_data dictionary.
#
#          :param name: The name for the dataset in the HDF5 format
#          :param data: The numpy array to be written to HDF5. The data write function
#                omsi_file_experiment.create_analysis used for writing of the data to file can
#                in principal also handel other primitive data types by explicitly converting them
#                to numpy. However, in this case the dtype is determined based on the numpy conversion
#                and correct behavior is not guaranteed. I.e., even single scalars should be stored as
#                a 1D numpy array here.
#          :param dtype: The data type to be used during writing. For standard numpy data types this is just
#                 the dtype  of the dataset, i.e., ['data'].dtype. Other allowed datatypes are:
#
#                 * For string:  omsi_format.str_type (omsi_format is located in omsi.dataformat.omsi_file )
#                 * To generate data links: ana_hdf5link   (analysis_data)
#
#        """
#        self.__parameter_list.append( analysis_data(name=name, data=data, dtype=dtype) )
#

#    def add_dependency_data(self , dependency ) :
#        """Add a new analysis parameter dataset to the list of data to be written to the HDF5 file
#
#          :param dependency: The data that should be added to __dependency_list member to be written to HDF5.
#          :type data: dependency_dict
#
#          :raises: ValueError in case data is not of dependency_dict object.
#        """
#        if isinstance( dependency , dependency_dict ) :
#            self.__dependency_list.append( dependency )
#        else :
#             raise ValueError( "Invalid input for add_data_dependency function" )

    def write_analysis_data(self, analysis_group=None):
        """
        This function is used to write the actual analysis data to file. If not implemented, then the
        omsi_file_analysis API's default behavior is used instead.

        :param analysis_group: The h5py.Group object where the analysis is stored. May be None on cores that
            do not perform any writing but which need to participate in communication, e.g., to collect data
            for writing.

        """
        # The default implementation does roughly the following
        # from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        # for ana_data in self..get_all_analysis_data():
        #        omsi_file_analysis.__write_omsi_analysis_data__(analysis_group, ana_data)
        raise NotImplementedError

    def add_custom_data_to_omsi_file(self, analysis_group):
        """
        This function can be optionally overwritten to implement a custom data write
        function for the analysis to be used by the omsi_file API.

        Note, this function should be used only to add additional data to the analysis
        group. The data that is written by default is still written by
        the `omsi_file_experiment.create_analysis()` function, i.e., the following data is
        written by default: i) analysis_identifier ,ii) get_analysis_type, iii)__data_list,
        iv) parameters, v) runinfo . Since the `omsi_file.experiment.create_analysis()`
        functions takes care of setting up the basic structure of the analysis storage
        (included the subgroubs for storing parameters and data dependencies) this setup can generally
        be assumed to exist before this function is called. This function is called
        automatically at the end omsi_file.experiment.create_analysis() (i.e, actually
        `omsi_file_analysis.__populate_analysis__(..)` so that this function typically does not need to
        be called explicitly.

        :param analysis_group: The h5py.Group object where the analysis is stored.

        """
        pass

    def read_from_omsi_file(self,
                            analysis_object,
                            load_data=True,
                            load_parameters=True,
                            load_runtime_data=True,
                            dependencies_omsi_format=True,
                            ignore_type_conflict=False):
        """
        This function can be optionally overwritten to implement a custom data read.

        The default implementation tries to reconstruct the original data as far
        as possible, however, in particular in case that a custom add_custom_data_to_omsi_file
        function has been implemented, the default implementation may not be sufficient.
        The default implementation reconstructs: i) analysis_identifier and reads all
        custom data into iii)__data_list. Note, an error will be raised in case that
        the analysis type specified in the HDF5 file does not match the analysis type
        specified by get_analysis_type()

        :param analysis_object: The omsi_file_analysis object associated with the hdf5 data group with \
            the analysis data_list
        :param load_data: Should the analysis data be loaded from file (default) or just stored as h5py data objects
        :param load_parameters: Should parameters be loaded from file (default) or just stored as h5py data objects.
        :param load_runtime_data: Should runtime data be loaded from file (default) or just stored as h5py data objects
        :param dependencies_omsi_format: Should dependencies be loaded as omsi_file API objects (default)
            or just as h5py objects.
        :param ignore_type_conflict: Set to True to allow the analysis to be loaded into the
               current analysis object even if the type indicated in the file does not match the
               class. Default value is False. This behavior can be useful when different analysis
               have compatible data structures or when we want to load the data in to a generic
               analysis container, e.g, analysis_generic.

        :returns bool: Boolean indicating whether the data was read successfully

        :raise: TypeError : A type error will be raised in case that the analysis type specified \
        by the file does not match the analysis type provided by self.get_analysis_type()


        """
        log_helper.debug(__name__, "Restoring the analysis from file", root=self.mpi_root, comm=self.mpi_comm)
        if not ignore_type_conflict:
            if str(analysis_object.get_analysis_type()[0]) != str(self.get_analysis_type()):
                error_message = "The type of the analysis specified in the omsi data file " + \
                                "does not match the analysis type of the object " +\
                                str(analysis_object.get_analysis_type()[0]) + \
                                " != " + str(self.get_analysis_type())
                raise TypeError(error_message)

        identifier = analysis_object.get_analysis_identifier()
        if identifier is not None:
            self.analysis_identifier = identifier[0]
        else:
            warnings.warn("The analysis identifier could not be read from the omsi file")

        self.__data_list = analysis_object.get_all_analysis_data(load_data=load_data)
        parameter_list = analysis_object.get_all_parameter_data(load_data=load_parameters,
                                                                exclude_dependencies=False)
        parameters_values = {param['name']: param['data'] for param in parameter_list}
        self.clear_parameter_data()
        self.set_parameter_values(**parameters_values)
        self.run_info = analysis_object.get_all_runinfo_data(load_data=load_runtime_data)
        self.omsi_analysis_storage.append(analysis_object)
        self.update_analysis = False
        return True

    def get_analysis_identifier(self):
        """Return the name of the analysis used as key when searching for a particular analysis"""
        return self.analysis_identifier

    def set_analysis_identifier(self,
                                identifier):
        """
        Set the name of the analysis to identifer

        Side Effects: This function modifies self.analysis_identifier

        :param identifier: The new analysis identifier string to be used (should be unique)
        :type identifier: str

        """
        try:
            log_helper.debug(__name__, "Setting analysis identifier to: " +
                             unicode(identifier) + " for " + unicode(self))
        except:
            pass
        self.analysis_identifier = identifier

    def analysis_identifier_defined(self):
        """
        Check whether the analysis identifier is defined by the user, i.e., set to value different than undefined
        :return: bool
        """
        return self.get_analysis_identifier() != "undefined"


