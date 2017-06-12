"""
This module defines the basic format for storing mass spectrometry imaging data,
metadata, and analysis in HDF5 in compliance with OpenMSI file format.
"""
import h5py
import numpy as np


class omsi_format_common(object):

    """
    Specification of common attributes, and names for the file format.

    :var str_type: `str_type  = h5py.new_vlen(str)` : Datatype used for storing strings in hdf5
    :var type_attribute: Name of the optional type attribute indicating which omsi_file_* class
                        should be used to interact with a given group.
    """

    def __init__(self):
        super(omsi_format_common, self).__init__()

    try:
        str_type = h5py.special_dtype(vlen=unicode)
        str_type_unicode = True
    except NotImplementedError:
        str_type = h5py.special_dtype(vlen=str)
        str_type_unicode = False
    type_attribute = "omsi_type"
    timestamp_attribute = "timestamp"
    version_attribute = "version"
    current_version = "0.1"


class omsi_format_file(omsi_format_common):
    """
    Specification of main-file related specific name conventions
    """
    def __init__(self):
        super(omsi_format_file, self).__init__()
    current_version = "0.1"


class omsi_format_experiment(omsi_format_common):

    """
    Specification of file format specific name conventions

    :var exp_groupname: `entry_` : The base name for a group containing data about an experiment
    :var exp_identifier_name: `experiment_identifier` : The identifier dataset for an experiment
    """

    def __init__(self):
        super(omsi_format_experiment, self).__init__()

    exp_groupname = "entry_"
    exp_identifier_name = "experiment_identifier"
    current_version = "0.1"

class omsi_format_metadata_collection(omsi_format_common):
    """
    Specification of the basic format for a general-purpose metadata storage

    :var metadata_collection_groupname_default: Default name for the group where
        the collection of metadata is stored.
    :var description_value_attribute: The attribute to be associated with each
        metadata value describing the content in a human-readable form
    :var unit_value_attribute: The attribute to be associated with each metatdata
        value describing the unit
    :var ontology_value_attribute: Optional ontology associated with a metadata value
    :var is_json_dict_attribute: Name of the attribute used to define whether an object is json
    """
    def __init__(self):
        super(omsi_format_metadata_collection, self).__init__()

    metadata_collection_groupname_default = 'metadata'
    description_value_attribute = 'description'
    unit_value_attribute = 'unit'
    ontology_value_attribute = 'ontology'
    is_json_dict_attribute = 'isjson'
    current_version = "0.1"

class omsi_format_methods(omsi_format_metadata_collection):

    """
    Specification of the basic format for storing method-related information

    :var methods_groupname: `methods` : The group storing all the information about the method
    :var methods_old_groupname: `method` : The group object was refactored to methods. To ensure that old
        files can still be read, this variable was added and is checked as well if needed.
    :var methods_name: `name` : The dataset with the name of the method
    """

    def __init__(self):
        super(omsi_format_methods, self).__init__()

    methods_groupname = "methods"
    methods_old_groupname = "sample"
    methods_name = "name"
    current_version = "0.3"


class omsi_format_instrument(omsi_format_metadata_collection):

    """
    Specification for storing instrument related information

    :var instrument_groupname: `instrument`  : Group with information about the intrument used
    :var instrument_mz_name: `mz` : Name of the dataset for the intrument's mz data values
    :var instrument_name: `name` : The dataset with the name of the instrument
    """

    def __init__(self):
        super(omsi_format_instrument, self).__init__()

    instrument_groupname = "instrument"
    instrument_mz_name = "mz"
    instrument_name = "name"
    current_version = "0.2"


class omsi_format_analysis(omsi_format_common):

    """
    Specification for storing analysis related data.

    :var analysis_groupname: `analysis_` : Group with additional analysis results
    :var analysis_identifier: `analysis_identifier` : Identifier for the analysis to enable look-up by analysis id
    :var analysis_type: `analysis_type` : Dataset used to store the analysis type descriptor string
    :var analysis_parameter_group: Group for storing analysis parameters. Dependent parameters are stored
                                   separately using a omsi_format_dependencies group.
    :var analysis_runinfo_group: Group for storing run information, e.g., where was the analysis run, how long
                                 did it take etc.
    """

    def __init__(self):
        super(omsi_format_analysis, self).__init__()

    analysis_groupname = "analysis_"
    analysis_identifier = "analysis_identifier"
    analysis_type = "analysis_type"
    analysis_class = "analysis_class"
    analysis_class_pickle_np_dtype = np.dtype('uint8')
    analysis_parameter_group = "parameter"
    analysis_parameter_help_attr = 'help'
    analysis_runinfo_group = "runinfo"
    current_version = "0.2"


class omsi_format_dependencydata(omsi_format_common):

    """
    Specification for the storage of a single dependency.

    This type of group does not have specific name to allow the user to specify a specific link_name
    to ease retrieval of the data.

    :var dependency_parameter: `parameter_name` : Name of string dataset used to store the name of
                the dependent parameter
    :var dependency_selection: `selection` : Name of the string dataset used to store a selection string if needed.
    :var dependency_mainname:  `main_name` : Name of the string dataset used to store the description of the
                link to the object that this depends on.
    :var dependency_datasetname: 'data_name` : Name fo the string dataset used to store the name of dataset
                within the mainname highlevel object.
    """

    def __init__(self):
        super(omsi_format_dependencydata, self).__init__()

    dependency_parameter = "parameter_name"
    dependency_selection = "selection"
    dependency_mainname = "main_name"
    dependency_datasetname = "data_name"
    dependency_typename = 'dependency_type'
    dependency_parameter_help_attr = 'help'
    current_version = "0.3"


class omsi_format_dependencies(omsi_format_common):

    """
    Specification for the management of a collection of dependencies.

    :var dependencies_groupname: `dependency` : Name of the group the dependencies are stored in.
    """

    def __init__(self):
        super(omsi_format_dependencies, self).__init__()

    dependencies_groupname = "dependency"
    current_version = "0.1"


class omsi_format_data(omsi_format_common):

    """
    Specification for storing raw data information.

    :var data_groupname: The base name for the hdf5 group containing the imaging data
    :var dataset_name: The base name for storing raw data. In the case of MSI data, this is \
                        the 3D data cube stored as 3D (full_cube), 2D (partial_cube) or \
                        1D (partial_spectra) dataset, depending on the format_type.
    :var data_dependency_group: Optional group for storing data dependencies
    """

    def __init__(self):
        super(omsi_format_data, self).__init__()

    data_groupname = "data_"
    dataset_name = "data_"
    current_version = "0.1"


class omsi_format_msidata(omsi_format_data):
    """
    Specification of the basic format for storing an MSI dataset consisting of a complete 3D cube
    (or a 3D cube completed with 0s for missing data)

    :var format_types: Data layout types supported for storing MSI data.
    :var mzdata_name: Global mz axis for the MSI data cube.
    :var format: Dataset in HDF5 with the format_type descriptor.
    """

    def __init__(self):
        super(omsi_format_msidata, self).__init__()

    format_types = {'full_cube': 1, 'partial_cube': 2, 'partial_spectra': 3}
    mzdata_name = "mz"
    format_name = "format"
    current_version = "0.1"


class omsi_format_msidata_partial_cube(omsi_format_msidata):

    """
    Specification of the basic format for storing an MSI datasets that define a partial
    cube with full spectra

    :var xy_index_name: 2D dataset indicating for each spectrum its start location in the the main dataset
    :var inv_xy_index_name: 2D dataset with n rows and 2 columns indicating for each spectrum i the (x,y)
                    pixel index the spectrum belongs to.  This index is stored for convenience purposes but
                    is not actually needed for data access.
    :var shape_name: Simple [3] indicating the true image size in x,y,mz
    """

    def __init__(self):
        super(omsi_format_msidata_partial_cube, self).__init__()

    xy_index_name = "xyindex"
    inv_xy_index_name = "inv_xy_index"
    shape_name = "shape"


class omsi_format_msidata_partial_spectra(omsi_format_msidata_partial_cube):

    """
    Specification of the basic format for storing an MSI dataset of a full or
    partial cube with partial spectra

    :var mz_index_name: 1D dataset of the same size as the spectrum data, indicating the indices
                into the global m/z list
    :var xy_index_end_name: 2D dataset indicating for each spectrum its end location (index not included)
                in the the main dataset
    """

    def __init__(self):
        super(omsi_format_msidata_partial_spectra, self).__init__()

    mz_index_name = "mz_index"
    xy_index_end_name = "xyindexend"
