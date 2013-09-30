"""This module defines the basic format for storing mass spectrometry imaging data, 
   metadata, and analysis in HDF5 in compliance with OpenMSI file format.""" 
import h5py

class omsi_format_common:
    """Specification of common attributes, and names for the file format.
    
       :var str_type: `str_type  = h5py.new_vlen(str)` : Datatype used for storing strings in hdf5
       :var type_attribute: Name of the optional type attribute indicating which omsi_file_* class should be used to interact with a given group.
    """
    
    str_type = h5py.special_dtype(vlen=str) 
    type_attribute = "omsi_type"
    timestamp_attribute = "timestamp"
    version_attribute = "version"
    current_version = "0.1"

class omsi_format_experiment( omsi_format_common ) :
    """Specification of file format specific name conventions

       :var exp_groupname: `entry_` : The base name for a group containing data about an experiment 
       :var exp_identifier_name: `experiment_identifier` : The identifier dataset for an experiment
    """
    exp_groupname = "entry_"
    exp_identifier_name = "experiment_identifier"
    current_version = "0.1"

class omsi_format_sample( omsi_format_common ):
    """Specification of the basic format for storing sample-related information
    
       :var sample_groupname: `sample` : The group storing all the information about the sample
       :var sample_name: `name` : The dataset with the name of the sample
    """
    sample_groupname = "sample"
    sample_name = "name"
    current_version = "0.1"

class omsi_format_instrument( omsi_format_common ) : 
    """ Specification for storing instrument related information
    
       :var instrument_groupname: `instrument`  : Group with information about the intrument used
       :var instrument_mz_name: `mz` : Name of the dataset for the intrument's mz data values
       :var instrument_name: `name` : The dataset with the name of the instrument
    """

    instrument_groupname = "instrument"
    instrument_mz_name = "mz"
    instrument_name = "name"
    current_version = "0.1"
    
class omsi_format_analysis( omsi_format_common ) :
    """ Specification for storing analysis related data.
    
       :var analysis_groupname: `analysis_` : Group with additional analysis results
       :var analysis_identifier: `analysis_identifier` : Identifier for the analysis to enable look-up by analysis id
       :var analysis_type: `analysis_type` : Dataset used to store the analysis type descriptor string
    """
    analysis_groupname = "analysis_"
    analysis_identifier = "analysis_identifier"
    analysis_type = "analysis_type"
    analysis_parameter_group = "parameter"
    analysis_dependency_group = "dependency"
    current_version = "0.1"

class omsi_format_dependencydata( omsi_format_common ):
    
    dependency_parameter = "parameter_name"
    dependency_selection = "selection"
    dependency_mainname  = "main_name"
    current_version = "0.1"
    
    
class omsi_format_data( omsi_format_common ) :
    """ Specification for storing raw data information. 
    
        :var data_groupname: The base name for the hdf5 group containing the imaging data
        :var dataset_name: The base name for storing raw data. In the case of MSI data, this is \ 
                            the 3D data cube stored as 3D (full_cube), 2D (partial_cube) or \
                            1D (partial_spectra) dataset, depending on the format_type.
        :var data_dependency_group: Optional group for storing data dependcies 
    """
    data_groupname = "data_"  
    dataset_name= "data_"
    current_version = "0.1"
    data_dependency_group = "dependency"

class omsi_format_msidata( omsi_format_data ):
    """Specification of the basic format for storing an MSI dataset consisting of a complete 3D cube 
       (or a 3D cube completed with 0s for missing data)
    
       :var format_types: Data layout types supported for storing MSI data.
       :var mzdata_name: Global mz axis for the MSI data cube.
       :var format: Dataset in HDF5 with the format_type descriptor.
    """
    format_types = { 'full_cube':1, 'partial_cube':2  , 'partial_spectra':3 }
    mzdata_name = "mz"
    format_name = "format"
    current_version = "0.1"
    
class omsi_format_msidata_partial_cube( omsi_format_msidata ) :
    """Specification of the basic format for storing an MSI datasets that define a partial 
       cube with full spectra
       
       :var xy_index_name: 2D dataset indicating for each spectrum its start location in the the main dataset
       :var inv_xy_index_name: 2D dataset with n rows and 2 columns indicating for each spectrum i the (x,y) pixel index the spectrum belongs to. \
                          This index is stored for convenience purposesbut is not actually needed for data access.
       :var shape_name: Simple [3] indicating the true image size in x,y,mz
    """
    xy_index_name= "xyindex"
    inv_xy_index_name= "inv_xy_index"
    shape_name   = "shape"
    
class omsi_format_msidata_partial_spectra( omsi_format_msidata_partial_cube ) :
   """Specification of the basic format for storing an MSI dataset of a full or 
      partial cube with partial spectra
      
     :var mz_index_name: 1D dataset of the same size as the spectrum data, indicating the indicies into the global m/z list
     :var xy_index_end_name: 2D dataset indicating for each spectrum its end location (index not included) in the the main dataset
   """
   mz_index_name = "mz_index"
   xy_index_end_name= "xyindexend"
