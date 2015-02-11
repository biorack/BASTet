"""
Helper module with data structures for managing analysis-related data.
"""

import numpy as np


class omsi_analysis_data(dict):
    """Define a dataset for the analysis that should be written to the omsi HDF5 file
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
             - To generate data links: ana_hdf5link   (omsi_analysis_data)

        """
        super(omsi_analysis_data, self).__init__()
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
            raise KeyError("\'"+str(key)+'\' key not in default key set of omsi_analysis_data')

