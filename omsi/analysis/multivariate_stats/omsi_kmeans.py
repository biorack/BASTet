"""
Module for performing kmeans clustering for MSI data.
"""
import numpy as np

from omsi.analysis.base import analysis_base


class omsi_kmeans(analysis_base):
    """Class defining a basic nmf analysis for a 2D MSI data file or slice of the data"""

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""
        super(omsi_kmeans, self).__init__()
        self.analysis_identifier = name_key
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI matrix to be analyzed',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='numClusters',
                           help='The number of output clusters to compute.',
                           default=20,
                           dtype=dtypes['int'],
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='numRepeat',
                           help='The number of times to repeat k-means.',
                           dtype=dtypes['int'],
                           default=20,
                           group=groups['stop'])
        self.add_parameter(name='threshold',
                           help='Terminates the k-means algorithm if the change in distortion ' +
                                'since the last k-means iteration is less than thresh.',
                           dtype=dtypes['float'],
                           default=1.0e-5,
                           group=groups['stop'])
        self.add_parameter(name='normalize',
                           help='Normalize the data by dividing each spectrum by the amax',
                           dtype=dtypes['bool'],
                           default=False,
                           group=groups['settings'])
        self.add_parameter(name='clusterImages',
                           help='Transpose the matrix to cluster images instead of spectra',
                           dtype=dtypes['bool'],
                           default=False,
                           group=groups['settings'])
        self.add_parameter(name='mask',
                           help='A boolean mask indicating the elements to be clustered',
                           dtype=dtypes['ndarray'],
                           default=None,
                           required=False,
                           group=groups['settings'])
        self.data_names = ['clusters', 'centers']

    def execute_analysis(self):
        """
        Execute the kmeans clustering for the given msidata
        """
        from scipy.cluster.vq import kmeans, vq

        # Assign parameters to local variables for convenience
        current_msidata = self['msidata']
        current_num_clusters = self['numClusters']
        current_num_repeat = self['numRepeat']
        current_threshold = self['threshold']
        normalize_vals = self['normalize']
        cluster_images = self['clusterImages']
        current_mask = self['mask']

        # Copy the input data
        data = current_msidata[:]
        indata_shape = data.shape
        # Set masked data to a single value
        if current_mask is not None:
            if not cluster_images:
                data = data[current_mask, :]
            else:
                data = data[...,current_mask]


        # Normalize the data values to a range 0,1
        if normalize_vals:
            temp = np.amax(data, axis=-1) + 1  # The last data dimension is assumed to contain the spectra
            data = np.divide(data, temp[:, :, np.newaxis])

        # Determine the input shape
        num_bins = data.shape[-1]
        num_pixels = data.size / num_bins  # The last data dimension is assumed to contain the spectra

        # Reshape the data
        data = data.reshape(num_pixels, num_bins)
        if cluster_images:
            data = np.transpose(data)

        # Execute clustering
        centers, _ = kmeans(data,
                            current_num_clusters,
                            iter=current_num_repeat,
                            thresh=current_threshold)
        cluster, _ = vq(data, centers)

        # Fill the clusters
        if current_mask is not None:
            out_cluster = np.zeros(current_mask.shape, dtype='int')
            out_cluster[:] = -1
            out_cluster[current_mask] = cluster
            cluster = out_cluster
        else:
            # Reshape the clusters to the original shape
            if cluster_images:
                pass  # Nothing to be done
            else:
                cluster = cluster.reshape(indata_shape[:-1])

        return cluster, centers


if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_kmeans).main()
