"""
Module for performing non-negative matrix factorization (NMF) for MSI data.
"""
import numpy as np

from omsi.analysis.base import analysis_base


class omsi_nmf(analysis_base):
    """
    Class defining a basic nmf analysis.

    The function has primarily been tested we MSI datasets but should support
    arbitrary n-D arrays (n>=2). The last dimension of the input array must be the
    spectrum dimnensions.
    """

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""
        super(omsi_nmf, self).__init__()
        self.analysis_identifier = name_key
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The data matrix to be analyzed.',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='numComponents',
                           help='The number of output components to compute.',
                           default=20,
                           dtype=int,
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='timeOut',
                           help='Time in seconds after which to terminate.',
                           default=600,
                           dtype=int,
                           required=True,
                           group=groups['stop'])
        self.add_parameter(name='numIter',
                           help='Maximum number of iterations.',
                           dtype=int,
                           default=2000,
                           required=True,
                           group=groups['stop'])
        self.add_parameter(name='tolerance',
                           help='The tolerance for a relative stopping condition.',
                           dtype=float,
                           default=0.0001,
                           required=True,
                           group=groups['stop'])
        self.add_parameter(name='mask',
                           help='A n-1 dimensional boolean mask indicating the elements to be clustered',
                           dtype=dtypes['ndarray'],
                           default=None,
                           required=False,
                           group=groups['settings'])
        self.data_names = ['wo', 'ho']

    @classmethod
    def v_qslice(cls, analysis_object, z, viewer_option=0):
        """Implement support for qslice URL requests for the viewer"""
        from omsi.shared.data_selection import selection_string_to_object
        if viewer_option == 0:
            dataset = analysis_object['ho']
            z_select = selection_string_to_object(selection_string=z)
            data = dataset[:, :, z_select]
            return data
        elif viewer_option > 0:
            return super(omsi_nmf, cls).v_qslice(analysis_object, z, viewer_option-1)

    @classmethod
    def v_qspectrum(cls, analysis_object, x, y, viewer_option=0):
        """Implement support for qspectrum URL requests for the viewer"""
        from omsi.shared.data_selection import selection_string_to_object
        data = None
        custom_mz = None
        if viewer_option == 0:  # Loadings
            dataset = analysis_object['ho']
            x_select = selection_string_to_object(selection_string=x)
            y_select = selection_string_to_object(selection_string=y)
            if isinstance(x_select, list) and isinstance(y_select, list):
                # Load the full data if multiple lists are given for
                # the selection and let numpy handle the sub-selection
                data = dataset[:][x_select, y_select, :]
            else:
                data = dataset[x_select, y_select, :]
            custom_mz = None
        elif viewer_option > 0:
            return super(omsi_nmf, cls).v_qspectrum(analysis_object, x, y, viewer_option-1)

        return data, custom_mz

    @classmethod
    def v_qmz(cls, analysis_object, qslice_viewer_option=0, qspectrum_viewer_option=0):
        """Implement support for qmz URL requests for the viewer"""
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None
        ho_cube_shape = analysis_object['ho'].shape
        valuesX = range(0, ho_cube_shape[0])
        labelX = 'pixel index X'
        valuesY = range(0, ho_cube_shape[1])
        labelY = 'pixel index Y'
        if len(ho_cube_shape) > 3:
            valuesZ = range(0, ho_cube_shape[2])
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None

        # We do not need to handle the qslice_viewer_option separately here since there is only one option right now
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  # Loadings
            mz_spectra = np.arange(0, analysis_object['ho'].shape[2])
            label_spectra = "Component Index"
            mz_slice = None
            label_slice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_nmf, cls).v_qmz(analysis_object,
                                           qslice_viewer_option=qslice_viewer_option-1,
                                           qspectrum_viewer_option=qspectrum_viewer_option-1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option > 0:
            mz_spectra = np.arange(0, analysis_object['ho'].shape[2])
            label_spectra = "Component Index"
            temp_a, temp_b, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_nmf, cls).v_qmz(analysis_object,
                                           qslice_viewer_option=qslice_viewer_option-1,
                                           qspectrum_viewer_option=0)
            # NOTE: if qspectrum and qslice share the same axis, this call will not return the
            # copied data, i.e., we need to copy the qspectrum values to the qslice values.
            if mz_slice is None:
                mz_slice = temp_a
                label_slice = temp_b
        elif qspectrum_viewer_option > 0 and qslice_viewer_option == 0:
            mz_slice = np.arange(0, analysis_object['ho'].shape[2])
            label_slice = "Component Index"
            # Ignore the spatial components as we are implementing our own spatial slicing here
            mz_spectra, label_spectra, temp_a, temp_b, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_nmf, cls).v_qmz(analysis_object,
                                           qslice_viewer_option=0,
                                           qspectrum_viewer_option=qspectrum_viewer_option-1)

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_nmf, cls).v_qspectrum_viewer_options(analysis_object)
        spectrum_viewer_options = ["NMF Loadings"] + dependent_options
        return spectrum_viewer_options

    @classmethod
    def v_qslice_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_nmf, cls).v_qslice_viewer_options(analysis_object)
        slice_viewer_options = ["NMF Images"] + dependent_options
        return slice_viewer_options

    def execute_analysis(self):
        """
        Execute the nmf for the given msidata
        """
        from omsi.analysis.multivariate_stats.third_party.nmf import nmf

        # Assign parameters to local variables for convenience
        current_msidata = self['msidata']
        current_num_components = self['numComponents']
        current_time_out = self['timeOut']
        current_num_iter = self['numIter']
        current_tolerance = self['tolerance']
        current_mask = self['mask']

        # Copy the input data
        data = current_msidata[:]
        indata_shape = data.shape
        output_shape = (indata_shape[:-1] + (current_num_components,))

        # Mask the data if requested
        if current_mask is not None:
            current_mask = current_mask[:]
            data = data[current_mask, :]

        # Determine the input shape after masking
        num_bins = data.shape[-1]
        num_pixels = data.size / num_bins  # The last data dimension is assumed to contain the spectra

        # Reshape the data
        data = data.reshape(num_pixels, num_bins)
        data = data.transpose()

        # Execute nmf
        (wo_matrix, ho_matrix) = nmf(data,
                                     np.random.randn(num_bins, current_num_components),
                                     np.random.randn(current_num_components, num_pixels),
                                     current_tolerance,
                                     current_time_out,
                                     current_num_iter)

        # Reshape the ho matrix to be a 3D image cube
        ho_matrix = ho_matrix.transpose()
        if current_mask is not None:
            out_ho = np.zeros(output_shape, dtype=ho_matrix.dtype)
            out_ho[current_mask, :] = ho_matrix
            ho_matrix = out_ho
        else:
            ho_matrix = ho_matrix.reshape(output_shape)

        return wo_matrix, ho_matrix

if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_nmf).main()


# def main(argv=None):
#     """Then main function"""
#
#     import sys
#     from sys import argv, exit
#     import matplotlib.pyplot as plt
#     import matplotlib.gridspec as gridspec
#     from math import log
#     from omsi.dataformat.omsi_file import  omsi_file
#
#     if argv is None:
#         argv = sys.argv
#
#     # Check for correct usage
#     if len(argv) < 2:
#         print "USAGE: Call \"omsi_nmf [expIndex dataIndex]   \" "
#         print "\n"
#         print "This is a simple viewer for looking at OMSI HDF5 files."
#         print "The viewer takes the index of the experiment and the"
#         print "dataset to be used as optional input"
#         exit(0)
#
#     # Read the input arguments
#     omsi_output_file = argv[1]
#     experiment_index = 0
#     data_index = 0
#     if len(argv) == 4:
#         experiment_index = int(argv[2])
#         data_index = int(argv[3])
#
#     # Open the input HDF5 file
#     omsi_input_file = omsi_file(omsi_output_file, 'r')  # Open file in read only mode
#
#     # Get the experiment and dataset
#     exp = omsi_input_file.get_experiment(experiment_index)
#     data = exp.get_msidata(data_index)
#
#     # Execute the nmf
#     test_nmf = omsi_nmf()
#     print "Executing nmf analysis"
#     test_nmf.execute(msidata=data)
#     print "Getting nmf analysis results"
#     wo_dataset = test_nmf.get_analysis_data('wo')['data']
#     ho_dataset = test_nmf.get_analysis_data('ho')['data']
#     print ho_dataset
#     print "Plotting nmf analysis results"
#     shape_x = data.shape[0]
#     shape_y = data.shape[1]
#     ho1 = ho_dataset[0, :]
#     ho1 = np.ravel(ho1)
#     ho1 = ho1.reshape(shape_x, shape_y)
#
#     ho2 = ho_dataset[1, :]
#     ho2 = np.ravel(ho2)
#     ho2 = ho2.reshape(shape_x, shape_y)
#
#     ho3 = ho_dataset[2, :]
#     ho3 = np.ravel(ho3)
#     ho3 = ho3.reshape(shape_x, shape_y)
#
#     main_figure = plt.figure()
#     figure_grid_spec = gridspec.GridSpec(1, 4)
#     image_figure = main_figure.add_subplot(figure_grid_spec[0])
# #    image_figure.autoscale(True,'both',tight=True)
#     image_plot = image_figure.pcolor(log(ho1 + 1))
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[1])
# #    image_figure.autoscale(True,'both',tight=True)
#     image_plot = image_figure.pcolor(log(ho2 + 1))
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[2])
# #    image_figure.autoscale(True,'both',tight=True)
#     image_plot = image_figure.pcolor(log(ho3 + 1))
#
#
# #    do the three color
#     ho_dataset = ho_dataset.transpose()
#     ho_dataset = ho_dataset.reshape(shape_x, shape_y, 3)
#     temp = log(ho_dataset[:, :, 0] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho_dataset[:, :, 0] = temp
#
#     temp = log(ho_dataset[:, :, 1] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho_dataset[:, :, 1] = temp
#
#     temp = log(ho_dataset[:, :, 2] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho_dataset[:, :, 2] = temp
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[3])
#     image_figure.autoscale(True, 'both', tight=True)
#     image_plot = image_figure.imshow(ho_dataset)
#
#     plt.show()
#
#
# if __name__ == "__main__":
#     main()

