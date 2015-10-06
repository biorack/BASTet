"""
Global peak finder computing peaks and associated ion-images
for the full MSI data.
"""
from omsi.analysis.base import analysis_base
from omsi.shared.log import log_helper

class omsi_findpeaks_global(analysis_base):
    """
    Basic global peak detection analysis. The default implementation
    computes the peaks on the average spectrum and then computes the peak-cube data,
    i.e., the values for the detected peaks at each pixel.

    TODO: The current version assumes 2D data
    """

    def __init__(self,
                 name_key="undefined"):
        """Initialize the basic data members"""
        super(omsi_findpeaks_global, self).__init__()
        self.analysis_identifier = name_key
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI dataset to be analyzed',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='mzdata',
                           help='The m/z values for the spectra of the MSI dataset',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='integration_width',
                           help='The window over which peaks should be integrated',
                           dtype=float,
                           default=0.1,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='peakheight',
                           help='Peak height parameter',
                           dtype=int,
                           default=2,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='slwindow',
                           help='Sliding window parameter',
                           dtype=int,
                           default=100,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='smoothwidth',
                           help='Smooth width parameter',
                           dtype=int,
                           default=3,
                           group=groups['settings'],
                           required=True)
        self.data_names = ['peak_cube',
                           'peak_mz']

    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """Implement support for qslice URL requests for the viewer"""
        from omsi.shared.data_selection import selection_string_to_object
        if viewer_option == 0:
            dataset = analysis_object['peak_cube']
            try:
                z_select = selection_string_to_object(selection_string=z)
                data = dataset[:, :, z_select]
                return data
            except:
                log_helper.error(__name__, "Global peak selection failed. ")
                return None
        elif viewer_option >= 0:
            return super(omsi_findpeaks_global, cls).v_qslice(analysis_object, z, viewer_option-1)
        else:
            return None

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """Implement support for qspectrum URL requests for the viewer"""
        # Get the h5py dataset with the peak_cube data
        data = None
        custom_mz = None
        if viewer_option == 0:
            from omsi.shared.data_selection import check_selection_string, \
                selection_type, \
                selection_string_to_object
            x_select = selection_string_to_object(selection_string=x)
            y_select = selection_string_to_object(selection_string=x)
            dataset = analysis_object['peak_cube']
            if (check_selection_string(x) == selection_type['indexlist']) and \
                    (check_selection_string(y) == selection_type['indexlist']):
                # The peak-cube data is usually small enough. To handle the multiple list selection case
                # we here just load the full data cube and use numpy to do the sub-selection. Note, this
                # version would work for all selection types but we would like to avoid loading the
                # full data if we don't have to.
                data = dataset[:][x_select, y_select, :]
            else:
                data = dataset[x_select, y_select, :]
            # Return the spectra and indicate that no custom_mz data values (i.e. None) are needed
            return data, None
        elif viewer_option > 0:
            return super(omsi_findpeaks_global, cls).v_qspectrum(analysis_object,
                                                                 x,
                                                                 y,
                                                                 viewer_option-1)

        return data, custom_mz

    @classmethod
    def v_qmz(cls,
              analysis_object,
              qslice_viewer_option=0,
              qspectrum_viewer_option=0):
        """Implement support for qmz URL requests for the viewer"""
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None
        peak_cube_shape = analysis_object['peak_cube'].shape
        valuesX = range(0, peak_cube_shape[0])
        labelX = 'pixel index X'
        valuesY = range(0, peak_cube_shape[1])
        labelY = 'pixel index Y'
        if len(peak_cube_shape) > 3:
            valuesZ = range(0, peak_cube_shape[2])
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None

        # We do not need to handle the qslice_viewer_option separately here since there is only one option right now
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  # Loadings
            mz_spectra = analysis_object['peak_mz'][:]
            label_spectra = "m/z"
            mz_slice = None
            label_slice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_findpeaks_global, cls).v_qmz(analysis_object,
                                                        qslice_viewer_option-1,
                                                        qspectrum_viewer_option-1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option > 0:
            mz_spectra = analysis_object['peak_mz'][:]
            label_spectra = "m/z"
            temp_a, temp_b, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_findpeaks_global, cls).v_qmz(analysis_object,
                                                        qslice_viewer_option=qslice_viewer_option-1,
                                                        qspectrum_viewer_option=0)
            # NOTE: if qspectrum and qslice share the same axis, this call will not return the
            # copied data, i.e., we need to copy the qspectrum values to the qslice values.
            if mz_slice is None:
                mz_slice = temp_a
                label_slice = temp_b
        elif qspectrum_viewer_option > 0 and qslice_viewer_option == 0:
            mz_slice = analysis_object['peak_mz'][:]
            label_slice = "m/z"
            # Ignore the spatial coordinates. We have to use the shape of the qcube
            mz_spectra, label_spectra, temp_a, temp_b, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_findpeaks_global, cls).v_qmz(analysis_object,
                                                        qslice_viewer_option=0,
                                                        qspectrum_viewer_option=qspectrum_viewer_option-1)

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_global, cls).v_qspectrum_viewer_options(analysis_object)
        return ["Peak cube"] + dependent_options

    @classmethod
    def v_qslice_viewer_options(cls,
                                analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_global, cls).v_qslice_viewer_options(analysis_object)
        return ["Peak cube"] + dependent_options

    def execute_analysis(self):
        """
        Execute the global peak finding for the given msidata and mzdata.
        """
        # Make sure all imports are here
        from omsi.analysis.findpeaks.third_party.findpeaks import findpeaks
        import numpy as np

        # Copy parameters to local variables for convenience
        msidata = self['msidata']
        mzdata = self['mzdata'][:]  # Load the mz data into memory so that we can use numpy int arrays for slicing
        integration_width = self['integration_width']
        peakheight = self['peakheight']
        slwindow = self['slwindow']
        smoothwidth = self['smoothwidth']

        # Load the input data
        data = msidata[:]
        # Ensure the our MSI dataset has sufficient numbers of dimensions
        if len(msidata.shape) == 1:
            data = data[:][np.newaxis, np.newaxis, :]
        elif len(msidata.shape) == 2:
            data = data[:][np.newaxis, :]

        # Determine the data dimensions
        shape_x, shape_y, shape_z = data.shape
        # Linearize the spectral data
        processed_msidata = data.reshape(shape_x*shape_y, shape_z)
        # Compute the average spectrum
        processed_msidata = np.mean(processed_msidata, axis=0)
        # Find peaks in the average spectrum

        # REMOVED 20140429
        findpeaks_data = findpeaks(mzdata[:],
                                   processed_msidata,
                                   smoothwidth,
                                   slwindow,
                                   peakheight)
        processed_msidata = findpeaks_data.smoothListGaussian()

        # from the smoothed spectra subtract a sliding minima

        # REMOVED 20140429
        findpeaks_data = findpeaks(mzdata[:],
                                   processed_msidata,
                                   smoothwidth,
                                   slwindow,
                                   peakheight)
        slmin = [x for x in findpeaks_data.sliding_window_minimum()]
        processed_msidata = processed_msidata - slmin

        # find peaks in the smoothed, background subtracted spectra
        findpeaks_data = findpeaks(mzdata[:],
                                   processed_msidata,
                                   smoothwidth,
                                   slwindow,
                                   peakheight)
        [pkmax, pkmin] = findpeaks_data.peakdet()
        xval_peak = [x[0] for x in pkmax]
        yval_peak = [x[1] for x in pkmax]
        pks = np.asarray(pkmax)
        mz_peaks = mzdata[pks[:, 0].astype(int)]
        flat_data = data.reshape(shape_x*shape_y, shape_z)
        peak_cube = np.zeros((shape_x*shape_y, mz_peaks.shape[0]))
        for i in range(len(mz_peaks)):
            xx = np.where(np.abs(mzdata - mz_peaks[i]) < integration_width)
            peak_cube[:, i] = np.amax(flat_data[:, xx[0]], 1)
        peak_cube = peak_cube.reshape(shape_x, shape_y, len(mz_peaks))
        # integrate peaks +/- integration_width bins around each of the peaks found in the total spectra
        # TODO : THIS LOOP NEEDS TO BE CONVERTED TO A MAX INSTEAD OF A SUM
        # im = data[:,:,xp]
        # im = im*0.0;
        # xp = np.array(xp)
        # for i in range(-integration_width,integration_width):
        #     idx = xp+i
        #     try:
        #         im = im + data[:,:,idx]
        #     except IndexError:
        #         pass

        # Add the analysis results and parameters to the analysis data so that it can be accessed and written to file
        # We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        # handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        # to ensure a consistent behavior we convert the values directly here

        # Clean up the biggest chunks of memory we allocated when we potentially loaded the data from file
        del data
        del processed_msidata
        del flat_data

        # Save the analysis data to the __data_list so that the data can be
        # saved automatically by the omsi HDF5 file API
        return peak_cube, mz_peaks

if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_findpeaks_global).main()


# def main(argv=None):
#     """Then main function"""
#
#     import sys
#     from omsi.dataformat.omsi_file.main_file import omsi_file
#     import matplotlib.pyplot as plt
#     import matplotlib.gridspec as gridspec
#     import numpy as np
#
#     if argv is None:
#         argv = sys.argv
#
#     # Check for correct usage
#     if len(argv) < 2:
#         print "USAGE: Call \"omsi_findpeaks_global OMSI_FILE [expIndex data_index]   \" "
#         print "\n"
#         print "This is a simple test function to test the peak finding on a given omsi HDF5 file"
#         exit(0)
#
#     # Read the input arguments
#     omsi_input_filename = argv[1]
#     experiment_index = 0
#     data_index = 0
#     if len(argv) == 4:
#         experiment_index = int(argv[2])
#         data_index = int(argv[3])
#
#     # Open the input HDF5 file
#     omsi_input_file = omsi_file(omsi_input_filename, 'r')  # Open file in read only mode
#
#     # Get the experiment and dataset
#     experiment = omsi_input_file.get_experiment(experiment_index)
#     msidata = experiment.get_msidata(data_index)
#     mzdata = experiment.get_instrument_info().get_instrument_mz()
#
#     # Execute the peak finding
#     ofpg_object = omsi_findpeaks_global()
#     print "Executing peakfinding analysis"
#     ofpg_object.execute(msidata=msidata, mzdata=mzdata)
#     print "Getting peak finding analysis results"
#     peak_cube = ofpg_object['peak_cube']['data']
#     print peak_cube
#
#     # Plot the first three peak images
#     print "Plotting example peak finding analysis results"
#     shape_x = peak_cube.shape[0]
#     shape_y = peak_cube.shape[1]
#     ho1 = peak_cube[:, :, 0]
#     ho2 = peak_cube[:, :, 0]
#     ho3 = peak_cube[:, :, 0]
#
#     main_figure = plt.figure()
#     figure_grid_spec = gridspec.GridSpec(1, 4)
#     image_figure = main_figure.add_subplot(figure_grid_spec[0])
#     image_figure.autoscale(True, 'both', tight=True)
#     _ = image_figure.pcolor(np.log(ho1 + 1))
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[1])
#     image_figure.autoscale(True, 'both', tight=True)
#     _ = image_figure.pcolor(np.log(ho2 + 1))
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[2])
#     image_figure.autoscale(True, 'both', tight=True)
#     _ = image_figure.pcolor(np.log(ho3 + 1))
#
#     # do the three color
#     ho = np.zeros(shape=(shape_x, shape_y, 3))
#     temp = np.log(peak_cube[:, :, 0] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho[:, :, 0] = temp
#
#     temp = np.log(peak_cube[:, :, 1] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho[:, :, 1] = temp
#
#     temp = np.log(peak_cube[:, :, 2] + 1)
#     temp = temp - temp.min()
#     temp = temp / temp.max()
#     ho[:, :, 2] = temp
#
#     image_figure = main_figure.add_subplot(figure_grid_spec[3])
#     image_figure.autoscale(True, 'both', tight=True)
#     _ = image_figure.imshow(ho)
#
#     plt.show()

