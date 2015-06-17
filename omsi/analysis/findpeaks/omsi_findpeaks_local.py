"""
Local peak finding analysis module.
"""

from omsi.analysis.base import analysis_base


class omsi_findpeaks_local(analysis_base):
    """
    Class defining a basic gloabl peak finding. The default implementation computes the peaks on the average
    spectrum and then computes the peak-cube data, i.e., the values for the detected peaks at each pixel.

    TODO: The current version assumes 2D data
    """

    def __init__(self, name_key="undefined"):
        """
        Initialize the basic data members
        """
        super(omsi_findpeaks_local, self).__init__()
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
                           default=10,
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
        self.add_parameter(name='printStatus',
                           help='Print progress status during the analysis',
                           dtype=bool,
                           required=True,
                           group=groups['settings'],
                           default=False)
        self.data_names = ['peak_mz',
                           'peak_value',
                           'peak_arrayindex',
                           'indata_mz']

    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """Implement support for qslice URL requests for the viewer"""
        # Use the dependency data for slicing here. We do not have a native option to reconstruct
        # images from local peak finding data
        return super(omsi_findpeaks_local, cls).v_qslice(analysis_object,
                                                         z,
                                                         viewer_option)

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """Implement support for qspectrum URL requests for the viewer"""
        # Retrieve the h5py objects for the requried datasets from the local peak finding
        if viewer_option == 0:

            from omsi.shared.data_selection import check_selection_string, selection_type, selection_to_indexlist
            import numpy as np
            peak_mz = analysis_object['peak_mz']
            peak_values = analysis_object['peak_value']
            array_indices = analysis_object['peak_arrayindex'][:]
            indata_mz = analysis_object['indata_mz']
            # Determine the shape of the original raw data
            if (indata_mz is None) or (array_indices is None):
                return None, None
            num_x = array_indices[:, 0].max()
            num_y = array_indices[:, 1].max()
            num_mz = indata_mz.shape[0]
            num_spectra = array_indices.shape[0]
            # Determine the size of the selection and the set of selected items
            x_list = selection_to_indexlist(x, num_x)
            y_list = selection_to_indexlist(y, num_y)
            if (check_selection_string(x) == selection_type['indexlist']) and \
                    (check_selection_string(y) == selection_type['indexlist']):
                if len(x_list) == len(y_list):
                    items = [(x_list[i], y_list[i]) for i in xrange(0, len(x_list))]
                else:
                    return None, None
            else:
                items = [0]*(len(x_list)*len(y_list))
                index = 0
                for xi in x_list:
                    for yi in y_list:
                        items[index] = (xi, yi)
                        index += 1

            shape_x = len(items)
            shape_y = 1
            shape_z = num_mz
            # Initalize the data cube to be returned
            data = np.zeros((shape_x, shape_y, shape_z), dtype=peak_values.dtype)
            # Fill the non-zero locations for the data cube with data
            for ni in xrange(0, len(items)):
                current_index = (items[ni][0]*num_y + items[ni][1])
                current_dx = ni
                current_dy = 0
                start_index = array_indices[current_index][2]
                if current_index < num_spectra:
                    end_index = array_indices[(current_index+1)][2]
                else:
                    end_index = peak_values.size
                if start_index != end_index:
                    temp_values = peak_values[start_index: end_index]
                    temp_mz = peak_mz[start_index: end_index]
                    data[current_dx, current_dy, temp_mz] = temp_values
                else:
                    # The start and end index may be the same in case that
                    # no peaks for found for the given spectrum
                    # The data is already initialized to 0 so there is nothing to do here
                    pass

            if len(items) == 1:
                data = data.reshape((shape_x, shape_z))

            # Return the spectra and indicate that no customMZ data values (i.e. None) are needed
            return data, None

        elif viewer_option > 0:
            return super(omsi_findpeaks_local, cls).v_qspectrum(analysis_object,
                                                                x,
                                                                y,
                                                                viewer_option-1)
        else:
            return None, None

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
        # We do not have native option for qslice, so we rely on the input data in all cases
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  # Loadings
            mz_spectra = analysis_object['indata_mz'][:]
            label_spectra = "m/z"
            mz_slice = None
            label_slice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mz_spectra, label_spectra, mz_slice, label_slice = \
                super(omsi_findpeaks_local, cls).v_qmz(analysis_object,
                                                       qslice_viewer_option=qslice_viewer_option,
                                                       qspectrum_viewer_option=qspectrum_viewer_option-1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option >= 0:
            mz_spectra = analysis_object['indata_mz'][:]
            label_spectra = "m/z"
            temp_a, temp_b, mz_slice, label_slice = \
                super(omsi_findpeaks_local, cls).v_qmz(analysis_object,
                                                       0,
                                                       qspectrum_viewer_option)
            # NOTE: if qspectrum and qslice share the same axis, this call will
            # not return the copied data, i.e., we need to copy the
            # qspectrum values to the qslice values.
            if mz_slice is None:
                mz_slice = temp_a
                label_slice = temp_b

        return mz_spectra, label_spectra, mz_slice, label_slice

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_local, cls).v_qspectrum_viewer_options(analysis_object)
        re = ["Local peaks"] + dependent_options
        return re

    @classmethod
    def v_qslice_viewer_options(cls,
                                analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        return super(omsi_findpeaks_local, cls).v_qslice_viewer_options(analysis_object)

    def execute_analysis(self, msidata_subblock=None):
        """
        Execute the local peak finder for the given msidata.

        :param msidata_subblock: Optional input parameter used for parallel execution of the
            analysis only. If msidata_subblock is set, then the given subblock will be processed
            in SERIAL instead of processing self['msidata'] in PARALLEL (if available). This
            parameter is strictly optional and intended for internal use only.

        """
        # Make sure needed imports are available
        from omsi.analysis.findpeaks.third_party.findpeaks import findpeaks
        import omsi.shared.mpi_helper as mpi_helper
        import numpy as np

        # Assign parameters to local variables for convenience
        msidata =  self['msidata']
        if msidata_subblock is not None:
            msidata = msidata_subblock
        mzdata = self['mzdata']
        integration_width = self['integration_width']
        peakheight = self['peakheight']
        slwindow = self['slwindow']
        smoothwidth = self['smoothwidth']
        print_status = self['printStatus']
        if print_status:
            import sys

        #############################################################
        # Parallel execution using MPI
        #############################################################
        import omsi.shared.mpi_helper as mpi_helper
        if mpi_helper.get_size() > 1 and len(self['msidata'].shape) > 1:
            if msidata_subblock is None:
                split_axis = range(len(self['msidata'].shape)-1)
                # TODO Implement the data collection to also allow the use of the DYNAMIC scheduler
                scheduler = mpi_helper.parallel_over_axes(task_function=self.execute_analysis,
                                                          task_function_params={},
                                                          main_data=self['msidata'],
                                                          split_axes=split_axis,
                                                          main_data_param_name='msidata_subblock',
                                                          root=self.mpi_root,
                                                          schedule=mpi_helper.parallel_over_axes.SCHEDULES['STATIC'],
                                                          collect_output=True,
                                                          comm=self.mpi_comm)
                result = scheduler.run()

                # Return the output result as is if we are a "worker" rank
                if mpi_helper.get_rank() != self.mpi_root:
                    return result[0]
                else:
                    # Compile the results from all ranks if we are on the "main" rank
                    peak_mz = np.concatenate(tuple([ri[0] for ri in result[0]]), axis=-1)
                    peak_values = np.concatenate(tuple([ri[1] for ri in result[0]]), axis=-1)
                    peak_arrayindex = np.concatenate(tuple([ri[2] for ri in result[0]]), axis=0)
                    mzdata = result[0][0][3]
                    return peak_mz, peak_values, peak_arrayindex, mzdata


        #############################################################
        # Serial processing of the current data block
        #############################################################
        # Ensure the our MSI dataset has sufficient numbers of dimensions
        if len(msidata.shape) == 1:
            msidata = msidata[:][np.newaxis, np.newaxis, :]
        elif len(msidata.shape) == 2:
            msidata = msidata[:][np.newaxis, :]

        # Determine the data dimensions
        shape_x = msidata.shape[0]
        shape_y = msidata.shape[1]

        peak_mz = []      # The x values for all peaks, stored in a linear array
        peak_values = []    # The y values for all peaks, stored in a linear array
        # List describing for each pixel the start index where its peaks
        # are stored in the peaks_MZ and peaks_values array
        peak_arrayindex = np.zeros(shape=(shape_x*shape_y, 3), dtype='int64')
        current_index = long(0)
        pixel_index = 0
        for xi in xrange(0, shape_x):
            for yi in xrange(0, shape_y):
                if print_status:
                    sys.stdout.write("[" + str(int(100. * float(pixel_index)/float(shape_x*shape_y))) + "%]" + "\r")
                    sys.stdout.flush()

                # Load the spectrum
                y = msidata[xi, yi, :]
                # Find peaks in the spectrum
                peak_finder = findpeaks(mzdata[:],
                                        y,
                                        smoothwidth,
                                        slwindow,
                                        peakheight)
                y = peak_finder.smoothListGaussian()
                # from the smoothed spectra subtract a sliding minima
                peak_finder = findpeaks(mzdata[:],
                                        y,
                                        smoothwidth,
                                        slwindow,
                                        peakheight)
                slmin = [x for x in peak_finder.sliding_window_minimum()]
                y = y - slmin
                # find peaks in the smoothed, background subtracted spectra
                peak_finder = findpeaks(mzdata[:],
                                        y,
                                        smoothwidth,
                                        slwindow,
                                        peakheight)
                [pkmax, pkmin] = peak_finder.peakdet()
                xp = [x[0] for x in pkmax]
                yp = [x[1] for x in pkmax]
                peak_mz = peak_mz + xp
                peak_values = peak_values + yp
                peak_arrayindex[pixel_index, 0] = xi
                peak_arrayindex[pixel_index, 1] = yi
                peak_arrayindex[pixel_index, 2] = current_index
                pixel_index += 1
                current_index += len(yp)

        # Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
        # We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        # handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        # to ensure a consitent behavior we convert the values directly here

        # Save the analysis data to the __data_list so that the data can be
        # saved automatically by the omsi HDF5 file API
        return np.asarray(peak_mz), np.asarray(peak_values), peak_arrayindex, mzdata[:]

if __name__ == "__main__":
    from omsi.workflow.analysis_driver.omsi_cl_driver import omsi_cl_driver
    omsi_cl_driver(analysis_class=omsi_findpeaks_local).main()


# def main(argv=None):
#     """Then main function"""
#
#     import sys
#     from omsi.dataformat.omsi_file.main_file import omsi_file
#
#     if argv is None:
#         argv = sys.argv
#
#     # Check for correct usage
#     if len(argv) < 2:
#         print "USAGE: Call \"omsi_findpeaks_global OMSI_FILE [expperiment_index data_index]   \" "
#         print "\n"
#         print "This is a simple test function to test the peak finding on a given omsi HDF5 file"
#         exit(0)
#
#     # Read the input arguments
#     omsi_input_filename = argv[1]
#     expperiment_index = 0
#     data_index = 0
#     if len(argv) == 4:
#         expperiment_index = int(argv[2])
#         data_index = int(argv[3])
#
#     # Open the input HDF5 file
#     omsi_input_file = omsi_file(omsi_input_filename, 'r')  # Open file in read only mode
#
#     # Get the experiment and dataset
#     exp = omsi_input_file.get_experiment(expperiment_index)
#     data = exp.get_msidata(data_index)
#     mzdata = exp.get_instrument_info().get_instrument_mz()
#
#     # Execute the peak finding
#     test_omsi_fpl = omsi_findpeaks_local()
#     print "Executing peakfinding analysis"
#     test_omsi_fpl.execute(msidata=data, mzdata=mzdata)
#     print "Getting peak finding analysis results"
#     pmz = test_omsi_fpl['peak_mz']['data']
#     print pmz
#     pv = test_omsi_fpl['peak_value']['data']
#     print pv
#     pai = test_omsi_fpl['peak_arrayindex']['data']
#     print pai
#
#
# if __name__ == "__main__":
#     main()

