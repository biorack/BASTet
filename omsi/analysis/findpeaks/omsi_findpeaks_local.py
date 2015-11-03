"""
Local peak finding analysis module.
"""

from omsi.analysis.base import analysis_base
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.log import log_helper

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
                           dtype=dtypes['bool'],
                           required=True,
                           group=groups['settings'],
                           default=False)
        self.add_parameter(name='schedule',
                           help='Scheduling to be used for parallel MPI runs',
                           dtype=str,
                           required=False,
                           choices=mpi_helper.parallel_over_axes.SCHEDULES.values(),
                           group=groups['parallel'],
                           default=mpi_helper.parallel_over_axes.SCHEDULES['DYNAMIC'])
        self.add_parameter(name='collect',
                           help='Collect results to the MPI root rank when running in parallel',
                           dtype=dtypes['bool'],
                           required=False,
                           group=groups['parallel'],
                           default=True)
        self.data_names = ['peak_mz',
                           'peak_value',
                           'peak_arrayindex',
                           'indata_mz']
        # TODO Allow the precursor_mz to be stored when doing local peak finding on MS2 data. This should be a 1D float array called 'precursor_mz' with the precursor m/z value for each spectrum

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
            # Initialize the data cube to be returned
            data = np.zeros((shape_x, shape_y, shape_z), dtype=peak_values.dtype)
            # Fill the non-zero locations for the data cube with data
            for ni, ci in enumerate(items):
                try:
                    # Pixel indices may be out of order (e.g, when we use MPI) so we look up the pixel location
                    current_index = np.nonzero(np.logical_and(array_indices[0] == ci[0],
                                                              array_indices[1] == ci[1]))[0][0]
                except:
                    log_helper.warning(__name__, "Requested pixel not found: " + str(items[ni]))
                    continue
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
        array_indices = analysis_object['peak_arrayindex'][:]
        x_size = array_indices[:, 0].max()+1
        y_size = array_indices[:, 1].max()+1
        valuesX = range(0, x_size)
        labelX = 'pixel index X'
        valuesY = range(0, y_size)
        labelY = 'pixel index Y'
        if array_indices.shape[1] > 2:
            z_size = array_indices[:, 2].max()+1
            valuesZ = range(0, z_size)
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None

        # We do not have native option for qslice, so we rely on the input data in all cases
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  # Loadings
            mz_spectra = analysis_object['indata_mz'][:]
            label_spectra = "m/z"
            mz_slice = None
            label_slice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ= \
                super(omsi_findpeaks_local, cls).v_qmz(analysis_object,
                                                       qslice_viewer_option=qslice_viewer_option,
                                                       qspectrum_viewer_option=qspectrum_viewer_option-1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option >= 0:
            mz_spectra = analysis_object['indata_mz'][:]
            label_spectra = "m/z"
            temp_a, temp_b, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_findpeaks_local, cls).v_qmz(analysis_object,
                                                       0,
                                                       qspectrum_viewer_option)
            # NOTE: if qspectrum and qslice share the same axis, this call will
            # not return the copied data, i.e., we need to copy the
            # qspectrum values to the qslice values.
            if mz_slice is None:
                mz_slice = temp_a
                label_slice = temp_b

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

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

    def write_analysis_data(self, analysis_group=None):
        """
        This function is used to write the actual analysis data to file. If not implemented, then the
        omsi_file_analysis API's default behavior is used instead.

        :param analysis_group: The h5py.Group object where the analysis is stored.

        """
        # Check if a user attempts to do parallel I/O with collect being disabled
        if mpi_helper.get_size() > 1 and not self['collect']:
            # Check if any of the other ranks have data
            num_elements = self['peak_arrayindex'].shape[0] if len(self['peak_arrayindex'].shape) == 2 else 0
            result_sizes = mpi_helper.gather(num_elements, comm=self.mpi_comm, root=self.mpi_root)
            if mpi_helper.get_rank() == self.mpi_root:
                for element_size in result_sizes[1:]:
                    if element_size > 0:
                        raise ValueError('Parallel I/O with collect parameter set to false not supported')
        raise NotImplementedError
        """
        import numpy as np

        # Serial no MPI, single rank SERIAL with MPI, or we are on the mpi root rank where we have all the data
        # This is a purely serial write using the standard mechanism
        if not mpi_helper.MPI_AVAILABLE or \
                mpi_helper.get_size()==1 or \
                (self['collect'] and mpi_helper.get_rank() == self.mpi_root):  #
            raise NotImplementedError   # Just let the default implementation of omsi_file_analysis handle this
        else:
            # TODO Implement the data write when the data is distributed
            rank = mpi_helper.get_rank(comm=self.mpi_comm)
            #if rank > 0:    # FIXE REMOVE this if statement when all is done
            #    dat = np.ones(5,dtype=float)
            #    dat = dat * rank
            #    self['peak_mz'] = dat
            #    self['peak_value'] = dat
            #    dat2  = np.ones(5*3, dtype='i').reshape(5,3)
            #    dat2 *= rank
            #    self['peak_arrayindex'] = dat2
            rank_empty = len(self['peak_mz'].shape) == 0
            if rank ==  self.mpi_root:
                if analysis_group is None:
                    raise ValueError("Not file open at the MPI root for writing")
                for dataset_name in ['peak_mz', 'peak_value']: # , 'peak_arrayindex']:
                    if rank_empty:
                        num_values = np.array(0, dtype='i')
                    else:
                        num_values = np.array(self[dataset_name].size, dtype='i')
                    rank_sizes = np.zeros(self.mpi_comm.size, dtype='i')
                    self.mpi_comm.Gather(num_values, rank_sizes, root=self.mpi_root)
                    total_num_values = rank_sizes.sum()
                    rank_displ = [0] + np.cumsum(rank_sizes[0:-1]).tolist()
                    rank_sizes = rank_sizes.tolist()

                    my_dataset = self[dataset_name] if not rank_empty else np.empty(0, dtype='float')
                    send_buff = (my_dataset, num_values[()])
                    # print total_num_values
                    target_buff = np.zeros(total_num_values, dtype=float)
                    recv_buff = [target_buff, rank_sizes, rank_displ, mpi_helper.mpi_type_from_dtype(np.dtype(float))]
                    self.mpi_comm.Gatherv(send_buff, recv_buff, root=self.mpi_root)
                    analysis_group[dataset_name] = target_buff
                    analysis_group.file.flush()
                    # print analysis_group[dataset_name][:]
                    del target_buff
                # TODO Add collecting and saving of peak_arrayindex
                # TODO Write mzdata from the raw file (which we have on all cores, i.e, no gather needed)
            else:
                for dataset_name in ['peak_mz', 'peak_value']: #, 'peak_arrayindex']:
                    if rank_empty:
                        num_values = np.array(0, dtype='i')
                    else:
                        num_values = np.array(self[dataset_name].size, dtype='i')
                    self.mpi_comm.Gather(num_values, None, root=self.mpi_root)
                    my_dataset = self[dataset_name] if not rank_empty else np.empty(0, dtype='float')
                    if dataset_name == 'peak_arrayindex':
                        send_buff = (my_dataset, num_values[()]*3)
                    else:
                        send_buff = (my_dataset, num_values[()])
                    recv_buff = None
                    self.mpi_comm.Gatherv(send_buff, recv_buff, root=self.mpi_root)
        """

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
        import numpy as np

        # Assign parameters to local variables for convenience
        msidata = self['msidata']
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
        # We have more than a single core AND we have multiple spectra to process
        if mpi_helper.get_size() > 1 and len(self['msidata'].shape) > 1:
            # We were not asked to process a specific data subblock from a parallel process
            # but we need to initiate the parallel processing.
            if msidata_subblock is None:
                # Setup the parallel processing using mpi_helper.parallel_over_axes
                split_axis = range(len(self['msidata'].shape)-1)  # The axes along which we can split the data
                scheduler = mpi_helper.parallel_over_axes(task_function=self.execute_analysis,  # Execute this function
                                                          task_function_params={},              # No added parameters
                                                          main_data=msidata,                    # Process the msidata
                                                          split_axes=split_axis,                # Split along axes
                                                          main_data_param_name='msidata_subblock',  # data input param
                                                          root=self.mpi_root,                   # The root MPI task
                                                          schedule=self['schedule'],            # Parallel schedule
                                                          comm=self.mpi_comm)                   # MPI communicator
                # Execute the analysis in parallel
                result = scheduler.run()
                # Collect the output data to the root rank if requested
                if self['collect']:
                    result = scheduler.collect_data()

                # TODO Record runtime information data from the scheduler in our provenance data
                # self.run_info['SCHEDULER_blocks'] = scheduler.blocks
                # self.run_info['SCHEDULER_block_times'] = scheduler.block_times
                # self.run_info['SCHEDULER_run_time'] = scheduler.run_time
                # self.run_info['SCHEDULER_schedule'] = scheduler.schedule

                # Compile the data from the parallel execution
                # Case Table:
                #
                # collect + worker       2
                #           worker       2
                # collect + root         3
                #           root         1
                use_dynamic_schedule = (self['schedule'] == mpi_helper.parallel_over_axes.SCHEDULES['DYNAMIC'])
                # Case 1: root rank without collect data disabled
                if mpi_helper.get_rank() == self.mpi_root and not self['collect']:
                    # We did not process any data on the root if DYNAMIC scheduling was used
                    if use_dynamic_schedule:
                        return None, None, None, mzdata
                    # We processed a data block using dynamic scheduling
                    else:
                        return result[0][0]
                # Case 2: Compile the data on the worker
                elif mpi_helper.get_rank() != self.mpi_root: # and use_dynamic_schedule:
                    # Compile the results from all processing task (on workers) or from all workers (on the root)
                    peak_mz = np.concatenate(tuple([ri[0] for ri in result[0]]), axis=-1)
                    peak_values = np.concatenate(tuple([ri[1] for ri in result[0]]), axis=-1)
                    if len(result[1])>1: # Correct indices from the individual runs since they all start at 0
                        peak_arrayindex = np.asarray([ [b[0], b[1], 0] for b in result[1] ])
                        peak_arrayindex[:,2] = np.cumsum([0] + [ len(ri[0]) for ri in result[0] ])[:-1]
                    else:
                        peak_arrayindex = result[0][0][2]
                    mzdata = result[0][0][3]
                    return peak_mz, peak_values, peak_arrayindex, mzdata
                # Case 3: Compile collected data on the root
                elif mpi_helper.get_rank() == self.mpi_root: # and use_dynamic_schedule:
                    # Compile the results from all processing task (on workers) or from all workers (on the root)
                    peak_mz = np.concatenate(tuple([ri[0] for ri in result[0]]), axis=-1)
                    peak_values = np.concatenate(tuple([ri[1] for ri in result[0]]), axis=-1)
                    # Dynamic scheduling uses selections of (int,int,slice) while the static
                    # scheduling uses (slice, slice, slice), hence we need to compile the peak_arrayindex
                    # slightly differently depending on the scheduler used
                    if use_dynamic_schedule:
                        peak_arrayindex = np.asarray([ [b[0], b[1], 0] for b in result[1]])
                        peak_arrayindex[:,2] = np.cumsum([0] + [ len(ri[0]) for ri in result[0]])[:-1]
                    else:
                        peak_arrayindex = np.concatenate(tuple([ri[2] for ri in result[0]]), axis=0)
                        d = np.cumsum([0] + [len(ri[0]) for ri in result[0]])
                        d2 = np.cumsum([0] + [len(ri[2]) for ri in result[0]])
                        for di in range(len(d2)-1):
                            peak_arrayindex[d2[di]:d2[di+1], 2] += d[di]
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
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_findpeaks_local).main()

