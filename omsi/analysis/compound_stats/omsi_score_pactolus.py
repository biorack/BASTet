"""
MIDAS spectrum analysis
"""

# TODO  Add support for using centroided MSI data directly by removing 0's from the intensity and mz array
# TODO  We can further parallelize the callculations by splitting the compound list up as well (not just the spectra)

from omsi.analysis.base import analysis_base
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.log import log_helper
try:
    from pactolus import score_frag_dag
except ImportError:
    log_helper.error(__name__, "Import of Pactolus failed. The  omsi_score_pactolus module will not work.")
import os
import numpy as np
import time
import sys

class omsi_score_pactolus(analysis_base):
    """
    Class for executing Pactolus on a local peak finding dataset.
    """

    def __init__(self, name_key="undefined"):
        """
        Initialize the basic data members
        """
        super(omsi_score_pactolus, self).__init__()
        self.analysis_identifier = name_key
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='fpl_data',
                           help='The raw findpeaks local dataset to be analyzed. Either the finpeaks local ' +
                                'group in the file or the omsi_findpeaks_local analysis object',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='precursor_mz',
                           help='Floating point precursor mass over charge value',
                           dtype=dtypes['float'],
                           required=False,
                           default=-1,
                           group=groups['input'])
        self.add_parameter(name='metabolite_database',
                           help='The database of metabolites to be used.',
                           dtype=dtypes['str'],
                           required=False,
                           default="",
                           group=groups['settings'])
        self.add_parameter(name='trees',
                           help='1) Path to the sub-directory with all relevent .h5 pactolus fragmentation trees. or ' +
                                '2) path to a text-file with a list of names of the tree files, or 3) Path to the' +
                                '.npy file with the pactolus file lookup table',
                           dtype=dtypes['str'],
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='ms1_mass_tolerance',
                           help='float, max. diff. in Da of tree_parent_mass & MS1_precursor_mass',
                           dtype=dtypes['float'],   # dtypes['ndarray'],
                           default=[0, ],
                           group=groups['settings'])
        self.add_parameter(name='ms2_mass_tolerance',
                           help='float, max. mass in Da by which two MS2/MSn peaks can differ',
                           dtype=dtypes['float'],  # dtypes['ndarray']
                           default=0.01,
                           group=groups['settings'])
        self.add_parameter(name='max_depth',
                           help='Maximum depth of fragmentation pathways',
                           dtype=dtypes['int'],
                           choices=[1, 2, 3, 4, 5],
                           default=5,
                           group=groups['settings'])
        self.add_parameter(name='neutralizations',
                           help='List of floats, adjustments (in Da) added to data peaks in order to neutralize them.',
                           dtype=dtypes['ndarray'],
                           default=[0, 1, 2],
                           group=groups['settings'])
        # Parallel execution parameters
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
        self.data_names = ['pixel_index', 'score', 'id', 'name', 'mass', 'n_peaks', 'n_match']

    def execute_analysis(self, spectrum_indexes=None, file_lookup_table=None):
        """
        Execute the local peak finder for the given msidata.

        :param spectrum_indexes: List with a list of integer indicies of the subset of sepctra
            that should be processed by this MPI task.  If spectrum_indexes is set, then the given
            subblock will be processed in SERIAL instead of processing self['fpl_data'] in PARALLEL
            (if available). This parameter is strictly optional and intended for internal use only
            to facilitate the efficient parallel implementation.

        :param file_lookup_table: The Pactolus lookup table with the list of tree files and their mass.

        :returns: A series of numpy arrays  with the score data for each pixel and a 2D array
            of pixel indices describing for each spectrum the (x,y) pixel location in the image.

            ['pixel_index', 'score', 'id', 'name', 'mass', 'n_peaks', 'n_match']
                * 'pixel_index'  , int,  2D array of pixel indices describing for each spectrum \
                   the (x,y) pixel location in the imag
                * 'score',  float,  MIDAS score of row
                * 'id',     str,    database ID e.g. 'MetaCyC_7884'
                * 'name',   str,    database name, e.g. 'glycine'
                * 'mass',   float,  mass in Da of IDed compound
                * 'n_peaks', int,   number of peaks in data
                * 'n_match', int,   number of peaks in data matched

        """
        log_helper.debug(__name__, 'Reading inputs', comm=self.mpi_comm, root=self.mpi_root)
        # Get the data we need to process
        fpl_data = self['fpl_data']
        fpl_peak_mz = fpl_data['peak_mz']
        fpl_peak_value = fpl_data['peak_value']
        fpl_peak_arrayindex = fpl_data['peak_arrayindex']
        # Calculate the parent_mass
        precursor_mz = self['precursor_mz']
        if precursor_mz == -1:
            precursor_mz = self['fpl_data']['precursor_mz'][:]
        # Assign parameter settings to local variables for convenience
        metabolite_database = self['metabolite_database']
        ms1_mass_tol = self['ms1_mass_tolerance']
        ms2_mass_tol = self['ms2_mass_tolerance']
        neutralizations = self['neutralizations']
        max_depth = self['max_depth']

        # Make the numpy array with the list of tree files and their MS1 masses
        if file_lookup_table is None:
            # TODO: Possible further optimization by reading only on self.mpi_root and then sending the list to all
            log_helper.debug(__name__, 'Preparing file lookup table', comm=self.mpi_comm, root=self.mpi_root)
            if os.path.isfile(self['trees']):
                if self['trees'].endswith('.npy'):
                    file_lookup_table = np.load(self['trees'])
                else:
                    in_treefile = open(self['trees'], 'r')
                    tree_files = [line.rstrip('\n') for line in in_treefile]
                    in_treefile.close()
                    file_lookup_table = score_frag_dag.make_file_lookup_table_by_MS1_mass(tree_files=tree_files)
            elif os.path.isdir(self['trees']):
                file_lookup_table = score_frag_dag.make_file_lookup_table_by_MS1_mass(path=self['trees'])

        # Define the common pactolus paramters
        pactolus_parameters = {'file_lookup_table': file_lookup_table,
                               'ms1_mass_tol': ms1_mass_tol,
                               'ms2_mass_tol': ms2_mass_tol,
                               'neutralizations': neutralizations,
                               'max_depth': max_depth}

        # Get the peak_arrayindex with [[x,y, array_offset], ...] values describing the
        # index of the pixel in (x,y) and the offset in the peak_mz and peak_value array
        # where we can find the spectrum that we need to processes
        num_spectra = fpl_peak_arrayindex.shape[0]
        if spectrum_indexes is None:
            # Get the complete peak array index data
            spectrum_indexes = np.arange(0, num_spectra)
            enable_parallel = True
        else:
            if isinstance(spectrum_indexes, int):
                spectrum_indexes = np.asarray([spectrum_indexes, ])
            enable_parallel = False

        #############################################################
        # Parallel execution using MPI
        #############################################################
        # We have more than a single core AND we have multiple spectra to process
        if mpi_helper.get_size() > 1 and len(spectrum_indexes) > 1:
            # We were not asked to process a specific data subblock from a parallel process
            # but we need to initiate the parallel processing.
            if enable_parallel:
                log_helper.debug(__name__, 'Preparing parallel execution', comm=self.mpi_comm, root=self.mpi_root)
                # Setup the parallel processing using mpi_helper.parallel_over_axes
                split_axis = [0, ]
                scheduler = mpi_helper.parallel_over_axes(
                    task_function=self.execute_analysis,                    # Execute this function
                    task_function_params={'file_lookup_table': file_lookup_table},  # Reuse the file_lookup_table
                    main_data=spectrum_indexes,                             # Process the spectra independently
                    split_axes=split_axis,                                  # Split along axes
                    main_data_param_name='spectrum_indexes',                # data input param
                    root=self.mpi_root,                                     # The root MPI task
                    schedule=self['schedule'],                              # Parallel scheduling scheme
                    comm=self.mpi_comm)                                     # MPI communicator
                # Execute the analysis in parallel
                result = scheduler.run()
                # Collect the output data to the root rank if requested
                if self['collect']:
                    result = scheduler.collect_data()

                # Compile the data from the parallel execution
                pixel_index = np.zeros((0, 2), dtype='int')
                score = np.zeros((0,), dtype='f4')
                id_data = np.zeros((0,), dtype='a100')
                name = np.zeros((0,), dtype='a100')
                mass = np.zeros((0,), dtype='f4')
                n_peaks = np.zeros((0,), dtype='i4')
                n_match = np.zeros((0,), dtype='i4')

                use_dynamic_schedule = (self['schedule'] == mpi_helper.parallel_over_axes.SCHEDULES['DYNAMIC'])

                # TODO NEED to update since collect now returns a single list not a list of lists
                if not self['collect'] and (mpi_helper.get_rank() == self.mpi_root and use_dynamic_schedule):
                    # We did not process any data on the root process when using dynamic scheduling
                    # and we did not collect the data to the root either
                    pass
                # elif self['collect'] and mpi_helper.get_rank() == self.mpi_root:
                #    temp_data = [ri[0] for rt in result[0] for ri in rt]
                #    if len(temp_data) > 0:
                #        hit_table = np.concatenate(tuple(temp_data), axis=-1)
                #    temp_data = [ri[1] for rt in result[0] for ri in rt]
                #    if len(temp_data) > 0:
                #        pixel_index = np.concatenate(tuple(temp_data), axis=0) # axis=-1
                else:
                    log_helper.debug(__name__, 'Compiling output')
                    # Compile pixel_index
                    temp_data = [ri[0] for ri in result[0]]
                    if len(temp_data) > 0:
                        pixel_index = np.concatenate(tuple(temp_data), axis=0)
                    temp_data = [ri[1] for ri in result[0]]
                    # Compile scores
                    if len(temp_data) > 0:
                        score = np.concatenate(tuple(temp_data), axis=0)
                    # Compile id
                    temp_data = [ri[2] for ri in result[0]]
                    if len(temp_data) > 0:
                        id_data = np.concatenate(tuple(temp_data), axis=0)
                    # Compile name
                    temp_data = [ri[3] for ri in result[0]]
                    if len(temp_data) > 0:
                        name = np.concatenate(tuple(temp_data), axis=0)
                    # Compile mass
                    temp_data = [ri[4] for ri in result[0]]
                    if len(temp_data) > 0:
                        mass = np.concatenate(tuple(temp_data), axis=0)
                    # Compile n_peaks
                    temp_data = [ri[5] for ri in result[0]]
                    if len(temp_data) > 0:
                        n_peaks = np.concatenate(tuple(temp_data), axis=0)
                    # Compile n_match
                    temp_data = [ri[6] for ri in result[0]]
                    if len(temp_data) > 0:
                        n_match = np.concatenate(tuple(temp_data), axis=0)
                    log_helper.log_var(__name__, score=score)
                # Return the compiled output
                return pixel_index, score, id_data, name, mass, n_peaks, n_match

        #############################################################
        # Serial processing of the current data block
        #############################################################
        log_helper.debug(__name__, 'Processing spectra', comm=self.mpi_comm, root=self.mpi_root)
        # Initialize the output data structures
        # pixel_index = fpl_peak_arrayindex[spectrum_indexes, 0:2]
        # if len(pixel_index.shape) == 1:
        #    pixel_index = pixel_index[np.newaxis, :]
        hit_matrix = []

        # Iterate through all the pixel we were asked to process in serial
        for current_index, spectrum_index in enumerate(spectrum_indexes):
            # Determine the start and stop index for the m/z and intensity data of the current spectrum
            start = int(fpl_peak_arrayindex[spectrum_index, 2])
            stop = int(fpl_peak_arrayindex[(spectrum_index+1), 2]
                   if spectrum_index < (num_spectra-1)
                   else fpl_peak_value.size)
            spectrum_length = stop - start
            # Skip empty spectra
            if spectrum_length == 0:
                time_str = "rank : " + str(mpi_helper.get_rank()) + " : pixel_index : " + \
                           str(fpl_peak_arrayindex[spectrum_index, 0:2]) + " Spectrum not scored."
                log_helper.info(__name__, time_str, comm=self.mpi_comm, root=None)
                continue
            # Load the m/z and intensity values for the current spectrum
            current_peaks_list = np.zeros(shape=(spectrum_length, 2), dtype=float)
            current_peaks_list[:, 0] = fpl_peak_mz[start:stop]
            current_peaks_list[:, 1] = fpl_peak_value[start:stop]

            # Get the parent mass
            current_parent_mass = precursor_mz if len(precursor_mz) == 1 else precursor_mz[spectrum_index]

            start_time = time.time()
            # Call MIDAS to score the current spectrum against all compounds in the database
            current_hits = score_frag_dag.score_scan_list_against_trees(scan_list=[current_peaks_list, ],
                                                                        ms1_mz=[current_parent_mass, ],
                                                                        params=pactolus_parameters)
            end_time = time.time()
            execution_time = end_time - start_time
            time_str = "rank : " + str(mpi_helper.get_rank()) + " : pixel_index : " + \
                       str(fpl_peak_arrayindex[spectrum_index, 0:2]) + " : time in s : " + str(execution_time)
            time_str += " : num hits : " + str((current_hits > 0).sum())
            #log_helper.info(__name__, time_str, comm=self.mpi_comm, root=None)
            #sys.stdout.flush()
            print time_str
            sys.stdout.flush()

            # Save the hits for the current pixel
            hit_matrix.append(current_hits[0, :])

        # Index the results based on the given metabolite database
        score = []
        id_data = []
        name = []
        mass = []
        n_peaks = []
        n_match = []
        pixel_index = []
        if len(metabolite_database) > 0:  # We don't have an empty string
            for current_index, spectrum_index in enumerate(spectrum_indexes):
                non_zero_scores = np.where(hit_matrix[current_index] > 0)
                if non_zero_scores.size > 0:
                    current_hit_table = np.asarray(score_frag_dag.make_pactolus_hit_table(
                        pactolus_results=hit_matrix[current_index],
                        table_file=file_lookup_table,
                        original_db=metabolite_database))
                    for score_index in non_zero_scores:
                        pixel_index.append(fpl_peak_arrayindex[spectrum_index, 0:2])
                        score.append(current_hit_table['score'][score_index])
                        id_data.append(current_hit_table['id'][score_index])
                        name.append(current_hit_table['name'][score_index])
                        mass.append(current_hit_table['mass'][score_index])
                        n_peaks.append(current_hit_table['n_peaks'][score_index])
                        n_match.append(current_hit_table['n_match'][score_index])
        else:
            pixel_index = fpl_peak_arrayindex[spectrum_indexes, 0:2]
            score = np.asarray(hit_matrix)

        # Return the hit_table and the index of the pixel each hit_table applies to
        print "rank : " + str(mpi_helper.get_rank()) + " : scores " + str(score)
        sys.stdout.flush()
        return np.asarray(pixel_index), \
               np.asarray(score), \
               np.asarray(id_data), \
               np.asarray(name), \
               np.asarray(mass), \
               np.asarray(n_peaks), \
               np.asarray(n_match)

if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_score_pactolus).main()

