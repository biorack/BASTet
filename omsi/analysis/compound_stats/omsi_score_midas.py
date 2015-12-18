"""
MIDAS spectrum analysis
"""

# TODO  Add support for using centroided MSI data directly by removing 0's from the intensity and mz array
# TODO  We can further parallelize the callculations by splitting the compound list up as well (not just the spectra)

from omsi.analysis.base import analysis_base
import omsi.shared.mpi_helper as mpi_helper
from omsi.shared.log import log_helper
try:
    import MIDAS
except ImportError:
    log_helper.error(__name__, "Import of MIDAS failed. The omsi_score_midas module will not work.")
import os
import numpy as np
import time
import sys


class omsi_score_midas(analysis_base):
    """
    Class for executing midas on an MSI or local peak finding dataset.
    """

    def __init__(self, name_key="undefined"):
        """
        Initialize the basic data members
        """
        super(omsi_score_midas, self).__init__()
        self.analysis_identifier = name_key
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        # FIXME Check that the group and required settings for parameters are set appropriately
        self.add_parameter(name='fpl_data',
                           help='The raw findpeaks local dataset to be analyzed. Either the finpeaks local ' +
                                'group in the file or the omsi_findpeaks_local analysis object',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        # FIXME This is not in the config but should be read from the MS2 data if available
        self.add_parameter(name='precursor_mz',
                           help='Floating point precursor mass over charge value',
                           dtype=dtypes['float'],
                           required=True,  # FIXME make this false once we can read it from file
                           default=-1,
                           group=groups['input'])
        # FIXME This is not in the config file. Needed to compute the neutralized *mass* in Da of the parent molecule
        self.add_parameter(name='default_charge',
                           help='Positive (1) or negative (-1) default charge state of the parent molecule',
                           dtype=dtypes['int'],
                           choices=[1, -1],
                           default=1,
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='precursor_type',     # dCurrentPrecursor_type
                           help='Positive (1) or negative (-1) ion',
                           dtype=dtypes['int'],
                           choices=[1, -1],
                           default=1,
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='metabolite_database',
                           help='The database of metabolites to be used.',
                           dtype=dtypes['str'],
                           required=True,
                           default=os.path.join(os.path.dirname(__file__), 'third_party/MetaCyc.mdb'),
                           group=groups['settings'])
        self.add_parameter(name='parent_mass_windows',
                           help='List of mass windows to be open around parent ion',
                           dtype=dtypes['ndarray'],
                           default=[0, ],
                           group=groups['settings'])
        self.add_parameter(name='positive_ion_fragment_mass_windows',
                           help='List of ints of offsets to apply to frag m/zs',
                           dtype=dtypes['ndarray'],
                           default=[0, 1, 2],
                           group=groups['settings'])
        self.add_parameter(name='negative_ion_fragment_mass_windows',
                           help='List of ints of offsets to apply to frag m/zs',
                           dtype=dtypes['ndarray'],
                           default=[-2, -1, 0],
                           group=groups['settings'])
        self.add_parameter(name='mass_tolerance_parent_ion',
                           help='Parent mass tolerance in Da',
                           dtype=dtypes['float'],
                           default=0.01,
                           group=groups['settings'])
        self.add_parameter(name='mass_tolerance_fragment_ions',
                           help='Fragment mass tolerance in Da',
                           dtype=dtypes['float'],
                           default=0.01,
                           group=groups['settings'])
        self.add_parameter(name='break_rings',
                           help='Should we attempt to break rings in the molecule?',
                           dtype=dtypes['bool'],
                           default=True,
                           group=groups['settings'])
        self.add_parameter(name='fragmentation_depth',
                           help='Depth of fragmentation pathways',
                           dtype=dtypes['int'],
                           choices=[1, 2, 3, 4, 5],
                           default=3,
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
        self.data_names = ['hit_table',
                           'pixel_index']

    def execute_analysis(self, spectrum_indexes=None, compound_list=None):
        """
        Execute the local peak finder for the given msidata.

        :param spectrum_indexes: List with a list of integer indicies of the subset of sepctra
            that should be processed by this MPI task.  If spectrum_indexes is set, then the given
            subblock will be processed in SERIAL instead of processing self['fpl_data'] in PARALLEL
            (if available). This parameter is strictly optional and intended for internal use only
            to facilitate the efficient parallel implementation.

        :param compound_list: List of the compounds from the database file. This parameter is used
            to avoid having to read the compound database on every compute task that calls this function
            when running in parallel.  This  parameter is strictly optional and intended for internal
            use only to facilitate the efficient parallel implementation.

        :returns: A tuple with an array of hit_tables with the scores for each pixel and a 2D array
            of pixel indices describing for each spectrum the (x,y) pixel location in the image. The
            hit_table is an array of (#spectra x #compounds). The hit_table is a structured numpy
            array with the following columns:

                * 'score',  float,  MIDAS score of row
                * 'id',     str,    database ID e.g. 'MetaCyC_7884'
                * 'name',   str,    database name, e.g. 'glycine'
                * 'mass',   float,  mass in Da of IDed compound
                * 'n_peaks', int,   number of peaks in data
                * 'n_match', int,   number of peaks in data matched

        """
        # Assign parameter settings to local variables for convenience
        metabolite_database = self['metabolite_database']
        precursor_type = self['precursor_type']
        parent_mass_windows = self['parent_mass_windows']
        positive_ion_fragment_mass_windows = self['positive_ion_fragment_mass_windows']
        negative_ion_fragment_mass_windows = self['negative_ion_fragment_mass_windows']
        mass_tolerance_parent_ion = self['mass_tolerance_parent_ion']
        mass_tolerance_fragment_ions = self['mass_tolerance_fragment_ions']
        break_rings = self['break_rings']
        fragmentation_depth = self['fragmentation_depth']

        # Calculate the parent_mass
        precursor_mz = self['precursor_mz']              # FIXME  Get the precursor_mz from the MS2 data
        if precursor_mz == -1:
             precursor_mz = self['fpl_data']['precursor_mz'][:]
        default_charge = self['default_charge']          # FIXME  Is this an input or should we get this from file
        proton_mass = 1.00782503207 - 5.4857990946e-4
        parent_mass = precursor_mz - (default_charge * proton_mass)

        # Get the data we need to process
        fpl_data = self['fpl_data']
        fpl_peak_mz = fpl_data['peak_mz']
        fpl_peak_value = fpl_data['peak_value']
        fpl_peak_arrayindex = fpl_data['peak_arrayindex']

        # Get the compound list if we have not read it previously.
        if compound_list is None:
            # TODO: Possible further optimization by reading only on self.mpi_root and then sending the list to all
            compound_list = MIDAS.ReadCompoundFile(metabolite_database)

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
                # Setup the parallel processing using mpi_helper.parallel_over_axes
                split_axis = [0, ]
                scheduler = mpi_helper.parallel_over_axes(
                    task_function=self.execute_analysis,                    # Execute this function
                    task_function_params={'compound_list': compound_list},  # Reuse the compound_list
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
                hit_table = np.zeros((0, 0), dtype=MIDAS.scoring_C.HIT_TABLE_DTYPE)  # initialize hit_table as empty
                pixel_index = np.zeros((0, 2), dtype='int')
                use_dynamic_schedule = (self['schedule'] == mpi_helper.parallel_over_axes.SCHEDULES['DYNAMIC'])

                # TODO NEED to update since collect now returns a single list not a list of lists
                if not self['collect'] and (mpi_helper.get_rank() == self.mpi_root and use_dynamic_schedule):
                    # We did not process any data on the root process when using dynamic scheduling
                    # and we did not collect the data to the root either
                    pass
                #elif self['collect'] and mpi_helper.get_rank() == self.mpi_root:
                #    temp_data = [ri[0] for rt in result[0] for ri in rt]
                #    if len(temp_data) > 0:
                #        hit_table = np.concatenate(tuple(temp_data), axis=-1)
                #    temp_data = [ri[1] for rt in result[0] for ri in rt]
                #    if len(temp_data) > 0:
                #        pixel_index = np.concatenate(tuple(temp_data), axis=0) # axis=-1
                else:
                    temp_data = [ri[0] for ri in result[0]]
                    if len(temp_data) > 0:
                        hit_table = np.concatenate(tuple(temp_data), axis=-1)
                    temp_data = [ri[1] for ri in result[0]]
                    if len(temp_data) > 0:
                        pixel_index = np.concatenate(tuple(temp_data), axis=0)
                return hit_table, pixel_index

        #############################################################
        # Serial processing of the current data block
        #############################################################
        # Initialize the output data structures
        pixel_index = fpl_peak_arrayindex[spectrum_indexes, 0:2]
        if len(pixel_index.shape) == 1:
            pixel_index = pixel_index[np.newaxis, :]
        hit_table = None  # FIXME The initalization of the hit_table is only valid if we assume that all spectra have the same precursor m/z, which may not be the case
        # Iterate through all the pixel we were asked to process in serial
        for current_index, spectrum_index in enumerate(spectrum_indexes):
            # Determine the start and stop index for the m/z and intensity data of the current spectrum
            start = fpl_peak_arrayindex[spectrum_index, 2]
            stop = fpl_peak_arrayindex[(spectrum_index+1), 2] \
                if spectrum_index < (num_spectra-1) \
                else fpl_peak_value.size
            spectrum_length = stop - start
            # Skip empty spectra
            if spectrum_length == 0:
                time_str =  "rank : " + str(mpi_helper.get_rank()) + " : pixel_index : " + str(fpl_peak_arrayindex[spectrum_index, 0:2]) + " Spectrum not scored."
                print time_str
                continue
            # Load the m/z and intensity values for the current spectrum
            current_peaks_list = np.zeros(shape=(spectrum_length, 3), dtype=float)
            current_peaks_list[:, 0] = fpl_peak_mz[start:stop]
            current_peaks_list[:, 1] = fpl_peak_value[start:stop]

            # Get the parent mass
            current_parent_mass = parent_mass if len(parent_mass) == 1 else parent_mass[spectrum_index]

            start_time = time.time()
            # Call MIDAS to score the current spectrum against all compounds in the database
            current_hits = MIDAS.scoring_C.score_main(
                Compound_list=compound_list,
                bBreakRing=break_rings,
                dCurrentPrecursor_type=precursor_type,
                dCurrentParentMass=current_parent_mass,
                current_peaks_list=current_peaks_list,
                iParentMassWindow_list=parent_mass_windows,
                dMass_Tolerance_Parent_Ion=mass_tolerance_parent_ion,
                dMass_Tolerance_Fragment_Ions=mass_tolerance_fragment_ions,
                iFragmentation_Depth=fragmentation_depth,
                iPositive_Ion_Fragment_Mass_Windows_list=positive_ion_fragment_mass_windows,
                iNegative_Ion_Fragment_Mass_Windows_list=negative_ion_fragment_mass_windows,
                top_n=None)

            end_time = time.time()
            execution_time = end_time - start_time
            time_str =  "rank : " + str(mpi_helper.get_rank()) + " : pixel_index : " + str(fpl_peak_arrayindex[spectrum_index, 0:2]) + " : time in s : " + str(execution_time)
            time_str += " : num hits : " + str(current_hits.shape[0])
            print time_str
            sys.stdout.flush()

            # Initialize the hit_table if necessary
            if hit_table is None:
                # If our compound database does not contain any related compounds then just finish
                if current_hits.shape[0] == 0:
                    # Initialize the results as empty and finish as there is nothing to do
                    hit_table = np.zeros(shape=(pixel_index.shape[0], 0),
                                         dtype=MIDAS.scoring_C.HIT_TABLE_DTYPE)  # FIXME the number of hits may be different for different spectra if we have varying precursor m/z
                    continue
                # If our compound database contains at least one relevant compound then check all spectra
                else:
                    # Create the data structure to store all results
                    hit_table = np.zeros(shape=(pixel_index.shape[0], current_hits.shape[0]),
                                         dtype=current_hits.dtype)  # FIXME the number of hits may be different for different spectra if we have varying precursor m/z
            # Save the hits for the current pixel
            hit_table[current_index] = current_hits

        if hit_table is None:
            hit_table = np.zeros(shape=(pixel_index.shape[0], 0),
                                 dtype=MIDAS.scoring_C.HIT_TABLE_DTYPE)

        # Return the hit_table and the index of the pixel each hit_table applies to
        return hit_table, pixel_index

if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_score_midas).main()

