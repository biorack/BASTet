"""
Module with the TIC normalization analysis.
"""
from tempfile import TemporaryFile

import numpy as np

from omsi.analysis.base import analysis_base
from omsi.shared.log import log_helper


###############################################################
#  1) Basic integration of your analysis with omsi (Required) #
###############################################################
class omsi_tic_norm(analysis_base):

    """
    TIC Normalization analysis.
    """
    def __init__(self,
                 name_key="undefined"):
        """Initalize the basic data members"""

        super(omsi_tic_norm, self).__init__()
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI matrix to be analyzed',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='mzdata',
                           help='The m/z values for the spectra of the MSI dataset',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='maxCount',
                           help='Maximum count',
                           dtype=dtypes['int'],
                           required=True,
                           group=groups['settings'],
                           default=0)
        self.add_parameter(name='mzTol',
                           help='MZ tolerance for the infIons',
                           dtype=dtypes['float'],
                           group=groups['settings'],
                           default=0.1,
                           required=True)
        self.add_parameter(name='infIons',
                           help='Optional list of informative ions',
                           dtype=dtypes['ndarray'],
                           required=False,
                           group=groups['settings'],
                           default=None)

        self.data_names = ['norm_msidata', 'norm_mz']
        self.analysis_identifier = name_key

    def execute_analysis(self):
        """
        Normalize the data based on the total intensity of a spectrum or the
        intensities of a select set of ions.

        Calculations are performed using a memory map approach to avoid loading
        all data into memory. TIC normalization can as such be performed even
        on large files (assuming sufficient disk space).

        Keyword Arguments:

        :param msidata: The input MSI dataset
        :type msidata: h5py.dataset or numpu array (3D)
        :param mzdata: The mz axsis do the dataset
        :param maxCount: ...
        :param mzTol: ...
        :param infIons: List of informative ions

        """
        # Setting Default Values and retrieving parameter data
        nx, ny, nz = self['msidata'].shape
        try:
            ion_list = self['infIons']
        except KeyError:
            ion_list = None
        mzdata = self['mzdata'][:]

        # Compute the mz values to be used for normalization
        if ion_list is not None:
            idx_mz = []
            for ion in ion_list:
                temp = np.where(abs(mzdata-ion) <= self['mzTol'])
                idx_mz = np.concatenate([idx_mz, temp[0]])
            idx_mz = idx_mz.astype(int)
        else:
            idx_mz = range(nz)

        output_filename = TemporaryFile()  # Create a temporary file to compute the normalization out-of-core
        idx = np.linspace(0, nx, 20, dtype='int')  # Divide the data into chunks in x to process the data in blocks

        # Compute the normalization factors and column maxs first (computed one-block-at-a-time)
        tic_norm_factors = np.zeros(shape=(nx, ny), dtype='float')
        msi_spectrum_maxs = np.zeros(shape=(nx, ny), dtype=self['msidata'].dtype)
        for i in range(len(idx)-1):  # Process the data one-block-at-a-time
            current_data = self['msidata'][idx[i]:idx[i+1], :, idx_mz]
            tic_norm_factors[idx[i]:idx[i+1], :] = np.sum(current_data, 2)
            msi_spectrum_maxs[idx[i]:idx[i+1], :] = np.amax(current_data[:, :, idx_mz], 2)
            del current_data
        non_zero_tic = tic_norm_factors > 0
        mean_tic_norm = float(tic_norm_factors[non_zero_tic].mean())
        tic_norm_factors[non_zero_tic] = 1.0 / (tic_norm_factors[non_zero_tic] / mean_tic_norm)

        # Determine the output data format
        max_output_value = np.max(np.multiply(msi_spectrum_maxs, tic_norm_factors))
        outformat = 'float'
        if max_output_value <= np.iinfo(np.uint16).max:
            outformat = 'uint16'
        elif max_output_value <= np.iinfo(np.uint32).max:
            outformat = 'uint32'
        log_helper.debug(__name__, "Output format: " + str(outformat))

        # Normalize the data one-block-at-a-time
        for i in range(len(idx)-1):
            num_spectra_x = idx[i+1]-idx[i]
            # Load the data for the current block amd reshape it
            current_data = self['msidata'][idx[i]:idx[i+1], :, :].reshape(num_spectra_x*ny, nz)
            tip = tic_norm_factors[idx[i]:idx[i+1], :].reshape(num_spectra_x*ny)
            mip = msi_spectrum_maxs[idx[i]:idx[i+1], :].reshape(num_spectra_x*ny)
            idx_thresh = np.where(np.logical_and(mip >= self['maxCount'], tip > 0))
            current_out_norm = np.zeros(current_data.shape)
            if len(idx_thresh[0]) > 0:
                current_out_norm[idx_thresh[0], :] = \
                    np.around(np.multiply(current_data[idx_thresh[0], :],
                                          tip.astype(float)[idx_thresh[0]][:, np.newaxis]))
                # Write the data to the temporary file
                ofp = np.memmap(output_filename, dtype=outformat, mode='w+', shape=(nx, ny, nz))
                ofp[idx[i]:idx[i+1], :, :] = current_out_norm.reshape(num_spectra_x, ny, nz)
                del ofp

        # for i in range(len(idx)-1): # Process the data one-block-at-a-time
        #     print "Current block: " +  str((idx[i], idx[i+1]))
        #     current_data = self['msidata'][idx[i]:idx[i+1], :, :]   # Load the data for the current block
        #     num_spectra = idx[i+1]-idx[i]  # Number of spectra in the current data block
        #     current_data = current_data.reshape(num_spectra*ny, nz)  # Reshape the data to a 2D array
        #     tip = np.sum(current_data[:, idx_mz], 1)   # Compute the sum of intensities per spectrum
        #     mip = np.amax(current_data[:, idx_mz], 1)  # Compute the maximum value per spectrum
        #     idx_thresh = np.where((mip >= self['maxCount']) & (tip > 0))
        #     current_out_norm = np.zeros(current_data.shape)  # The normalized data for the current block
        #     if len(idx_thresh[0]) > 0:
        #         # tic normalize the data for the current block
        #         current_out_norm[idx_thresh[0], :] = np.divide(current_data[idx_thresh[0], :],
        #                                                        tip.astype(float)[idx_thresh[0]][:, np.newaxis])
        #         # Write the data to the temporary file
        #         ofp = np.memmap(output_filename, dtype='float', mode='w+', shape=(nx, ny, nz))
        #         ofp[idx[i]:idx[i+1], :, :] = current_out_norm.reshape(num_spectra, ny, nz)
        #         del ofp

        # Save the data
        self['norm_msidata'] = np.memmap(output_filename, dtype=outformat, mode='w+', shape=(nx, ny, nz))
        self['norm_mz'] = mzdata
        output_filename.close()

    def record_execute_analysis_outputs(self, analysis_output):
        """
        We are not returning any outputs here, but we are going to record them manually.
        :param analysis_output: The output of the execute_analysis(...) function.
        """
        pass

    ###############################################################
    #  2) Integrating your analysis with the OpenMSI              #
    #     web-based viewer (Recommended)                          #
    ###############################################################
    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer

           :param analysis_object: The omsi_file_analysis object for which slicing should be performed
           :param z: Selection string indicting which z values should be selected.
           :param viewer_option: If multiple default viewer behaviors are available for a given analysis
                  then this option is used to switch between them.

           :returns: numpy array with the data to be displayed in the image slice viewer. Slicing will
                     be performed typically like [:,:,zmin:zmax].

        """
        # Convert the z selection to a python selection
        from omsi.shared.data_selection import selection_string_to_object
        z_select = selection_string_to_object(z)  # Convert the selection string to a python selection

        num_custom_viewer_options = 1

        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_tic_norm, cls).v_qslice(analysis_object,
                                                      z,
                                                      viewer_option=viewer_option-num_custom_viewer_options)
        elif viewer_option == 0:
            dataset = analysis_object['norm_msidata']
            return dataset[:, :, z_select]
        return None

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer

           Developer Note: h5py currently supports only a single index list. If the user provides an index-list for both
                           x and y, then we need to construct the proper merged list and load the data manually, or if
                           the data is small enough, one can load the full data into a numpy array which supports
                           multiple lists in the selection.

           :param analysis_object: The omsi_file_analysis object for which slicing should be performed
           :param x: x selection string
           :param y: y selection string
           :param viewer_option: If multiple default viewer behaviors are available for a given analysis then
                                this option is used to switch between them.

           :returns: The following two elements are expected to be returned by this function :

                1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The mass (m/z) axis must be the last \
                  axis. For index selection x=1,y=1 a 1D array is usually expected. For indexList selections \
                  x=[0]&y=[1] usually a 2D array is expected. For ragne selections x=0:1&y=1:2 we one usually \
                 expects a 3D array/
                2) None in case that the spectra axis returned by v_qmz are valid for the returned spectrum. \
                 Otherwise, return a 1D numpy array with the m/z values for the spectrum (i.e., if custom m/z \
                 values are needed for interpretation of the returned spectrum).This may be needed, e.g., in \
                 cases where a per-spectrum peak analysis is performed and the peaks for each spectrum appear \
                 at different m/z values.
        """

        # Convert the x,y selection to a python selection
        from omsi.shared.data_selection import selection_string_to_object
        x_select = selection_string_to_object(x)  # Convert the selection string to a python selection
        y_select = selection_string_to_object(y)  # Convert the selection string to a python selection

        num_custom_viewer_options = 1
        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_tic_norm, cls).v_qspectrum(analysis_object,
                                                         x,
                                                         y,
                                                         viewer_option=viewer_option-num_custom_viewer_options)
        elif viewer_option == 0:
            dataset = analysis_object['norm_msidata']
            # h5py supports only one indexlist as part of a selection. We have to work-around this
            # by retrieving the spectra we want one-at-a-time if the user provides two index lists
            if isinstance(x_select, list) and isinstance(y_select, list):
                xsize = len(x_select)
                ysize = len(y_select)
                if xsize != ysize:
                    raise ValueError("x and y selection size do not match")
                data = np.zeros(shape=(xsize, dataset.shape[2]), dtype=dataset.dtype)
                for i in xrange(0, xsize):
                    data[i, :] = dataset[x_select[i], y_select[i], :]
            else:
                data = dataset[x_select, y_select, :]
            return data, None
        return None, None

    @classmethod
    def v_qmz(cls,
              analysis_object,
              qslice_viewer_option=0,
              qspectrum_viewer_option=0):
        """ Get the mz axes for the analysis

            :param analysis_object: The omsi_file_analysis object for which slicing should be performed
            :param qslice_viewer_option: If multiple default viewer behaviors are available for a given analysis then
                                        this option is used to switch between them for the qslice URL pattern.
            :param qspectrum_viewer_option: If multiple default viewer behaviors are available for a given analysis
                                     then this option is used to switch between them for the qspectrum URL pattern.

            :returns: The following four arrays are returned by the analysis:

                - mz_spectra : Array with the static mz values for the spectra.
                - label_spectra : Lable for the spectral mz axis
                - mz_slice : Array of the static mz values for the slices or None if identical to the mz_spectra.
                - label_slice : Lable for the slice mz axis or None if identical to label_spectra.
        """
        num_custom_slice_options = 1
        num_custom_spectrum_options = 0

        # Compute the output
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None

        norm_msidata_shape = analysis_object['norm_msidata'].shape
        valuesX = range(0, norm_msidata_shape[0])
        labelX = 'pixel index X'
        valuesY = range(0, norm_msidata_shape[1])
        labelY = 'pixel index Y'
        if len(norm_msidata_shape) > 3:
            valuesZ = range(0, norm_msidata_shape[2])
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None

        # Both viewer_options point to a data dependency
        if qspectrum_viewer_option >= num_custom_spectrum_options and \
                qslice_viewer_option >= num_custom_slice_options:
            super_qslice_option = qslice_viewer_option-num_custom_slice_options
            super_qspectrum_option = qspectrum_viewer_option-num_custom_spectrum_options
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ= \
                super(omsi_tic_norm, cls).v_qmz(analysis_object,
                                                qslice_viewer_option=super_qslice_option,
                                                qspectrum_viewer_option=super_qspectrum_option)

        # Implement the qmz pattern for all the custom qslice and qspectrum viewer options.
        if qslice_viewer_option == 0 and qspectrum_viewer_option >= num_custom_spectrum_options:
            super_qslice_option = 0
            super_qspectrum_option = qspectrum_viewer_option-num_custom_spectrum_options
            # Ignore the spatial dimensions. We need to use the spatial axes of our norm_msidata
            mz_spectra, label_spectra, mz_slice, label_slice, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_tic_norm, cls).v_qmz(analysis_object,
                                                qslice_viewer_option=super_qslice_option,
                                                qspectrum_viewer_option=super_qspectrum_option)
            label_slice = 'm/z'
            mz_slice = analysis_object['norm_mz'][:]
        elif qspectrum_viewer_option == 0 and qslice_viewer_option >= num_custom_slice_options:
            super_qslice_option = qslice_viewer_option-num_custom_slice_options
            super_qspectrum_option = 0
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_tic_norm, cls).v_qmz(analysis_object,
                                                qslice_viewer_option=super_qslice_option,
                                                qspectrum_viewer_option=super_qspectrum_option)
            label_spectra = 'm/z'
            mz_spectra = analysis_object['norm_mz'][:]
        else:  # qslice_viewer_option == 0 and qspectrum_viewer_option==0
            label_spectra = 'm/z'
            mz_spectra = analysis_object['norm_mz'][:]
            mz_slice = None
            label_slice = None

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """Get a list of strings describing the different default viewer options for the analysis for qspectrum.
           The default implementation tries to take care of handling the spectra retrieval for all the depencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.

            :param analysis_object: The omsi_file_analysis object for which slicing should be performed.  For most cases
                   this is not needed here as the support for slice operations is usually a static decission based
                   on the class type, however, in some cases additional checks may be needed (e.g., ensure that
                   the required data is available).

            :returns: List of strings indicating the different available viewer options. The list should be empty
                    if the analysis does not support qspectrum requests (i.e., v_qspectrum(...) is not available).
        """
        custom_options = ['Tic Norm']
        dependent_options = super(omsi_tic_norm, cls).v_qspectrum_viewer_options(analysis_object)
        spectrum_viewer_options = custom_options + dependent_options
        return spectrum_viewer_options

    @classmethod
    def v_qslice_viewer_options(cls,
                                analysis_object):
        """
        Get a list of strings describing the different default viewer options for the analysis for qslice.
        The default implementation tries to take care of handling the spectra retrieval for all the depencies
        but can naturally not decide how the qspectrum should be handled by a derived class. However, this
        implementation is often called at the end of custom implementations to also allow access to data from
        other dependencies.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.  For
            most cases this is not needed here as the support for slice operations is usually a static
            decision based on the class type, however, in some cases additional checks may be needed (e.g.,
            ensure that the required data is available).

        :returns: List of strings indicating the different available viewer options. The list should be empty if
                  the analysis does not support qslice requests (i.e., v_qslice(...) is not available).
        """
        custom_options = ['Tic Norm']
        dependent_options = super(omsi_tic_norm, cls).v_qslice_viewer_options(analysis_object)
        slice_viewer_options = custom_options + dependent_options
        return slice_viewer_options


############################################################
#  3) Making your analysis self-sufficient   (Recommended) #
############################################################
if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_tic_norm).main()








