from omsi.analysis.omsi_analysis_base import omsi_analysis_base
from omsi.analysis.omsi_analysis_data import omsi_analysis_data
from omsi.shared.omsi_dependency import *
import numpy as np
from tempfile import TemporaryFile


###############################################################
#  1) Basic integration of your analysis with omsi (Required) #
###############################################################
class omsi_tic_norm(omsi_analysis_base):

    """Template intended to help with the development of new analysis classes.
    
       Search for EDIT_ME to find locations that need to be changed.
       
       EDIT_ME: Replace this text with the appropriate documentation for the analysis.
    
    """
    def __init__(self, nameKey="undefined"):
        """Initalize the basic data members"""

        super(omsi_tic_norm, self).__init__()
        self.parameter_names = ['msidata', 'mzdata', 'maxCount', 'mzTol', 'infIons']
        self.data_names = ['norm_msidata', 'norm_mz']
        self.analysis_identifier = nameKey

    def execute_analysis(self):
        """
        Normalize the data based on the total intensity of a spectrum or the
        intentsity of a select set of ions:
           
        Keyword Arguments:

        :param msidata: The input MSI dataset
        :type msidata: h5py.dataset or numpu array (3D)
        :param mzdata: The mz axsis do the dataset
        :param maxCount: ...
        :param mzTol: ...
        :param infIons: List of informative ions
           
        """
        #Setting Default Values and retrieving parameter data
        nx, ny, nz = self['msidata'].shape
        if not self['maxCount']:
            self['maxCount'] = 0
        if not self['mzTol']:
            self['mzTol'] = 0.1
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

        outputFilename = TemporaryFile() # Create a temporary file to compute the normalization out-of-core
        idx = np.linspace(0, nx, 20).astype(int)  # Divide the data into chunks in x to process the data in blocks

        # Compute the normalization factors and column maxs first (computed one-block-at-a-time)
        tic_norm_factors = np.zeros(shape=(nx,ny), dtype='float')
        msi_spectrum_maxs = np.zeros(shape=(nx,ny), dtype=self['msidata'].dtype)
        for i in range(len(idx)-1): # Process the data one-block-at-a-time
            current_data = self['msidata'][idx[i]:idx[i+1], :, idx_mz]
            tic_norm_factors[idx[i]:idx[i+1], :] = np.sum(current_data, 2)
            msi_spectrum_maxs[idx[i]:idx[i+1], :] = np.amax(current_data[:,:, idx_mz], 2)
            del current_data
        non_zero_tic = tic_norm_factors > 0
        mean_tic_norm = float(tic_norm_factors[non_zero_tic].mean())
        tic_norm_factors[non_zero_tic] = 1.0 / (tic_norm_factors[non_zero_tic] / mean_tic_norm)

        # Normalize the data one-block-at-a-time
        for i in range(len(idx)-1):
            num_spectra_x = idx[i+1]-idx[i]
            current_data = self['msidata'][idx[i]:idx[i+1], :, :].reshape(num_spectra_x*ny, nz)   # Load the data for the current block
            tip = tic_norm_factors[idx[i]:idx[i+1], :].reshape(num_spectra_x*ny)
            mip = msi_spectrum_maxs[idx[i]:idx[i+1], :].reshape(num_spectra_x*ny)
            idx_thresh = np.where(np.logical_and(mip >= self['maxCount'], tip > 0))
            current_out_norm = np.zeros(current_data.shape)
            if len(idx_thresh[0]) > 0:
                  current_out_norm[idx_thresh[0], :] = np.multiply(current_data[idx_thresh[0], :],
                                                                   tip.astype(float)[idx_thresh[0]][:, np.newaxis])
                  # Write the data to the temporary file
                  ofp = np.memmap(outputFilename, dtype='float', mode='w+', shape=(nx, ny, nz))
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
        #         ofp = np.memmap(outputFilename, dtype='float', mode='w+', shape=(nx, ny, nz))
        #         ofp[idx[i]:idx[i+1], :, :] = current_out_norm.reshape(num_spectra, ny, nz)
        #         del ofp

        # Save the data
        self['norm_msidata'] = np.memmap(outputFilename, dtype='float', mode='w+', shape=(nx, ny, nz))
        self['norm_mz'] = mzdata
        outputFilename.close()

    ###############################################################
    #  2) Integrating your analysis with the OpenMSI              #
    #     web-based viewer (Recommended)                          #
    ###############################################################
    @classmethod
    def v_qslice(cls, anaObj, z, viewerOption=0):
        """Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer

           :param anaObj: The omsi_file_analysis object for which slicing should be performed
           :param z: Selection string indicting which z values should be selected.
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis
                  then this option is used to switch between them.

           :returns: numpy array with the data to be displayed in the image slice viewer. Slicing will
                     be performed typically like [:,:,zmin:zmax].

        """
        #Convert the z selection to a python selection
        from omsi.shared.omsi_data_selection import selection_string_to_object
        zselect = selection_string_to_object(z)  # Convert the selection string to a python selection

        numCustomViewerOptions = 1

        #Expose the qslice viewer functionality of any data dependencies
        if viewerOption >= numCustomViewerOptions:
            return super(omsi_tic_norm, cls).v_qslice(anaObj,
                                                      z,
                                                      viewerOption=viewerOption-numCustomViewerOptions)
        elif viewerOption == 0:
            dataset = anaObj['norm_msidata']
            return dataset[:, :, zselect]
        return None

    @classmethod
    def v_qspectrum(cls, anaObj, x, y, viewerOption=0):
        """Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer

           Developer Note: h5py currently supports only a single index list. If the user provides an index-list for both
                           x and y, then we need to construct the proper merged list and load the data manually, or if
                           the data is small enough, one can load the full data into a numpy array which supports
                           mulitple lists in the selection.

           :param anaObj: The omsi_file_analysis object for which slicing should be performed
           :param x: x selection string
           :param y: y selection string
           :param viewerOption: If multiple default viewer behaviors are available for a given analysis then
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

        #Convert the x,y selection to a python selection
        from omsi.shared.omsi_data_selection import selection_string_to_object
        xselect = selection_string_to_object(x)  # Convert the selection string to a python selection
        yselect = selection_string_to_object(y)  # Convert the selection string to a python selection

        numCustomViewerOptions = 1
        #Expose the qslice viewer functionality of any data dependencies
        if viewerOption >= numCustomViewerOptions:
            return super(omsi_tic_norm, cls).v_qspectrum(anaObj, x, y, viewerOption=viewerOption-numCustomViewerOptions)
        elif viewerOption == 0:
            dataset = anaObj['norm_msidata']
            # h5py supports only one indexlist as part of a selection. We have to work-around this
            # by retrieving the spectra we want one-at-a-time if the user provides two index lists
            if isinstance(xselect, list) and isinstance(yselect, list):
                xsize = len(xselect)
                ysize = len(yselect)
                if xsize != ysize:
                    raise ValueError("x and y selection size do not match")
                data = np.zeros(shape=(xsize,dataset.shape[2]), dtype=dataset.dtype)
                for i in xrange(0, xsize):
                    data[i, :] = dataset[xselect[i], yselect[i], :]
            else:
                data = dataset[xselect, yselect, :]
            return data, None
        return None, None

    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0):
        """ Get the mz axes for the analysis

            :param anaObj: The omsi_file_analysis object for which slicing should be performed
            :param qslice_viewerOption: If multiple default viewer behaviors are available for a given analysis then
                                        this option is used to switch between them for the qslice URL pattern.
            :param qspectrum_viewerOption: If multiple default viewer behaviors are available for a given analysis
                                     then this option is used to switch between them for the qspectrum URL pattern.

            :returns: The following four arrays are returned by the analysis:

                - mzSpectra : Array with the static mz values for the spectra.
                - labelSpectra : Lable for the spectral mz axis
                - mzSlice : Array of the static mz values for the slices or None if identical to the mzSpectra.
                - labelSlice : Lable for the slice mz axis or None if identical to labelSpectra.
        """
        numCustomSliceViewerOptions = 1
        numCustomSpectrumViewerOptions = 0

        #Compute the output
        mzSpectra = None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        #Both viewerOptions point to a data dependency
        if qspectrum_viewerOption >= numCustomSpectrumViewerOptions and \
                qslice_viewerOption >= numCustomSliceViewerOptions:
            super_qslice_option = qslice_viewerOption-numCustomSliceViewerOptions
            super_qspectrum_option = qspectrum_viewerOption-numCustomSpectrumViewerOptions
            mzSpectra, labelSpectra, mzSlice, labelSlice = \
                super(omsi_tic_norm, cls).v_qmz(anaObj,
                                                qslice_viewerOption=super_qslice_option,
                                                qspectrum_viewerOption=super_qspectrum_option)

        #Implement the qmz pattern for all the custom qslice and qspectrum viewer options.
        if qslice_viewerOption == 0 and qspectrum_viewerOption >= numCustomSpectrumViewerOptions:
            super_qslice_option = 0
            super_qspectrum_option = qspectrum_viewerOption-numCustomSpectrumViewerOptions
            mzSpectra, labelSpectra, mzSlice, labelSlice = \
                super(omsi_tic_norm, cls).v_qmz(anaObj,
                                                qslice_viewerOption=super_qslice_option,
                                                qspectrum_viewerOption=super_qspectrum_option)
            labelSlice = 'm/z'
            mzSlice = anaObj['norm_mz'][:]
        elif qspectrum_viewerOption == 0 and qslice_viewerOption >= numCustomSliceViewerOptions:
            super_qslice_option = qslice_viewerOption-numCustomSliceViewerOptions
            super_qspectrum_option = 0
            mzSpectra, labelSpectra, mzSlice, labelSlice = \
                super(omsi_tic_norm, cls).v_qmz(anaObj,
                                                qslice_viewerOption=super_qslice_option,
                                                qspectrum_viewerOption=super_qspectrum_option)
            labelSpectra = 'm/z'
            mzSpectra = anaObj['norm_mz'][:]
        else:  # qslice_viewerOption == 0 and qspectrum_viewerOption==0
            labelSpectra = 'm/z'
            mzSpectra = anaObj['norm_mz'][:]
            mzSlice = None
            labelSlice = None

        return mzSpectra, labelSpectra, mzSlice, labelSlice

    @classmethod
    def v_qspectrum_viewerOptions(cls, anaObj):
        """Get a list of strings describing the different default viewer options for the analysis for qspectrum.
           The default implementation tries to take care of handling the spectra retrieval for all the depencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.

            :param anaObj: The omsi_file_analysis object for which slicing should be performed.  For most cases
                   this is not needed here as the support for slice operations is usually a static decission based
                   on the class type, however, in some cases additional checks may be needed (e.g., ensure that
                   the required data is available).

            :returns: List of strings indicating the different available viewer options. The list should be empty
                    if the analysis does not support qspectrum requests (i.e., v_qspectrum(...) is not available).
        """
        customOptions = ['Tic Norm']
        dependentOptions = super(omsi_tic_norm, cls).v_qspectrum_viewerOptions(anaObj)
        re = customOptions + dependentOptions
        return re

    @classmethod
    def v_qslice_viewerOptions(cls, anaObj):
        """Get a list of strings describing the different default viewer options for the analysis for qslice.
           The default implementation tries to take care of handling the spectra retrieval for all the depencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.

            :param anaObj: The omsi_file_analysis object for which slicing should be performed.  For most cases this
                   is not needed here as the support for slice operations is usually a static decission based on the
                   class type, however, in some cases additional checks may be needed (e.g., ensure that the
                   required data is available).

            :returns: List of strings indicating the different available viewer options. The list should be empty if
                      the analysis does not support qslice requests (i.e., v_qslice(...) is not available).
        """
        customOptions = ['Tic Norm']
        dependentOptions = super(omsi_tic_norm, cls).v_qslice_viewerOptions(anaObj)
        re = customOptions + dependentOptions
        return re


############################################################
#  3) Making your analysis self-sufficient   (Recommended) #
############################################################
def main(argv=None):
    """EDIT_ME : Optional
    
       Implement this function to enable a user to use your module also as a stand-alone script.
       Remember, you should always call execute(...) to run your analysis and NOT execute_analysis(...)
       
    """
    #Get the input arguments
    import sys
    from sys import argv, exit
    if argv is None:
        argv = sys.argv


if __name__ == "__main__":
    main()








