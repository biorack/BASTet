"""
Module for performing CX factorization on MSI data.
"""
import numpy as np

from omsi.analysis.base import analysis_base



###############################################################
#  1) Basic integration of your analysis with omsi (Required) #
###############################################################
class omsi_cx(analysis_base):
    """
    Class used to implement CX factorization on MSI data.
    """

    # This internal dict is used to avoid errors due to misinterpretation of the usage of dimensions
    dimension_index = {'imageDim': 0, 'pixelDim': 1}

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""
        super(omsi_cx, self).__init__()
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI matrix to be analyzed',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='rank',
                           help='Rank approximation to be used',
                           dtype=dtypes['int'],
                           default=10,
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='objectiveDim',
                           help='Objective dimension',
                           required=True,
                           default=self.dimension_index['imageDim'],
                           dtype=dtypes['int'],
                           group=groups['settings'])
        self.add_parameter(name='pixelMask',
                           help='A boolean or integer mask indicating the spectra to be used',
                           dtype=dtypes['ndarray'],
                           default=None,
                           required=False,
                           group=groups['settings'])
        self.add_parameter(name='maskIndex',
                           help='If an integer mask is used, which index should we do CX for',
                           dtype=dtypes['int'],
                           default=None,
                           required=False,
                           group=groups['settings'])
        self.data_names = ['infIndices', 'levScores']
        self.analysis_identifier = name_key

    def execute_analysis(self):
        """ EDIT_ME:

            Replace this text with the appropriate documentation for the analysis.
            Describe what your analysis does and how a user can use it. Note, a user will
            call the function execute(...) which takes care of storing parameters, collecting
            execution data etc., so that you only need to implement your analysis, the rest
            is taken care of by analysis_base. omsi uses Sphynx syntax for the
            documentation.

            Keyword Arguments:

            :param mydata: ...
            :type mydata: ...

        """
        # getting the values into local variables
        msidata = np.copy(self['msidata'][:])  # Copy all MSI data
        original_shape = msidata.shape
        rank = self['rank']
        objective_dimensions = self['objectiveDim']
        mask_index = self['maskIndex']
        current_pixel_mask = self['pixelMask']
        if mask_index is not None and current_pixel_mask is not None:
            temp_mask = np.zeros(shape=current_pixel_mask.shape, dtype=bool)
            temp_mask[current_pixel_mask == mask_index] = True
            current_pixel_mask = temp_mask

        # Mask the data if requested
        if current_pixel_mask is not None:
            msidata = msidata[current_pixel_mask, :]

        # Determine the number of pixels and bins
        num_bins = msidata.shape[-1]
        num_pixels = msidata.size / num_bins  # The last data dimension is assumed to contain the spectra

        # Reshape the data to a 2D matrix for CX
        msidata = msidata.reshape(num_pixels, num_bins).transpose()

        # Compute the CX decomposition
        leverage_scores = self.comp_lev_exact(msidata, rank, objective_dimensions)
        informative_indices = leverage_scores.argsort()[::-1]

        # Fill in the values that we removed with the mask and reshape the data for the output
        if current_pixel_mask is not None:
            if self['objectiveDim'] == self.dimension_index['pixelDim']:
                out_leverage_scores = np.zeros(shape=original_shape[0:-1], dtype=leverage_scores.dtype)
                out_leverage_scores[:] = -1
                out_leverage_scores[current_pixel_mask] = leverage_scores
                out_informative_indices = np.zeros(shape=original_shape[0:-1], dtype=informative_indices.dtype)
                out_informative_indices[:] = -1
                out_informative_indices[current_pixel_mask] = informative_indices
            elif self['objectiveDim'] == self.dimension_index['imageDim']:
                out_leverage_scores = leverage_scores
                out_informative_indices = informative_indices
        # Reshape the data for the output
        else:
            # If the leverage scores are computed for pixels, then, convert back to image space
            if self['objectiveDim'] == self.dimension_index['pixelDim']:
                leverage_scores = leverage_scores.reshape(original_shape[0:-1])
                informative_indices = informative_indices.reshape(original_shape[0:-1])
            out_leverage_scores = leverage_scores
            out_informative_indices = informative_indices

        # Safe the output results
        return out_informative_indices, out_leverage_scores

    @classmethod
    def comp_lev_exact(cls,
                       A,
                       k,
                       axis):
        """
        This function computes the column or row leverage scores of the input matrix.


        :param A: n-by-d matrix
        :param k: rank parameter, k <= min(n,d)
        :param axis: 0: compute row leverage scores; 1: compute column leverage scores.

        :returns: 1D array of leverage scores. If axis = 0, the length of lev is n.  otherwise, the length of lev is d.
        """
        U, D, V = np.linalg.svd(A, full_matrices=False)

        if axis == 0:
            lev = np.sum(U[:, :k]**2, axis=1)
        else:
            lev = np.sum(V[:k, :]**2, axis=0)

        return lev

    ###############################################################
    #  2) Integrating your analysis with the OpenMSI              #
    #     web-based viewer (Recommended)                          #
    ###############################################################
    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """
        Get 3D analysis dataset for which z-slices should be extracted for presentation in the OMSI viewer

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

        """EDIT_ME Specify the number of custom  viewer_options you are going to provide for qslice"""
        current_objective_dimension = analysis_object['objectiveDim'][0]
        if current_objective_dimension == cls.dimension_index['imageDim']:
            num_custom_viewer_options = 1
        else:
            num_custom_viewer_options = 2

        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_cx, cls).v_qslice(analysis_object,
                                                z,
                                                viewer_option=viewer_option-num_custom_viewer_options)
        elif viewer_option == 0 and current_objective_dimension == cls.dimension_index['imageDim']:
            informative_indices = analysis_object['infIndices'][z_select]
            cx_analysis_object = omsi_cx()
            cx_analysis_object.read_from_omsi_file(analysis_object=analysis_object,
                                                   load_data=False,
                                                   load_parameters=False,
                                                   load_runtime_data=False)
            return cx_analysis_object['msidata'][:, :, informative_indices]
        elif viewer_option == 0 and current_objective_dimension == cls.dimension_index['pixelDim']:
            return analysis_object['levScores'][:]
        elif viewer_option == 1 and current_objective_dimension == cls.dimension_index['pixelDim']:
            return analysis_object['infIndices'][:]

        """EDIT_ME

           Define your custom qslice viewer options. Here you need to handle all the different
           behaviors that are custom to your analysis. Below a simple example.

           if viewer_option == 0 :
               dataset = anaObj[ 'my_output_data' ] #This is e.g, an output dataset of your analysis
               return dataset[ : , :, z_select ]
           elif viewer_option == 1 :
               ...
        """
        return None

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """
        Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer

        Developer Note: h5py currently supports only a single index list. If the user provides an index-list for both
                       x and y, then we need to construct the proper merged list and load the data manually, or if
                       the data is small enough, one can load the full data into a numpy array which supports
                       mulitple lists in the selection.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param x: x selection string
        :param y: y selection string
        :param viewer_option: If multiple default viewer behaviors are available for a given analysis then
            this option is used to switch between them.

        :returns: The following two elemnts are expected to be returned by this function :

            1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The mass (m/z) axis must be the last \
                axis. For index selection x=1,y=1 a 1D array is usually expected. For indexList selections x=[0]&y=[1] \
                usually a 2D array is expected. For ragne selections x=0:1&y=1:2 we one usually expects a 3D array.
            2) None in case that the spectra axis returned by v_qmz are valid for the returned spectrum. Otherwise, \
                return a 1D numpy array with the m/z values for the spectrum (i.e., if custom m/z values are needed \
                for interpretation of the returned spectrum).This may be needed, e.g., in cases where a per-spectrum \
                peak analysis is performed and the peaks for each spectrum appear at different m/z values.
        """

        # Convert the x,y selection to a python selection
        from omsi.shared.data_selection import selection_string_to_object
        x_select = selection_string_to_object(x)  # Convert the selection string to a python selection
        y_select = selection_string_to_object(y)  # Convert the selection string to a python selection

        """EDIT_ME Specify the number of custom viewer_options you are going to provide for qslice"""
        num_custom_viewer_options = 0
        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_cx, cls).v_qspectrum(analysis_object,
                                                   x,
                                                   y,
                                                   viewer_option=num_custom_viewer_options)

        """EDIT_ME

           Define your custom qspectrum viewer options. Here you need to handle all the different
           behaviors that are custom to your analysis. Note, this function is expected to return
           two object: i) The data for the spectrum and ii) the m/z axis information for the spectrum
           or None, in case that the m/z data is identical to what the v_qmz function returns.
           Below a simple example.

           if viewer_option == 0 :
               dataset = anaObj[ 'my_output_data' ] #This is e.g, an output dataset of your analysis
               data = dataset[ x_select , y_select, : ]
               return data, None
           elif viewer_option == 1 :
               ...
        """
        return None, None

    @classmethod
    def v_qmz(cls,
              analysis_object,
              qslice_viewer_option=0,
              qspectrum_viewer_option=0):
        """
        Get the mz axes for the analysis

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param qslice_viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them for the qslice URL pattern.
        :param qspectrum_viewer_option: If multiple default viewer behaviors are available for a
            given analysis then this option is used to switch between them for the qspectrum URL pattern.

        :returns: The following four arrays are returned by the analysis:

            - mz_spectra : Array with the static mz values for the spectra.
            - label_spectra : Lable for the spectral mz axis
            - mz_slice : Array of the static mz values for the slices or None if identical to the mz_spectra.
            - label_slice : Lable for the slice mz axis or None if identical to label_spectra.
        """

        """EDIT_ME: Define the number of custom viewer options for qslice and qspectrum."""
        current_objective_dimension = analysis_object['objectiveDim'][0]
        if current_objective_dimension == cls.dimension_index['imageDim']:
            num_custom_slice_options = 1
        else:
            num_custom_slice_options = 2

        num_custom_spectrum_options = 0

        # Compute the output
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None
        cx_analysis_object = omsi_cx()
        cx_analysis_object.read_from_omsi_file(analysis_object=analysis_object,
                                               load_data=False,
                                               load_parameters=False,
                                               load_runtime_data=False)
        cx_msidata_shape = cx_analysis_object['msidata'].shape
        valuesX = range(0, cx_msidata_shape[0])
        labelX = 'pixel index X'
        valuesY = range(0, cx_msidata_shape[1])
        labelY = 'pixel index Y'
        if len(cx_msidata_shape) > 3:
            valuesZ = range(0, cx_msidata_shape[2])
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None

        # Both viewer_options point to a data dependency
        if qspectrum_viewer_option >= num_custom_spectrum_options and \
                qslice_viewer_option >= num_custom_slice_options:
            """EDIT_ME Replace the omsi_cx class name with your class name"""
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_cx, cls).v_qmz(analysis_object,
                                          qslice_viewer_option=qslice_viewer_option-num_custom_slice_options,
                                          qspectrum_viewer_option=qspectrum_viewer_option-num_custom_spectrum_options)

        # Implement the qmz pattern for all the custom qslice and qspectrum viewer options.
        if qslice_viewer_option == 0 and current_objective_dimension == cls.dimension_index['imageDim']:
            # Ignore the spatial dimensions as this is a custom qslice option
            mz_spectra, label_spectra, mz_slice, label_slice, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_cx, cls).v_qmz(analysis_object,
                                          qslice_viewer_option=0,
                                          qspectrum_viewer_option=qspectrum_viewer_option-num_custom_spectrum_options)
            label_slice = 'Informative Images'
            mz_slice = np.arange(analysis_object['infIndices'].shape[0])

        if (qslice_viewer_option == 0 or qslice_viewer_option == 1) and \
                current_objective_dimension == cls.dimension_index['pixelDim']:
            # Ignore the spatial dimensions as this is a custom qslice option
            mz_spectra, label_spectra, mz_slice, label_slice, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_cx, cls).v_qmz(analysis_object,
                                          qslice_viewer_option=0,
                                          qspectrum_viewer_option=qspectrum_viewer_option-num_custom_spectrum_options)
            if qslice_viewer_option == 0:
                label_slice = 'Pixel Leverage Scores'
                mz_slice = np.arange(1)
            elif qslice_viewer_option == 1:
                label_slice = 'Pixel Rank'
                mz_slice = np.arange(1)

        #
        # elif viewer_option == 0 and current_objective_dimension == cls.dimension_index['imageDim']:
        #     informative_indices = analysis_object['infIndices'][z_select]
        #     cx_analysis_object = omsi_cx()
        #     cx_analysis_object.read_from_omsi_file(analysis_object=analysis_object,
        #                                            load_data=False,
        #                                            load_parameters=False,
        #                                            load_runtime_data=False)
        #     return cx_analysis_object['msidata'][:, :, informative_indices]
        # elif viewer_option == 0 and current_objective_dimension == cls.dimension_index['pixelDim']:
        #     return analysis_object['levScores'][:]
        # elif viewer_option == 1 and current_objective_dimension == cls.dimension_index['pixelDim']:
        #     return analysis_object['infIndices'][:]

        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """
        Get a list of strings describing the different default viewer options for the analysis for qspectrum.
        The default implementation tries to take care of handling the spectra retrieval for all the depencies
        but can naturally not decide how the qspectrum should be handled by a derived class. However, this
        implementation is often called at the end of custom implementations to also allow access to data from
        other dependencies.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.  For most cases
            this is not needed here as the support for slice operations is usually a static decission based on the
            class type, however, in some cases additional checks may be needed (e.g., ensure that the
            required data is available).

        :returns: List of strings indicating the different available viewer options. The list should be empty
            if the analysis does not support qspectrum requests (i.e., v_qspectrum(...) is not available).
        """

        """EDIT_ME Define a list of custom viewer_options are supported. E.g:

           custom_options = ['Peak cube']
        """
        custom_options = []
        dependent_options = super(omsi_cx, cls).v_qspectrum_viewer_options(analysis_object)
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
            most cases this is not needed here as the support for slice operations is usually a static decision
            based on the class type, however, in some cases additional checks may be needed (e.g., ensure that
            the required data is available).

        :returns: List of strings indicating the different available viewer options. The list should be empty
            if the analysis does not support qslice requests (i.e., v_qslice(...) is not available).
        """

        # Define a list of custom viewer_options are supported.
        if analysis_object['objectiveDim'][0] == cls.dimension_index['imageDim']:
            custom_options = ['Informative Images']
        else:
            custom_options = ['Pixel Leverage Scores', 'Pixel Rank']

        dependent_options = super(omsi_cx, cls).v_qslice_viewer_options(analysis_object)
        slice_viewer_options = custom_options + dependent_options
        return slice_viewer_options


############################################################
#  3) Making your analysis self-sufficient   (Recommended) #
############################################################
if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    cl_analysis_driver(analysis_class=omsi_cx).main()


"""Test script:

from omsi.dataformat.omsi_file import *
from omsi.analysis.multivariate_stats.omsi_cx import *
inputFile = 'test_cx.h5'
f = omsi_file( inputFile , 'a' )
e = f.get_experiment(0)
a = e.get_analysis(0)
d = a['peak_cube']
ocxi = omsi_cx(nameKey='testCX_Images')
ocxi.execute( msidata=d , rank=10, objectiveDim=omsi_cx.dimension_index['imageDim'] )
ocxp = omsi_cx(nameKey='testCX_Pixel')
ocxp.execute( msidata=d , rank=10, objectiveDim=omsi_cx.dimension_index['pixelDim'] )
e.create_analysis( ocxi , flushIO=True)
e.create_analysis( ocxp , flushIO=True)
f.close_file()
"""







