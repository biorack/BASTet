"""
File defining classes and methods for randomized PCA analysis of OMSI images via scikit
"""
import numpy as np
import time
from sklearn import decomposition
from omsi.analysis.base import analysis_base

###############################################################
#  1) Basic integration of your analysis with omsi (Required) #
###############################################################
class omsi_rpca(analysis_base):
    """
    Class representing Randomized PCA analysis
    REF: http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.RandomizedPCA.html
    REF: http://arxiv.org/pdf/0909.4061v2.pdf
    """

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""

        super(omsi_rpca, self).__init__()
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI matrix to be analyzed',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])
        self.add_parameter(name='numComponents',
                           help='The number of output components to compute.',
                           default=20,
                           dtype=int,
                           required=True,
                           group=groups['settings'])
        self.add_parameter(name='iteratedPower',
                           help="""Number of iterations for the power method.  3 by default.
                                   See sklearn.decomposition.RandomizedPCA()""",
                           default=3,
                           dtype=int,
                           required=True,
                           group=groups['settings'])

        # should "NoneType" be added to permissible dtypes?
        # how is the dtype parameters that could have multiple possible types specified?
        self.add_parameter(name='randomState',
                           help="""Pseudo Random Number generator seed control. If None, use the numpy.random singleton.
                                See sklearn.decomposition.RandomizedPCA()""",
                           default=1,
                           dtype=['int'],  # should be ['int', 'NoneType']
                           required=True,
                           group=groups['settings'])


        self.data_names = ['newImageCube', 'components', 'explainedVariance', 'mean', 'analysisTime']
        self.analysis_identifier = name_key

    def execute_analysis(self):
        """
        function to do randomized PCA on openMSI data via sklearn.decomposition.RandomizedPCA()
        The "copy" parameter of RandomizedPCA() is not supported (is always False)
        The "whiten" parameter of RandomizedPCA() is not currently supported (assumed False)
        """
        # extract input parameters

        start = time.time()

        msidata = self['msidata']
        n_components = self['numComponents']
        iterated_power = self['iteratedPower']
        random_state = self['randomState']

        nx, ny, nmz = msidata.shape

        # Randomized PCA
        # # reshape msidata from 3D (x by y by mz) to 2D (xy by mz)
        flatdata = np.array(
                            [np.array(msidata[:, :, i]).flatten()
                             for i in range(nmz)
                            ]).T

        # # do randomized PCA
        pca = decomposition.RandomizedPCA(n_components=n_components,
                                          iterated_power=iterated_power,
                                          random_state=random_state)
        pca.fit(flatdata, n_components)

        # # make new image and reshape to expected size
        newImageCubeFlat = pca.transform(flatdata)
        newImageCube = newImageCubeFlat.reshape(nx, ny, n_components)

        # # return other pca data (redundant with returning full pca object but not sure which way is best right now)
        components = pca.components_
        explainedVariance = pca.explained_variance_ratio_
        mean = pca.mean_

        # # return analysis time
        stop = time.time()
        analysisTime = stop-start

        # # return results
        return newImageCube, components, explainedVariance, mean, analysisTime

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
        :param viewer_option: If multiple default viewer behaviors are available for
            a given analysis then this option is used to switch between them.

        :returns: numpy array with the data to be displayed in the image slice viewer.
            Slicing will be performed typically like [:,:,zmin:zmax].

        """

        # Convert the z selection to a python selection
        from omsi.shared.data_selection import selection_string_to_object
        zselect = selection_string_to_object(z)  # Convert the selection string to a python selection

        """EDIT_ME Specify the number of custom viewer_options you are going to provide for qslice"""
        num_custom_viewer_options = 0

        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_rpca, cls).v_qslice(analysis_object,
                                                               z,
                                                               viewer_option=viewer_option-num_custom_viewer_options)

        """
        EDIT_ME

        Define your custom qslice viewer options. Here you need to handle all the different
        behaviors that are custom to your analysis. Below a simple example.

        if viewer_option == 0 :
           dataset = anaObj[ 'my_output_data' ] #This is e.g, an output dataset of your analysis
           return dataset[ : , :, zselect ]
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
        :param viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them.

        :returns: The following two elemnts are expected to be returned by this function :

            1) 1D, 2D or 3D numpy array of the requested spectra. NOTE: The mass (m/z) axis must be the last axis. \
               For index selection x=1,y=1 a 1D array is usually expected. For indexList selections x=[0]&y=[1] \
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

        """
        EDIT_ME

        Specify the number of custom viewer_options you are going to provide for qslice
        """
        num_custom_viewer_options = 0

        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            """
            EDIT_ME

            Replace omsi_rpca with your classname
            """
            return super(omsi_rpca, cls).v_qspectrum(analysis_object,
                                                                  x,
                                                                  y,
                                                                  viewer_option=viewer_option-num_custom_viewer_options)

        """
        EDIT_ME

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
    def v_qmz(cls, analysis_object, qslice_viewer_option=0, qspectrum_viewer_option=0):
        """
        Get the mz axes for the analysis

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param qslice_viewer_option: If multiple default viewer behaviors are available for a
            given analysis then this option is used to switch between them for the qslice URL pattern.
        :param qspectrum_viewer_option: If multiple default viewer behaviors are available for a
            given analysis then this option is used to switch between them for the qspectrum URL pattern.

        :returns: The following four arrays are returned by the analysis:

            - mz_spectra : Array with the static mz values for the spectra.
            - label_spectra : Lable for the spectral mz axis
            - mz_slice : Array of the static mz values for the slices or None if identical to the mz_spectra.
            - label_slice : Lable for the slice mz axis or None if identical to label_spectra.
        """

        """
        EDIT_ME

        Define the number of custom viewer options for qslice and qspectrum.
        """
        num_custom_slice_viewer_options = 0
        num_custom_spectrum_viewer_options = 0

        # Compute the output
        mz_spectra = None
        label_spectra = None
        mz_slice = None
        label_slice = None
        valuesX = None
        labelX = None
        valuesY = None
        labelY = None
        valuesZ = None
        labelZ = None
        # Both viewer_options point to a data dependency
        if qspectrum_viewer_option >= num_custom_spectrum_viewer_options \
                and qslice_viewer_option >= num_custom_slice_viewer_options:
            """EDIT_ME Replace the omsi_rpca class name with your class name"""
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_rpca, cls)\
                    .v_qmz(analysis_object,
                           qslice_viewer_option=qslice_viewer_option-num_custom_slice_viewer_options,
                           qspectrum_viewer_option=qspectrum_viewer_option-num_custom_spectrum_viewer_options)

        """
        EDIT_ME

        Implement the qmz pattern for all the custom qslice and qspectrum viewer options. E.g:

        if qspectrum_viewer_option == 0 and qslice_viewer_option==0:
            mz_spectra =  anaObj[ 'peak_mz' ][:]
            label_spectra = "m/z"
            mz_slice  = None
            label_slice = None
        """
        return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls, analysis_object):
        """Get a list of strings describing the different default viewer options for the analysis for qspectrum.
           The default implementation tries to take care of handling the spectra retrieval for all the dependencies
           but can naturally not decide how the qspectrum should be handled by a derived class. However, this
           implementation is often called at the end of custom implementations to also allow access to data from
           other dependencies.

            :param analysis_object: The omsi_file_analysis object for which slicing should be performed.
                For most cases this is not needed here as the support for slice operations is usually a
                static decission based on the class type, however, in some cases additional checks
                may be needed (e.g., ensure that the required data is available).

            :returns: List of strings indicating the different available viewer options. The list should
                be empty if the analysis does not support qspectrum requests
                (i.e., v_qspectrum(...) is not available).
        """

        """
        EDIT_ME

        Define a list of custom viewer_options are supported. E.g:

        custom_options = ['Peak cube']
        """
        custom_options = []

        """
        EDIT_ME

        Change the omsi_rpca class-name to your class. If you did a
        replace all, then this should be done already.
        """
        dependent_options = super(omsi_rpca, cls).v_qspectrum_viewer_options(analysis_object)
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

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.
            For most cases this is not needed here as the support for slice operations is usually a
            static decision based on the class type, however, in some cases additional checks may be
            needed (e.g., ensure that the required data is available).

        :returns: List of strings indicating the different available viewer options. The list should be
            empty if the analysis does not support qslice requests (i.e., v_qslice(...) is not available).
        """

        """
        EDIT_ME

        Define a list of custom viewer_options are supported. E.g:

        custom_options = ['Peak cube']
        """
        custom_options = []

        """
        EDIT_ME

        Change the omsi_rpca class-name to your class.  If you did
        a replace all, then this should be done already.
        """
        dependent_options = super(omsi_rpca, cls).v_qslice_viewer_options(analysis_object)
        slice_viewer_options = custom_options + dependent_options
        return slice_viewer_options


############################################################
#  3) Making your analysis self-sufficient   (Recommended) #
############################################################
if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    """
    EDIT_ME

    Simply replace the omsi_rpca class name with your class name
    """
    cl_analysis_driver(analysis_class=omsi_rpca).main()


