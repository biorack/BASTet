from omsi.analysis.omsi_analysis_base import omsi_analysis_base
import numpy as np

###############################################################
#  1) Basic integration of your analysis with omsi (Required) #
###############################################################
class omsi_mz_upsample(omsi_analysis_base):
    """
    Class representing upsampling the m/z dimension of MSI data by *interpolation*.
    """

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""

        super(omsi_mz_upsample, self).__init__()
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()

        self.add_parameter(name='msidata',
                           help='The MSI matrix to be analyzed',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])

        self.add_parameter(name='mzdata',
                           help='The current m/z values for the spectra of the MSI dataset',
                           dtype=dtypes['ndarray'],
                           required=True,
                           group=groups['input'])

        self.add_parameter(name='new_mzdata',
                           help='The desired new m/z values for the spectra of the MSI dataset',
                           dtype=dtypes['ndarray'],
                           required=False,
                           default=np.array([]),
                           group=groups['input'])

        self.add_parameter(name='new_spacing',
                           help='The desired spacing in Daltons between new equally-spaced m/z datapoints',
                           dtype=float,
                           required=False,
                           default=0,
                           group=groups['input'])

        self.add_parameter(name='method',
                           help='method used by griddata for interpolation, either "linear", "nearest", or "cubic"',
                           dtype=str,
                           required=True,
                           default='linear',
                           group=groups['settings'])

        self.data_names = ['new_msidata', 'new_mzdata']
        self.analysis_identifier = name_key

    def execute_analysis(self):
        """
        Input:   msidata     the original msidata cube [ndarray with shape X by Y by M]
                    mzdata      the original mz data vector  [ndarray with shape M]
                    new_mzdata  the desired mz data vector after re-binning [ndarray with shape N, optional]
                                    incompatible with new spacing
                                    this must be used for non-uniform rebinning
                    new_spacing the desired spacing in Daltons after re-binning [float, optional]
                                    incompatible with new_mzdata.  If this option is picked, the
                                    new_mzdata is calculated as np.arange(min(mzdata), max(mzdata), new_spacing)
                    method      the method used by griddata for interpolation
                                        [{'linear', 'nearest', 'cubic'}, optional]
                                        ***NOTE: 'cubic' will be unusably slow for most datasets

        Output:     new_msidata a new datacube [ndarray with shape X by Y by N]
                    new_mzdata  the new m/z data vector [ndarray with shape N]"""


        #import required function
        from scipy.interpolate import griddata

        #unpack variables
        msidata = self['msidata']
        mzdata = self['mzdata']
        new_mzdata = self['new_mzdata']
        new_spacing = self['new_spacing']
        method = self['method']

        # decide whether to use new_mzdata or new_spacing
        if new_mzdata.shape == (0,) and new_spacing == 0:
            print "at least one of new_mzdata and new_spacing must be specified"
            raise Exception
        elif new_mzdata.shape != (0,) and new_spacing != 0:
            print "Cannot define both new_mzdata and new_spacing"
            raise Exception
        elif new_mzdata.shape == (0,):       # if using new_spacing, form new_mzdata vector
            minmz = np.min(mzdata)
            maxmz = np.max(mzdata)
            new_mzdata = np.arange(minmz, maxmz, new_spacing)

        # ensure new m/z spacing is tighter than old m/z spacing everwhere
        # if not, then this function will downsample by *interpolation* rather than histogramming which will result in
        # reduction of data quality
        mindiff_old = np.min(np.diff(mzdata))
        maxdiff_new = np.max(np.diff(new_mzdata))
        if mindiff_old < maxdiff_new:
            print """This function is meant only for upsampling, not downsampling, but your new m/z axis
                        is more widely spaced than the existing one at at least some m/z regions"""
            raise Exception

        # do upsampling
        ## !this could be improved by multithreading or multiprocessing!

        ## initialize new datacube
        nx = msidata.shape[0]
        ny = msidata.shape[1]
        nmz = len(new_mzdata)

        new_msidata = np.zeros([nx, ny, nmz])

        ## interpolation loop
        for ix in range(nx):
            for iy in range(ny):
                new_msidata[ix, iy, :] = \
                    griddata(points=mzdata, values=msidata[ix, iy, :], xi=new_mzdata, method=method)

        #return variables
        return new_msidata, new_mzdata


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
        from omsi.shared.omsi_data_selection import selection_string_to_object
        zselect = selection_string_to_object(z)  # Convert the selection string to a python selection

        """EDIT_ME Specify the number of custom viewer_options you are going to provide for qslice"""
        num_custom_viewer_options = 0

        # Expose the qslice viewer functionality of any data dependencies
        if viewer_option >= num_custom_viewer_options:
            return super(omsi_mz_upsample, cls).v_qslice(analysis_object,
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
        from omsi.shared.omsi_data_selection import selection_string_to_object
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

            Replace omsi_mz_upsample with your classname
            """
            return super(omsi_mz_upsample, cls).v_qspectrum(analysis_object,
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
            """EDIT_ME Replace the omsi_mz_upsample class name with your class name"""
            mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_mz_upsample, cls)\
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

        Change the omsi_mz_upsample class-name to your class. If you did a
        replace all, then this should be done already.
        """
        dependent_options = super(omsi_mz_upsample, cls).v_qspectrum_viewer_options(analysis_object)
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

        Change the omsi_mz_upsample class-name to your class.  If you did
        a replace all, then this should be done already.
        """
        dependent_options = super(omsi_mz_upsample, cls).v_qslice_viewer_options(analysis_object)
        slice_viewer_options = custom_options + dependent_options
        return slice_viewer_options


############################################################
#  3) Making your analysis self-sufficient   (Recommended) #
############################################################
if __name__ == "__main__":
    from omsi.workflow.driver.cl_analysis_driver import cl_analysis_driver
    """
    EDIT_ME

    Simply replace the omsi_mz_upsample class name with your class name
    """
    cl_analysis_driver(analysis_class=omsi_mz_upsample).main()


