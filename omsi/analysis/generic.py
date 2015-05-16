"""
Generic analysis class used to represent analyses of unknown type, e.g., when loading
a custom user-defined analysis from file for which the indicate class may not be
available with the local installation. In this case we want to at least be able
to load and investigate the data.
"""
from omsi.analysis.base import analysis_base


class analysis_generic(analysis_base):
    """
    This analysis class is used if the specific anlaysis type is unknown, e.g., when loading
    custom user-defined analysis data that may have not be available in the standard
    omsi package used.
    """
    def __init__(self, name_key="undefined"):
        """
        Initialize the basic data members

        :param name_key: The name for the analysis
        """
        super(analysis_generic, self).__init__()
        self.analysis_identifier = name_key
        self.real_analysis_type = None  # This is the analysis type indicated in the HDF5 file

    def execute_analysis(self):
        """
        Nothing to do here.
        """
        raise NotImplementedError("We cannot run this analysis. analysis_generic cannot run an analysis.")

    @classmethod
    def v_qslice(cls,
                 analysis_object,
                 z,
                 viewer_option=0):
        """
        Implement support for qslice URL requests for the viewer
        """
        return None

    @classmethod
    def v_qspectrum(cls,
                    analysis_object,
                    x,
                    y,
                    viewer_option=0):
        """
        Implement support for qspectrum URL requests for the viewer
        """
        super(analysis_generic, cls).v_qspectrum(analysis_object, x, y, viewer_option)

    @classmethod
    def v_qmz(cls,
              analysis_object,
              qslice_viewer_option=0,
              qspectrum_viewer_option=0):
        """
        Implement support for qmz URL requests for the viewer
        """
        super(analysis_generic, cls).v_qmz(analysis_object, qslice_viewer_option, qspectrum_viewer_option)

    @classmethod
    def v_qspectrum_viewer_options(cls,
                                   analysis_object):
        """
        Define which viewer_options are supported for qspectrum URL's
        """
        return super(analysis_generic, cls).v_qspectrum_viewer_options(analysis_object)

    @classmethod
    def v_qslice_viewer_options(cls,
                                analysis_object):
        """
        Define which viewer_options are supported for qspectrum URL's
        """
        return super(analysis_generic, cls).v_qslice_viewer_options(analysis_object)

    @classmethod
    def get_analysis_type(cls):
        """
        Return a string indicating the type of analysis performed
        """
        return "generic"

    def read_from_omsi_file(self,
                            analysis_object,
                            load_data=True,
                            load_parameters=True,
                            load_runtime_data=True,
                            dependencies_omsi_format=True,
                            ignore_type_conflict=False):
        """
        See `omsi.analysis.analysis_base.read_from_omsi_file(...)` for details.
        The function is overwritten here mainly to initialize the self.real_analysis_type
        instance variable but otherwise uses the default behavior.

        """
        output_val = super(analysis_generic, self).read_from_omsi_file(
            analysis_object=analysis_object,
            load_data=load_data,
            load_parameters=load_parameters,
            load_runtime_data=load_runtime_data,
            dependencies_omsi_format=dependencies_omsi_format,
            ignore_type_conflict=ignore_type_conflict)
        self.real_analysis_type = unicode(analysis_object.get_analysis_type()[:])
        return output_val

    def get_real_analysis_type(self):
        """
        This class is designed to handle generic (including unkown) types of analysis.
        In cases, e.g., were this class is used to store analysis data from an HDF5
        file we may have an actual analysis type available even if we do not have
        a special analysis class may not be available in the current installation
        """
        return self.real_analysis_type

