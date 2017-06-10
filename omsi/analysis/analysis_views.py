"""
Helper module with functions and classes for interfacing with different
analysis algorithms. Many of these functions are used to ease interaction
with the analysis module in a generic fashion, without having to explicitly
know about all the different available modules, e.g., we can just look
up modules by name and interact with them directly.
"""
# from omsi.dataformat.omsi_file import *
# from omsi.shared.omsi_data_selection import *
from omsi.shared.data_selection import transform_and_reduce_data
from omsi.analysis import *


class analysis_views(object):
    """
    Helper class for interfacing different analysis algorithms with the web-based viewer
    """

    def __init__(self):
        """Nothing to do here."""
        pass

    @classmethod
    def get_slice(cls,
                  analysis_object,
                  z,
                  operations=None,
                  viewer_option=0):
        """
        Get 3D analysis dataset for which z-slices should be extracted for presentation
        in the OMSI viewer

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param z: Selection string indicting which z values should be selected.
        :param operations: JSON string with list of dictionaries or a python
            list of dictionaries. Each dict specifies a single data
            transformation or data reduction that are applied in order.
            See omsi.shared.omsi_data_selection.transform_and_reduce_data(...)
            for details.
        :param viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them.

        :returns: numpy array with the data to be displayed in the image slice viewer. Slicing
            will be performed typically like [:,:,zmin:zmax].
        """
        analysis_type = analysis_object.get_analysis_type()[0]
        data = cls.analysis_name_to_class(analysis_type).v_qslice(analysis_object, z, viewer_option)
        if data is None:
            return None

        # We expect a 3D dataset (x,y,m/z) so if only a single 2D slice is returned then we reshape the data first
        if len(data.shape) == 2:
            data = data.reshape((data.shape[0], data.shape[1], 1))

        if operations:
            data = transform_and_reduce_data(data=data,
                                             operations=operations,
                                             http_error=True)

        return data

    @classmethod
    def get_spectra(cls,
                    analysis_object,
                    x,
                    y,
                    operations=None,
                    viewer_option=0):
        """
        Get from which 3D analysis spectra in x/y should be extracted for presentation in the OMSI viewer

        **Developer Note:**
        h5py currently supports only a single index list. If the user provides an index-list for both
        x and y, then we need to construct the proper merged list and load the data manually, or if
        the data is small enough, one can load the full data into a numpy array which supports
        multiple lists in the selection.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param x: x selection string
        :param y: y selection string
        :param operations: JSON string with list of dictionaries or a python
            list of dictionaries. Each dict specifies a single data
            transformation or data reduction that are applied in order.
            See omsi.shared.omsi_data_selection.transform_and_reduce_data(...)
            for details.
        :param viewer_option: If multiple default viewer behaviors are available for a given
            analysis then this option is used to switch between them.

        :returns: 2D or 3D numpy array of the requested spectra. The mass (m/z) axis must be the last axis.

        """
        analysis_type = str(analysis_object.get_analysis_type()[0])
        data, custom_mz = cls.analysis_name_to_class(analysis_type).v_qspectrum(analysis_object, x, y, viewer_option)
        if data is None:
            return None, None

        if operations:
            data = transform_and_reduce_data(data=data,
                                             operations=operations,
                                             http_error=True)

        return data, custom_mz

    @classmethod
    def get_axes(cls,
                 analysis_object,
                 qslice_viewer_option=0,
                 qspectrum_viewer_option=0):
        """
        Get the mz axes for the analysis

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed
        :param qslice_viewer_option: If multiple default viewer behaviors are available for
            a given analysis then this option is used to switch between them for the qslice URL pattern.
        :param qspectrum_viewer_option: If multiple default viewer behaviors are available
            for a given analysis then this option is used to switch between them for the
            qspectrum URL pattern.

        :returns: The following four arrays are returned by the analysis:

            - mz_spectra : Array with the static mz values for the spectra.
            - label_spectra : Lable for the spectral mz axis
            - mz_slice : Array of the static mz values for the slices or None if identical to the mz_spectra.
            - label_slice : Lable for the slice mz axis or None if identical to label_spectra.
            - values_x: The values for the x axis of the image (or None)
            - label_x: Label for the x axis of the image
            - values_y: The values for the y axis of the image (or None)
            - label_y: Label for the y axis of the image
            - values_z: The values for the z axis of the image (or None)
            - label_z: Label for the z axis of the image

        """
        analysis_type = str(analysis_object.get_analysis_type()[0])
        return cls.analysis_name_to_class(analysis_type).v_qmz(analysis_object,
                                                               qslice_viewer_option=qslice_viewer_option,
                                                               qspectrum_viewer_option=qspectrum_viewer_option)
        # return mz_spectra, label_spectra, mz_slice, label_slice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def supports_slice(cls,
                       analysis_object):
        """
        Get whether a default slice selection behavior is defined for the analysis.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed

        :return: Boolean indicating whether get_slice(...) is defined for the analysis object.
        """
        supports_slice = len(cls.get_qslice_viewer_options(analysis_object)) > 0
        return supports_slice

    @classmethod
    def supports_spectra(cls,
                         analysis_object):
        """
        Get wheter a default spectra selection behavior is defined for the analysis.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.

        :return: Boolean indicating whether get_spectra(...) is defined for the analysis object.
        """
        supports_spectra = len(cls.get_qspectrum_viewer_options(analysis_object)) > 0
        return supports_spectra

    @classmethod
    def get_qspectrum_viewer_options(cls,
                                     analysis_object):
        """
        Get a list of strings describing the different default viewer options for qspectrum.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.

        :returns: Array of strings indicating the different available viewer options.
            The array may be empty if now viewer_options are available, i.e., get_slice
            and get_spectrum are undefined for the given analysis.
        """
        # viewer_options = []
        analysis_type = str(analysis_object.get_analysis_type()[0])
        viewer_options = cls.analysis_name_to_class(analysis_type).v_qspectrum_viewer_options(analysis_object)
        return viewer_options

    @classmethod
    def get_qslice_viewer_options(cls,
                                  analysis_object):
        """
        Get a list of strings describing the different default viewer options for qslice.

        :param analysis_object: The omsi_file_analysis object for which slicing should be performed.

        :returns: Array of strings indicating the different available viewer options.
            The array may be empty if now viewer_options are available, i.e., get_slice
            and get_spectrum are undefined for the given analysis.
        """
        analysis_type = str(analysis_object.get_analysis_type()[0])
        viewer_options = cls.analysis_name_to_class(analysis_type).v_qslice_viewer_options(analysis_object)
        return viewer_options

    @classmethod
    def analysis_name_to_class(cls,
                               class_name):
        """
        Convert the given string indicating the class to a python class.

        :param class_name: Name of the analysis class. This may be a fully qualified
               name, e.g., `omsi.analysis.multivariate_stat.omsi_nmf` or a name
               relative to the omis.analysis module, e.g, `multivariate_stat.omsi_nmf`.

        :raises: NameError in case that the class cannot be restored.
        """
        import types
        import sys
        try:
            # Remove the omsi.analysis module name from the beginning if present
            if class_name.startswith('omsi.analysis.'):
                class_name = class_name.split(".")[-1]
            if class_name == "generic":
                from omsi.analysis.generic import analysis_generic
                class_object = analysis_generic
            else:
                # Get the class that corresponds to the given name
                class_object = getattr(sys.modules[__name__], class_name)
            if isinstance(class_object, (types.ClassType, types.TypeType)):
                return class_object
            else:
                raise AttributeError("Class could not be restored.")
        except AttributeError:
            raise NameError(class_name+" doesn't exist or is not a class.")

    @classmethod
    def available_analysis(cls):
        """
        Get all available analysis, i.e., all analysis that are subclasses of
        analysis_base.

        :return: Dictionary where the dict-keys are the full qualified name of the
                module and the values are the analysis class corresponding to that
                module.
        """
        return {sub.__module__: sub for sub in analysis_base.__subclasses__()}

    @classmethod
    def available_analysis_descriptions(cls):
        """
        Get all available analysis, i.e., all analysis that are subclasses of
        analysis_base. For each analysis compile the list of input parameters,
        outputs, the corresponding class etc.

        :return: Dictionary where the dict-keys are the full qualified name of the
                module and the values are dicts with class, list of analysis paremeter names,
                list of analysis outputs.
        """
        return {sub.__module__: {'class': sub,
                                 'parameters': sub().get_parameter_names(),
                                 'outputs': sub().get_analysis_data_names(),
                                 'help': sub.__doc__ + "\n\n" + sub.execute_analysis.__doc__}
                for sub in analysis_base.__subclasses__()}

