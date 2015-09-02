.. _custom-filereader:

Integrating new file formats
============================

Developing a file reader
------------------------

In order to develop a new file reader we need to implement a corresponding class for the file format that inherits \
form ``file_reader_base`` (for formats that always contain a single region) or \
``file_reader_base_with_regions`` (for formats that support multiple imaging regions in a single file). Both
base classes are availble in the ``omsi.dataformat.file_reader_base`` module. The developer then needs
to implement the following functions:


.. code-block:: python

    from omsi.dataformat.file_reader_base import file_reader_base

    class formatname_file(file_reader_base):

        # 1. Implement the init function which must accept basename and requires_slicing as inputs
        def __init__(self, basename=None, requires_slicing=True):
            """
            basename : Name of the file/folder to be opened
            requires_slicing: Boolean indicating whether the user requires array slicing via
                the __getitem__ function to work or not. This is an optimization, because many MSI
                data formats do not easily support arbitrary slicing of data but rather only
                iteration over spectra.
            """
            super(formatname_file, self).__init__(basename, requires_slicing)  # 1.1 Call super __init__
            self.data_type = 'uint16'   # 1.2 Define the data type used
            self.shape = [0, 0, 0]      # 1.3 Define the shape of the dataset
            self.mz = 0                 # 1.4 Define the m/z axis

        # 2. Implement __getitem__
        def __getitem__(self, key):
            # Implement array-based slicing for the format so that the data can be read
            # via [x,y,z]

        # 3. Implement close_file
        def close_file(self):
            # Function called to close any open files when deleting the object

        # 4. Implement is_valid_dataset
        @classmethod
        def is_valid_dataset(cls, name):
            # Given the name of a file or directory, determine whether the given
            # data defines a valid dataset for the current format


For file formats that support multiple regions the implementation is aside from a few additions he same. The main difference are:

    1) Inherit from ``file_reader_base_with_regions`` instead of ``file_reader_base``
    2) Set the ``self.region_dicts`` and ``self.select_region`` attributes in the ``__init__`` function
    3) Implement the ``set_region_selection`` function


.. code-block:: python

    from omsi.dataformat.file_reader_base import file_reader_base_with_regions

    class formatname_file(file_reader_base):

        # 1. Implement the init function which must accept basename and requires_slicing as inputs
        def __init__(self, basename=None, requires_slicing=True):
            super(formatname_file, self).__init__(basename, requires_slicing)  # 1.1 Call super __init__
            self.data_type = 'uint16'    # 1.2 Define the data type used
            self.shape = [0, 0, 0]       # 1.3 Define the shape of the dataset
            self.mz = 0                  # 1.4 Define the m/z axis
            self.region_dicts = []       # 1.5 Define a list of dicts where each dict describes
                                         #     the parameters of a region. In the simplest case
                                         #     of rectangular regions, a dict would specify the
                                         #     `origin` and `extends` of a region.
            self.select_region = None    # 1.6 Define the index of the selected region.
                                         #     Set to None to treat dataset as a whole (i.e.,
                                         #     merge all regions). If set to a region index,
                                         #     then the __getitem__ function is expected to
                                         #     behave as if the file consisted of just the
                                         #     selected region and the self.shape parameter
                                         #     must be set accordingly

        # 2-4: Implement the other functions of file_reader_base as described above.
        #      Note, __getitem__ must consider the value of self.select_region and
        #      treat any data requests as if they referred to the selected region
        #      only. Depending on the data format this may require transformation
        #      of the selection keys to locate the appropriate data

        # 5. Implement the set_region_selection function to allow a user to select a region
        def set_region_selection(self, region_index=None):
        """Define which region should be selected for local data reads.

           :param region_index: The index of the region that should be read. The shape of the
                    data will be adjusted accordingly. Set to None to select all regions and
                    treat the data as a single full 3D image.
        """
        # 5.1 Select all data
        if region_index is None:
            self.select_region = None           # 5.1.1 Set the region selection to None
            self.shape = self.full_shape        # 5.1.2 Define shape of the complete data
        # 5.2 Select a particular region
        elif region_index < self.get_number_of_regions():
            self.select_region = region_index   # 5.2.1 Set region index
            self.shape = ...                    # 5.2.2 Define the 3D shape of the region



Integrating the file reader with OpenMSI
----------------------------------------

Integrating a new file reader with OpenMSI is simple:

    1) Add the file reader module to the ``omsi.dataformat`` format module and
    2) Add the name of your module to the ``__all__`` variable in ``omis.dataformat.__init__.py``

Once these steps are complete, the ``omsi.dataformat.file_reader_base`` module will automatically detect the new format and make it available as part of the file conversion script ``omsi.tools.convertToOMSI``. To check which formats are registered, simply do:

.. code-block:: python

    >>> from omsi.dataformat import *
    >>> formats = file_reader_base.file_reader_base.available_formats()
    >>> print formats
    {'bruckerflex_file': <class 'omsi.dataformat.bruckerflex_file.bruckerflex_file'>,
     'img_file': <class 'omsi.dataformat.img_file.img_file'>}


Using this basic feature makes it possible to easily iterate over all available formats and the consistent interface described by the ``file_format_base`` module allows us to use all the formats in a consistent manner (avoiding special cases). E.g., if we want to read a file of an unknown format we can simply:

.. code-block:: python

    from omsi.dataformat import *
    formats = file_reader_base.file_reader_base.available_formats()
    filename = 'my_unknown_file'
    filereader = None
    formatname = None
    for fname, fclass in formats.items():
        if fclass.is_valid_dataset(filename):
            filereader = fclass
            formatname = fname
    if filereader is not None:
        print "Using "+str(formatname)+" to read the file"
        openfile = filereader(basename=filename, requires_slicing=True)



