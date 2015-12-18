"""
Module for managing OpenMSI HDF5 data files.
"""

import time

import h5py
import numpy as np

from omsi.dataformat.omsi_file.common import omsi_file_common
from omsi.dataformat.omsi_file.format import omsi_format_file
from omsi.dataformat.omsi_file.experiment import omsi_experiment_manager
try:
    from omsi.shared import mpi_helper
    if mpi_helper.MPI_AVAILABLE:
        from mpi4py import MPI
except:
    pass


class omsi_file(omsi_experiment_manager,
                omsi_file_common):
    """
    API for creating and managing a single OpenMSI data file.

    **Use of supe()r**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common` .
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the `__init__` function calls
    `super(...).__init__(manager_group)` with a single  parameter indicating the
    parent group.

    **Inherited Instance Variables**

    :ivar managed_group: The group that is managed by this object
    :ivar name: Name of the managed group

    """
    @classmethod
    def __create__(cls,
                   filename,
                   mode,
                   create_omsi_file_object=True,
                   **kwargs):
        import os
        if os.path.exists(filename):
            raise ValueError('The input filename already exists.')
        output_file = h5py.File(filename, mode=mode, **kwargs)
        return cls.__populate_file__(input_file=output_file,
                                     create_omsi_file_object=create_omsi_file_object)

    @classmethod
    def __populate_file__(cls,
                          input_file,
                          create_omsi_file_object=True):
        """
        Populate a new file with the basic required format data

        :param input_file: The h5py.File input file
        :param create_omsi_file_object: Boolean indicating whether we should create a new omsi_file object
                           for the given input_file after it has been populated. (Default is True).
        :return: omsi_file if create_omsi_file_object is True otherwise the input_file is returned.
        """
        main_group = input_file['/']
        main_group.attrs[omsi_format_file.type_attribute] = "omsi_file"
        main_group.attrs[omsi_format_file.version_attribute] = omsi_format_file.current_version
        main_group.attrs[omsi_format_file.timestamp_attribute] = str(time.ctime())
        if create_omsi_file_object:
            return omsi_file(main_group.file)
        else:
            return input_file

    @classmethod
    def is_valid_dataset(cls, name):
        """
        Perform basic checks for the given filename, whether it is a valid OMSI file.

        :param name: Name of the file to be checked.

        :return: Boolean indicating whether the file is valid
        """

        try:
            checkfile = omsi_file(filename=name, mode='r')
            valid = checkfile.get_num_experiments() > 0
            checkfile.close_file()
            del checkfile
            return valid
        except:
            return False

    def __init__(self, filename, mode='a', **kwargs):
        """
        Open the given file or create it if does not exit.

        The creation of the object may fail  if the file does not exist, and
        the selected mode is 'r' or 'r+'.

        Keyword arguments:

        :param filename: string indicating the name+path of the OpenMSI data file. Alternatively
                        this may also be an h5py.File instance (or an h5py.Group, h5py.Dataset
                        instance from which we can get the file)
        :param mode: read/write mode. One of : \n
                    r = readonly, file must exist. \n
                    r+ = read/write, file must exist. \n
                    w = Create file, truncate if exists. \n
                    w- = create file, fail if exists. \n
                    a = read/write if exists, create otherwise (default)
        :param **kargs: Other keyword arguments to be used for opening the file using h5py.
                    See the h5py.File documentation for details. For example to use parallel HDF5,
                    the following additional parameters can be given driver='mpio', comm:MPI.COMM_WORLD.
        """
        if isinstance(filename, h5py.File):
            self.hdf_filename = filename.filename
            self.hdf_file = filename
        elif isinstance(filename, h5py.Group) or isinstance(filename, h5py.Dataset):
            self.hdf_file = filename.file
            self.hdf_filename = self.hdf_file.filename
        elif isinstance(filename, basestring):
            import os
            self.hdf_filename = filename  # Name of the HDF5 file
            if os.path.exists(filename):
                self.hdf_file = h5py.File(filename, mode=mode, **kwargs)
            else:
                self.hdf_file = self.__create__(filename=filename,
                                                mode=mode,
                                                create_omsi_file_object=False,
                                                **kwargs)
                # self.hdf_file = h5py.File(filename, mode=mode, **kwargs)
                # _ = self.__populate_file__(self.hdf_file,
                #                           create_omsi_file_object=False)
        else:
            raise ValueError("filename must be h5py.File, h5py.Group, h5py.Dataset, or a valid filename")

        super(omsi_file, self).__init__(self.hdf_file['/'])
        # The super call takes care of the following initializations
        # self.managed_group = self.hdf_file['/']
        # self.name = self.hdf_file.name

    def get_h5py_file(self):
        """
        Get the h5py object for the omsi file

        :returns: h5py redernce to the HDF5 file
        """
        return self.hdf_file

    def get_filename(self):
        """
        Get the name of the omsi file

        :returns: String indicating the filename (possibly including the full path,
                  depending on how the object has been initalized)
        """
        return self.hdf_filename

    #####################################################
    #  File related functions
    ####################################################
    def close_file(self):
        """
        Close the msi data file
        """

        self.hdf_file.flush()  # Making sure all data has been written
        self.hdf_file.close()  # Close the file

    def flush(self):
        """
        Flush all I/O
        """
        self.hdf_file.flush()

    def write_xdmf_header(self, xdmf_filename):
        """
        Write XDMF header file for the current HDF5 datafile

        :param xdmf_filename: The name of the xdmf XML header file to be created for the HDF5 file.
        """
        import os
        from omsi.dataformat.omsi_file.format import omsi_format_experiment, omsi_format_msidata
        # Open the output header file
        xdmf = open(xdmf_filename, 'w')
        # Write the starting header data
        xdmf.write("<?xml version=\"1.0\" ?>\n")
        xdmf.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
        xdmf.write("<Xdmf Version=\"2.0\">\n")
        xdmf.write("<Domain>\n")
        # For each experiment (/entry_#) and dataset (/entry_#/data_#) write an according entry
        # poly.write("%d %f %f" % ( i , vertices[i,0] , vertices[i,1]  )  )
        numexp = self.get_num_experiments()
        for exp_index in xrange(0, numexp):
            exp = self.get_experiment(exp_index)
            numdat = exp.get_num_msidata()
            # For all datasets available for the experiment
            for dataset_index in xrange(0, numdat):
                data = exp.get_msidata(dataset_index)
                size_x = data.shape[0]
                size_y = data.shape[1]
                size_z = data.shape[2]
                xdmf.write("<Grid Name=\"ImageGrid\" GridType=\"Uniform\">\n")
                xdmf.write("  <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n" %
                           (size_x, size_y, size_z))
                xdmf.write("  <Geometry GeometryType=\"Origin_DxDyDz\">\n")
                xdmf.write(
                    "    <DataItem Format=\"XML\" Dimensions=\"1 3\">\n")
                xdmf.write("      0.0\n")
                xdmf.write("      0.0\n")
                xdmf.write("      0.0\n")
                xdmf.write("    </DataItem>\n")
                xdmf.write(
                    "    <DataItem Format=\"XML\" Dimensions=\"1 3\">\n")
                xdmf.write("      %f\n" %
                           (float(data.shape[0]) / float(data.shape[2])))
                xdmf.write("      1.0\n")
                xdmf.write("      1.0\n")
                xdmf.write("    </DataItem>\n")
                xdmf.write("  </Geometry>\n")
                xdmf.write(
                    "  <Attribute Name=\"Spectrum\" AttributeType=\"Scalar\" Center=\"Node\">\n")
                xdmf.write(
                    "    <DataItem Dimensions=\"%d %d %d\" NumberType=\"UInt\" Format=\"HDF\">\n" %
                    (size_x, size_y, size_z))
                xdmf.write("    " + os.path.basename(self.hdf_filename) + ":/" +
                           omsi_format_experiment.exp_groupname + str(exp_index) + "/" +
                           omsi_format_msidata.data_groupname + str(dataset_index) + "\n")
                xdmf.write("    </DataItem>\n")
                xdmf.write("  </Attribute>\n")
                xdmf.write(" </Grid>\n")

        # Write ending of xdmf file
        xdmf.write("</Domain> \n")
        xdmf.write("</Xdmf> \n")
