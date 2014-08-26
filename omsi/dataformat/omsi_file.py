"""
This module defines the API for reading and writing HDF5 files for storing
mass spectrometry imaging data, metadata, and analysis results according to
the format defined by the OpenMSI project

TODO: omsi_file_msidata
    -- The partial spectra case is untested and needs attention
    -- The __getitem__ function for the partial_spectra case is not fully implemented yet
    -- The estimates in def __best_dataset__(self,keys) are fairly crude at this point
    -- The __setitem__ function should be implemented for the different cases as well (AHHHH).
    -- Note, for the partial_spectra and partial_cube case the datasets needs to be dynamically expandable.

"""
# TODO Rename omsi_file_methods to omsi_file_methods and rename the group from method to methods in a backward compatible way

import h5py
import numpy as np
import math
import time
from omsi.dataformat.omsi_format import *
try:
    from mpi4py import MPI
except:
    pass


##########################################################################
##########################################################################
##   The main class for omsi data files                                  ##
##########################################################################
##########################################################################
class omsi_file(object):

    """API for creating and managing a single OpemMSI data file.
    """
    @classmethod
    def is_managed(cls, in_object):
        """ Check whether the given objet is managed by any omsi API class.

            :param in_object: The object to be checked
            :type in_object: Any omsi_file API object or h5py.Dataset or h5py.Group or h5py.File object.
        """
        managed_object = omsi_file.get_omsi_object(in_object)
        if isinstance(managed_object, omsi_file) or \
           isinstance(managed_object, omsi_file_experiment) or \
           isinstance(managed_object, omsi_file_methods) or\
           isinstance(managed_object, omsi_file_instrument) or \
           isinstance(managed_object, omsi_file_analysis) or \
           isinstance(managed_object, omsi_file_msidata) or \
           isinstance(managed_object, omsi_file_dependencies) or \
           isinstance(managed_object, omsi_file_dependencydata):
            return True
        else:
            return False

    @classmethod
    def get_h5py_object(cls, omsi_object, resolve_dependencies=False):
        """ This static method is a convenience function used to retrieve the corresponding h5py
            interface object for any omsi file API object.

            :param omsi_object: omsi file API input object for which the corresponding h5py.Group, h5py.File, or
                            h5py.Dataset object should be retrieved. The omsi_object input may itself also be a h5py.Group,
                            h5py.File, or h5py.Dataset, in which case omsi_object itself is returned by the function.
            :param resolve_dependencies: Set to True if omsi_file_dependencydata objects should be resolved to retrieve
                            the dependent object the dependency is pointing to. Dependencies are resolved recursively,
                            i.e., if a dependency points to another dependency then that one will be resolved as well.
                            Default value is False, i.e., the omis_file_dependency object itself is returned.
            :returns:  h5py.Group, h5py.File, or h5py.Dataset corresponding to the given omsi_object.

            :raises ValueError: A ValueError is raised in case that an unsupported omsi_object object is given, i.e., the
                                input object is not a omsi_file API object nor a h5py Group, File, or Dataset object.
        """
        if isinstance(omsi_object, omsi_file):
            h5pyobject = omsi_object.get_h5py_file()
        elif isinstance(omsi_object, omsi_file_experiment):
            h5pyobject = omsi_object.get_h5py_experimentgroup()
        elif isinstance(omsi_object, omsi_file_methods):
            h5pyobject = omsi_object.get_h5py_methodgroup()
        elif isinstance(omsi_object, omsi_file_instrument):
            h5pyobject = omsi_object.get_h5py_instrumentgroup()
        elif isinstance(omsi_object, omsi_file_analysis):
            h5pyobject = omsi_object.get_h5py_analysisgroup()
        elif isinstance(omsi_object, omsi_file_dependencydata):
            if resolve_dependencies:
                h5pyobject = omsi_file.get_h5py_object(
                    omsi_object.get_dependency_omsiobject(), resolve_dependencies)
            else:
                h5pyobject = omsi_object.get_h5py_dependencygroup()
        elif isinstance(omsi_object, omsi_file_dependencies):
            h5pyobject = omsi_object.get_h5py_dependenciesgroup()
        elif isinstance(omsi_object, omsi_file_msidata):
            h5pyobject = omsi_object.get_h5py_datagroup()
        elif isinstance(omsi_object, h5py.File) or \
                isinstance(omsi_object, h5py.Group) or \
                isinstance(omsi_object, h5py.Dataset):
            h5pyobject = omsi_object
        else:
            raise ValueError(
                "Unsupported input object given to omsi_file.get_h5py_object.")
        return h5pyobject

    @classmethod
    def get_omsi_object(cls, h5py_object, resolve_dependencies=False):
        """This static method is convenience function used to retrieve the corresponding interface class for a
           given h5py group object.

            :param h5py_object: h5py object for which the corresponding omsi_file API object should be generated.
            :param resolve_dependencies: Set to True if omsi_file_dependencydata objects should be resolved to retrieve
                            the dependent object the dependency is pointing to. Dependencies are resolved recursively,
                            i.e., if a dependency points to another dependency then that one will be resolved as well.
                            Default value is False, i.e., the omis_file_dependency object itself is returned.

            :returns: None in case no corresponding object was found. Otherwise an instance of:

                * omsi_file : If the given object is a h5py.File object
                * omsi_file_experiment : If the given object is an experiment groupt
                * omsi_file_methods : If the given object is a method group
                * omsi_file_instrument : If the given object is an instrument group
                * omsi_file_analysis : If the given object is an analysis group
                * omsi_file_msidata : If the given object is a MSI data group
                * omsi_file_dependencydata : If the fiven object is a dependency group
                * The input h5py_object: If the given objet is a h5py.Dataset
                * None: In case that an unknown type is given. E.g. a user may have \
                        given an unmanaged Group object that does not have a \
                        corresponding omsi file API object.
        """

        # If the input object is already an omsi API object then return it as
        # is
        if isinstance(h5py_object, omsi_file) or \
           isinstance(h5py_object, omsi_file_experiment) or \
           isinstance(h5py_object, omsi_file_methods) or\
           isinstance(h5py_object, omsi_file_instrument) or \
           isinstance(h5py_object, omsi_file_analysis) or \
           isinstance(h5py_object, omsi_file_msidata) or \
           isinstance(h5py_object, omsi_file_dependencies) or \
           isinstance(h5py_object, omsi_file_dependencydata):
            return h5py_object

        if isinstance(h5py_object, h5py.File):
            return omsi_file(h5py_object)
        elif isinstance(h5py_object, h5py.Group):
            # Check if the group has an explicit type attribute
            try:
                # Try to determine the type of the group based on the
                # attributes
                type_attribute = h5py_object.attrs[
                    omsi_format_common.type_attribute]
                if type_attribute == "omsi_file_experiment":
                    return omsi_file_experiment(h5py_object)
                elif type_attribute == "omsi_file_methods" or type_attribute == "omsi_file_methods":
                    return omsi_file_methods(h5py_object)
                elif type_attribute == "omsi_file_instrument":
                    return omsi_file_instrument(h5py_object)
                elif type_attribute == "omsi_file_analysis":
                    return omsi_file_analysis(h5py_object)
                elif type_attribute == "omsi_file_msidata":
                    return omsi_file_msidata(h5py_object)
                elif type_attribute == "omsi_file":
                    return omsi_file(h5py_object)
                elif type_attribute == "omsi_file_dependencydata":
                    omsiobject = omsi_file_dependencydata(h5py_object)
                    if resolve_dependencies:
                        return omsi_file.get_omsi_object(
                            omsi_file.get_h5py_object(
                                omsiobject.get_dependency_omsiobject(),
                                resolve_dependencies))
                    else:
                        return omsiobject
                elif type_attribute == "omsi_file_dependencies":
                    return omsi_file_dependencies(h5py_object)
                else:
                    return None
            except:
                # If the attribute is missing, then try to determin the type
                # based on he name of group
                groupname = h5py_object.name.split("/")[-1]
                parentgroupname = h5py_object.parent.name.split("/")[-1]
                if groupname.startswith(omsi_format_experiment.exp_groupname):
                    return omsi_file_experiment(h5py_object)
                elif groupname.startswith(omsi_format_methods.methods_groupname) or \
                        groupname.startswith(omsi_format_methods.methods_old_groupname):
                    return omsi_file_methods(h5py_object)
                elif groupname.startswith(omsi_format_instrument.instrument_groupname):
                    return omsi_file_instrument(h5py_object)
                elif groupname.startswith(omsi_format_analysis.analysis_groupname):
                    return omsi_file_analysis(h5py_object)
                elif groupname.startswith(omsi_format_data.data_groupname):
                    return omsi_file_msidata(h5py_object)
                elif groupname.startswith(omsi_format_dependencies.dependencies_groupname):
                    return omsi_file_dependencies(h5py_object)
                elif parentgroupname.startswith(omsi_format_dependencies.dependencies_groupname):
                    omsiobject = omsi_file_dependencydata(h5py_object)
                    if resolve_dependencies:
                        return omsi_file.get_omsi_object(
                            omsi_file.get_h5py_object(
                                omsiobject.get_dependency_omsiobject(),
                                resolve_dependencies))
                    else:
                        return omsiobject
                elif groupname == "":  # We are at the root group
                    return omsi_file(h5py_object.file)
                else:
                    return None
        elif isinstance(h5py_object, h5py.Dataset):
            return h5py_object
        else:
            return None

    @classmethod
    def is_valid_dataset(cls, name):
        try:
            checkfile = omsi_file(filename=name, mode='r')
            valid = checkfile.get_num_exp() > 0
            checkfile.close_file()
            del checkfile
            return valid
        except:
            return False

    def __init__(self, filename, mode='a', **kwargs):
        """Open the given file or create it if does not exit.

           The creation of the object may fail  if the file does not exist, and
           the selected mode is 'r' or 'r+'.

           Keyword arguments:

           :param filename: string indicating the name+path of the OpenMSI data file. Alternatively
                            this may also be an h5py.File instance.
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
        else:
            self.hdf_filename = filename  # Name of the HDF5 file
            # This is a public attribute.
            self.hdf_file = h5py.File(filename, mode=mode, **kwargs)
        self.name = self.hdf_file.name

    def __getitem__(self, key):
        """Support direct read interaction with the h5py file"""
        return self.hdf_file[key]

    def __setitem__(self, key, value):
        """Support direct write interaction with the h5py file"""
        self.hdf_file[key] = value

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.hdf_file.items()

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def get_h5py_file(self):
        """Get the h5py object for the omsi file

           :returns: h5py redernce to the HDF5 file
        """
        return self.hdf_file

    def get_filename(self):
        """Get the name of the omsi file

           :returns: String indicating the filename (possibly including the full path,
                     depending on how the object has been initalized)
        """
        return self.hdf_filename

    #####################################################
    #  Access the different experiments, associated
    #  properties etc.
    ####################################################
    @staticmethod
    def get_exp_path(exp_index=None):
        """Based on the index of experiment return the full path to the hdf5
           group containing the data for an experiment.

           :param exp_index: The index of the experiment.

           :returns: String indicating the path to the experiment.
        """

        expindexstr = ""
        if exp_index is not None:
            expindexstr = str(exp_index)
        return "/" + omsi_format_experiment.exp_groupname + expindexstr

    def get_num_exp(self):
        """Get the number of experiments in this file.

           :returns: Integer indicating the number of experiments.
        """

        return omsi_file.get_num_items(self.hdf_file.get("/"), omsi_format_experiment.exp_groupname)

    def get_exp(self, exp_index):
        """Get the omsi_format_experiment object for the experiment with the given index

           :param exp_index: The index of the requested experiment
           :type exp_index: uint

           :returns: h5py reference to the experiment with the given index. Returns None in case \
                     the experiment does not exist.
        """

        if exp_index < self.get_num_exp():
            return omsi_file_experiment(self.hdf_file[unicode(omsi_format_experiment.exp_groupname + str(exp_index))])
        else:
            return None

    def get_exp_by_identifier(self, exp_identifier_string):
        """Get the omsi_format_experiment object for the experiment with the given identifier.

           :param exp_identifier_string: The string used to identify the analysis
           :type exp_identifier_string: string

           :returns: Returns h5py object of the experiment group or None in case the experiment is not found."""

        # Iterate through all groups of the root folder
        for it in self.hdf_file.get("/").items():
            if it[0].startswith(omsi_format_experiment.exp_groupname):
                cur_exp_id = omsi_file_experiment(
                    self.hdf_file[it[0]]).get_exp_identifier()
                if cur_exp_id is not None:
                    if cur_exp_id[0] == exp_identifier_string:
                        return omsi_file_experiment(self.hdf_file[it[0]])
        return None

    #####################################################
    #  Functions used for creation of the standard
    #  data of the OMSI file format
    ####################################################
    def create_exp(self, exp_identifier=None, flush_io=True):
        """Create a new group in the file for a new experiment and return the omsi_file_experiment
           object for the new experiment.

           :param exp_identifier: The string used to identify the analysis
           :type exp_identifier: string or None (default)
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that
                           all data has been written to file.

           :returns: omsi_file_experiment object for the newly created group for the experiment
        """
        numexp = self.get_num_exp()
        expindex = numexp
        expname = omsi_format_experiment.exp_groupname + str(expindex)
        # Using require_group ensures that the group is not overwritten if it
        # already exists
        exp = self.hdf_file.get("/").require_group(expname)
        exp.attrs[omsi_format_common.type_attribute] = "omsi_file_experiment"
        exp.attrs[omsi_format_common.version_attribute] = omsi_format_experiment.current_version
        exp.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        exp_identifier_dataset = exp.require_dataset(name=unicode(
            omsi_format_experiment.exp_identifier_name), shape=(1,), dtype=omsi_format_common.str_type)
        if exp_identifier is not None:
            if omsi_format_common.str_type_unicode:
                exp_identifier_dataset[0] = exp_identifier
            else:
                exp_identifier_dataset[0] = str(exp_identifier)
        else:
            exp_identifier_dataset[0] = "undefined"

        if flush_io:
            self.hdf_file.flush()

        return omsi_file_experiment(exp)

    #####################################################
    #  File related functions
    ####################################################
    def close_file(self):
        """Close the msi data file"""

        self.hdf_file.flush()  # Making sure all data has been written
        self.hdf_file.close()  # Close the file

    def flush(self):
        """Flush all I/O"""
        self.hdf_file.flush()

    def write_xdmf_header(self, xdmf_filename):
        """Write XDMF header file for the current HDF5 datafile

           :param xdmf_filename: The name of the xdmf XML header file to be created for the HDF5 file.
        """

        import os
        # Open the output header file
        xdmf = open(xdmf_filename, 'w')
        # Write the starting header data
        xdmf.write("<?xml version=\"1.0\" ?>\n")
        xdmf.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
        xdmf.write("<Xdmf Version=\"2.0\">\n")
        xdmf.write("<Domain>\n")
        # For each experiment (/entry_#) and dataset (/entry_#/data_#) write an according entry
        #poly.write("%d %f %f" % ( i , vertices[i,0] , vertices[i,1]  )  )
        numexp = self.get_num_exp()
        for ei in xrange(0, numexp):
            exp = self.get_exp(ei)
            numdat = exp.get_num_msidata()
            # For all datasets available for the experiment
            for di in xrange(0, numdat):
                data = exp.get_msidata(di)
                nx = data.shape[0]
                ny = data.shape[1]
                nz = data.shape[2]
                xdmf.write("<Grid Name=\"ImageGrid\" GridType=\"Uniform\">\n")
                xdmf.write(
                    "  <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n" %
                    (nx, ny, nz))
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
                    (nx, ny, nz))
                xdmf.write("    " + os.path.basename(self.hdf_filename) + ":/" +
                           omsi_format_experiment.exp_groupname + str(ei) + "/" +
                           omsi_format_msidata.data_groupname + str(di) + "\n")
                xdmf.write("    </DataItem>\n")
                xdmf.write("  </Attribute>\n")
                xdmf.write(" </Grid>\n")

        # Write ending of xdmf file
        xdmf.write("</Domain> \n")
        xdmf.write("</Xdmf> \n")

    #####################################################
    #  Private helper functions and class methods
    ####################################################
    @classmethod
    def get_num_items(cls, file_group, basename=""):
        """Get the number of object with the given basename at the given path

           :param file_group: The h5py object to be examined
           :param basename: The name that should be searched for.

           :returns: Number of objexts with the given basename at the given path
        """

        numitems = 0
        # Iterate through all groups of the root folder
        for it in file_group.items():
            if it[0].startswith(basename):
                numitems += 1
        return numitems

    @classmethod
    def unit_test(cls, test_filename="test.h5"):
        """This simple unit test function creates a new OpenMSI HDF5
           data file and tries to populate it with dummy data

          :param test_filename: The name of the test HDF5 file to be created.
        """
        from omsi.analysis.omsi_analysis_generic import omsi_analysis_generic
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data

        testfile = omsi_file(test_filename)
        print "1. Creating new test experiment"
        exp = testfile.create_exp()
        if exp is None:
            print "FAILED TEST"
            exit()

        print "Retrieving experiement wiht index 0"
        print testfile.get_exp(0)

        print "Setting experiment identifier"
        exp_identifier_text = "Experiment 1"
        exp.set_exp_identifier(exp_identifier_text)

        print "Retrieving the experiment identifier"
        print exp.get_exp_identifier()

        print "Retrieving experiment by identifier"
        print testfile.get_exp_by_identifier(exp_identifier_text)

        print "Adding a full cube dataset to the experiment"
        tempshape = tuple([10, 10, 1000])
        tempnumelements = tempshape[0] * tempshape[1] * tempshape[2]
        specchunking = tuple([1, 1, 1000])
        slicchunking = tuple([10, 10, 1])
        data_dataset, mz_dataset, data_group = exp.create_msidata_full_cube(
            data_shape=tempshape, chunks=specchunking)

        print "Assigning data to dataset"
        t = omsi_file_msidata(data_group)
        t[:] = np.arange(tempnumelements).reshape(tempshape)
        #data_dataset[:] = np.arange( tempnumelements ).reshape( tempshape )

        print "Assigning data to mz_data"
        mz_dataset[:] = np.arange(tempshape[2])

        print "Creating optimized copy of the data"
        testdataobj = omsi_file_msidata(data_group)
        testdataobj.create_optimized_chunking(
            chunks=slicchunking, compression=None, compression_opts=None, copy_data=True, print_status=True)

        print "Retrieving the dataset"
        dataset = exp.get_msidata(0)
        print dataset

        print "Testing data slicing operations"
        print dataset[5, 5, :]
        print dataset[:, :, 1]
        print dataset[3:4, 1:4, 5:50]
        print dataset[[1, 3], 8:9, :]
        print dataset[:, :, :]

        print "Adding a partial cube dataset to the experiment"
        mask = np.zeros((tempshape[0], tempshape[1]), dtype='bool')
        mask[1:8, 5:10] = True
        data_dataset1, mz_dataset1, xy_index_dataset1, inv_xy_index_dataset1, datagroup1 = exp.create_msidata_partial_cube(
            data_shape=tempshape,
            mask=mask,
            chunks=specchunking)
        testdataobj1 = omsi_file_msidata(datagroup1)
        testdataobj1[1:8, 5:10, :] = np.arange(7 * 5 * 1000).reshape((7, 5, 1000))
        mz_dataset[:] = np.arange(tempshape[2])
        print "Creating optimized copy of the data"
        testdataobj1.create_optimized_chunking(
            chunks=slicchunking, compression=None, compression_opts=None, copy_data=True, print_status=True)
        print "Retrieving the dataset"
        dataset = exp.get_msidata(1)
        print dataset
        print "Testing data slicing operations"
        print dataset[1, 1, :]
        print dataset[:, :, 1]
        print dataset[1:2, 1:2, 5:50]
        print dataset[[1, 3], 8:9, :]
        print dataset[:, :, :]

        print "Checking for number of experiments"
        print testfile.get_num_exp()

        print "Checking for number of msidata for the experiment"
        print exp.get_num_msidata()

        print "Adding method information"
        methodinfo = exp.create_method_info()

        print "Getting the method information group"
        print exp.get_method_info()

        print "Getting the method name"
        methodname = methodinfo.get_method_name()
        print methodname[0]
        print "Modifing method name"
        methodinfo.set_method_name("Mouse")
        print methodinfo.get_method_name()[0]

        print "Adding instrument information"
        mzdata = np.arange(10, dtype='float32')
        instrumentinfo = exp.create_instrument_info("undefined", mzdata)

        print "Getting the instrument information group"
        print exp.get_instrument_info()

        print "Getting the instrument name"
        instrumentname = instrumentinfo.get_instrument_name()
        if instrumentname is not None:
            print instrumentname[0]
        else:
            print "FAILED TEST"
        print "Modifing instrument name"
        instrumentinfo.set_instrument_name("LBNL_OpenMSI Python Instrument")
        print instrumentinfo.get_instrument_name()[0]

        print "Getting the instrument mz data"
        print instrumentinfo.get_instrument_mz()

        print "Creating derived analysis"
        testanaidname = "Peak Finding 123"
        testana = omsi_analysis_generic(nameKey=testanaidname)
        print "Creating dummy analysis data"
        testanadata = omsi_analysis_data()
        testanadata['name'] = 'peakcube'
        testanadata['data'] = np.zeros(shape=(5, 5, 5), dtype='float32')
        testanadata['dtype'] = 'float32'
        # print "Adding the dummy analysis data to the analys object"
        #testana.add_analysis_data( testanadata )
        print "Saving the analysis data to file"
        analysis, ai = exp.create_analysis(testana)
        print "Getting the analysis identifier"
        tempanaid = analysis.get_analysis_identifier()
        if tempanaid is not None:
            print tempanaid[0]
        else:
            print "FAILED TEST"
        print "Getting the analysis type"
        tempanatype = analysis.get_analysis_type()
        if tempanatype is not None:
            print tempanatype[0]
        else:
            print "FAILED TEST"
        print "Getting the analysis by identifier"
        print exp.get_analysis_by_identifier(testanaidname)

        print "Closing the file"
        testfile.close_file()


##########################################################################
##########################################################################
#   The main class for experiments within omsi data files
##########################################################################
##########################################################################
class omsi_file_experiment(object):

    """Class for managing experiment specific data"""

    def __init__(self, exp_group):
        """Initalize the experiment object given the h5py object of the experiment group

           :param exp_group: The h5py object with the experiment group of the omsi hdf5 file.
        """
        self.experiment = exp_group
        self.name = self.experiment.name

    def __getitem__(self, key):
        """Support direct read interaction with the experiment h5py group"""
        return self.experiment[key]

    def __setitem__(self, key, value):
        """Support direct write interaction with the experiment h5py group"""
        self.experiment[key] = value

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.experiment.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.experiment.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the experiment group was created in the HDF5 file.

           :returns: Python timestamp string generated using time.ctime().
                     None may be returned in case that the timestamp does not exists
                     or cannot be retrieved from the file for some reason.

        """
        try:
            return self.experiment.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    ##########################################################
    #  Access the data directly associated with the experiment
    ##########################################################
    def get_h5py_experimentgroup(self):
        """Get the h5py object with the experiment data.

           :returns: The hpy reference for the experiment with which this instance of
                     omsi_file_experiment has been initalized.
        """
        return self.experiment

    def get_exp_identifier(self):
        """Get the HDF5 dataset with the identifier description for the experiment.

           :returns: h5py object of the experiment identifer or None in case not present"""

        try:
            return self.experiment[
                unicode(omsi_format_experiment.exp_identifier_name)]
        except:
            return None

    def get_num_msidata(self):
        """Get the number of raw mass spectrometry images stored for a given experiment

           :returns: Integer indicating the number of msi datasets available for the experiment.
        """

        return omsi_file.get_num_items(self.experiment, omsi_format_msidata.data_groupname)

    def get_msidata(self, data_index, fill_space=True, fill_spectra=True, fill_value=0, preload_mz=True, preload_xy_index=True):
        """Get the dataset with the given index for the given experiment.

           For more detailed information about the use of the fill_space and fill_spectra and preload_mz and
           preload_xy_index options, see the init function of omsi.dataformat.omsi_file_msidata.

           :param data_index: Index of the dataset.
           :type data_index: unsigned int
           :param fill_space: Define whether the data should be padded in space (filled with 0's)
                            when accessing the data using [..] operator so that the data behaves like a 3D cube.
           :param fill_spectra: Define whether the spectra should completed by adding 0's so that
                            all spectra retrived via the [..] opeator so that always spectra of the full
                            length are returned.
           :param preload_mz: Should the data for the mz axis be kept in memory or loaded on the fly when needed.
           :param preload_xy_index: Should the xy index (if available) be preloaderd into memory or should the
                            required data be loaded on the fly when needed.
           :param fill_value: The integer value that should be used to fill in missing data values.

           :returns: omsi_file_msidata object for the given data_index or None in case the data
                     with given index does not exist or the access failed for any
                     other reason.
        """
        try:
            return omsi_file_msidata(
                self.experiment[unicode(omsi_format_msidata.data_groupname + str(data_index))],
                fill_space=fill_space,
                fill_spectra=fill_spectra,
                fill_value=fill_value,
                preload_mz=preload_mz,
                preload_xy_index=preload_xy_index)
        except:
            return None

    def get_msidata_by_name(self, data_name):
        """Get the h5py data object for the the msidata with the given name.

           :param data_name: The name of the dataset
           :type data_name: string

           :returns: h5py object of the dataset or None in case the dataset is not found."""
        # Iterate through all items of the experiment
        for it in self.experiment.items():
            if it[0] == data_name:
                return self.experiment[it[0]]
        return None

    def set_exp_identifier(self, identifier):
        """Overwrite the current identfier string for the experiment with the given string

           :param identifier: The new experiment identifier string.
        """

        # Get the name of the intrument
        expid = self.get_exp_identifier()
        # Create the dataset for the id name if it does not exist
        if expid is None:
            expid = self.experiment.require_dataset(name=unicode(
                omsi_format_experiment.exp_identifier_name), shape=(1,), dtype=omsi_format_common.str_type)

        if omsi_format_common.str_type_unicode:
            expid[0] = identifier
        else:
            expid[0] = str(identifier)

    ###########################################################
    # Get sub-group object associated with the experiment
    ###########################################################
    def has_method_info(self):
        """
        Check whether the experiment has a method info object.

        :returns: Bool indicating whether the experiment has a method info object.
        """
        return (self.get_method_info() is not None)

    def has_instrument_info(self):
        """
        Check whether the experiment has an instrument info object.

        :returns: Bool indicating whether the experiment has an instrument info object.
        """
        return (self.get_instrument_info() is not None)


    def get_method_info(self):
        """Get the omsi_file_methods object with the method information.

           :returns: omsi_file_methods object for the requested method info. The function returns \
                     None in case no method information was found for the experiment
        """
        try:
            return omsi_file_methods(self.experiment[unicode(omsi_format_methods.methods_groupname)])
        except:
            try:
                return omsi_file_methods(self.experiment[unicode(omsi_format_methods.methods_old_groupname)])
            except:
                return None

    def get_instrument_info(self):
        """Get the HDF5 group opbject with the instrument information.

           :returns:  omsi_file_instrument object for the requested instrument info. The function returns \
                      None in case no instrument information was found for the experiment
        """
        instrument_info = None
        try:
            return omsi_file_instrument(
                self.experiment[unicode(omsi_format_instrument.instrument_groupname)])
        except:
            return None

    def get_num_analysis(self):
        """Get the number of raw mass spectrometry images stored for a given experiment

           :returns: Integer indicating the number of analyses available for the experiment.
        """
        return omsi_file.get_num_items(self.experiment, omsi_format_analysis.analysis_groupname)

    def get_analysis_identifiers(self):
        """Get a list of all identifiers for all analysis stored for the experiment

           :returns: List of strings of analysis identifiers.
        """

        re = []
        # Iterate through all groups of the root folder
        for it in self.experiment.items():
            if it[0].startswith(omsi_format_analysis.analysis_groupname):
                cur_ana_id = omsi_file_analysis(
                    self.experiment[it[0]]).get_analysis_identifier()
                if cur_ana_id is not None:
                    re.append(cur_ana_id[0])
                else:
                    re.append("")

        return re

    def get_analysis(self, analysis_index):
        """Get the omsi_format_analysis analysis object for the experiment with
           the given index.

           :param analysis_index: The index of the analysis
           :type analysis_index: Unsigned integer

           :returns: omsi_file_analysis object for the requested analysis. The function \
                          returns None in case the analysis object was not found.
           """
        try:
            return omsi_file_analysis(
                self.experiment[unicode(omsi_format_analysis.analysis_groupname + str(analysis_index))])
        except:
            return None

    def get_analysis_by_identifier(self, analysis_identifier_string):
        """Get the omsi_format_analysis analysis object for the the analysis with the given identifier.

           :param analysis_identifier_string: The string used as identifier for the analysis.
           :type analysis_identifier_string: string

           :returns: h5py obejct of the analysis or None in case the analysis is not found."""

        # Iterate through all groups of the root folder
        for it in self.experiment.items():
            if it[0].startswith(omsi_format_analysis.analysis_groupname):
                cur_ana_id = omsi_file_analysis(
                    self.experiment[it[0]]).get_analysis_identifier()
                if cur_ana_id is not None:
                    if cur_ana_id[0] == analysis_identifier_string:
                        return omsi_file_analysis(self.experiment[it[0]])
        return None

    ###########################################################
    # Add new data to the experiment
    ###########################################################
    def __create_msidata_group__(self):

        nummsi = self.get_num_msidata()
        dataindex = nummsi
        data_group = self.experiment.require_group(
            omsi_format_msidata.data_groupname + str(dataindex))
        data_group.attrs[
            omsi_format_common.type_attribute] = "omsi_file_msidata"
        data_group.attrs[
            omsi_format_common.version_attribute] = omsi_format_msidata.current_version
        data_group.attrs[
            omsi_format_common.timestamp_attribute] = str(time.ctime())
        return data_group

    def create_msidata_full_cube(self, data_shape, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, flush_io=True):
        """Create a new mass spectrometry imaging dataset for the given experiment written as a full 3D cube.

           :param data_shape: Shape of the dataset. Eg. shape=(10,10,10) creates a 3D dataset with
                             10 entries per dimension
           :param data_type:  numpy style datatype to be used for the dataset.
           :param mzdata_type: numpy style datatype to be used for the mz data array.
           :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk
                           sizes to be used in x,y, and m/z explicitly.
           :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                              Can also use an integer in range(10) indicating gzip.
           :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

           :returns: The following two empty (but approbriately sized) h5py datasets are returned in order
                     to be filled with data:

                * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

           :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the
        # group
        data_group = self.__create_msidata_group__()
        data_dataset, mz_dataset = omsi_file_msidata.__create_msidata_full_cube__(
            data_group=data_group, data_shape=data_shape,
            data_type=data_type, mzdata_type=mzdata_type,
            chunks=chunks, compression=compression, compression_opts=compression_opts)
        if flush_io:
            self.experiment.file.flush()
        return data_dataset, mz_dataset, data_group

    def create_msidata_partial_cube(self, data_shape, mask, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, flush_io=True):
        """Create a new mass spectrometry imaging dataset for the given experiment written as a partial 3D cube
           of complete spectra.

           :param data_shape: Shape of the dataset. Eg. shape=(10,10,10) creates a 3D dataset with
                              10 entries per dimension
           :param mask: 2D boolean NumPy array used as mask to indicate which (x,y) locations have
                        spectra associated with them.
           :param data_type:  numpy style datatype to be used for the dataset.
           :param mzdata_type: numpy style datatype to be used for the mz data array.
           :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk sizes
                           to be used in x,y, and m/z explicitly.
           :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                              Can also use an integer in range(10) indicating gzip.
           :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has
                           been written to file

           :returns: The following two empty (but approbriately sized) h5py datasets are returned in order to
                     be filled with data:

                * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

           :returns: The following already complete dataset

                * ``xy_index_dataset`` : This dataset indicates for each xy location to which index in ``data_dataset`` the \
                                         location corresponds to. This dataset is needed to identify where spectra need to be written to.

           :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the
        # group
        data_group = self.__create_msidata_group__()
        data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset = omsi_file_msidata.__create_msidata_partial_cube__(
            data_group=data_group,
            data_shape=data_shape,
            mask=mask,
            data_type=data_type,
            mzdata_type=mzdata_type,
            chunks=chunks,
            compression=compression,
            compression_opts=compression_opts)
        if flush_io:
            self.experiment.file.flush()
        return data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset, data_group

    def create_msidata_partial_spectra(self, spectra_length, len_global_mz, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, flush_io=True):
        """Create a new mass spectrometry imaging dataset for the given experiment written as a partial 3D
           cube of partial spectra.

           :param spectra_length: 2D boolean NumPy array used indicating for each (x,y) locations the length
                    of the corresponding partial spectrum.
           :param len_global_mz: The total number of m/z values in the global m/z axis for the full 3D cube
           :param data_type: The dtype for the MSI dataset
           :param mzdata_type: The dtype for the mz dataset
           :param mzdata_type: numpy style datatype to be used for the mz data array.
           :param chunks:  Specify whether chunkning should be used (True,False), or specify the chunk sizes to
                           be used in x,y, and m/z explicitly.
           :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                              Can also use an integer in range(10) indicating gzip.
           :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data
                            has been written to file

           :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to be
                filled with data:

                * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * ``mz_index_dataset`` : The
                * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

           :returns: The following already complete dataset

                * ``xy_index_dataset`` : This dataset indicates for each xy location at which index in ``data_dataset``
                                the corresponding spectrum starts. This dataset is needed to identify where spectra
                                need to be written to.
                * ``xy_index_end_dataset`` : This dataset indicates for each xy location at which index in
                                ``data_dataset`` the corresponding spectrum ends (exclusing the given value).
                                This dataset is needed to identify where spectra need to be written to.

           :returns: ``data_group`` : The h5py object with the group in the HDF5 file where the data should be stored.
        """
        # Create the HDF5 group and initalize the OMSI object for managing the
        # group
        data_group = self.__create_msidata_group__()
        data_dataset, mz_index_dataset, mz_dataset, xy_index_start_dataset, xy_index_end_dataset = omsi_file_msidata.__create_msidata_partial_spectra__(
            data_group=data_group,
            spectra_length=spectra_length, len_global_mz=len_global_mz,
            data_type=data_type, mzdata_type=mzdata_type,
            chunks=chunks, compression=compression, compression_opts=compression_opts)
        if flush_io:
            self.experiment.file.flush()
        return data_dataset, mz_index_dataset, mz_dataset, xy_index_start_dataset, xy_index_end_dataset, data_group

    def create_method_info(self, method_name=None, flush_io=True):
        """Add information about the method imaged to the experiment

           :param method_name: Optional name of the method
           :type method_name: string, None
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

           :returns: h5py object of the newly created method group.
        """
        method_group = self.experiment.require_group(omsi_format_methods.methods_groupname)
        method_group.attrs[omsi_format_common.type_attribute] = "omsi_file_methods"
        method_group.attrs[omsi_format_common.version_attribute] = omsi_format_methods.current_version
        method_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            self.experiment.file.flush()
        return omsi_file_methods.__create_method_info__(method_group=method_group, method_name=method_name)

    def create_instrument_info(self, instrument_name=None, mzdata=None, flush_io=True):
        """Add information about the instrument used for creating the images for this experiment.

           :param instrument_name: The name of the instrument
           :type instrument_name: string, None
           :param mzdata: Numpy array of the mz data values of the instrument
           :type mzdata: numpy array or None
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

           :returns: The function returns the h5py HDF5 handler to the instrument info group created for the experiment.

        """
        # Create the group for instrument specific data
        instrument_group = self.experiment.require_group(
            omsi_format_instrument.instrument_groupname)
        instrument_group.attrs[
            omsi_format_common.type_attribute] = "omsi_file_instrument"
        instrument_group.attrs[
            omsi_format_common.version_attribute] = omsi_format_instrument.current_version
        instrument_group.attrs[
            omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            self.experiment.file.flush()
        return omsi_file_instrument.__create_instrument_info__(instrument_group=instrument_group, instrument_name=instrument_name, mzdata=mzdata)

    def create_analysis(self, analysis, flush_io=True):
        """Add a new group for storing derived analysis results for the current experiment

          Create the analysis group and use omsi_file_analysis.__create_analysis__(...) to populate the
          group with the approbriate data.

          :param analysis: Instance of omsi.analysis.omsi_analysis_base defining the analysis
          :type analysis: omsi.analysis.omsi_analysis_base:
          :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

          :returns: The omsi_file_analysis object for the newly created analysis group and the integer index of the analysis
         """
        from omsi.analysis.omsi_analysis_base import omsi_analysis_base
        if not isinstance(analysis, omsi_analysis_base):

            errormessage = "Could not write the analysis. The input object was of type " + \
                str(type(analysis)) + " not of omsi_analysis_base as expected."
            raise NameError(errormessage)

        # Create a new group for the analysis data
        numana = self.get_num_analysis()
        anaindex = numana
        analysis_group = self.experiment.require_group(
            omsi_format_analysis.analysis_groupname + str(anaindex))
        analysis_group.attrs[
            omsi_format_common.type_attribute] = "omsi_file_analysis"
        analysis_group.attrs[
            omsi_format_common.version_attribute] = omsi_format_analysis.current_version
        analysis_group.attrs[
            omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            self.experiment.file.flush()
        return omsi_file_analysis.__create_analysis__(analysis_group, analysis), anaindex


##########################################################################
##########################################################################
#   The main class for method information within omsi data files
##########################################################################
##########################################################################
class omsi_file_methods(object):

    """Class for managing method specific data"""

    @classmethod
    def __create_method_info__(cls, method_group, method_name=None):
        """Add information about the method  imaged to the experiment

           NOTE: This is a private helper function used to populate the given group for the \
           method data. Use the corresponding omsi_file_experiment.create_method_info(...) \
           to create a new method group and populate it with data.

           :param method_group: h5py group object that should be populated with the method data.
           :param method_name: Optional name of the method
           :type method_name: string, None

           :returns: h5py object of the newly created method group.
        """
        methodnamedata = method_group.require_dataset(name=unicode(
            omsi_format_methods.methods_name), shape=(1,), dtype=omsi_format_common.str_type)
        if method_name is None:
            if len(methodnamedata[0]) == 0:
                methodnamedata[0] = "undefined"
        else:
            if omsi_format_common.str_type_unicode:
                methodnamedata[0] = method_name
            else:
                methodnamedata[0] = str(method_name)
        return omsi_file_methods(method_group)

    def __init__(self, methods_group):
        """Initialize the method object given the h5py object of the method group

           :param methods_group: The h5py object with the method group of the omsi hdf5 file.
        """
        self.method = methods_group
        self.name = self.method.name

    def __getitem__(self, key):
        """Support direct read interaction with the method h5py group"""
        return self.method[key]

    def __setitem__(self, key, value):
        """Support direct write interaction with the method h5py group"""
        # If the version of h5py does not support automatic storage of strings and we have a string then do-it-yourself
        if (isinstance(value, unicode) or isinstance(value, str)) and not omsi_format_common.str_type_unicode:
            d = self.method.require_dataset(name=unicode(key),
                                                shape=(1,),
                                                dtype=omsi_format_common.str_type)
            d[0] = str(value)
        else:
            self.method[key] = value

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.method.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.method.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the method group was created in the HDF5 file.

           :returns: Python timestamp string generated using time.ctime().
                     None may be returned in case that the timestamp does not exists
                     or cannot be retrieved from the file for some reason.

        """
        try:
            return self.method.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def get_h5py_methodgroup(self):
        """Get the h5py object for the HDF5 group of the method data

           :returns: h5py reference for the method group.
        """
        return self.method

    def get_method_name(self):
        """Get the HDF5 dataset with the name of the method.

           To retrieve the name string use get_method_name()[...]

           :returns: h5py object where the method name is stored.
                     Returns None in case no method name is found.
        """
        method_name = None
        if self.method is None:
            return None
        try:
            method_name = self.method[unicode(omsi_format_methods.methods_name)]
        except:
            method_name = None
        return method_name

    def set_method_name(self, name_string):
        """Overwrite the name string for the method with the given name string

           :param name_string: The new method name.
           :type name_string: string
        """

        # Get the name of the intrument
        name = self.get_method_name()
        # Create the dataset for the instrument name if it does not exist
        if name is None:
            name = self.method.require_dataset(name=unicode(
                omsi_format_methods.methods_name), shape=(1,), dtype=omsi_format_common.str_type)

        if omsi_format_common.str_type_unicode:
            name[0] = name_string
        else:
            name[0] = str(name_string)


##########################################################################
##########################################################################
#   The main class for instrument information within omsi data files
##########################################################################
##########################################################################
class omsi_file_instrument(object):

    """Class for managing instrument specific data"""

    @classmethod
    def __create_instrument_info__(cls, instrument_group, instrument_name=None, mzdata=None):
        """Populate the empty instrument group with the given data.

           NOTE: This is a private helper function used to populate the instrument group with data. \
           Use the corresponding omsi_file_experiment.create_instrument_info(...) function to  \
           generate a new instrument information data in HDF5

           :param instrument_group: The h5py group to which the instrument group data should be written to.
           :param instrument_name: The name of the instrument used.
           :param mzdata: The mz data for the instrument.
        """
        # Name of the instrument
        instrumentnamedata = instrument_group.require_dataset(name=unicode(
            omsi_format_instrument.instrument_name), shape=(1,), dtype=omsi_format_common.str_type)
        if instrument_name is None:
            if len(instrumentnamedata[0]) == 0:
                instrumentnamedata[0] = "undefined"
        else:
            if omsi_format_common.str_type_unicode:
                instrumentnamedata[0] = instrument_name
            else:
                instrumentnamedata[0] = str(instrument_name)
        # MZ data for the instrument
        if mzdata is not None:
            instrumentmzdata = instrument_group.require_dataset(
                name=omsi_format_instrument.instrument_mz_name, shape=mzdata.shape, dtype=mzdata.dtype)
            instrumentmzdata[:] = mzdata[:]

        return omsi_file_instrument(instrument_group)

    def __init__(self, instrument_group):
        """Initalize the instrument object given the h5py object of the instrument group

          :param instrument_group: The h5py object with the instrument group of the omsi hdf5 file.
        """
        self.instrument = instrument_group
        self.name = self.instrument.name

    def __getitem__(self, key):
        """Support direct read interaction with the instrument h5py group"""
        return self.instrument[key]

    def __setitem__(self, key, value):
        """Support direct write interaction with the instrument h5py group"""
        # If the version of h5py does not support automatic storage of strings and we have a string then do-it-yourself
        if (isinstance(value, unicode) or isinstance(value, str)) and not omsi_format_common.str_type_unicode:
            d = self.instrument.require_dataset(name=unicode(key),
                                                shape=(1,),
                                                dtype=omsi_format_common.str_type)
            d[0] = str(value)
        else:
            self.instrument[key] = value

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.instrument.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.instrument.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the instrument group was created in the HDF5 file.

            :returns: Python timestamp string generated using time.ctime().
                      None may be returned in case that the timestamp does not exists
                      or cannot be retrieved from the file for some reason.

            """
        try:
            return self.instrument.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def get_h5py_instrumentgroup(self):
        """Get the h5py object with the instrument group.

           :returns: h5py object of the instrument group.
         """
        return self.instrument

    def get_instrument_name(self):
        """Get the HDF5 dataset with the name of the instrument.

           To get the string of the instrument name use:
           get_instrument_name()[...]

           :returns: h5py object to the dataset with the instrument name.
                     Returns None in case no method name is found.
        """
        instrument_name = None
        if self.instrument is None:
            return None
        try:
            instrument_name = self.instrument[
                unicode(omsi_format_instrument.instrument_name)]
        except:
            instrument_name = None
        return instrument_name

    def get_instrument_mz(self):
        """Get the HDF5 dataset with the mz data for the instrument.

           To get the numpy array of the full mz data use:
            get_instrument_mz()[:]

           :returns: Returns the h5py object with the instrument mz data.
                     Returns None in case no mz data was found for the instrument.
        """
        if self.instrument is None:
            return None
        try:
            return self.instrument[
                unicode(omsi_format_instrument.instrument_mz_name)]
        except:
            return None

    def set_instrument_name(self, name):
        """Overwrite the current identfier string for the experiment with the given string.

           :param name: The new instrument name.
           :type name: string.
        """

        # Get the name of the intrument
        namedataset = self.get_instrument_name()
        # Create the dataset for the id name if it does not exist
        if namedataset is None:
            namedataset = self.instrument.require_dataset(name=unicode(
                omsi_format_instrument.instrument_name), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            namedataset[0] = name
        else:
            namedataset[0] = str(name)

##########################################################################
##########################################################################
#   The main class for analysis data within omsi data files                    #
##########################################################################
##########################################################################


class omsi_file_analysis(object):

    """Class for managing analysis specific data in omsi hdf5 files"""

    @classmethod
    def __create_analysis__(cls, analysis_group, analysis):
        """Populate the given h5py group with the analysis data.

            NOTE: This is a private helper function. Use the corresponding create_analysis function \
            of omsi_file_experiment to create a completely new analysis.

            :param analysis_group: h5py group in which the analysis data should be stored.
            :param analysis: Instance of omsi.analysis.omsi_analysis_base defining the analysis
            :type analysis: omsi.analysis.omsi_analysis_base:

            :returns: The omsi_file_analysis object for the newly created analysis group. The analysis data is \
                      automatically written to file by this function so no addition work is required.

        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data

        # Write the analysis name
        analysisidentifierdata = analysis_group.require_dataset(name=unicode(
            omsi_format_analysis.analysis_identifier), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysisidentifierdata[0] = analysis.get_analysis_identifier()
        else:
            analysisidentifierdata[0] = str(analysis.get_analysis_identifier())

        # Write the analysis type
        analysistypedata = analysis_group.require_dataset(name=unicode(
            omsi_format_analysis.analysis_type), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            analysistypedata[0] = analysis.get_analysis_type()
        else:
            analysistypedata[0] = str(analysis.get_analysis_type())

        # Write the analysis data
        for a in analysis.get_all_analysis_data():

            cls.__write_omsi_analysis_data__(analysis_group, a)

        # Write all the parameters
        parameter_group = analysis_group.require_group(
            omsi_format_analysis.analysis_parameter_group)
        for p in analysis.get_all_parameter_data():

            cls.__write_omsi_analysis_data__(parameter_group, p)

        # Write all the runtime execution information
        runinfo_group = analysis_group.require_group(
            omsi_format_analysis.analysis_runinfo_group)
        for k, v in analysis.get_all_run_info().items():
            # Generate an omsi_analysis_data object in order to use the
            # __write_omsi_analysis_data function to write the data
            if isinstance(v, unicode) or isinstance(v, str):
                anadata = omsi_analysis_data(name=unicode(k),
                                             data=v,
                                             dtype=omsi_format_common.str_type)
            else:
                dat = np.asarray(v)
                if len(dat.shape) == 0:
                    dat = np.asarray([v])
                anadata = omsi_analysis_data(name=unicode(k),
                                             data=dat,
                                             dtype=dat.dtype)
            cls.__write_omsi_analysis_data__(runinfo_group, anadata)

        # Write all data dependencies
        dependencies_group = analysis_group.require_group(
            omsi_format_dependencies.dependencies_groupname)
        dependencies_group.attrs[
            omsi_format_common.type_attribute] = "omsi_file_dependencies"
        dependencies_group.attrs[
            omsi_format_common.version_attribute] = omsi_format_dependencies.current_version
        dependencies_group.attrs[
            omsi_format_common.timestamp_attribute] = str(time.ctime())
        omsidep = omsi_file_dependencies.__create_dependencies__(
            dependencies_group=dependencies_group,
            dependencies_data_list=analysis.get_all_dependency_data())

        # Execute the custom data write for the analysis
        analysis.write_to_omsi_file(analysis_group)

        return omsi_file_analysis(analysis_group)

    @classmethod
    def __write_omsi_analysis_data__(cls, data_group, ana_data):
        """Private helper function used to write the data defined by a omsi_analysis_data object to HDF5.

            :param data_group: The h5py data group to which the data should be written to.
            :param ana_data: The omsi_analysis_data object with the description of the data to be written.
            :type ana_data: omsi.analysis.omsi_analysis_data
        """

        # Create link in HDF5 to an existing dataset within the file
        if isinstance(ana_data['dtype'], int):
            if ana_data['dtype'] == ana_data.ana_hdf5link:
                linkobject = data_group.file.get(ana_data['data'])
                data_group[ana_data['name']] = linkobject
                omsiobj = omsi_file.get_omsi_object(linkobject)
                try:
                    # Check if we already have a type attribute
                    a = data_group[ana_data['name']].attrs[
                        omsi_format_common.type_attribute]
                except:
                    # Generate the type attribute from scratch
                    if omsiobj is not None:
                        omsiobjtype = omsiobj.__class__.__name__
                    else:
                        omsiobjtype = ""
                    data_group[ana_data['name']].attrs[
                        omsi_format_common.type_attribute] = omsiobjtype
        # Create a new string-type dataset
        elif (ana_data['dtype'] == omsi_format_common.str_type) or (ana_data['dtype'] == h5py.special_dtype(vlen=str)):
            tempdata = data_group.require_dataset(name=unicode(ana_data['name']),
                                                  shape=(1,),
                                                  dtype=omsi_format_common.str_type)
            if len(unicode(ana_data['data'])) > 0:
                if omsi_format_common.str_type_unicode:
                    tempdata[0] = ana_data['data']
                else:
                    tempdata[0] = str(ana_data['data'])
            else:
                print "WARNGING: " + ana_data['name'] + " dataset generated but not written. The given dataset was empty."
        # Create a new dataset to store the current numpy-type dataset
        elif 'numpy' in str(type(ana_data['data'])):
            # Decide whether we want to enable chunking for the current
            # analysis dataset
            chunks = None
            if ana_data['data'].size > 1000:
                chunks = True
            # Write the current analysis dataset
            tempdata = data_group.require_dataset(
                name=ana_data['name'], shape=ana_data['data'].shape, dtype=ana_data['dtype'], chunks=chunks)
            if ana_data['data'].size > 0:
                tempdata[:] = ana_data['data']
            else:
                print "WARNGING: " + ana_data['name'] + " dataset generated but not written. The given dataset was empty."
        # Unkown dtype. Attempt to convert the dataset to numpy and write it to
        # file.
        else:
            print "WARNING: " + str(ana_data['name']) + ": The data specified by the analysis object is not in numpy format. Attempting to convert the data to numpy"
            try:
                dat = np.asarray(ana_data['data'])
                if len(dat.shape) == 0:
                    dat = np.asarray([ana_data['data']])
                tempdata = data_group.require_dataset(
                    name=ana_data['name'], shape=dat.shape, dtype=str(dat.dtype))
                if dat.size > 0:
                    tempdata[:] = dat
                else:
                    print "WARNGING: " + ana_data['name'] + " dataset generated but not written. The given dataset was empty."
            except:
                print "ERROR: " + str(ana_data['name']) + ": The data specified by the analysis could not be converted to numpy for writing to HDF5"

    def __init__(self, analysis_group):
        """Initalize the analysis object given the h5py object of the analysis group.

           :param analysis_group: The h5py object with the analysis group of the omsi hdf5 file.

        """
        self.analysis = analysis_group
        self.parameter = self.analysis[
            unicode(omsi_format_analysis.analysis_parameter_group)]
        try:
            self.dependencies = omsi_file_dependencies(
                self.analysis[unicode(omsi_format_dependencies.dependencies_groupname)])
        except:
            print "WARNGING: omsi_file_analysis.__init__ : No dependencies group found."
            self.dependencies = None
        self.analysis_omsi_object = None
        self.name = self.analysis.name

    def __getitem__(self, key):
        """Support direct read interaction with the analysis h5py group"""
        try:
            return self.analysis[key]
        except:
            try:
                return self.parameter[key]
            except:
                try:
                    return self.dependencies[key]
                except:
                    return None
        # No errors have occured but no matching object has been found either
        return None

    def __setitem__(self, key, value):
        """Support direct write interaction with the analysis h5py group"""
        try:
            self.analysis[key] = value
        except:
            try:
                self.parameter[key] = value
            except:
                try:
                    if key in self.dependencies.items():
                        raise KeyError(
                            "Assignment to dependcies is not permitted via this mechanism")
                    else:
                        raise KeyError(
                            "Requested assignment operation failed")
                except:
                    raise KeyError("Requested assignment operation failed")

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.analysis.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.analysis.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the analysis group was created in the HDF5 file.

            :returns: Python timestamp string generated using time.ctime().
                      None may be returned in case that the timestamp does not exists
                      or cannot be retrieved from the file for some reason.

            """
        try:
            return self.analysis.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def get_h5py_analysisgroup(self):
        """Retrun the h5py object with the analysis data.

           The returned object can be used to read data directly from the HDF5 file.
           Write operations to the analysis group can be performed only if the
           associated omsi_file was opened with write permissions.

           :returns: h5py object for the analysis group.
        """
        return self.analysis

    def get_analysis_identifier(self):
        """Get the identifier name of the analysis.

           Use get_analysis_identifier()[...] to retrive the identifier string.

           :returns: h5py object for the dataset with the identifier string.
                     Returns None, in case no identifer exisits. This should
                     not be the case for a valid OpenMSI file.
        """
        if self.analysis is None:
            return None

        try:
            return self.analysis[unicode(omsi_format_analysis.analysis_identifier)]
        except:
            return None

    def get_analysis_type(self):
        """Get the type for the analysis.

           Use get_analysis_type()[...] tor retrieve the type string.

           :returns: h5py object with the dataset of the analysis string. Returns,
                     None in case no analysis type exists. This should not be the
                     case in a valid omsi file.
        """
        if self.analysis is None:
            return None
        try:
            return self.analysis[unicode(omsi_format_analysis.analysis_type)]
        except:
            return None

    def get_analysis_data_names(self):
        """This function returns all dataset names (and groups) that are custom
           to the analysis, i.e., that are not part of the omsi file standard.

           :returns: List of analysis-specific dataset names.
        """
        re = []
        if self.analysis is not None:
            for it in self.analysis.items():
                if it[0] != omsi_format_analysis.analysis_identifier and \
                        it[0] != omsi_format_analysis.analysis_type and \
                        it[0] != omsi_format_analysis.analysis_parameter_group and \
                        it[0] != omsi_format_dependencies.dependencies_groupname and \
                        it[0] != omsi_format_analysis.analysis_runinfo_group:
                    re.append(it[0])
        return re

    def get_analysis_data_shapes_and_types(self):
        """This function returns two dictionaries with all dataset names (and groups)
           that are custom to the analysis, i.e., that are not part of the omsi
           file standard, and idenifies the shape of the analysis data objects.

           :returns: Dictonary indicating for each analysis-specific dataset its name (key) and
                     shape (value). And a second dictionariy indicating the name (key) and dtype
                     of the dataset.
        """
        re = {}
        ret = {}
        if self.analysis is not None:
            for it in self.analysis.items():
                if it[0] != omsi_format_analysis.analysis_identifier and \
                        it[0] != omsi_format_analysis.analysis_type and \
                        it[0] != omsi_format_analysis.analysis_parameter_group and \
                        it[0] != omsi_format_dependencies.dependencies_groupname and \
                        it[0] != omsi_format_analysis.analysis_runinfo_group:
                    try:
                        re[it[0]] = self.analysis[unicode(it[0])].shape
                    except:
                        re[it[0]] = (0,)
                    try:
                        ret[it[0]] = self.analysis[unicode(it[0])].dtype
                    except:
                        ret[it[0]] = "Unkown type"
        return re, ret

    def get_all_analysis_data(self, load_data=False):
        """Get all analysis data associated with the analysis.

           :param load_data: load_data: Should the data be loaded or just the h5py objects be
                             stored in the dictionary.

           :returns: List of omsi_analysis_data objects with the names and h5py or numpy objects.
                     Access using [index]['name'] and [index]['data'].
        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data
        re = []
        if self.analysis is not None:
            for it in self.analysis.items():
                if it[0] != omsi_format_analysis.analysis_identifier and \
                        it[0] != omsi_format_analysis.analysis_type and \
                        it[0] != omsi_format_analysis.analysis_parameter_group and \
                        it[0] != omsi_format_dependencies.dependencies_groupname and \
                        it[0] != omsi_format_analysis.analysis_runinfo_group:
                    re.append(omsi_analysis_data())
                    re[-1]['name'] = str(it[0])
                    if load_data:
                        re[-1]['data'] = self.analysis[unicode(it[0])][:]
                    else:
                        re[-1]['data'] = self.analysis[unicode(it[0])]
                    re[-1]['dtype'] = str(re[-1]['data'].dtype)
        return re

    def get_all_parameter_data(self, load_data=False):
        """Get all parameter data associated with the analysis.

           :param load_data: Should the data be loaded or just the h5py objects be stored in the dictionary.

           :returns: List of omsi_analysis_data objects with names and h5py or numpy object. Access using
                     [index]['name'] and [index]['data'].
        """
        from omsi.analysis.omsi_analysis_data import omsi_analysis_data
        re = []
        if self.parameter is not None:
            for it in self.parameter.items():
                re.append(omsi_analysis_data())
                re[-1]['name'] = str(it[0])
                if load_data:
                    re[-1]['data'] = self.parameter[unicode(it[0])][:]
                else:
                    re[-1]['data'] = self.parameter[unicode(it[0])]
                re[-1]['dtype'] = str(re[-1]['data'].dtype)
        return re

    def get_all_dependency_data(self, omsi_dependency_format=True):
        """Get all direct dependencies associdated with the analysis.

           :param omsi_dependency_format: Should the dependcies be retrieved as omsi_analysis_dependency
                    object (True) or as an omsi_file_dependencydata object (False).

           :returns: List omsi_dependency objects containing either omsi file API objects or h5py objects
                    for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def get_all_dependency_data_recursive(self, omsi_dependency_format=True):
        """Get all direct and indirect dependencies associdated with the analysis.

           NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

           :param omsi_dependency_format: Should the dependencies be retrieved as omsi_dependency object
                        (True) or as an omsi_file_dependencydata object (False)

           :returns: List omsi_analysis_data objects containing either omsi file API interface objects or
                        h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data_recursive(omsi_dependency_format=omsi_dependency_format)
        else:
            return []


##########################################################################
##########################################################################
#   The main class for managing collections of dependencies                    #
##########################################################################
##########################################################################
class omsi_file_dependencies(object):

    """Class for managing collections of dependencies"""

    @classmethod
    def __create_dependencies__(cls, dependencies_group, dependencies_data_list=None):
        """Populate the given h5py group with the dependencies data.

            NOTE: This is a private helper function. Use the corresponding create_dependencies function \
            of omsi_file_analysis and omsi_file_msidata to create a completely new collection of dependencies.

            :param dependencies_group: h5py.Group object that should be used to store the dependencies data.
            :type dependencies_group: h5py.Group
            :param dependencies_data_list: List of dependencies to be stored
            :type dependencies_data_list: List of omsi.shared.omsi_dependency_data objects to be stored.
        """
        omsiobj = omsi_file_dependencies(dependencies_group)
        if dependencies_data_list:
            for d in dependencies_data_list:
                omsiobj.add_dependency(dependency_data=d)
        return omsiobj

    def add_dependency(self, dependency_data):
        """Add a new dependency to the collection.

           :param dependency_data: The analysis dependency specification.
           :type dependency_data: omsi.shared.omsi_dependency_data
        """
        depname = dependency_data['link_name']
        # Check if the group alreay exists and raise an error if it does
        try:
            temp = self.dependencies[depname]
            raise KeyError(
                "Another dependencies with the same index already exists.")
        except:
            pass
        new_dep_group = self.dependencies.require_group(depname)
        new_dep_group.attrs[
            omsi_format_common.type_attribute] = "omsi_file_dependencydata"
        new_dep_group.attrs[
            omsi_format_common.version_attribute] = omsi_format_dependencydata.current_version
        new_dep_group.attrs[
            omsi_format_common.timestamp_attribute] = str(time.ctime())
        depobj = omsi_file_dependencydata.__create_dependency__(
            data_group=new_dep_group, dependency_data=dependency_data)
        return depobj

    def __init__(self, dependencies_group):

        self.dependencies = dependencies_group
        self.name = self.dependencies.name

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.dependencies.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.dependencies.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the dependency group was created in the HDF5 file.

            :returns: Python timestamp string generated using time.ctime().
                      None may be returned in case that the timestamp does not exists
                      or cannot be retrieved from the file for some reason.

            """
        try:
            return self.dependencies.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def get_h5py_dependenciesgroup(self):
        """Retrun the h5py object with the dependency data.

           The returned object can be used to read data directly from the HDF5 file.
           Write operations to the analysis group can be performed only if the
           associated omsi_file was opened with write permissions.

           :returns: h5py object for the dependency group.
        """
        return self.dependencies

    def __getitem__(self, key):
        """Retrieve the h5py dataset for any of the dependency datasets or the omsi file API object.
           If the key is a numpy selection (slicing) then the data is retrieved from the link the
           dependency is pointing so that we can interact with the dependency as if it were the
           linked object.
        """
        return self.dependencies[key]

    def get_omsi_file_dependencydata(self, name):
        """Retrieve the omsi_file_dependencydata object for the dependency with the given name"""
        return omsi_file.get_omsi_object(self.dependencies[name])

    def get_dependency_omsiobject(self, name, recursive=True):
        """Get the omsi file API object corresponding to the object the dependency is pointing to.

           :param name: Name of the dependency opbject to be loaded .
           :param recursive: Should the dependency be resolved recursively, i.e., if the dependeny points
                    to another dependencies. Default=True.

           :returns: An omsi file API object (e.g., omsi_file_analysis or omsi_file_msidata) if the link points
                    to a group or the h5py.Dataset the link is pointing to.
        """
        return self.get_omsi_file_dependencydata(name).get_dependency_omsiobject(recursive=recursive)

    def get_all_dependency_data(self, omsi_dependency_format=True):
        """Get all direct dependencies associdated with the analysis.

           :param omsi_dependency_format: Should the dependcies be retrieved as omsi_analysis_dependency object
                    (True) or as an omsi_file_dependencydata object (False).

           :returns: List omsi_dependency objects containing either omsi file API objects or h5py objects for
                    the dependcies. Access using [index]['name'] and [index]['data'].
        """
        re = []
        for it in self.dependencies.items():
            # try :
            omsi_object = omsi_file_dependencydata(
                self.dependencies[unicode(it[0])])
            if omsi_dependency_format:
                re.append(omsi_object.get_omsi_dependency())
            else:
                re.append(omsi_object)
            # except :
            #    pass
        return re

    def get_all_dependency_data_recursive(self, omsi_dependency_format=True):
        """Get all direct and indirect dependencies associdated with the analysis.

           NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

           :param omsi_dependency_format: Should the dependencies be retrieved as omsi_dependency object (True)
                    or as an omsi_file_dependencydata object (False)

           :returns: List omsi_analysis_data objects containing either omsi file API interface objects or
                    h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        from omsi.shared.omsi_dependency import omsi_dependency

        re = []
        for it in self.dependencies.items():
            omsi_obj = omsi_file.get_omsi_object(self.dependencies[it[0]])
            # If we can have recusive dependencies
            if isinstance(omsi_obj, omsi_file_analysis):
                it_depend = omsi_obj.get_all_dependency_data_recursive(
                    omsi_dependency_format=omsi_dependency_format)
                re = re + it_depend
            try:
                omsi_object = omsi_file_dependencydata(
                    self.dependencies[unicode(it[0])])
                if omsi_dependency_format:
                    re.append(omsi_object.get_omsi_dependency())
                else:
                    re.append(omsi_object)
            except:
                import sys
                print "WARNING: Error occured in omsi_file_dependencies::get_all_dependency_data_recursive(...):  " + \
                      unicode(it[0]) + "   :" + str(sys.exc_info())
        return re

    def get_all_dependency_data_graph(self, include_omsi_dependency=False, include_omsi_file_dependencydata=False, recursive=True, level=0, name_key='name', prev_nodes=None, prev_links=None, parent_index=None, metadata_generator=None, metadata_generator_kwargs=None):
        """ Get all direct and indirect dependencies associdated with the analysis in form of a graph describing
            all nodes and links in the provenance hierarchy.

            NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

            :param include_omsi_dependency: Should the omsi_dependency object be included in the entries
                        in the nodes dict?
            :param include_omsi_file_dependencydata: Should the omsi_file_dependencydata object be included in
                        the entries in the nodes dict?
            :param recursive: Should we trace dependencies recursively to construct the full graph, or only the
                        direct dependencies. Default true (ie., trace recursively)
            :param name_key: Which key should be used in the dicts to indicate the name of the object?
                        Default value is 'name'
            :param level: Integer used to indicated the recursion level. Default value is 0.
            :param prev_nodes: List of nodes that have been previously generated. Note, this list will be
                        modified by the call. Note, each node is represented by a dict which is expected to
                        have at least the following keys defined, path, name_key, level (name_key refers to the key
                        defined by the input parameter name_key.
            :param prev_links: Previouly established links in the list of nodes. Note, this list will be
                        modified by the call.
            :param parent_index: Index of the parent node in the list of prev_nodes for this call.
            :param metadata_generator: Optional parameter. Pass in a function that generates additional metadata
                        about a given omsi API object. Note, the key's level and path and name (i.e., name_key) are
                        already set by this function. The metadata_generator may overwrite these key's, however,
                        the path has to be unique as it is used to identify duplicate nodes. Overwriting the path
                        with a non-unique value, hence, will lead to errors (missing entries) when generating the graph.
                        Note, the metadata_generator function needs to support the following keyword arguments:

                              * inDict : The dictionary to which the metadata should be added to.
                              * obj : The omsi file API object for which metadata should be generated
                              * name : A qualifying name for the object
                              * name_key : The key to be used for storing the name
            :param metadata_generator_kwargs: Dictionary of additional keyword arguments that should be passed to
                        the metadata_generator function.

            :returns: Dictionary containing two lists. 1) nodes : List of dictionaries, describing the elements
                      in the dependency graph. 2) links : List of tuples with the links in the graph. Each
                      tuple consists of two integer indices for the nodes list. For each node the following
                      entries are given:

                    * omsi_dependency: Optional key used to store the corresponding omsi_dependency object.
                            Used only of include_omsi_dependency is True.
                    * omsi_file_dependencydata: Optional key used to store the corresponding
                            omsi_file_dependencydata object. Used only of include_omsi_file_dependencydata is True.
                    * name : Name of the dependency. The actual key is sepecified by name_key
                    * level : The recursion level at which the object occurs.
                    * ... : Any other key/value pairs from the omsi_dependency dict.

        """
        from omsi.shared.omsi_dependency import omsi_dependency
        import os
        if prev_nodes:
            nodes = prev_nodes
        else:
            nodes = []
        if prev_links:
            links = prev_links
        else:
            links = []
        if metadata_generator_kwargs:
            metagenkwargs = metadata_generator_kwargs
        else:
            metagenkwargs = {}


        def find_node_index(path_string):
            """Internal helper funtion used to check whether a node already exists in the graph

                :param path_string: The path to object used for comparison.
            """
            if prev_nodes:
                for i in range(len(prev_nodes)):
                    if nodes[i]['path'] == path_string:
                        return i
            return None

        # Iterate through all dependencies
        for it in self.dependencies.items():
            # print "______________________________"
            # print it
            # print "-----A----"+str(links)
            try:
                # 1) Check if a node exists that represents the current
                # dependency
                dep_obj = omsi_file_dependencydata(
                    self.dependencies[unicode(it[0])])
                omsi_obj = dep_obj.get_dependency_omsiobject(recursive=True)
                #omsi_obj = omsi_file.get_omsi_object( self.dependencies[ it[0] ] )
                curr_path = omsi_obj.name
                curr_name = os.path.basename(curr_path)
                curr_index = find_node_index(omsi_obj.name)
                # 2) If no node exists for the current dependency then add a
                # new node
                if not curr_index:
                    # 2.1)Add a new node with the basic metadata
                    nodes.append({})
                    curr_index = len(nodes) - 1
                    curr_node = nodes[curr_index]
                    curr_node['level'] = level
                    curr_node[name_key] = curr_name
                    curr_node['path'] = curr_path
                    # print curr_path
                    if include_omsi_dependency:
                        curr_node[
                            'omsi_dependency'] = dep_obj.get_omsi_dependency()
                    if include_omsi_file_dependencydata:
                        curr_node['omsi_file_dependencydata'] = dep_obj
                    # 2.2) Expand the metadata dict for the current object with
                    # any user-defined metadata
                    if metadata_generator:
                        metadata_generator(inDict=curr_node, obj=omsi_obj, name=os.path.basename(
                            omsi_obj.name), nameKey=name_key, **metagenkwargs)

                # print "PARENT:"+str(parent_index)
                # 3) Add a link from the current node to its parent
                if parent_index is not None:
                    templink = {'source': parent_index, 'target': curr_index}
                    links.append(templink)
                    # print "LINK: "+str(templink)

                # 4) If we can have recusive dependencies that we should trace,
                # then trace them back and add them to the nodes and links
                # lists
                if recursive:
                    if isinstance(omsi_obj, omsi_file_analysis) or isinstance(omsi_obj, omsi_file_msidata):
                        if omsi_obj.dependencies:
                            omsi_obj.dependencies.get_all_dependency_data_graph(
                                include_omsi_dependency=include_omsi_dependency,
                                include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                                level=level +
                                1,
                                name_key=name_key,
                                prev_nodes=nodes,
                                prev_links=prev_links,
                                parent_index=curr_index,
                                metadata_generator=metadata_generator,
                                metadata_generator_kwargs=metagenkwargs
                            )
                # print "-----B----"+str(links)
            except:
                import sys
                print "WARNING: Error occured in omsi_file_dependencies::get_all_dependency_data_graph(...):" + str(sys.exc_info())

        return nodes, links


##########################################################################
##########################################################################
#   The main class for accessing dependency datasets                           #
##########################################################################
##########################################################################
class omsi_file_dependencydata(object):

    """Class for managing data groups used for storing data dependencies"""

    @classmethod
    def __create_dependency__(cls, data_group, dependency_data):
        """Create a new dependency group within the given datagroup and initalize all the required datasets based on the given dependency information.

           :param data_group: The h5py group opbject to which the dependency dataset should be added
           :type data_group: h5py.Group
           :param dependency_data: The analysis dependency specification.
           :type dependency_data: omsi.shared.omsi_dependency_data

        """
        dep_group = data_group
        # Save the name of the parameter
        param_name_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_parameter), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            param_name_data[0] = dependency_data['param_name']
        else:
            param_name_data[0] = str(dependency_data['param_name'])
        # Save the selection
        selection_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_selection), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['selection'] is not None:
            from omsi.shared.omsi_data_selection import check_selection_string
            # This should always be True since omsi_dependency checks for this
            # but we need to be sure.
            if check_selection_string(dependency_data['selection']):
                if omsi_format_common.str_type_unicode:
                    selection_data[0] = dependency_data['selection']
                else:
                    selection_data[0] = str(dependency_data['selection'])
            else:
                selection_data[0] = ""
                errormessage = "Invalid selection string given for data dependency : " + \
                    str(dependency_data['selection'])
                raise ValueError(errormessage)
        else:
            selection_data[0] = ""
        # Save the main omsi object
        mainname_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_mainname), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            mainname_data[0] = unicode(dependency_data['omsi_object'].name)
        else:
            mainname_data[0] = str(dependency_data['omsi_object'].name)
        # Save the additional dataset name
        dataset_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_datasetname), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['dataname']:
            if omsi_format_common.str_type_unicode:
                dataset_data[0] = dependency_data['dataname']
            else:
                dataset_data[0] = str(dependency_data['dataname'])
        else:
            dataset_data[0] = u''

        return omsi_file_dependencydata(dep_group)

    def __init__(self, dependency_group):
        """Create a new omsi_file_dependencydata object for the given h5py.Group

            :param dependency_group: h5py.Group object with the dependency data
            :type dependency_group: h5py.Group with a corresponding omsi type
        """
        self.dependency = dependency_group
        self.name = self.dependency.name

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self.dependency.items()

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self.dependency.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the dependency group was created in the HDF5 file.

            :returns: Python timestamp string generated using time.ctime().
                      None may be returned in case that the timestamp does not exists
                      or cannot be retrieved from the file for some reason.

            """
        try:
            return self.dependency.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def get_h5py_dependencygroup(self):
        """Retrun the h5py object with the dependency data.

           The returned object can be used to read data directly from the HDF5 file.
           Write operations to the analysis group can be performed only if the
           associated omsi_file was opened with write permissions.

           :returns: h5py object for the dependency group.
        """
        return self.dependency

    def __getitem__(self, key):
        """Retrieve the h5py dataset for any of the dependency datasetes or the omsi file API object.
           If the key is a numpy selection (slicing) then the data is retrieved from the link the
           dependency is pointing so that we can interact with the dependency as if it were the
           linked object.
        """
        if isinstance(key, str) or isinstance(key, unicode):
            if key == omsi_format_dependencydata.dependency_mainname or \
               key == omsi_format_dependencydata.dependency_datasetname:
                return omsi_file.get_omsi_object(self.dependency[key])
            else:
                return self.dependency[key][0]
        else:
            omsi_object_name = self.dependency[
                unicode(omsi_format_dependencydata.dependency_mainname)][0]
            h5py_object = self.dependency.file[unicode(omsi_object_name)]
            omsi_object = omsi_file.get_omsi_object(h5py_object)
            dataset_name = self.dependency[
                unicode(omsi_format_dependencydata.dependency_datasetname)][0]
            if len(dataset_name) > 0:
                return omsi_object[dataset_name][key]
            else:
                return omsi_object[key]

    def get_dependency_type(self, recursive=True):
        """Indicated the type of the object the dependency is pointing to.

          :param recursive: Should dependencies be resolved recursively, i.e., if the dependeny points to another dependencies. Default=True.

          :returns: String indicating the class of the omsi file API class that is suited to manage the dependency link or the name of the corresponding h5py class.
        """
        omsi_object = self.get_dependency_omsiobject(recursive=recursive)
        return omsi_object.__class__.__name__

    def get_selection_string(self):
        """String indicating the applied selection. This is an empty string in case no selection was applied.

           :returns: Selection string. See the omsi.shared.omsi_data_selection for helper functions to deal with selection strings.
        """
        return self.dependency[unicode(omsi_format_dependencydata.dependency_selection)][0]

    def get_parameter_name(self):
        """Get the string indicating the name of the dependend parameter of the analysis.

           :returns: String of the parameter name that has the dependency.
        """
        return self.dependency[unicode(omsi_format_dependencydata.dependency_parameter)][0]

    def get_dataset_name(self):
        """Get the string indicating the name of dataset. This may be empty as it is only used if the dependency points to an objec within a managed omsi API object.

            :returns: String indicating the name of the optional dataset.
        """
        try:
            return self.dependency[unicode(omsi_format_dependencydata.dependency_datasetname)][0]
        except:
            return ""

    def get_dependency_omsiobject(self, recursive=True):
        """Get the omsi file API object corresponding to the object the dependency is pointing to.

           :param recursive: Should dependencies be resolved recursively, i.e., if the dependeny points to another dependencies. Default=True.

           :returns: An omsi file API object (e.g., omsi_file_analysis or omsi_file_msidata) if the link points to a group or the h5py.Dataset the link is pointing to.
        """
        omsi_object_name = self.dependency[
            unicode(omsi_format_dependencydata.dependency_mainname)][0]
        h5py_object = self.dependency.file[unicode(omsi_object_name)]
        omsi_object = omsi_file.get_omsi_object(
            h5py_object, resolve_dependencies=recursive)
        return omsi_object

    def get_link_name(self):
        """Get the name of the dependency link

            :returns: String indicating the name of the dependency link.
        """
        return self.dependency.name.split("/")[-1]

    def get_omsi_dependency(self):
        """Get the dependency information as an omsi.shared.omsi_dependency object (as defined in the omsi.shared.omsi_dependency module)

           :returns: omsi_dependency object with all the dependency data.
        """
        from omsi.shared.omsi_dependency import omsi_dependency
        re = omsi_dependency()
        try:
            re['param_name'] = self.get_parameter_name()
        except:
            re['param_name'] = None
        try:
            re['link_name'] = self.get_link_name()
        except:
            re['link_name'] = None
        try:
            re['selection'] = self.get_selection_string()
        except:
            re['selection'] = None
        try:
            re['omsi_object'] = self.get_dependency_omsiobject()
        except:
            re['omsi_object'] = None
        try:
            re['dataname'] = self.get_dataset_name()
        except:
            re['dataname'] = None

        return re


class omsi_file_msidata(object):

    """Interface for interacting with mass spectrometry imaging datasets stored in omis HDF5 files.
       The interface allows users to interact with the data as if it where a 3D cube even if data
       is missing. Full spectra may be missing in cases where only a region of interest in space
       has been imaged. Spectra may further be pre-processed so that each spectrum has only information
       about its peaks so that each spectrum has it's own mz-axis.

       To load data ue standard array syntax, e.g., [1,1,:] can be used to retrieve the spectrums at
       location (1,1).

       Current limitations:

            * The estimates in def __best_dataset__(self,keys) are fairly crude at this point
            * The __getitem__ function for the partial_spectra case is not implemented yet.
            * The __setitem__ function for the partial spectra case is not implemented yet (Note, it \
              should also support dynamic expansion of the cube by adding previously missing spectra).
            * For the partial cube case, assignement using __setitem__ function is only supported to \
              valid spectra, i.e., spectra that were specified as occupied during the intital creation process.

       Public object variables:

            :ivar shape: Define the full 3D shape of the dataset (i.e., even if the data is stored in sparse manner)
            :ivar dtype: The numpy datatyp of the main MSI data. This is the same as dataset.dtype
            :ivar name: The name of the corresponding groupt in the HDF5 file. Used to generate hard-links to the group.
            :ivar format_type: Define according to which standard the data is stored in the file
            :ivar datasets: List of h5py objects containing possibly multiple different version of the same MSI data (spectra). There may be multiple versions stored with different layouts in order to optimize the selection process.
            :ivar mz: dataset with the global mz axis information. If prelaod_mz is set in the constructor, then this is a numpy dataset with the preloaded data. Otherwise, this is the h5py dataset pointing to the data on disk.
            :ivar xy_index: None if format_type is 'full_cube'. Otherwise, this is the 2D array indicating for each x/y location the index of the spectrum in dataset. If prelaod_xy_index is set in the constructor, then this is a numpy dataset with the preloaded data. Otherwise, this is the h5py dataset pointing to the data on disk. Negative (-1) entries indicate that no spectrum has been recored for the given pixel.
            :ivar inv_xy_index:  2D dataset with n rows and 2 columns indicating for each spectrum i the (x,y) pixel index the spectrum belongs to. \
                                 This index is stored for convenience purposesbut is not actually needed for data access.
            :ivar mz_index: None if format_type is not 'partial_spectra'. Otherwise this is a dataset of the same size as the spectra data stored in dataset. Each entry indicates an index into the mz dataset to determine the mz_data value for a spectrum. This means mz[ mx_index ] gives the true mz value.
            :ivar xy_index_end: None if format_type is not 'partial_spectra'. Otherwise this is a 2D array indicating for each x/y location the index where the given spectrum ends in the dataset. If prelaod_xy_index is set in the constructor, then this is a numpy dataset with the preloaded data. Otherwise, this is the h5py dataset pointing to the data on disk. Negative (-1) entries indicate that no spectrum has been recored for the given pixel.

       Private object variables:

            :ivar _data_group: Store the pointer to the HDF5 group with all the data
            :ivar _fill_xy: Define whether the data should be reconstructed as a full image cube Set using the set_fill_space function(..)
            :ivar _fill_mz: Define whether spectra should be remapped onto a global m/z axis. Set using the set_fill_spectra function(..)

    """

    def __init__(self, data_group, fill_space=True, fill_spectra=True, fill_value=0, preload_mz=False, preload_xy_index=False):
        """ Initialize the omsi_msidata object.

            The fill options are provided to enable a more convenient access to the data independent of how the
            data is stored in the file. If the fill options are enabled, then the user can interact with the
            data as if it where a 3D cube while missing is data is filled in by the given fill value.

            The prelaod options provided here refer to generally smaller parts of the data for which it may be
            more efficient to load the data and keep it around rather than doing repeated reads. If the object is
            used only for a single read and destroyed afterwards, then disabling the preload options may give a
            slight advantage but in most cases enabling the preload should be Ok (default).

            :param data_group: The h5py object for the group with the omsi_msidata.
            :param fill_space: Define whether the data should be padded in space (filled with 0's) when accessing the data using [..]  \
                               operator so that the data behaves like a 3D cube.
            :param fill_spectra: Define whether the spectra should completed by adding 0's so that all spectra retrived via the [..] \
                               opeator so that always spectra of the full length are returned. This option is provided to ease extension
                               of the class to cases where only partial spectra are stored in the file but is not used at this point.
            :param preload_mz: Should the data for the mz axis be kept in memory or loaded on the fly when needed.
            :param preload_xy_index: Should the xy index (if available) be preloaderd into memory or should the required data be loaded on the fly when needed.
            :param fill_value: The integer value that should be used to fill in missing data values.
        """
        # Define the 3D shape of the dataset (i.e., even if the data is stored
        # otherwise)
        self.shape = []
        self.dtype = None  # Which datatyp does the main MSI data have
        # The name of the data group. This is replicated here to ease
        # interaction with the interface.
        self.name = None
        # Store the pointer to the HDF5 group with all the data
        self._data_group = None
        # Define according to which standard the data is stored in the file
        self.format_type = None
        # Store list of pointers to all copies of the HDF5 dataset
        self.datasets = []
        self.mz = None
        self.xy_index = None
        self.inv_xy_index = None
        self.xy_index_end = None
        self.mz_index = None
        # Define whether the data should be reconstructed as a full image cube
        self._fill_xy = fill_space
        # Define whether spectra should be remapped onto a global m/z axis
        self._fill_mz = fill_spectra
        self.is_valid = False
        if data_group is not None:
            self._data_group = data_group
            self.name = self._data_group.name
            self.format_type = omsi_format_msidata.format_types[
                str(self._data_group .get(unicode(omsi_format_msidata.format_name))[0])]
            try:
                self.dependencies = omsi_file_dependencies(
                    self._data_group[unicode(omsi_format_dependencies.dependencies_groupname)])
            except:
                self.dependencies = None
        else:
            raise ValueError(
                "Data group must be h5py.Group object. None give.")
        if self.format_type is not None:
            # Initalize the dataset
            self.datasets = [self._data_group[unicode(x[0])] for x in self._data_group.items(
            ) if x[0].startswith(omsi_format_msidata.dataset_name)]
            self.dtype = self.datasets[0].dtype
            # Initalize the mz data
            if preload_mz:
                self.mz = self._data_group[
                    unicode(omsi_format_msidata.mzdata_name)][:]
            else:
                self.mz = self._data_group[
                    unicode(omsi_format_msidata.mzdata_name)]
            # Initalize the data shape
            if self.format_type == omsi_format_msidata.format_types['full_cube']:
                self.shape = self.datasets[0].shape
                self.is_valid = True
            else:
                self.shape = self._data_group[
                    unicode(omsi_format_msidata_partial_cube.shape_name)][:]

            # We are done initalizaing all variables for the 'full_cube' format
            # Initlize the xy index, inv_xy_index and mz-index for the other
            # formats
            if self.format_type != omsi_format_msidata.format_types['full_cube']:
                if preload_xy_index:
                    self.xy_index = self._data_group[
                        unicode(omsi_format_msidata_partial_cube.xy_index_name)][:]
                    self.inv_xy_index = self._data_group[
                        unicode(omsi_format_msidata_partial_cube.inv_xy_index_name)][:]
                    if self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                        self.xy_index_end = self._data_group[
                            unicode(omsi_format_msidata_partial_spectra.xy_index_end_name)][:]
                        self.mz_index = self._data_group[
                            unicode(omsi_format_msidata_partial_spectra.mz_index_name)]
                        # We have successfully loaded all data for the
                        # partial_spectra case
                        self.is_valid = True
                    elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
                        # We have successfully loaded all data for the
                        # partial_cube case
                        self.is_valid = True
                else:
                    self.xy_index = self._data_group[
                        unicode(omsi_format_msidata_partial_cube.xy_index_name)]
                    self.inv_xy_index = self._data_group[
                        unicode(omsi_format_msidata_partial_cube.inv_xy_index_name)]
                    if self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                        self.xy_index_end = self._data_group[
                            unicode(omsi_format_msidata_partial_spectra.xy_index_end_name)]
                        self.mz_index = self._data_group[
                            unicode(omsi_format_msidata_partial_spectra.mz_index_name)]
                        # We have successfully loaded all data for the
                        # partial_spectra case
                        self.is_valid = True
                    elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
                        # We have successfully loaded all data for the
                        # partial_cube case
                        self.is_valid = True

    def get_h5py_datasets(self, index=0):
        """Get the h5py dataset object for the given dataset.

           :param index: The index of the dataset.
           :returns: h5py object for the requested dataset.
           :raises: and Index error is generated in case an invalid index is given.
        """
        return self.datasets[index]

    def get_h5py_mzdata(self):
        """Get the h5py object for the mz datasets.

           :returns: h5py object of the requested mz dataset.
        """
        return self._data_group[unicode(omsi_format_msidata.mzdata_name)]

    def get_h5py_datagroup(self):
        """Get the h5py group object associated with this instance of omsi_file_msidata.

           :returns: h5py object of the requested data group (or None if the object is invalid)
        """
        return self._data_group

    def items(self):
        """Get the list of items associdated with the h5py.Group object managed by this object"""
        return self._data_group.items()

    def __eq__(self, value):
        """Check whether the two objects have the same h5py name"""
        try:
            equal = (value.name == self.name)
            return equal
        except:
            return False

    def __ne__(self, value):
        """Check whether the two objects have different h5py names"""
        return not self.__eq__(value)

    def get_version(self):
        """Get the omsi version for the representation of this object in the HDF5 file"""
        try:
            return self._data_group.attrs[omsi_format_common.version_attribute]
        except:
            return None

    def get_timestamp(self):
        """Get the timestamp when the dataset group was created in the HDF5 file.

            :returns: Python timestamp string generated using time.ctime().
                      None may be returned in case that the timestamp does not exists
                      or cannot be retrieved from the file for some reason.

            """
        try:
            return self._data_group.attrs[omsi_format_common.timestamp_attribute]
        except:
            return None

    def add_dependency(self, dependency, flush_io=True):
        """Create a new dependency for this dataset

           :param dependency" omsi.shared.omsi_dependency object describing the data dependency
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

           :returns: omsi_file_dependencydata object with the dependency data or None in case that an error occured and the dependency has not been generated.
        """
        if not self.dependencies:
            self.dependencies = self.__create_dependencies__(
                self._data_group, dependencies_data_list=[])

        dep = self.dependencies.add_dependency(dependency)
        if flush_io:
            self._data_group.file.flush()
        return dep

    def create_method_info(self, method_name=None, flush_io=True):
        """Add information about the method imaged to the experiment

           :param method_name: Optional name of the method
           :type method_name: str, None
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 buffers are flushed so that all data has been written to file

           :returns: h5py object of the newly created method group.
        """
        method_group = self._data_group.require_group(omsi_format_methods.methods_groupname)
        method_group.attrs[omsi_format_common.type_attribute] = "omsi_file_methods"
        method_group.attrs[omsi_format_common.version_attribute] = omsi_format_methods.current_version
        method_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            self._data_group.file.flush()
        return omsi_file_methods.__create_method_info__(method_group=method_group, method_name=method_name)

    def create_instrument_info(self, instrument_name=None, mzdata=None, flush_io=True):
        """Add information about the instrument used for creating the images for this experiment.

           :param instrument_name: The name of the instrument
           :type instrument_name: string, None
           :param mzdata: Numpy array of the mz data values of the instrument
           :type mzdata: numpy array or None
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has been written to file

           :returns: The function returns the h5py HDF5 handler to the instrument info group created for the experiment.

        """
        # Create the group for instrument specific data
        instrument_group = self._data_group.require_group(omsi_format_instrument.instrument_groupname)
        instrument_group.attrs[omsi_format_common.type_attribute] = "omsi_file_instrument"
        instrument_group.attrs[omsi_format_common.version_attribute] = omsi_format_instrument.current_version
        instrument_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        if flush_io:
            self._data_group.file.flush()
        return omsi_file_instrument.__create_instrument_info__(instrument_group=instrument_group, instrument_name=instrument_name, mzdata=mzdata)

    def get_method_info(self, check_parent=True):
        """
        Get the omsi_file_methods object with the method information.

        :param check_parent: If no method group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the method. (default=True)

        :returns: omsi_file_methods object for the requested method info. The function returns
                  None in case no method information was found for the experiment

        """
        try:
            return omsi_file_methods(self._data_group[unicode(omsi_format_methods.methods_groupname)])
        except:
            try:
                return omsi_file_methods(self._data_group[unicode(omsi_format_methods.methods_old_groupname)])
            except:
                if check_parent:
                    try:
                        return omsi_file.get_omsi_object(self._data_group.parent).get_method_info()
                    except:
                        pass
        return None

    def get_instrument_info(self, check_parent=True):
        """
        Get the HDF5 group opbject with the instrument information.

        :param check_parent: If no method group is available for this dataset should we check
                             whether the parent object (i.e., the experiment group containing the dataset)
                             has information about the method. (default=True)

        :returns:  omsi_file_instrument object for the requested instrument info. The function returns \
                   None in case no instrument information was found for the experiment
        """
        try:
            return omsi_file_instrument(self._data_group[unicode(omsi_format_instrument.instrument_groupname)])
        except:
            # Check whether the parent group has information about the instrument
            if check_parent:
                try:
                    return omsi_file.get_omsi_object(self._data_group.parent).get_instrument_info()
                except:
                    pass
        return None

    def get_all_dependency_data(self, omsi_dependency_format=True):
        """Get all direct dependencies associdated with the analysis.

           :param omsi_dependency_format: Should the dependcies be retrieved as omsi_analysis_dependency object (True) or as an omsi_file_dependencydata object (False).

           :returns: List omsi_dependency objects containing either omsi file API objects or h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def get_all_dependency_data_recursive(self, omsi_dependency_format=True):
        """Get all direct and indirect dependencies associdated with the analysis.

           NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

           :param omsi_dependency_format: Should the dependencies be retrieved as omsi_dependency object (True) or as an omsi_file_dependencydata object (False)

           :returns: List omsi_analysis_data objects containing either omsi file API interface objects or h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data_recursive(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def has_method_info(self, check_parent=False):
        """Check whether custom method information is available for this dataset.

           :param check_parent: If no method group is available for this dataset should we check \
                                whether the parent object (i.e., the experiment group containing the dataset)\
                                has information about the method. (default=False)
           :returns: Boolean indicating whether method info is available.
        """
        return (self.get_method_info(check_parent=check_parent) is not None)

    def has_instrument_info(self, check_parent=False):
        """Check whether custom instrument information is available for this dataset.

           :param check_parent: If no instrument group is available for this dataset should we check \
                                whether the parent object (i.e., the experiment group containing the dataset) \
                                has information about the instrument. (default=False)

            :returns: Boolean indicating whether instrument info is available.
        """
        return (self.get_instrument_info(check_parent=check_parent) is not None)

    def has_dependencies(self):
        """Check whether any dependencies exisits for this datasets"""
        if self.dependencies:
            return True
        else:
            return False

    def __setitem__(self, key, value):
        """The __getitem__ function is used in python to implement the [..] operator for setting data values.
           This function allows the user to write the data as if it were a 3D cube using array notation [:,:,:]
           even if the data may be stored in a different fashion in the file. If less than three
           selections are specified then the remaining dimensions are assumed to be ":", i.e., all.

           :param key: Three elements of type index, list or slice defining the data selection
        """
        # The object is not fully initalized
        if self._data_group is None:
            raise ValueError("The msidata object has not been initalized.")

        # Complete the input selection if it is only partially specified. In this way we can
        # assume in the following code that we always have three key selection
        # parameters
        if not isinstance(key, tuple):
            key = (key, slice(None), slice(None))
        elif len(key) == 1:
            key = (key, slice(None), slice(None))
        elif len(key) == 2:
            key = (key[0], key[1], slice(None))
        elif len(key) != 3:
            raise ValueError("Invalid selection")

        # Check the data format and call the approbriate getitem function
        if self.format_type == omsi_format_msidata.format_types['full_cube']:
            return self.__setitem_fullcube__(key, value)
        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
            return self.__setitem_partialcube__(key, value)
        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
            return self.__setitem_partialspectra__(key, value)

    def __setitem_fullcube__(self, key, value):
        """Private helper function used in case that the data is stored as a full cube

           :param key: Three elements of type index, list or slice defining the data selection
        """
        # Update the data in all available version of the dataset
        for d in self.datasets:
            d[key] = value

    def __setitem_partialcube__(self, key, value):
        """Private helper function used in case that the data is stored as a partial cube

           :param key: Three elements of type index, list or slice defining the data selection
        """
        # Check which elements need to be updated
        indexlist = self.xy_index[key[0], key[1]]
        # Check whether only valid spectra are selected
        validselectsize = (indexlist >= 0).sum()
        try:
            totalsize = indexlist.size
            indexlist = indexlist.reshape(indexlist.size)
        except:
            totalsize = 1
        if validselectsize != totalsize:
            raise KeyError(
                "Assignment to uninitalized spectra is currently not supported")

        # Update the data in all available version of the dataset
        for d in self.datasets:
            d[indexlist.tolist(), key[2]] = value.reshape(
                totalsize, self.__num_elements__(key[2]))

    def __setitem_partialspectra__(self, key, value):
        """Private helper function used in case that the data is stored as a partial cube of partial spectra.

           :param key: Three elements of type index, list or slice defining the data selection
        """
        raise NotImplementedError(
            "Assignement of partial spectra datasets is currently not implemented")

    def __getitem__(self, key):
        """The __getitem__ function is used in python to implement the [..] operator. This function
           allows the user to access the data as if it were a 3D cube using array notation [:,:,:]
           even if the data may be stored in a different fashion in the file. If less than three
           selections are specified then the remaining dimensions are assumed to be ":", i.e., all.

           :param key: Three elements of type index, list or slice defining the data selection
        """

        # The object is not fully initalized
        if self._data_group is None:
            raise ValueError("The msidata object has not been initalized.")

        # Complete the input selection if it is only partially specified. In this way we can
        # assume in the following code that we always have three key selection
        # parameters
        if not isinstance(key, tuple):
            key = (key, slice(None), slice(None))
        elif len(key) == 1:
            key = (key, slice(None), slice(None))
        elif len(key) == 2:
            key = (key[0], key[1], slice(None))
        elif len(key) != 3:
            raise ValueError("Invalid selection")

        # Check the data format and call the approbriate getitem function
        if self.format_type == omsi_format_msidata.format_types['full_cube']:
            return self.__getitem_fullcube__(key)
        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:
            return self.__getitem_partialcube__(key)
        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
            return self.__getitem_partialspectra__(key)

    def __getitem_fullcube__(self, key):
        """Private helper function used in case that the data is stored as a full cube

           :param key: Three elements of type index, list or slice defining the data selection
        """
        # Get the dataset that is best suited for the selection
        d = self.__best_dataset__(key)
        # Access data using h5py
        return d[key]

    def __getitem_partialcube__(self, key):
        """Private helper function used in case that the data is stored as a partial cube.

          :param key: Three elements of type index, list or slice defining the data selection
        """
        # Get the dataset that is best suited for the selection
        d = self.__best_dataset__(key)
        # Compute the list elements to be loaded in (x,y)
        indexlist = self.xy_index[key[0], key[1]]
        if indexlist.size > 1:
            cleanindexlist = indexlist[indexlist >= 0]
        else:
            indexlist = np.asarray([indexlist])
            if indexlist[0] < 0:
                cleanindexlist = np.empty(0, dtype=int)
            else:
                cleanindexlist = np.asarray([indexlist], dtype=int)
        cleanindexlist = cleanindexlist.reshape(cleanindexlist.size)
        # Load the data
        if cleanindexlist.size > 0:
            data = d[cleanindexlist, key[2]]
        else:
            data = np.empty(0, dtype=d.dtype)

        # Check if we need to complete the data with 0's. If the selection defines a range in
        # x and y (i.e we select a full rectangular region in space), then reconstruct the
        # full region by filling any missing data with 0's.
        if self._fill_xy:
            # Compute how many elements we need to retrieve in the mz-axis
            mzsize = self.__num_elements__(key[2])
            renshape = indexlist.shape + (mzsize, )
            filldata = np.zeros(renshape, dtype=d.dtype)
            filldata[indexlist >= 0] = data.reshape(
                (cleanindexlist.size, self.__num_elements__(key[2])))
            return filldata
            # if len( indexlist.shape )==2 :
            #xsize  = indexlist.shape[0]
            #ysize  = indexlist.shape[1]
            #filldata = np.zeros( shape=(xsize ,  ysize , mzsize) , dtype=d.dtype )
            # if not isinstance( key[0] , list ) and not isinstance( key[1] , list  ) :
            #    filldata[ self.inv_xy_index[cleanindexlist,0]-self.__offset__(key[0]) , self.inv_xy_index[cleanindexlist,1]-self.__offset__(key[1])  ] = data.reshape((cleanindexlist.size, self.__num_elements__(key[2]) ))
            #    return filldata
            # if isinstance( key[0] , list ) :
        else:
            # Otherwise return a list of numpy arrays with the requested data
            # completed with None objects for missing data
            filldata = [None] * indexlist.size
            for i in xrange(0, indexlist.size):
                if indexlist[i] >= 0:
                    filldata[i] = data[cleanindexlist[i], :]
            return filldata

    def __getitem_partialspectra__(self, key):
        """Private helper function used in case that the data is stored as full or partial cube with partial spectra.

           #:param key: Three elements of type index, list or slice defining the data selection
        """
        raise NotImplementedError(
            "Selection for partial spectra datasets is currently not implemented. See __getitem_partialspectra__(...)")

        # Get the dataset that is best suited for the selection
        #d = self.__best_dataset__(key)
        # Determine the data values to be loaded in xy
        #indexList_start = self.xy_index[ key[0] , key[1] ]
        #indexList_end   = self.xy_index_end[ key[0] , key[1] ]
        #index_select = (indexList_start>=0)
        # if not isinstance( indexList_start , int ) :
            #cleanIndexList_start = indexList_start[ index_select ]
            #cleanIndexList_end = indexList_end[ index_select ]
        # else :
            # if indexList_start < 0 :
                #cleanIndexList_start = np.empty(0, dtype=int)
                #cleanIndexList_end   = np.empty(0, dtype=int)
            # else :
                #cleanIndexList_start = np.asarray( [ indexList_start ] , dtype=int )
                #cleanIndexList_end   = np.asarray( [ indexList_end ]   , dtype=int )

        # Load the full specta for each x/y location
        #mz_data = [[]] * cleanIndexList.shape[0]
        # for i in xrange(0,cleanIndexList.shape[0]) :
            #mz_data[i] = self.mz_index[ cleanIndexList_start[i]:cleanIndexList_end[i] ]

        # Determine which parts of the m/z data we actually need
        # if isinstance(key[2] , int ) :
            #mz_select = [ i==key[2] for i in mz_data  ]
        # if isinstance(key[2] , slice ) :
            # if key[2].step>1:
                #raise NotImplementedError("Selections with a stepping are currently not supported by the interface for the m/z dimension and particle_spectra data")
            # Select full specta
            # if (key[2].start is None or key[2].start == 0) and \
                #(key[2].end is None or key[2].end == mz_indexList.shape[1]) and \
                #(key[2].step is None or key[2].step == 1 ) :
                    #mz_select = [ np.ones(dtype='bool', shape=i.shape ) for i in mz_data ]
            # Select an mz-range of the spectra
            # else :
                #mz_select = [ np.logical_and( i>=key[2].start , i<key[2].end ) for i in mz_data  ]
                ##mz_select = np.where( np.logical_and( mz_indexList>=key[2].start , mz_indexList<key[2].end )  , mz_indexList , -1 )
        # elif isinstance(key[2] , list ) :
            # Treat the list the same way as an index selection if we only have one entry in the list
            # if len( key[2] ) == 1 :
                #mz_select = [ i==key[2][0] for i in mz_data  ]
            # Check if the list defines a continues selection
            # for i in xrange(1,len(key[2])) :
                # if key[2][i] != (key[2][i-1]+1) :
                    #raise NotImplementedError("Selections using discontinoues lists are currently not supported by the interface for the m/z dimension and particle_spectra data.")
                    # Implementing discontinoues lists for selection will require a different treatment also when filling the data
            # Treat the list as a continues slice-based selection.
            #mz_select = [ np.logical_and( i>=key[2][0] , i<key[2][-1] ) for i in mz_data  ]

        # Load the requested spectra
        # if self._fill_mz and self._fill_xy and (len( indexList_start.shape )==2)  :
            # Remap both the spatial and m/z coordinates so that the data behaves completely like a 3D cube
            # Compute how many elements we need to retrieve in the mz-axis
            #mzsize = self.__num_elements__(key[2])
            #ysize  = indexList_start.shape[0]
            #xsize  = indexList_start.shape[1]
            #filldata = np.zeros( shape=(xsize, ysize, mzsize) , dtype = d.dtype )
            #mzi = 0
            # for xi in xrange(0,xsize) :
                # for yi in xrange(0,ysize) :
                    # if indexList_start[xi,yi] >= 0 :
                        #filldata[xi,yi,mz_select[mzi]] = d[ indexList_start[i]:indexList_end[i] ][ mz_select[mzi] ]
                        #mzi = mzi+1
        # elif self._fill_mz: #This includes the case where both fill_mz and fill_xy are set but fill_xy is not needed.
            #mzsize = self.__num_elements__(key[2])
            #filldata = [np.zeros( mzsize) ] * len(mz_select)
            # for i in xrange(0, len(mz_select) ):
                #filldata[i][mz_select[i]] =  d[ cleanIndexList_start[i]:cleanIndexList_end[i] ][ mz_select[i] ]
            # return filldata
        # elif self._fill_xy: #Only fill_xy is set but not the fill_mz
            # Remap only the spatial coordinates
            #raise NotImplementedError("fill_xy without fill_mz set is currently not supported by the API")
        # else : #This is the case where neither fill_mz nor fill_xy are set
            # Just load the data and return both the data and mz indicies
            #filldata     =  [[]] * cleanIndexList.shape[0]
            #fillmzdata   =  [[]] * cleanIndexList.shape[0]
            # for i in xrange(0,cleanIndexList.shape[0]) :
                #filldata[i]   = d[ cleanIndexList_start[i]:cleanIndexList_end[i] ][mz_select[i]]
                #fillmzdata[i] = mz_data[i][mz_select[i]]
            # return filldata, fillmzdata

    def set_fill_space(self, fill_space):
        """Define whether spatial selection should be filled with 0's to retrieve full image slices"""
        self._fill_xy = fill_space

    def set_fill_spectra(self, fill_spectra):
        """Define whether spectra should be filled with 0's to map them to the global mz axis when retrieved"""
        self._fill_mz = fill_spectra

    @classmethod
    def __create_dependencies__(cls, data_group, dependencies_data_list=None):
        """Create a group for storing data dependencies

           :param data_group: The h5py.Group for which the dataset should be initalized.
           :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies. Default is empty list []

           :returns: omsi_file_dependencies object created by the function.
        """
        dependencies_group = data_group.require_group(
            omsi_format_dependencies.dependencies_groupname)
        dependencies_group.attrs[omsi_format_common.type_attribute] = "omsi_file_dependencies"
        dependencies_group.attrs[omsi_format_common.version_attribute] = omsi_format_dependencies.current_version
        dependencies_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        omsidep = omsi_file_dependencies.__create_dependencies__(
            dependencies_group=dependencies_group,
            dependencies_data_list=dependencies_data_list)
        return omsidep

    @classmethod
    def __create_msidata_full_cube__(cls, data_group, data_shape, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, dependencies_data_list=None):
        """ Create a new msi data group with all necessay datasets for storing a full 3D cube.

            NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
            corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

            Required input parameters

            :param data_group: The h5py.Group for which the dataset should be initalized.
            :param data_shape: The 3D shape of the MSI dataset in x,y, and m/z
            :param data_type: The dtype for the MSI dataset
            :param mzdata_type: The dtype for the mz dataset

            Optional layout optimization parameters (these refer to the main MSI dataset only)

            :param chunks:  Specify whether chunking should be used (True,False), or specify the chunk sizes to be used explicitly.
            :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'. \
                              Can also use an integer in range(10) indicating gzip.
            :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc.. \
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).

            Other optional parameters

            :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies. Default is empty list []

            :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to be filled with data:

                * data_dataset : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * mz_dataset : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

        """
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(
            omsi_format_msidata.format_name), shape=(1,), dtype=omsi_format_common.str_type)
        format_dataset[0] = 'full_cube'
        # mz data
        mz_dataset = data_group.require_dataset(
            name=omsi_format_msidata.mzdata_name, shape=(data_shape[2], ), dtype=mzdata_type)
        # Main MSI dataset
        if compression is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=data_shape, dtype=data_type, chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=data_shape, dtype=data_type, chunks=chunks, compression=compression)
        else:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=data_shape,
                dtype=data_type, chunks=chunks, compression=compression, compression_opts=compression_opts)
        # Create any dependencies
        omsi_file_msidata.__create_dependencies__(
            data_group, dependencies_data_list)
        # Return the datasets that need to be written
        return data_dataset, mz_dataset

    @classmethod
    def __create_msidata_partial_cube__(cls, data_group, data_shape, mask, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, dependencies_data_list=None):
        """ Create a new msi data group with all necessay datasets for storing a full 3D cube

            NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
            corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

            Required input parameters

            :param data_group: The hHDF5 group for which the dataset should be initalized.
            :param data_shape: The 3D shape of the MSI dataset in x,y, and m/z
            :param mask: 2D boolean NumPy array used as mask to indicate which (x,y) locations have spectra associated with them.
            :param data_type: The dtype for the MSI dataset
            :param mzdata_type: The dtype for the mz dataset

            Optional layout optimization parameters (these refer to the main MSI dataset only)

            :param chunks:  Specify whether chunking should be used (True,False), or specify the chunk sizes to be used explicitly.
            :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'. \
                              Can also use an integer in range(10) indicating gzip.
            :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc.. \
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).

            Other optional parameters

            :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies. Default is empty list []

            :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order to be filled with data:

                * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

            :returns: The following already complete dataset

                * ``xy_index_dataset`` : This dataset indicates for each xy location to which index in ``data_dataset`` the \
                                         location corresponds to. This dataset is needed to identify where spectra need to be written to.
                * ``inv_xy_index_dataset`` : This datasets indicates for each spectrum its xy locating. The dataset is of the from \
                                         n x 2 with inv_xy_index_dataset[:,0] being the x indices and inv_xy_index_dataset[:,1] \
                                         being the y indicies.

        """
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(
            omsi_format_msidata.format_name), shape=(1,), dtype=omsi_format_common.str_type)
        format_dataset[0] = 'partial_cube'
        # mz data
        mz_dataset = data_group.require_dataset(
            name=omsi_format_msidata.mzdata_name, shape=(data_shape[2], ), dtype=mzdata_type)

        # Main MSI dataset
        xyshape = data_shape[0] * data_shape[1]
        numspectra = mask.sum()
        mzshape = data_shape[2]
        if isinstance(chunks, tuple):
            chunks = (chunks[0] * chunks[1], chunks[2])
        if compression is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=(numspectra, mzshape), dtype=data_type, chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=(numspectra, mzshape), dtype=data_type, chunks=chunks, compression=compression)
        else:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"), shape=(numspectra, mzshape),
                dtype=data_type, chunks=chunks, compression=compression, compression_opts=compression_opts)

        # Determine the minimal dtype for the xy index
        xy_index_dtype = 'int16'
        if xyshape < np.iinfo('int16').max:
            xy_index_dtype = 'int16'
        elif xyshape < np.iinfo('int32').max:
            xy_index_dtype = 'int32'
        elif xyshape < np.iinfo('int64').max:
            xy_index_dtype = 'int64'

        # Create the xy_index_dataset
        xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.xy_index_name, shape=(data_shape[0], data_shape[1]), dtype=xy_index_dtype)
        xy_index_dataset[:, :] = -1
        xy_index_dataset[mask] = np.arange(0, numspectra)

        # Create the inverse inv_xy_index dataset
        inv_xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.inv_xy_index_name, shape=(numspectra, 2), dtype=xy_index_dtype)
        for xi in xrange(0, data_shape[0]):
            for yi in xrange(0, data_shape[1]):
                if xy_index_dataset[xi, yi] >= 0:
                    inv_xy_index_dataset[xy_index_dataset[xi, yi], :] = np.asarray([xi, yi])

        # Create the shape dataset
        shape_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.shape_name, shape=(3, ), dtype='uint32')
        shape_dataset[:] = np.asarray(data_shape, dtype='uint32')

        # Create any dependencies
        omsi_file_msidata.__create_dependencies__(
            data_group, dependencies_data_list)

        # Return the datasets that need to be written
        return data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset

    @classmethod
    def __create_msidata_partial_spectra__(cls, data_group, spectra_length, len_global_mz, data_type='f', mzdata_type='f', chunks=None, compression=None, compression_opts=None, dependencies_data_list=None):
        """Create a new msi data group with all necessay datasets for storing a full or partial cube of partial spectra

            NOTE: This is a private helper function used to initalize the content of the dataset group in HDF5. Use the
            corresponding create_msidata functions in omsi_file_experiment to create a new MSI dataset in HDF5.

            Required input parameters

            :param data_group: The hDF5 group for which the dataset should be initalized.
            :param spectra_length: 2D boolean NumPy array used indicating for each (x,y) locations the length of
                                   the corresponding partial spectrum.
            :param len_global_mz: The total number of m/z values in the global m/z axis for the full 3D cube
            :param data_type: The dtype for the MSI dataset
            :param mzdata_type: The dtype for the mz dataset

            Optional layout optimization parameters (these refer to the main MSI dataset only)

            :param chunks:  Specify whether chunking should be used (True,False), or specify the
                            chunk sizes to be used explicitly.
            :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                                Can also use an integer in range(10) indicating gzip.
            :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                                For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                                a number between zero and nine (inclusive) to indicate the tradeoff between speed
                                and compression ratio (zero is fastest, nine is best ratio).

            Other optional parameters

            :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies.
                                           Default is empty list []

            :returns: The following two empty (but approbriatelu sized) h5py datasets are returned in order
                      to be filled with data:

                * ``data_dataset`` : The primary h5py dataset for the MSI data with shape data_shape and dtype data_type.
                * ``mz_index_dataset`` : The
                * ``mz_dataset`` : The h5py dataset for the mz axis data with shape [data_shape[2]] and dtype mzdata_type.

            :returns: The following already complete dataset

                * ``xy_index_dataset`` : This dataset indicates for each xy location at which index in ``data_dataset``
                                        the corresponding spectrum starts. This dataset is needed to identify where
                                        spectra need to be written to.
                * ``xy_index_end_dataset`` : This dataset indicates for each xy location at which index in
                                ``data_dataset`` the corresponding spectrum ends (exclusing the given value).
                                This dataset is needed to identify where spectra need to be written to.

        """
        # Data format
        format_dataset = data_group.require_dataset(name=unicode(
            omsi_format_msidata.format_name), shape=(1,), dtype=omsi_format_common.str_type)
        format_dataset[0] = 'partial_spectra'
        # mz data
        mz_dataset = data_group.require_dataset(
            name=omsi_format_msidata.mzdata_name, shape=(len_global_mz, ), dtype=mzdata_type)

        # Main MSI dataset
        mask = (spectra_length > 0)
        numspectra = mask.sum()
        total_len_spectra = spectra_length[mask].sum()
        if isinstance(chunks, tuple):
            chunks = (chunks[0] * chunks[1] * chunks[2], )
        if compression is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks)
        elif compression_opts is None:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks,
                compression=compression)
        else:
            data_dataset = data_group.require_dataset(
                name=(omsi_format_msidata.dataset_name + "0"),
                shape=(total_len_spectra, ),
                dtype=data_type,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        # Determine the minimal dtype for the mz index
        mz_index_dtype = 'uint16'
        if len_global_mz < np.iinfo('uint16').max:
            mz_index_dtype = 'uint16'
        elif len_global_mz < np.iinfo('uint32').max:
            mz_index_dtype = 'uint32'
        elif len_global_mz < np.iinfo('uint64').max:
            mz_index_dtype = 'uint64'

        # Create the mz-index dataset
        if compression is None:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks)
        elif compression_opts is None:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks,
                compression=compression)
        else:
            mz_index_dataset = data_group.require_dataset(
                name=omsi_format_msidata_partial_spectra.mz_index_name,
                shape=(total_len_spectra, ),
                dtype=mz_index_dtype,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        # Determine the minimal dtype for the xy index
        xy_index_dtype = 'uint16'
        xyshape = spectra_length.shape[0] * spectra_length.shape[1]
        if xyshape < np.iinfo('uint16').max:
            xy_index_dtype = 'uint16'
        elif xyshape < np.iinfo('uint32').max:
            xy_index_dtype = 'uint32'
        elif xyshape < np.iinfo('uint64').max:
            xy_index_dtype = 'uint64'

        # Create the xy index start dataset
        xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.xy_index_name,
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        xy_index_dataset[:] = -1
        xy_index_dataset[mask] = np.insert(
            np.cumsum(spectra_length[mask]), 0, 0)[0:-1]

        # Create the xy index end dataset
        xy_index_end_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_spectra.xy_index_end_name,
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        xy_index_end_dataset[:] = -1
        xy_index_end_dataset[mask] = np.insert(
            np.cumsum(spectra_length[mask]), 0, 0)[1:]

        # Create the inverse inv_xy_index dataset
        inv_xy_index_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.inv_xy_index_name,
            shape=(numspectra, 2),
            dtype=xy_index_dtype)
        temp_count_xy_index = np.zeros(
            shape=spectra_length.shape,
            dtype=xy_index_dtype)
        temp_count_xy_index[:] = -1
        temp_count_xy_index[mask] = np.insert(
            np.cumsum(np.ones(shape=spectra_length.shape, dtype='uint16')[mask]), 0, 0)[0:-1]
        for xi in xrange(0, spectra_length.shape[0]):
            for yi in xrange(0, spectra_length.shape[1]):
                if temp_count_xy_index[xi, yi] >= 0:
                    inv_xy_index_dataset[temp_count_xy_index[xi, yi], :] = np.asarray([xi, yi])

        # Create the shape dataset
        shape_dataset = data_group.require_dataset(
            name=omsi_format_msidata_partial_cube.shape_name,
            shape=(3, ),
            dtype='uint32')
        shape_dataset[:] = np.asarray(mask.shape, dtype='uint32')

        # Create any dependencies
        omsi_file_msidata.__create_dependencies__(
            data_group, dependencies_data_list)

        # Return the datasets that need to be written
        return data_dataset, mz_index_dataset, mz_dataset, xy_index_dataset, xy_index_end_dataset

    def create_optimized_chunking(self, chunks=None, compression=None, compression_opts=None, copy_data=True, print_status=False, flush_io=True):
        """Helper function to allow one to create optimized copies of the dataset with different internal data
           layouts to speed up selections. The function expects that the original data has already been written
           to the data group. The function takes

           :param chunks:  Specify whether chunking should be used (True,False), or specify the chunk sizes
                           to be used explicitly.
           :param compression: h5py compression option. Compression strategy.  Legal values are 'gzip', 'szip',  'lzf'.
                              Can also use an integer in range(10) indicating gzip.
           :param compression_opts: h5py compression settings.  This is an integer for gzip, 2-tuple for szip, etc..
                              For gzip (H5 deflate filter) this is the aggression paramter. The agression parameter is
                              a number between zero and nine (inclusive) to indicate the tradeoff between speed
                              and compression ratio (zero is fastest, nine is best ratio).
           :param copy_data: Should the MSI data be copied by this function to the new dataset or not. If False, then
                             it is up to the user of the function to copy the appropriate data into the returned h5py
                             dataset (not recommended but may be useful for performance optimization).
           :param print_status: Should the function print the status of the conversion process to the command line?
           :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all data has
                           been written to file

           :returns: h5py dataset with the new copy of the data
        """
        if len(self.datasets) == 0:
            raise ValueError(
                "No datasets are currently stored for the dataset group that could be replicated")

        # Get the donor dataset. Try tp get the main dataset first. If it is
        # not found take the first one from the list.
        try:
            d = self._data_group .get(
                unicode(omsi_format_msidata.dataset_name + "0"))
        except:
            d = self.datasets[0]

        # Check and adjust the chunking
        if isinstance(chunks, tuple):
            if self.format_type == omsi_format_msidata.format_types['partial_cube']:
                chunks = (chunks[0] * chunks[1], chunks[2])
                # Correct the chunking if it is larger than the actual data
                if chunks[0] > d.shape[0]:
                    chunks = (d.shape[0], chunks[1])
                if chunks[1] > d.shape[1]:
                    chunks = (chunks[0], d.shape[1])
            elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:
                chunks = (chunks[0] * chunks[1] * chunks[2], )
                # Correct the chunking if the chunk is larger than the actual
                # data
                if chunks[0] > d.shape[0]:
                    chunks = (d.shape[0], )

        # Generate the new dataset in the HDF5 file
        new_data_name = omsi_format_msidata.dataset_name + \
            str(len(self.datasets))
        if compression is None:
            data_dataset = self._data_group.require_dataset(
                name=new_data_name, shape=d.shape, dtype=self.dtype, chunks=chunks)
        elif compression_opts is None:
            data_dataset = self._data_group.require_dataset(
                name=new_data_name, shape=d.shape, dtype=self.dtype, chunks=chunks, compression=compression)
        else:
            data_dataset = self._data_group.require_dataset(
                name=new_data_name,
                shape=d.shape,
                dtype=self.dtype,
                chunks=chunks,
                compression=compression,
                compression_opts=compression_opts)

        if copy_data:
            self.copy_msidataset(
                source=d, destination=data_dataset, print_status=print_status)

        if flush_io:
            self._data_group.file.flush()

        # Add the dataset to the list of datasets
        self.datasets.append(data_dataset)

        # Return the dataset. This is needed in case the caller decided to set
        # copy_data to False and copy the data themselfs
        return data_dataset

    def copy_msidataset(self, source, destination, print_status=False):
        """Helper function used to copy a source msi dataset one chunk at a time to the destination dataset.
           The data copy is done one destination chunk at a time to achieve chunk-aligned write

           :param source: The source h5py dataset
           :param destination: The desitnation dataset.
           :param print_status: Should the function print the status of the conversion process to the command line?

        """
        # Write the data from the donor to the target dataset
        if print_status:
            import sys
        if self.format_type == omsi_format_msidata.format_types['full_cube']:

            chunks = destination.chunks
            numchunksx = int(math.ceil(float(source.shape[0]) / float(chunks[0])))
            numchunksy = int(math.ceil(float(source.shape[1]) / float(chunks[1])))
            numchunksz = int(math.ceil(float(source.shape[2]) / float(chunks[2])))
            numchunks = numchunksx * numchunksy * numchunksz
            itertest = 0
            for xt in xrange(0, numchunksx):
                xstart = xt * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                for yt in xrange(0, numchunksy):
                    ystart = yt * chunks[1]
                    yend = min(ystart + chunks[1], source.shape[1])
                    for zt in xrange(0, numchunksz):
                        zstart = zt * chunks[2]
                        zend = min(zstart + chunks[2], source.shape[2])
                        destination[xstart:xend, ystart:yend, zstart:zend] = source[xstart:xend, ystart:yend, zstart:zend]
                        itertest += 1
                        if print_status:
                            sys.stdout.write(
                                "[" + str(int(100. * float(itertest) / float(numchunks))) + "%]" + "\r")
                            sys.stdout.flush()

        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:

            chunks = destination.chunks
            numchunksx = int(
                math.ceil(float(source.shape[0]) / float(chunks[0])))
            numchunksy = int(
                math.ceil(float(source.shape[1]) / float(chunks[1])))
            numchunks = numchunksx * numchunksy
            itertest = 0
            for xt in xrange(0, numchunksx):
                xstart = xt * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                for yt in xrange(0, numchunksy):
                    ystart = yt * chunks[1]
                    yend = min(ystart + chunks[1], source.shape[1])
                    destination[xstart:xend, ystart:
                                yend] = source[xstart:xend, ystart:yend]
                    itertest += 1
                    if print_status:
                        sys.stdout.write(
                            "[" + str(int(100. * float(itertest) / float(numchunks))) + "%]" + "\r")
                        sys.stdout.flush()

        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:

            chunks = destination.chunks
            numchunksx = int(
                math.ceil(float(source.shape[0]) / float(destination.chunks[0])))
            numchunks = numchunksx
            itertest = 0
            for xt in xrange(0, numchunksx):
                xstart = xt * chunks[0]
                xend = min(xstart + chunks[0], source.shape[0])
                destination[xstart:xend] = source[xstart:xend]
                itertest += 1
                if print_status:
                    sys.stdout.write(
                        "[" + str(int(100. * float(itertest) / float(numchunks))) + "%]" + "\r")
                    sys.stdout.flush()

    def __best_dataset__(self, keys, print_info=False):
        """Compute the index of the dataset that is best suited for executing the given selection

           :param keys: List of three keys indicting which elements should be selected from the MSI dataset.

           :returns: Integer indicating the index of the dataset to be used for the given selection.
        """
        # If there is only one version of the dataset then use that one
        if len(self.datasets) == 1:
            return self.datasets[0]
        # Get the basic properties of the selection
        sizex = self.__num_elements__(keys[0])
        sizey = self.__num_elements__(keys[1])
        sizez = self.__num_elements__(keys[2])
        offsetx = self.__offset__(keys[0])
        offsety = self.__offset__(keys[1])
        offsetz = self.__offset__(keys[2])

        # Try to determine how many chunks need to be touched in order to fullfill the given selection
        # Currently this is just an estimate assuming that we have spatially consecutive selections,
        # i.e., if we have a stepping >1 in a slice or lists of non-neighboring elements, then this
        # estimate will be wrong. This may result in a suggestion of a sub-optimal dataset with
        # respect to performance.
        # For the particle_cube case the estimate is also not accurate for spatial selections since
        # both x and y are linearized in this case. However, the estimate should still distinguish
        # properly between selection in space vs. spectra.
        # For the partial spectra case, currently no chunking optimization is
        # currently available.
        suggestion = self.datasets[0]
        if self.format_type == omsi_format_msidata.format_types['full_cube']:

            currenttouch = -1
            currentload = -1
            for d in self.datasets:
                chunking = d.chunks
                chunkload = chunking[0] * chunking[1] * chunking[2]
                # Account for the missalignment of the selection with chunk
                # boundaries
                ox = offsetx - \
                    (math.floor(float(offsetx) / chunking[0]) * chunking[0])
                oy = offsety - \
                    (math.floor(float(offsety) / chunking[1]) * chunking[1])
                oz = offsetz - \
                    (math.floor(float(offsetz) / chunking[2]) * chunking[2])
                # Determine how many chunks need to be touched, assuming that
                # we have continues selection.
                touchx = math.ceil(float(sizex + ox) / chunking[0])
                touchy = math.ceil(float(sizey + oy) / chunking[1])
                touchz = math.ceil(float(sizez + oz) / chunking[2])
                touchtotal = touchx * touchy * touchz
                # Determine the amount of data that need to be loaded
                loadtotal = touchtotal * chunkload
                # If less chunks have to be touched for this dataset then use
                # that one instead
                if loadtotal < currentload or currentload < 0:
                    suggestion = d
                    currenttouch = touchtotal
                    currentload = loadtotal

        elif self.format_type == omsi_format_msidata.format_types['partial_cube']:

            currenttouch = -1
            currentload = -1
            for d in self.datasets:
                chunking = d.chunks
                chunkload = chunking[0] * chunking[1]
                # Account for the missalignment of the selection with chunk
                # boundaries
                ox = offsetx - \
                    (math.floor(float(offsetx) / chunking[0]) * chunking[0])
                oy = offsety - \
                    (math.floor(float(offsety) / chunking[0]) * chunking[0])
                oz = offsetz - \
                    (math.floor(float(offsetz) / chunking[1]) * chunking[1])
                # Determine how many chunks need to be touched, assuming that
                # we have continues selection.
                touchx = math.ceil(float(sizex + ox) / chunking[0])
                touchy = math.ceil(float(sizey + oy) / chunking[0])
                touchz = math.ceil(float(sizez + oz) / chunking[1])
                touchtotal = touchx * touchy * touchz
                # Determine the amount of data that need to be loaded
                loadtotal = touchtotal * chunkload
                # If less chunks have to be touched for this dataset then use
                # that one instead
                if loadtotal < currentload or currentload < 0:
                    suggestion = d
                    currenttouch = touchtotal
                    currentload = loadtotal

        elif self.format_type == omsi_format_msidata.format_types['partial_spectra']:

            # ToDo currently no optimization is implemented for the partial
            # spectra case
            pass

        if print_info:
            print "Selection: " + str(keys)
            print "Suggest: " + str(suggestion.chunks)
            print "Alternatives: " + str([d.chunks for d in self.datasets])

        return suggestion

    @staticmethod
    def __offset__(key):
        """Determine the start offset of the given single key

           :param key: List, slice or interger indicating a single selection

           :returns: Integer indicating the lower bound of the selection defined by the key.
        """
        if isinstance(key, int):
            return key
        elif isinstance(key, list):
            return min(key)
        elif isinstance(key, slice):
            start = key.start
            if start is None:
                start = 0
            return start
        else:
            return 0

    def __num_elements__(self, key):
        """Compute the number of elements selected by a given single selection key

           :param key: List, slice or interger indicating a single selection

           :returns: The number of elements selected by the given key.
        """
        if isinstance(key, int):
            return 1
        elif isinstance(key, list) or isinstance(key, np.ndarray):
            return len(key)
        elif isinstance(key, slice):
            start = key.start
            if start is None:
                start = 0
            end = key.stop
            if end is None:
                end = self.shape[2]
            step = key.step
            if key.step is None:
                step = 1
            return math.ceil(float(end - start) / float(step))
        else:
            raise ValueError("Unexpected key value " + str(key) + " " +str(type(key)))
