"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile
import sys
import h5py
import numpy as np
from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.dataformat.omsi_file.msidata import omsi_file_msidata


class test_omsi_file_msidata(unittest.TestCase):

    def setUp(self):
        self.named_temporary_file = tempfile.NamedTemporaryFile()
        self.test_filename = self.named_temporary_file.name
        self.testfile = omsi_file(self.test_filename)
        self.exp = self.testfile.create_experiment()

    def tearDown(self):
        # Clean up the test suite
        del self.testfile
        del self.test_filename
        del self.named_temporary_file

    def test_create_msidata_partial_cube(self):
        # Test the creation and interaction with a partial cube dataset
        tempshape = tuple([10, 10, 1000])
        tempnumelements = tempshape[0] * tempshape[1] * tempshape[2]
        specchunking = tuple([1, 1, 1000])
        slicchunking = tuple([10, 10, 1])
        mask = np.zeros((tempshape[0], tempshape[1]), dtype='bool')
        mask[1:8, 5:10] = True
        data_dataset, mz_dataset, xy_index_dataset, inv_xy_index_dataset, datagroup = \
            self.exp.create_msidata_partial_cube(data_shape=tempshape,
                                            mask=mask,
                                            chunks=specchunking)
        self.assertIsInstance(data_dataset, h5py.Dataset, msg='Check type of partial msi dataset')
        self.assertIsInstance(mz_dataset, h5py.Dataset, msg='Check type of partial mz dataset')
        self.assertIsInstance(xy_index_dataset, h5py.Dataset, msg='Check type of partial xy_index')
        self.assertIsInstance(inv_xy_index_dataset, h5py.Dataset, msg='Check type of partial mz dataset')
        self.assertIsInstance(datagroup, h5py.Group, msg='Check type of data group')

        # Create the omsi_file_msidata object
        test_omsi_file_msidata_object = omsi_file_msidata(datagroup)
        self.assertIsInstance(test_omsi_file_msidata_object, omsi_file_msidata)

        # Writing the data into the object
        temp_data_msidata = np.arange(7 * 5 * 1000).reshape((7, 5, 1000))
        test_omsi_file_msidata_object[1:8, 5:10, :] = temp_data_msidata
        self.assertTrue(np.all(test_omsi_file_msidata_object[1:8, 5:10, :] == temp_data_msidata))

        # Write the mz data
        temp_data_mz = np.arange(tempshape[2])
        mz_dataset[:] = temp_data_mz
        self.assertTrue(np.all(temp_data_mz == mz_dataset[:]),
                        msg="Ensure that all mz data has been written properly")

        # Create the optimized chunking
        test_omsi_file_msidata_object.create_optimized_chunking(chunks=slicchunking,
                                                                compression=None,
                                                                compression_opts=None,
                                                                copy_data=True,
                                                                print_status=True)
        self.assertEquals(len(test_omsi_file_msidata_object.datasets), 2,
                          msg='Check if we have two datasets after replication')

        # Test retrieving the msi dataset
        self.assertIsInstance(self.exp.get_msidata(0), omsi_file_msidata, msg="Check retrieval of msi data")

        # Test slicing against the data
        self.assertTrue(np.all(test_omsi_file_msidata_object[1, 5, :] == temp_data_msidata[0, 0, :] ),
                        msg='Slicing with [1, 5, :] failed.')
        # TODO Add more tests for other slicing options
        # print dataset[:, :, 1]
        # print dataset[1:2, 1:2, 5:50]
        # print dataset[[1, 3], 8:9, :]
        # print dataset[:, :, :]

        # Test that the number of msi datasets is 1
        self.assertEquals(self.exp.get_num_msidata(), 1)


if __name__ == '__main__':
    unittest.main()
