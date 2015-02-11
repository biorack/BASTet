"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile
import numpy as np
from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.dataformat.omsi_file.instrument import omsi_file_instrument


class test_omsi_file_instrument(unittest.TestCase):

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

    def test_create_instrument_info(self):
        # Testing the creation of instrument information
        instrument_info = self.exp.create_instrument_info('undefined', np.arange(10))
        self.assertIsInstance(instrument_info, omsi_file_instrument)
        instrument_info = self.exp.get_instrument_info()
        self.assertIsInstance(instrument_info, omsi_file_instrument)
        instrument_info.set_instrument_name('Test')
        self.assertEquals(instrument_info.get_instrument_name()[:], 'Test')
        self.assertTrue(np.all(instrument_info.get_instrument_mz()[:] == np.arange(10)))



if __name__ == '__main__':
    unittest.main()
