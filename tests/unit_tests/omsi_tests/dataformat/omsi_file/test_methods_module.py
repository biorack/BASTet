"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile
from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.dataformat.omsi_file.methods import omsi_file_methods


class test_omsi_file_methods(unittest.TestCase):

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

    def test_create_method_info(self):
        # Testing the creation of instrument information
        method_info = self.exp.create_method_info('undefined')
        self.assertIsInstance(method_info, omsi_file_methods)
        method_info = self.exp.get_method_info()
        self.assertIsInstance(method_info, omsi_file_methods)
        method_info.set_method_name('Test')
        self.assertEquals(method_info.get_method_name()[:], 'Test')


if __name__ == '__main__':
    unittest.main()
