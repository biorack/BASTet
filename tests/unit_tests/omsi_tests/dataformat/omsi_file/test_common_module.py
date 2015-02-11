"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile
from omsi.dataformat.omsi_file.main_file import omsi_file


class test_omsi_file_common(unittest.TestCase):

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

if __name__ == '__main__':
    unittest.main()
