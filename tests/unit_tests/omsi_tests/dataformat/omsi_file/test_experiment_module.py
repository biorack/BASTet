"""
Test the experiment module
"""
import unittest
import tempfile
import sys
from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.dataformat.omsi_file.experiment import omsi_file_experiment


class test_omsi_file_experiment(unittest.TestCase):

    def setUp(self):
        self.named_temporary_file = tempfile.NamedTemporaryFile()
        self.test_filename = self.named_temporary_file.name
        self.testfile = omsi_file(self.test_filename)

    def tearDown(self):
        # Clean up the test suite
        del self.testfile
        del self.test_filename
        del self.named_temporary_file

    def test_create_and_get_experiment(self):
        # Test the creation of a new experiment
        exp = self.testfile.create_experiment()
        self.assertIsNotNone(exp,
                             msg='Creation of experiment returned None')
        self.assertIsInstance(exp, omsi_file_experiment,
                              msg='Creation of experiment failed')
        self.assertIsInstance(self.testfile.get_experiment(0), omsi_file_experiment,
                              msg='Getting the experiment failed')

        # Setting and retrieving the the experiment identifier
        exp_identifier_text = "Experiment 1"
        try:
            exp.set_experiment_identifier(exp_identifier_text)
        except:
            self.fail("Setting the experiment identifier failed. " + unicode(sys.exc_info()))
        self.assertEquals(exp.get_experiment_identifier()[:], exp_identifier_text)

        # Checking the number of experiments
        self.assertEquals(self.testfile.get_num_experiments(), 1)

if __name__ == '__main__':
    unittest.main()
