"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile

from omsi.dataformat.omsi_file.main_file import omsi_file
from omsi.analysis.generic import analysis_generic
from omsi.datastructures.analysis_data import analysis_data
from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
import numpy as np


class test_omsi_file_analysis(unittest.TestCase):

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

    def test_create_omsi_file(self):
        # Check whether the omsi_file has been created correctly
        self.assertIsInstance(self.testfile, omsi_file)

    def test_create_analysis(self):
        # Testing the creation of derived analysis
        testanaidname = "Peak Finding 123"
        def my_funct(a):
            return np.sum(a)
        testana = analysis_generic.from_function(my_funct, name_key=testanaidname)
        _ = testana.execute(a=np.arange(10))
        analysis, _ = self.exp.create_analysis(testana)
        self.assertIsInstance(analysis, omsi_file_analysis)
        tempanaid = analysis.get_analysis_identifier()[:]
        self.assertEquals(tempanaid, testanaidname)
        tempanatype = analysis.get_analysis_type()[:]
        self.assertEquals(tempanatype, testana.get_analysis_type())
        self.assertIsInstance(self.exp.get_analysis_by_identifier(testanaidname), omsi_file_analysis)


    """
    print "Creating derived analysis"
    testanaidname = "Peak Finding 123"
    testana = analysis_generic(name_key=testanaidname)
    print "Creating dummy analysis data"
    testanadata = analysis_data()
    testanadata['name'] = 'peakcube'
    testanadata['data'] = np.zeros(shape=(5, 5, 5), dtype='float32')
    testanadata['dtype'] = 'float32'
    # print "Adding the dummy analysis data to the analys object"
    # testana.add_analysis_data( testanadata )
    print "Saving the analysis data to file"
    analysis, _ = exp.create_analysis(testana)
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

"""

if __name__ == '__main__':
    unittest.main()
