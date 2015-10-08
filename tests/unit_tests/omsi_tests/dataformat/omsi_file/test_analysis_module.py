"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile

from omsi.dataformat.omsi_file.main_file import omsi_file


class test_omsi_file_analysis(unittest.TestCase):

    def setUp(self):
        self.named_temporary_file = tempfile.NamedTemporaryFile()
        self.test_filename = self.named_temporary_file.name
        self.testfile = omsi_file(self.test_filename)

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
        pass
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
