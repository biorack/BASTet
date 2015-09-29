"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
import tempfile
from omsi.analysis.generic import analysis_generic
from omsi.dataformat.omsi_file.main_file import omsi_file
import numpy as np


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

    def test_wrap_function(self):
        def f(a):
            return np.sum(a)
        g = analysis_generic.from_function(f)
        res = g.execute(a=np.arange(10))
        self.assertEquals(res, 45)

    def test_wrap_function_and_save(self):
        def f(a):
            return np.sum(a)
        g = analysis_generic.from_function(f)
        g.execute(a=np.arange(10))
        f = omsi_file(self.test_filename, 'a')
        e = f.create_experiment()
        a = e.create_analysis(g)
        f.flush()
        self.assertNotEquals(a, None)

    def test_warp_function_save_and_recreate(self):
        def f(a):
            return np.sum(a)
        g = analysis_generic.from_function(f)
        res1 = g.execute(a=np.arange(10))
        f = omsi_file(self.test_filename, 'a')
        e = f.create_experiment()
        e.create_analysis(g)
        f.flush()
        del f
        f = omsi_file(self.test_filename, 'a')
        e = f.get_experiment(0)
        a = e.get_analysis(0)
        g2 = a.restore_analysis()
        res2 = g2.execute()
        self.assertEquals(res1, res2)


if __name__ == '__main__':
    unittest.main()





