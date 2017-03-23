"""
Test the omsi.shared.omsi_dependency module
"""
import unittest
import numpy as np

import omsi.shared.data_selection as data_selection


# TODO Add tests for most other funcstions

class test_omsi_data_selection(unittest.TestCase):

    def setUp(self):
        # Create a base object of dependency_dict that we keep modifying
        # from omsi.dataformat.omsi_file.main_file import omsi_file
        # import tempfile
        pass

    def tearDown(self):
        # Clean up the test suite
        pass

    def test_perform_reduction(self):
        # Test that all reduciton operations run without error
        test_data = np.arange(1000).reshape(10,10,10)
        secondary_data = {}
        min_dim = None
        http_error = False
        failed_reductions = []
        custom_args = {'percentile': {'q': 0.95},
                       'append': {'values': test_data}}
        for reduction in data_selection.reduction_allowed_functions:
            try:
                kwargs = custom_args[reduction] if reduction in custom_args else {}
                tdata = test_data
                if reduction is 'diag':
                    tdata =  np.arange(100).reshape(10,10)
                elif reduction is 'bincount':
                    tdata = test_data.flatten()
                elif reduction is 'diag_indices':
                    tdata = 4
                data_selection.perform_reduction(data=tdata,
                                                 reduction=reduction,
                                                 secondary_data=secondary_data,
                                                 min_dim=min_dim,
                                                 http_error=http_error,
                                                 **kwargs)
            except:
                failed_reductions.append(reduction)
        if len(failed_reductions) > 0:
            self.fail("Data reduction failed: " + str(failed_reductions))

    def test_check_selection_string(self):
        # test check selection string for multiple inputs
        data = {'indexlist': {'data': '[1,2,3]', 'return': 'indexlist'},
                'all': {'data': ':', 'return': 'all'},
                'range': {'data': '3:10', 'return': 'range'},
                'index': {'data': '1', 'return': 'index'},
                'invalid' : {'data': '?', 'return': 'invalid'},
                'multiaxis' : {'data': '1:2|3:4', 'return': 'multiaxis'}}
        failed_test = []
        for key, case in data.iteritems():
            re = data_selection.check_selection_string(case['data'])
            if re !=  data_selection.selection_type[case['return']]:
                failed_test.append([re, case])
        self.assertTrue(len(failed_test)==0, "Failed check selection string for: " + str(failed_test))



if __name__ == '__main__':
    unittest.main()
