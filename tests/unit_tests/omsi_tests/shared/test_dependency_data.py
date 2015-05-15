"""
Test the omsi.shared.omsi_dependency module
"""
import unittest
from omsi.shared.dependency_data import dependency_dict

# TODO create a test openmsi data file during setup
# TODO Add tests for get_data()
# TODO Add tests for __getitem__ with slicing into getdata()

class test_omsi_dependency(unittest.TestCase):

    def setUp(self):
        # Create a base object of dependency_dict that we keep modifying
        # from omsi.dataformat.omsi_file.main_file import omsi_file
        # import tempfile
        pass

    def tearDown(self):
        # Clean up the test suite
        pass

    def test__init___(self):
        # Test the creation of an empty dependency_dict object.
        dependency_object = dependency_dict()
        self.assertIsInstance(dependency_object, dependency_dict)

    def test__setitem__dependency_type(self):
        # Test setting the dependency type
        dependency_object = dependency_dict()
        my_dependency_type = dependency_dict.dependency_types['co_modality']
        dependency_object['dependency_type'] = my_dependency_type
        self.assertEquals(dependency_object['dependency_type'], dependency_dict.dependency_types['co_modality'])

    def test__setitem__illegal_dependency_type(self):
        # Test setting the dependency type
        dependency_object = dependency_dict()
        my_dependency_type = 'illegal_dependency_type'
        with self.assertRaises(ValueError):
            dependency_object['dependency_type'] = my_dependency_type

    def test__setitem__param_name(self):
        # Test setting the parameter name
        dependency_object = dependency_dict()
        my_parameter = 'my_name'
        dependency_object['param_name'] = my_parameter
        self.assertEquals(dependency_object['param_name'], my_parameter)

    def test__setitem__link_name(self):
        # Test setting the _data which is not permitted
        dependency_object = dependency_dict()
        my_parameter = 'my_name'
        dependency_object['link_name'] = my_parameter
        self.assertEquals(dependency_object['link_name'], my_parameter)

    def test__setitem__illegal_param_name(self):
        # Test setting param_name to a non-string object
        dependency_object = dependency_dict()
        my_parameter = 5
        with self.assertRaises(ValueError):
            dependency_object['param_name']=  my_parameter

    def test__setitem__illegal_link_name(self):
        # Test setting link_name to an illegal object
        dependency_object = dependency_dict()
        my_parameter = 5
        with self.assertRaises(ValueError):
            dependency_object['link_name']=  my_parameter
        self.assertIsNone(dependency_object['_data'])

    def test__setitem__dataname(self):
        # Test setting the dataname parameter
        dependency_object = dependency_dict()
        my_parameter = 'my_name'
        dependency_object['dataname'] = my_parameter
        self.assertEquals(dependency_object['dataname'], my_parameter)
        self.assertIsNone(dependency_object['_data'])

    def test__setitem__illegal_dataname(self):
        # Test setting dataname to an illegal non-string value
        dependency_object = dependency_dict()
        my_parameter = 5
        with self.assertRaises(ValueError):
            dependency_object['dataname'] = my_parameter

    def test__setitem___data_key(self):
        # Test setting the _data which is not permitted
        dependency_object = dependency_dict()
        with self.assertRaises(KeyError):
            dependency_object['_data'] = 5

    def test__setitem__illegal_key(self):
        # Test setting an illegal key
        dependency_object = dependency_dict()
        with self.assertRaises(KeyError):
            dependency_object['illegal_key'] = 5

    def test__setitem__selection_string_index(self):
        # Test setting an index-based selection string
        from omsi.shared.omsi_data_selection import selection_to_string
        dependency_object = dependency_dict()
        selection_string = selection_to_string(5)
        dependency_object['selection'] = selection_string
        self.assertEquals(dependency_object['selection'], selection_string)

    def test__setitem__selection_string_list(self):
        # Test setting an list-based selection string
        from omsi.shared.omsi_data_selection import selection_to_string
        dependency_object = dependency_dict()
        selection_string = selection_to_string([0,1,2,3])
        dependency_object['selection'] = selection_string
        self.assertEquals(dependency_object['selection'], selection_string)

    def test__setitem__selection_string_slice(self):
        # Test setting an slice-based selection string
        from omsi.shared.omsi_data_selection import selection_to_string
        dependency_object = dependency_dict()
        selection_string = selection_to_string(slice(0,10,1))
        dependency_object['selection'] = selection_string
        self.assertEquals(dependency_object['selection'], selection_string)

    def test__setitem__selection_string_multi_axis(self):
        # Test setting an index-based selection string
        from omsi.shared.omsi_data_selection import selection_to_string
        dependency_object = dependency_dict()
        selection_string = selection_to_string((5,[0,1,3],slice(0,1,4)))
        dependency_object['selection'] = selection_string
        self.assertEquals(dependency_object['selection'], selection_string)

    def test__getitem__selection(self):
        # Test whether getting the selection value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['selection']
        except KeyError:
            self.fail("Getting the selection failed.")

    def test__getitem__link_name(self):
        # Test whether getting the link_name value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['link_name']
        except KeyError:
            self.fail("Getting the link_name failed.")

    def test__getitem__dataname(self):
        # Test whether getting the dataname value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['dataname']
        except KeyError:
            self.fail("Getting the dataname failed.")

    def test__getitem__param_name(self):
        # Test whether getting the selection value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['selection']
        except KeyError:
            self.fail("Getting the selection failed.")

    def test__getitem__omsi_object(self):
        # Test whether getting the selection value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['selection']
        except KeyError:
            self.fail("Getting the selection failed.")

    def test__getitem___data(self):
        # Test whether getting the _data value works
        dependency_object = dependency_dict()
        try:
            a = dependency_object['_data']
        except KeyError:
            self.fail("Getting the _data failed.")


if __name__ == '__main__':
    unittest.main()
