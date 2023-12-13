import unittest
from unittest.mock import patch, MagicMock
from cratermaker import to_config, set_properties, check_properties, create_catalogue, validate_and_convert_location, normalize_coords  
from cratermaker.utils.custom_types import FloatLike

mock_properties = [
    "name",       "property1",     "property2"
]
mock_values = [
    ("foo",     1.00,     2.00),
    ("bar",     3.00,     4.00),
]  
class TestSetProperties(unittest.TestCase):
    def setUp(self):
        # Set up any initial configurations or objects you need before each test
        self.mock_object = MagicMock()
        return
        
    def test_kwargs(self):
        set_properties(self.mock_object, name="baz", property1=5.0, property2=6.0)
        self.assertEqual(self.mock_object.name, "baz")
        self.assertEqual(self.mock_object.property1, 5.0)
        self.assertEqual(self.mock_object.property2, 6.0)
       
        # Reset the properties and try to cause an exception when a name property is not passed 
        with self.assertRaises(ValueError): 
            del self.mock_object.name  # Explicitly delete the name attribute
            set_properties(self.mock_object, property1=5.0, property2=6.0) # name omitted
        return    
    
    def test_catalogue(self):
        # Test a regular catalogue with multiple entries (must be selected by  name)
        catalogue = create_catalogue(mock_properties, mock_values)
        for k,v in catalogue.items():
            set_properties(self.mock_object, catalogue=catalogue, name=k)
            self.assertEqual(self.mock_object.name, k)
            self.assertEqual(self.mock_object.property1, v['property1'])
            self.assertEqual(self.mock_object.property2, v['property2'])
            
        # Test a single entry catalogue (name can be omitted)
        name = next(iter(catalogue))
        set_properties(self.mock_object, catalogue={name:catalogue[name]})
    
        # Create an invalid catalogue structure
        invalid_catalogue = {"invalid_key": "invalid_value"}
        
        # Test if ValueError is raised for invalid catalogue structure or other problems
        with self.assertRaises(ValueError):  
            set_properties(self.mock_object, catalogue="bad_value") # Not even a dict
        with self.assertRaises(ValueError):  
            set_properties(self.mock_object, catalogue=invalid_catalogue, name="invalid_key") # Not a nested dict
        with self.assertRaises(ValueError):  
            set_properties(self.mock_object, catalogue=catalogue) # Name omitted
        
        return

    @patch('builtins.open', create=True)  # Adjust 'your_module' accordingly
    def test_json_file_setting(self, mock_open):
        # Mock reading from a JSON file
        mock_open.return_value.__enter__.return_value.read.return_value = '{"qux":{"property1":7.0,"property2":8.0},"quux":{"property1":9.0,"property2":10.0}}'
        set_properties(self.mock_object, filename='dummy.json', name='qux')
        self.assertEqual(self.mock_object.name, 'qux')
        self.assertEqual(self.mock_object.property1, 7.0)
        self.assertEqual(self.mock_object.property2, 8.0)
        set_properties(self.mock_object, filename='dummy.json', name='quux')
        self.assertEqual(self.mock_object.name, 'quux')
        self.assertEqual(self.mock_object.property1, 9.0)
        self.assertEqual(self.mock_object.property2, 10.0)
        return 
        
    def test_kwarg_override(self):
        property2_new = 42.0
        catalogue = create_catalogue(mock_properties, mock_values)
        name = next(iter(catalogue))
        set_properties(self.mock_object, catalogue=catalogue, name=name, property2=property2_new)
        self.assertEqual(self.mock_object.name, name)
        self.assertEqual(self.mock_object.property1, catalogue[name]['property1'])
        self.assertEqual(self.mock_object.property2, property2_new)
        
        # test that overrides are ignored when None is passed:
        set_properties(self.mock_object, catalogue=catalogue, name=name, property2=None) 
        self.assertEqual(self.mock_object.property2, catalogue[name]['property2'])
        return

if __name__ == '__main__':
    unittest.main()
