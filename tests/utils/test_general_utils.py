import unittest
from unittest.mock import patch, mock_open
from cratermaker.utils.general_utils import _set_properties, parameter, _create_catalogue

mock_properties = [
    "name",       "prop",     "param"
]
mock_values = [
    ("foo",     1.00,     2.00),
    ("bar",     3.00,     4.00),
]  
class Dummy:
    def __init__(self):
        self._name = None
        self._prop = None
        self._param = None
        self._user_defined = set()

    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        self._name = value

    @property
    def prop(self):
        return self._prop

    @prop.setter
    def prop(self, value):
        self._prop = value

    @parameter
    def param(self):
        return self._param

    @param.setter
    def param(self, value):
        self._param = value

    @property
    def catalogue_key(self):
        return "name"

class TestGeneralUtils(unittest.TestCase):
    def setUp(self):
        # Set up any initial configurations or objects you need before each test
        self.dummy = Dummy()
        return
        
    def test_kwargs(self):
        _set_properties(self.dummy, name="baz", prop=5.0, param=6.0)
        self.assertEqual(self.dummy.name, "baz")
        self.assertEqual(self.dummy.prop, 5.0)
        self.assertEqual(self.dummy.param, 6.0)
       
        return    
    
    def test_catalogue(self):
        # Test a regular catalogue with multiple entries (must be selected by  name)
        catalogue = _create_catalogue(mock_properties, mock_values)
        for k,v in catalogue.items():
            _set_properties(self.dummy, catalogue=catalogue, name=k)
            self.assertEqual(self.dummy.name, k)
            self.assertEqual(self.dummy.prop, v['prop'])
            self.assertEqual(self.dummy.param, v['param'])
            
        # Test a single entry catalogue (name can be omitted)
        name = next(iter(catalogue))
        _set_properties(self.dummy, catalogue={name:catalogue[name]})
    
        # Create an invalid catalogue structure
        invalid_catalogue = {"invalid_key": "invalid_value"}
        
        # Test if ValueError is raised for invalid catalogue structure or other problems
        with self.assertRaises(ValueError):  
            _set_properties(self.dummy, catalogue="bad_value") # Not even a dict
        with self.assertRaises(ValueError):  
            _set_properties(self.dummy, catalogue=invalid_catalogue, key="invalid_key") # Not a nested dict
        
        return


    def test_kwarg_override(self):
        param_new = 42.0
        catalogue = _create_catalogue(mock_properties, mock_values)
        name = next(iter(catalogue))
        _set_properties(self.dummy, catalogue=catalogue, key=name, param=param_new)
        self.assertEqual(self.dummy.name, name)
        self.assertEqual(self.dummy.prop, catalogue[name]['prop'])
        self.assertEqual(self.dummy.param, param_new)
        
        # test that overrides are ignored when None is passed:
        _set_properties(self.dummy, catalogue=catalogue, key=name, param=None) 
        self.assertEqual(self.dummy.param, catalogue[name]['param'])
        return

if __name__ == '__main__':
    unittest.main()
