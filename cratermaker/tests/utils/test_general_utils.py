import unittest
from unittest.mock import patch, MagicMock
from cratermaker import general_utils as gu  # Adjust the import according to your project's structure

class TestSetProperties(unittest.TestCase):

    def setUp(self):
        # Set up any initial configurations or objects you need before each test
        self.mock_object = MagicMock()
        
   

    #  def test_direct_kwargs_override(self):
    #      # Test that direct kwargs override catalogue values
    #      catalogue = {'property1': 'catalogue_value', 'property2': 'catalogue_value'}
    #      gu.set_properties(self.mock_object, catalogue=catalogue, property1='kwarg_value')
    #      self.assertEqual(self.mock_object.property1, 'kwarg_value')
    #      self.assertEqual(self.mock_object.property2, 'catalogue_value')

    #  @patch('builtins.open', create=True)  # Adjust 'your_module' accordingly
    #  def test_json_file_setting(self, mock_open):
    #      # Mock reading from a JSON file
    #      mock_open.return_value.__enter__.return_value.read.return_value = '{"property1": "file_value"}'
    #      gu.set_properties(self.mock_object, filename='dummy.json')
    #      self.assertEqual(self.mock_object.property1, 'file_value')


if __name__ == '__main__':
    unittest.main()
