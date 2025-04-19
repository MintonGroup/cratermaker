import unittest
from cratermaker import Target, get_scaling_model, available_scaling_models
import numpy as np

Scale = get_scaling_model("richardson2009")
class TestScale(unittest.TestCase):
   
    def test_scale_from_catalogue(self):
        # Test creating a scale from the catalogue
        scale = Scale(material_name="Soft Rock")
        self.assertEqual(scale.K1, 0.20)
        self.assertEqual(scale.mu, 0.55)
        self.assertEqual(scale.Ybar, 7.60e6)
        self.assertEqual(scale.target_density, 2250.0)
        self.assertIsInstance(scale.K1, np.float64)
        self.assertIsInstance(scale.mu, np.float64)
        self.assertIsInstance(scale.Ybar, np.float64)
        self.assertIsInstance(scale.target_density, np.float64)
        return

    def test_scale_override_catalogue(self):
        # Test overriding a property from the catalogue
        scale = Scale(material_name="Soft Rock", Ybar=0.0)
        self.assertEqual(scale.Ybar, 0.0)
        return

    def test_custom_scale(self):
        # Test creating a custom scale
        scale = Scale(material_name="Flubber", K1=3.8, mu=0.1, Ybar=1e7, target_density=2000.0)
        self.assertEqual(scale.K1, 3.8)
        self.assertEqual(scale.mu, 0.1)
        self.assertEqual(scale.Ybar, 1e7)
        self.assertEqual(scale.target_density, 2000.0)
        return

    def test_incomplete_custom_scale(self):
        # Test incomplete custom scale creation
        with self.assertRaises(ValueError):
            Scale(material_name="Flubber", target_density=2000.0)
        with self.assertRaises(ValueError):
            Scale(material_name="Blorp", mu=0.1, Ybar=1e7, target_density=2000.0)
        return
    
    def test_scale_override_catalogue(self):
        # Test overriding a property from the catalogue
        target = Target(name="Mars")
        scale = Scale(target=target,material_name="Sand")
        self.assertEqual(scale.material_name, "Sand")    

if __name__ == '__main__':
    unittest.main()