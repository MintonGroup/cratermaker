import unittest
from cratermaker import Material
import numpy as np

class TestMaterial(unittest.TestCase):
   
    def test_material_from_catalogue(self):
        # Test creating a material from the catalogue
        material = Material("Soft Rock")
        self.assertEqual(material.K1, 0.20)
        self.assertEqual(material.mu, 0.55)
        self.assertEqual(material.Ybar, 7.60e6)
        self.assertEqual(material.density, 2250.0)
        self.assertIsInstance(material.K1, np.float64)
        self.assertIsInstance(material.mu, np.float64)
        self.assertIsInstance(material.Ybar, np.float64)
        self.assertIsInstance(material.density, np.float64)
        return

    def test_material_override_catalogue(self):
        # Test overriding a property from the catalogue
        material = Material(name="Soft Rock", Ybar=0.0)
        self.assertEqual(material.Ybar, 0.0)
        return

    def test_custom_material(self):
        # Test creating a custom material
        material = Material(name="Flubber", K1=3.8, mu=0.1, Ybar=1e7, density=2000.0)
        self.assertEqual(material.K1, 3.8)
        self.assertEqual(material.mu, 0.1)
        self.assertEqual(material.Ybar, 1e7)
        self.assertEqual(material.density, 2000.0)
        return

    def test_incomplete_custom_material(self):
        # Test incomplete custom material creation
        with self.assertRaises(ValueError):
            Material("Flubber", density=2000.0)
            Material(K1=3.8, mu=0.1, Ybar=1e7, density=2000.0)
        return

if __name__ == '__main__':
    unittest.main()