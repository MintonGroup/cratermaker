import unittest
from cratermaker import Target
from cratermaker import general_utils as gu
import numpy as np
gEarth = np.float64(9.80665)

class TestTarget(unittest.TestCase):

    def test_target_from_catalogue(self):
        # Test creating a target from the built-in catalogue
        planets = ["Mercury", "Venus", "Earth", "Moon", "Mars", "Ceres", "Vesta" ]
        for planet in  planets:
            target = Target(name=planet)
            v = target.catalogue[planet]
            self.assertEqual(target.radius, v['radius'])
            self.assertEqual(target.gravity, v['gravity'])
            self.assertEqual(target.material_name, v['material_name'])
            self.assertEqual(target.mean_impact_velocity, v['mean_impact_velocity'])
            self.assertIsInstance(target.radius, np.float64)
            self.assertIsInstance(target.gravity, np.float64)
            self.assertIsInstance(target.mean_impact_velocity, np.float64)
            
        # Test creating a target from a custom catalogue
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity"
        ]
        body_values = [
            ("Arrakis", 5995.0e3,  0.900 * gEarth, "Sand", 18000.0),
            ("Salusa Secundus", 7200.e3, 1.20 * gEarth, "Hard Rock", 29100.0),
        ]          
        catalogue = gu.create_catalogue(body_properties, body_values) 
        for planet,v in catalogue.items():
            target = Target(name=planet, catalogue=catalogue)
            self.assertEqual(target.radius, v['radius'])
            self.assertEqual(target.gravity, v['gravity'])
            self.assertEqual(target.material_name, v['material_name'])
            self.assertEqual(target.mean_impact_velocity, v['mean_impact_velocity'])
            self.assertIsInstance(target.radius, np.float64)
            self.assertIsInstance(target.gravity, np.float64)
            self.assertIsInstance(target.mean_impact_velocity, np.float64)        

    def test_target_override_catalogue(self):
        # Test overriding a property from the catalogue
        target = Target(name="Mars",material_name="Sand")
        self.assertEqual(target.material.name, "Sand")

    def test_custom_target(self):
        # Test creating a custom target
        target = Target(name="Arrakis", radius=5995.0e3, gravity=0.9*gEarth, material_name="Sand", mean_impact_velocity=18000.0)
        self.assertEqual(target.radius, 5995.0e3)
        self.assertEqual(target.gravity, 0.9*gEarth)
        self.assertEqual(target.material_name, "Sand")
        self.assertEqual(target.mean_impact_velocity, 18000.0)

    def test_incomplete_custom_target(self):
        # Test incomplete custom target creation
        with self.assertRaises(ValueError):
            Target("Arrakis", mean_impact_velocity=18000.0)
            Target(radius=5995.0, gravity=0.9*gEarth, material_name="Sand", mean_impact_velocity=18000.0)

if __name__ == '__main__':
    unittest.main()