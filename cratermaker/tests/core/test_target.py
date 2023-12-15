import unittest
from cratermaker import Target
from cratermaker import create_catalogue
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
            self.assertEqual(target.transition_scale_type, v['transition_scale_type'])
            self.assertIsInstance(target.radius, np.float64)
            self.assertIsInstance(target.gravity, np.float64)
            self.assertIsInstance(target.material_name, str)
            self.assertIsInstance(target.transition_scale_type, str)
            
        # Test creating a target from a custom catalogue
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity", "transition_scale_type"
        ]
        body_values = [
            ("Arrakis", 5995.0e3,  0.900 * gEarth, "Sand", 18000.0, "silicate"),
            ("Salusa Secundus", 7200.e3, 1.20 * gEarth, "Hard Rock", 29100.0, "silicate"),
        ]          
        catalogue = create_catalogue(body_properties, body_values) 
        for planet,v in catalogue.items():
            target = Target(name=planet, catalogue=catalogue)
            self.assertEqual(target.radius, v['radius'])
            self.assertEqual(target.gravity, v['gravity'])
            self.assertEqual(target.material_name, v['material_name'])
            self.assertEqual(target.mean_impact_velocity, v['mean_impact_velocity'])
            self.assertEqual(target.transition_scale_type, v['transition_scale_type'])
            self.assertIsInstance(target.radius, np.float64)
            self.assertIsInstance(target.gravity, np.float64)
            self.assertIsInstance(target.mean_impact_velocity, np.float64)        

    def test_target_override_catalogue(self):
        # Test overriding a property from the catalogue
        target = Target(name="Mars",material_name="Sand")
        self.assertEqual(target.material.name, "Sand")

    def test_custom_target(self):
        # Test creating a custom target
        target = Target(name="Arrakis", radius=5995.0e3, gravity=0.9*gEarth, material_name="Sand", mean_impact_velocity=18000.0, transition_scale_type="silicate")
        self.assertEqual(target.radius, 5995.0e3)
        self.assertEqual(target.gravity, 0.9*gEarth)
        self.assertEqual(target.material_name, "Sand")
        self.assertEqual(target.mean_impact_velocity, 18000.0)

    def test_invalid_target(self):
        # Test incomplete custom target creation or invalid arguments
        with self.assertRaises(ValueError):
            Target("Arrakis", mean_impact_velocity=18000.0)
        with self.assertRaises(ValueError):
            Target(radius=5995.0, gravity=0.9*gEarth, material_name="Sand", mean_impact_velocity=18000.0, transition_scale_type="silicate")
        with self.assertRaises(ValueError):
            Target("Moon",transition_scale_type="flubber")

if __name__ == '__main__':
    unittest.main()