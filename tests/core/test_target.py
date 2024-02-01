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
            "name",    "radius",   "gravity",      "material_name", "transition_scale_type"
        ]
        body_values = [
            ("Arrakis", 5995.0e3,  0.900 * gEarth, "Sand", "silicate"),
            ("Salusa Secundus", 7200.e3, 1.20 * gEarth, "Hard Rock", "silicate"),
        ]          
        catalogue = create_catalogue(body_properties, body_values) 
        for planet,v in catalogue.items():
            target = Target(name=planet, catalogue=catalogue)
            self.assertEqual(target.radius, v['radius'])
            self.assertEqual(target.gravity, v['gravity'])
            self.assertEqual(target.material_name, v['material_name'])
            self.assertEqual(target.transition_scale_type, v['transition_scale_type'])
            self.assertIsInstance(target.radius, np.float64)
            self.assertIsInstance(target.gravity, np.float64)


    def test_custom_target(self):
        # Test creating a custom target
        target = Target(name="Arrakis", radius=5995.0e3, gravity=0.9*gEarth, material_name="Sand", transition_scale_type="silicate")
        self.assertEqual(target.radius, 5995.0e3)
        self.assertEqual(target.gravity, 0.9*gEarth)
        self.assertEqual(target.material_name, "Sand")

    def test_invalid_target(self):
        # Test incomplete custom target creation or invalid arguments
        with self.assertRaises(ValueError):
            Target("Arrakis")
        
        # Missing name
        with self.assertRaises(TypeError):
            Target(radius=5995.0, gravity=0.9*gEarth, material_name="Sand", transition_scale_type="silicate")
        with self.assertRaises(ValueError):
            Target("Moon",transition_scale_type="flubber")

if __name__ == '__main__':
    unittest.main()