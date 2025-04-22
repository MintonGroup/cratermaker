import unittest
from cratermaker import Crater, Target
import numpy as np
from numpy.random import default_rng

class TestCrater(unittest.TestCase):
    
    def test_crater_initialization_with_diameter(self):
        diameter = 1000
        crater = Crater(diameter=diameter)
        self.assertEqual(crater.final_diameter, diameter)
        self.assertEqual(crater.final_radius, diameter / 2)
        self.assertIsNotNone(crater.transient_diameter)
        self.assertIsInstance(crater.transient_diameter, np.float64)        

    def test_crater_initialization_with_radius(self):
        radius = 500
        crater = Crater(radius=radius)
        self.assertEqual(crater.final_radius, radius)
        self.assertEqual(crater.final_diameter, radius * 2)
        
    def test_crater_initialization_with_transient_diameter(self):
        transient_diameter = 800
        crater = Crater(transient_diameter=transient_diameter)
        self.assertEqual(crater.transient_diameter, transient_diameter)
        self.assertGreater(crater.final_diameter,transient_diameter)

    def test_crater_initialization_with_transient_radius(self):
        transient_radius = 400
        crater = Crater(transient_radius=transient_radius)
        self.assertEqual(crater.transient_radius, transient_radius)
        # Add more assertions as needed

    def test_invalid_negative_values(self):
        with self.assertRaises(ValueError):
            Crater(diameter=-1000)

    def test_invalid_combinations(self):
        with self.assertRaises(ValueError):
            Crater(diameter=1000, radius=500)

    def test_non_default_target_and_rng(self):
        rng = default_rng()
        target = Target("Mars")
        crater = Crater(diameter=1000,target=target,rng=rng)
        self.assertEqual(crater.final_diameter, 1000)

    def test_invalid_target_or_rng_type(self):
        with self.assertRaises(ValueError):
            Crater(diameter=1000, target="invalid_target")
    
if __name__ == '__main__':
    unittest.main()
