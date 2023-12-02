import unittest
import cratermaker
from cratermaker import Projectile, Target
import numpy as np
from numpy.random import default_rng

class TestProjectile(unittest.TestCase):
    
    def test_projectile_initialization_with_diameter(self):
        diameter = 1000
        projectile = Projectile(diameter=diameter)
        self.assertEqual(projectile.diameter, diameter)
        self.assertEqual(projectile.radius, diameter / 2)

    def test_projectile_initialization_with_radius(self):
        radius = 500
        projectile = Projectile(radius=radius)
        self.assertEqual(projectile.radius, radius)
        self.assertEqual(projectile.diameter, radius * 2)
        
    def test_invalid_negative_values(self):
        with self.assertRaises(ValueError):
            Projectile(diameter=-1000)

    def test_invalid_combinations(self):
        with self.assertRaises(ValueError):
            Projectile(diameter=1000, radius=500)

    def test_non_default_target_and_rng(self):
        rng = default_rng()
        target = Target("Mars")
        projectile = Projectile(diameter=1000,target=target,rng=rng)
        self.assertEqual(projectile.diameter, 1000)

    def test_invalid_target_or_rng_type(self):
        with self.assertRaises(TypeError):
            Projectile(diameter=1000, target="invalid_target")
    
if __name__ == '__main__':
    unittest.main()