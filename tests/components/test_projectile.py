import unittest

from cratermaker import Projectile


class TestProjectile(unittest.TestCase):
    def test_sample_true_requires_mean_velocity(self):
        with self.assertRaises(ValueError):
            Projectile.maker("generic", sample=True)

    def test_sample_false_requires_velocity(self):
        with self.assertRaises(ValueError):
            Projectile.maker("generic", sample=False)

    def test_sample_false_defaults(self):
        projectile = Projectile.maker("generic", sample=False, velocity=15e3)
        self.assertEqual(projectile.angle, 90.0)
        self.assertEqual(projectile.direction, 0.0)
        self.assertEqual(projectile.density, 1000.0)

    def test_sample_false_custom_angle_direction(self):
        projectile = Projectile.maker(
            "generic", sample=False, velocity=15e3, angle=45.0, direction=180.0
        )
        self.assertEqual(projectile.angle, 45.0)
        self.assertEqual(projectile.direction, 180.0)
        self.assertEqual(projectile.density, 1000.0)

    def test_density_default_and_override(self):
        default = Projectile.maker("generic", sample=False, velocity=10e3)
        self.assertEqual(default.density, 1000.0)
        custom = Projectile.maker(
            "generic", sample=False, velocity=10e3, density=3000.0
        )
        self.assertEqual(custom.density, 3000.0)


if __name__ == "__main__":
    unittest.main()
