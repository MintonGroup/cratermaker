import unittest

import numpy as np

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
        newvals = projectile.new_projectile()
        self.assertEqual(newvals["projectile_velocity"], 15e3)
        self.assertEqual(newvals["projectile_angle"], 90.0)
        self.assertEqual(newvals["projectile_direction"], 0.0)
        self.assertEqual(newvals["projectile_density"], 1000.0)

    def test_sample_false_custom_angle_direction(self):
        projectile = Projectile.maker(
            "generic", sample=False, velocity=15e3, angle=45.0, direction=180.0
        )
        self.assertEqual(projectile.angle, 45.0)
        self.assertEqual(projectile.direction, 180.0)
        self.assertEqual(projectile.density, 1000.0)
        newvals = projectile.new_projectile()
        self.assertEqual(newvals["projectile_velocity"], 15e3)
        self.assertEqual(newvals["projectile_angle"], 45.0)
        self.assertEqual(newvals["projectile_direction"], 180.0)
        self.assertEqual(newvals["projectile_density"], 1000.0)

    def test_density_default_and_override(self):
        default = Projectile.maker("generic", sample=False, velocity=10e3)
        self.assertEqual(default.density, 1000.0)
        projectile = Projectile.maker(
            "generic", sample=False, velocity=10e3, density=3000.0
        )
        self.assertEqual(projectile.density, 3000.0)

        newvals = projectile.new_projectile()
        self.assertEqual(newvals["projectile_velocity"], 10e3)
        self.assertEqual(newvals["projectile_angle"], 90.0)
        self.assertEqual(newvals["projectile_direction"], 0.0)
        self.assertEqual(newvals["projectile_density"], 3000.0)

    def test_sampled_vals(self):
        popnames = ["generic", "asteroids", "comets"]
        for popname in popnames:
            projectile = Projectile.maker(popname, sample=True, mean_velocity=20e3)

            velocities = []
            angles = []
            directions = []
            Nsamples = 1000
            for _ in range(Nsamples):
                newvals = projectile.new_projectile()
                velocities.append(newvals["projectile_velocity"])
                angles.append(newvals["projectile_angle"])
                directions.append(newvals["projectile_direction"])
            velocities = np.array(velocities)
            angles = np.array(angles)
            directions = np.array(directions)

            vmean = np.mean(velocities)
            anglemean = np.mean(angles)
            directionmean = np.mean(directions)

            self.assertAlmostEqual(vmean / projectile.mean_velocity, 1.0, delta=0.1)
            self.assertAlmostEqual(anglemean / 45, 1.0, delta=0.1)
            self.assertAlmostEqual(directionmean / 180, 1.0, delta=0.1)


if __name__ == "__main__":
    unittest.main()
