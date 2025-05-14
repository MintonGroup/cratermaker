import unittest

import numpy as np

from cratermaker import Projectile, Target
from cratermaker.utils.general_utils import _create_catalogue


class TestProjectile(unittest.TestCase):
    def test_sample_true_requires_mean_velocity(self):
        with self.assertRaises(ValueError):
            Projectile.maker("generic", sample=True)

    def test_sample_false_requires_velocity(self):
        with self.assertRaises(ValueError):
            Projectile.maker("generic", sample=False)

    def test_sample_false_defaults(self):
        projectile = Projectile.maker("generic", sample=False, velocity=15e3)
        self.assertEqual(projectile.velocity, 15e3)
        self.assertEqual(projectile.angle, 90.0)
        self.assertEqual(projectile.direction, 0.0)
        self.assertEqual(projectile.density, 1000.0)
        self.assertEqual(projectile.location, (0.0, 0.0))
        projectile = projectile.new_projectile()
        self.assertEqual(projectile.velocity, 15e3)
        self.assertEqual(projectile.angle, 90.0)
        self.assertEqual(projectile.direction, 0.0)
        self.assertEqual(projectile.density, 1000.0)
        self.assertEqual(projectile.location, (0.0, 0.0))

    def test_sample_false_custom_angle_direction(self):
        projectile = Projectile.maker(
            "generic",
            sample=False,
            velocity=15e3,
            angle=45.0,
            direction=180.0,
            location=(120.0, -50.0),
        )
        self.assertEqual(projectile.angle, 45.0)
        self.assertEqual(projectile.direction, 180.0)
        self.assertEqual(projectile.density, 1000.0)
        self.assertEqual(projectile.location, (120.0, -50.0))
        projectile = projectile.new_projectile()
        self.assertEqual(projectile.velocity, 15e3)
        self.assertEqual(projectile.angle, 45.0)
        self.assertEqual(projectile.direction, 180.0)
        self.assertEqual(projectile.density, 1000.0)
        self.assertEqual(projectile.location, (120.0, -50.0))

    def test_density_default_and_override(self):
        default = Projectile.maker("generic", sample=False, velocity=10e3)
        self.assertEqual(default.density, 1000.0)
        projectile = Projectile.maker(
            "generic", sample=False, velocity=10e3, density=3000.0
        )
        self.assertEqual(projectile.density, 3000.0)

        projectile = projectile.new_projectile()
        self.assertEqual(projectile.velocity, 10e3)
        self.assertEqual(projectile.angle, 90.0)
        self.assertEqual(projectile.direction, 0.0)
        self.assertEqual(projectile.density, 3000.0)

    def test_sampled_vals(self):
        popnames = ["generic", "asteroids", "comets"]
        for popname in popnames:
            projectile = Projectile.maker(popname, sample=True, mean_velocity=20e3)

            velocities = []
            angles = []
            directions = []
            lons = []
            lats = []
            Nsamples = 1000
            for _ in range(Nsamples):
                projectile = projectile.new_projectile()
                velocities.append(projectile.velocity)
                angles.append(projectile.angle)
                directions.append(projectile.direction)
                lons.append(projectile.location[0])
                lats.append(projectile.location[1])
            velocities = np.array(velocities)
            angles = np.array(angles)
            directions = np.array(directions)
            lons = np.array(lons)
            lats = np.array(lats)

            vmean = np.mean(velocities)
            anglemean = np.mean(angles)
            directionmean = np.mean(directions)

            self.assertAlmostEqual(vmean / projectile.mean_velocity, 1.0, delta=0.1)
            self.assertAlmostEqual(anglemean / 45, 1.0, delta=0.1)
            self.assertAlmostEqual(directionmean / 180, 1.0, delta=0.1)

    def test_defaults(self):
        inner = ["Mercury", "Venus", "Earth", "Moon", "Vesta", "Ceres"]
        outer = ["Europa", "Umbriel", "Titan", "Triton", "Pluto", "Charon"]
        unknown = ["Arrakis", "Salusa Secundus"]
        body_properties = [
            "name",
            "radius",
            "mass",
            "material",
            "transition_scale_type",
        ]
        body_values = [
            ("Arrakis", 5995.0e3, 4.98e24, "Sand", "silicate"),
            ("Salusa Secundus", 7200.0e3, 8.62e24, "Hard Rock", "silicate"),
        ]
        catalogue = _create_catalogue(body_properties, body_values)

        for target in inner:
            projectile = Projectile.maker(target=target)
            self.assertEqual(projectile.target.name, target)
            self.assertEqual(projectile.population, "asteroids")
        for target in outer:
            projectile = Projectile.maker(target=target)
            self.assertEqual(projectile.target.name, target)
            self.assertEqual(projectile.population, "comets")
        for planet in unknown:
            target = Target.maker(target=planet, catalogue=catalogue)
            projectile = Projectile.maker(target=target)
            self.assertEqual(projectile.target.name, planet)
            self.assertEqual(projectile.population, "generic")


if __name__ == "__main__":
    unittest.main()
