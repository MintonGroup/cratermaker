import unittest

from numpy.random import default_rng

from cratermaker import Crater, Target


class TestCrater(unittest.TestCase):
    def test_crater_initialization_with_diameter(self):
        final_diameter = 1000
        crater = Crater.maker(final_diameter=final_diameter)
        self.assertEqual(crater.final_diameter, final_diameter)
        self.assertEqual(crater.final_radius, final_diameter / 2)
        self.assertIsNotNone(crater.transient_diameter)
        self.assertIsInstance(crater.transient_diameter, float)

    def test_crater_initialization_with_radius(self):
        final_radius = 500
        crater = Crater.maker(final_radius=final_radius)
        self.assertEqual(crater.final_radius, final_radius)
        self.assertEqual(crater.final_diameter, final_radius * 2)

    def test_crater_initialization_with_transient_diameter(self):
        transient_diameter = 800
        crater = Crater.maker(transient_diameter=transient_diameter)
        self.assertEqual(crater.transient_diameter, transient_diameter)
        self.assertGreater(crater.final_diameter, transient_diameter)

    def test_crater_initialization_with_transient_radius(self):
        transient_radius = 400
        crater = Crater.maker(transient_radius=transient_radius)
        self.assertEqual(crater.transient_radius, transient_radius)

    def test_invalid_negative_crater_values(self):
        with self.assertRaises(ValueError):
            Crater.maker(final_diameter=-1000)

    def test_invalid_crater_combinations(self):
        with self.assertRaises(ValueError):
            Crater.maker(final_diameter=1000, final_radius=500)

    def test_non_default_target_and_rng(self):
        rng = default_rng()
        target = Target("Mars")
        crater = Crater.maker(final_diameter=1000, target=target, rng=rng)
        self.assertEqual(crater.final_diameter, 1000)

    def test_invalid_target_or_rng_type(self):
        with self.assertRaises(ValueError):
            Crater.maker(final_diameter=1000, target="invalid_target")

    def test_projectile_initialization_with_diameter(self):
        diameter = 1000
        crater = Crater.maker(projectile_diameter=diameter, projectile_velocity=20e3)
        self.assertEqual(crater.projectile_diameter, diameter)
        self.assertEqual(crater.projectile_radius, diameter / 2)

    def test_projectile_initialization_with_radius(self):
        radius = 500
        crater = Crater.maker(projectile_radius=radius, projectile_velocity=20e3)
        self.assertEqual(crater.projectile_radius, radius)
        self.assertEqual(crater.projectile_diameter, radius * 2)

    def test_invalid_negative_projectile_values(self):
        with self.assertRaises(ValueError):
            Crater.maker(projectile_diameter=-1000, projectile_velocity=20e3)

    def test_invalid_projectile_combinations(self):
        with self.assertRaises(ValueError):
            Crater.maker(
                projectile_diameter=1000,
                projectile_radius=500,
                projectile_velocity=20e3,
            )
        with self.assertRaises(ValueError):
            Crater.maker(projectile_diameter=1000, final_diameter=10000)

    def test_non_default_target_and_rng_crater(self):
        rng = default_rng()
        crater = Crater.maker(
            projectile_diameter=1000, rng=rng, projectile_mean_velocity=10e3
        )
        self.assertEqual(crater.projectile_diameter, 1000)

    def test_invalid_target_or_rng_type_projectile(self):
        with self.assertRaises(ValueError):
            Crater.maker(
                projectile_diameter=1000,
                target="invalid_target",
                projectile_mean_velocity=10e3,
            )

    def test_argument_override(self):
        crater = Crater.maker(final_diameter=1000, location=(0.0, 0.0))
        self.assertEqual(crater.location, (0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
