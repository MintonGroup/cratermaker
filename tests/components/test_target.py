import unittest

from cratermaker import Target
from cratermaker.utils.general_utils import _create_catalogue


class TestTarget(unittest.TestCase):
    def test_target_from_catalogue(self):
        # Test creating a target from the built-in catalogue
        planets = ["Mercury", "Venus", "Earth", "Moon", "Mars", "Ceres", "Vesta"]
        for planet in planets:
            target = Target(name=planet)
            v = Target._catalogue[planet]
            self.assertEqual(target.radius, v["radius"])
            self.assertEqual(target.mass, v["mass"])
            self.assertEqual(target.material, v["material"])
            self.assertEqual(target.transition_scale_type, v["transition_scale_type"])
            self.assertIsInstance(target.radius, float)
            self.assertIsInstance(target.mass, float)
            self.assertIsInstance(target.material, str)
            self.assertIsInstance(target.transition_scale_type, str)

        # Test creating a target from a custom catalogue
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
        for planet, v in catalogue.items():
            target = Target(name=planet, catalogue=catalogue)
            self.assertEqual(target.radius, v["radius"])
            self.assertEqual(target.mass, v["mass"])
            self.assertEqual(target.material, v["material"])
            self.assertEqual(target.transition_scale_type, v["transition_scale_type"])
            self.assertIsInstance(target.radius, float)
            self.assertIsInstance(target.mass, float)

    def test_custom_target(self):
        # Test creating a custom target
        target = Target(
            name="Arrakis",
            radius=5995.0e3,
            mass=4.98e24,
            material="Sand",
            transition_scale_type="silicate",
        )
        self.assertEqual(target.radius, 5995.0e3)
        self.assertEqual(target.mass, 4.98e24)
        self.assertEqual(target.material, "Sand")

    def test_invalid_target(self):
        # Test incomplete custom target creation or invalid arguments
        with self.assertRaises(ValueError):
            Target("Arrakis")

        # Missing name
        with self.assertRaises(TypeError):
            Target(
                radius=5995.0,
                mass=4.98e24,
                material="Sand",
                transition_scale_type="silicate",
            )
        with self.assertRaises(ValueError):
            Target("Moon", transition_scale_type="flubber")


if __name__ == "__main__":
    unittest.main()
