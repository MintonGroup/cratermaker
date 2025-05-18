import unittest

import numpy as np

from cratermaker.components.production import Production

production_models = Production.available()


class TestProduction(unittest.TestCase):
    def test_production_N_to_time(self):
        for model_name in production_models:
            production = Production.maker(production=model_name)
            for age_orig in np.linspace(0, 4500, num=10):
                D = np.logspace(-1, 3, num=1000)
                N = production.function(diameter=D, age=age_orig, check_valid_age=False)
                age_new = production.function_inverse(
                    cumulative_number_density=N, diameter=D
                )
                np.testing.assert_array_almost_equal(age_orig, age_new, decimal=2)

    def test_sample_arguments(self):
        for model_name in production_models:
            production = Production.maker(production=model_name)

            # Test valid arguments
            try:
                production.sample(age=500.0, diameter_range=(10.0, 100.0), area=1e6)
            except ValueError:
                self.fail(
                    "sample() raised ValueError unexpectedly with valid arguments!"
                )

            # Test providing both age and cumulative_number
            with self.assertRaises(ValueError):
                production.sample(
                    age=500, diameter_number=(100, 10), diameter_range=(10, 100)
                )

            # Test providing both age_end and diameter_number_end
            with self.assertRaises(ValueError):
                production.sample(
                    age=1500,
                    age_end=1000,
                    diameter_number_end=(100, 10),
                    diameter_range=(10, 100),
                )

            # Test missing diameter_range
            with self.assertRaises(ValueError):
                production.sample(age=500)

            # Test diameter_range with invalid length
            with self.assertRaises(ValueError):
                production.sample(age=500, diameter_range=(10,), area=1000)

            # Test diameter_range with non-positive minimum diameter
            with self.assertRaises(ValueError):
                production.sample(age=500, diameter_range=(-1, 100), area=1000)

            # Test diameter_range with max diameter less than min diameter
            with self.assertRaises(ValueError):
                production.sample(age=500, diameter_range=(100, 10), area=1000)

            # Test invalid diameter_number length
            with self.assertRaises(ValueError):
                production.sample(diameter_number=(100,), diameter_range=(10, 100))

            # Test negative value in diameter_number
            with self.assertRaises(ValueError):
                production.sample(diameter_number=(-100, 10), diameter_range=(10, 100))

            # Test negative area
            with self.assertRaises(ValueError):
                production.sample(age=500, diameter_range=(10, 100), area=-1000)

            # Test non-scalar age
            with self.assertRaises(ValueError):
                production.sample(age=[500, 1000], diameter_range=(10, 100))

            # This checks if a ValueError is raised for an invalid length of diameter_number_end.
            with self.assertRaises(ValueError):
                production.sample(diameter_number_end=(100,), diameter_range=(10, 100))

            # Testing when a negative value is provided in diameter_number_end.
            with self.assertRaises(ValueError):
                production.sample(
                    diameter_number_end=(-100, 10), diameter_range=(10, 100)
                )

    def test_small_sample(self):
        neukum = Production.maker(production="neukum", version="Moon")
        # Test that the number of craters is zero when the age and area are ridiculously tiny
        diameter, age = neukum.sample(
            age=1e-6, diameter_range=(300e3, 1000e3), area=1e-6
        )
        self.assertEqual(diameter.size, 0)
        self.assertEqual(age.size, 0)

    def testmake_production(self):
        # Test that Production.maker returns a valid production model
        for model_name in production_models:
            production = Production.maker(production=model_name)
            self.assertIn(model_name, production_models)
            self.assertIsInstance(production, Production)


if __name__ == "__main__":
    unittest.main()
