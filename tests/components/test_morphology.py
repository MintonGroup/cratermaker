import tempfile
import unittest

import numpy as np

from cratermaker import Crater, Simulation
from cratermaker.components.morphology import Morphology
from cratermaker.components.surface import Surface
from cratermaker.components.target import Target

morphology_models = Morphology.available()


class TestMorphology(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 5

        self.dummy_crater = Crater.maker(
            location=(0.0, 0.0),
            final_diameter=100000.0,
        )

    def test_model_registration(self):
        models = Morphology.available()
        self.assertIn("simplemoon", models)

    def test_model_instantiation(self):
        for model_name in morphology_models:
            morphology = Morphology.maker(model_name)
            self.assertIsInstance(morphology, object)
            self.assertEqual(morphology.name, model_name)
            morphology.crater = self.dummy_crater
            self.assertIs(morphology.crater, self.dummy_crater)

    def test_form_crater_executes(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(
                simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel
            )
            for model_name in morphology_models:
                morphology = Morphology.maker(model_name, surface=surface)
                morphology.emplace(self.dummy_crater)

    def test_make_morphology(self):
        # Test the make_morphology function
        for model_name in morphology_models:
            morphology = Morphology.maker(moprhology=model_name)
            self.assertIsInstance(morphology, Morphology)

    def test_finite_profile_values(self):
        for model_name in morphology_models:
            morphology = Morphology.maker(model_name)
            crater_radius_values = [1.0, 1e3, 15e3, 50e3, 500e3, 3000e3]
            rvals = np.linspace(0, 10, 1000)
            for final_radius in crater_radius_values:
                crater = Crater.maker(final_radius=final_radius)
                crater_shape = morphology.crater_profile(crater, rvals * final_radius)
                self.assertTrue(
                    np.all(np.isfinite(crater_shape)),
                    f"Crater profile for {model_name} contains NaN or Inf values.",
                )
                ejecta_shape = morphology.ejecta_profile(crater, rvals * final_radius)
                self.assertTrue(
                    np.all(np.isfinite(ejecta_shape)),
                    f"Ejecta profile for {model_name} contains NaN or Inf values.",
                )

    def test_crater_depth_surface(self):
        from cratermaker.components.morphology.simplemoon import SimpleMoonCrater
        # Tests that the surface elevations are expected

        final_diameter_list = [100e3, 200e3, 500e3, 1000e3]
        delta_vals = [0.4, 0.3, 0.3, 0.2]

        gridargs = {
            "icosphere": {"gridlevel": 6},
            "arbitrary_resolution": {
                "pix": 10e3,
            },
            "hireslocal": {
                "pix": 1.0e3,
                "local_location": (0, 0),
                "local_radius": 100e3,
            },
        }

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            for name, args in gridargs.items():
                sim = Simulation(simdir=simdir, surface=name, **args)
                for final_diameter, delta in zip(final_diameter_list, delta_vals):
                    sim.surface.reset()
                    # verify that the surface is flat
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.min(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.min(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.max(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.max(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )

                    crater = SimpleMoonCrater.maker(
                        final_diameter=final_diameter, location=(0, 0)
                    )
                    sim.emplace(crater)

                    # Verify that the crater depth and rim heights are close to the expected values
                    self.assertAlmostEqual(
                        -sim.surface.node_elevation.min() / crater.floor_depth,
                        1.0,
                        delta=delta,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        -sim.surface.face_elevation.min() / crater.floor_depth,
                        1.0,
                        delta=delta,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.max() / crater.rim_height,
                        1.0,
                        delta=2 * delta,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.max() / crater.rim_height,
                        1.0,
                        delta=2 * delta,
                        msg=f"Failed for {name} with diameter {final_diameter}",
                    )


if __name__ == "__main__":
    unittest.main()
