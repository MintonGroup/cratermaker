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
            diameter=100000.0,
        )

    def test_model_registration(self):
        models = Morphology.available()
        self.assertIn("basicmoon", models)

    def test_model_instantiation(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True)
            for model_name in morphology_models:
                morphology = Morphology.maker(model_name, surface=surface)
                self.assertIsInstance(morphology, object)
                self.assertEqual(morphology.component_name, model_name)

    def test_form_crater_executes(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(model_name, surface=surface)
                morphology.emplace(self.dummy_crater)

    def test_make_morphology(self):
        # Test the make_morphology function
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(morphology=model_name, surface=surface)
                self.assertIsInstance(morphology, Morphology)

    def test_finite_profile_values(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(morphology=model_name, surface=surface)
                crater_radius_values = [1.0, 1e3, 15e3, 50e3, 500e3, 3000e3]
                rvals = np.linspace(0, 10, 1000)
                bvals = np.zeros_like(rvals)
                for radius in crater_radius_values:
                    crater = Crater.maker(radius=radius)
                    crater_shape = morphology.crater_profile(crater, rvals * radius, bvals)
                    self.assertTrue(
                        np.all(np.isfinite(crater_shape)),
                        f"Crater profile for {model_name} contains NaN or Inf values.",
                    )
                    ejecta_shape = morphology.ejecta_profile(crater, rvals * radius, bvals)
                    self.assertTrue(
                        np.all(np.isfinite(ejecta_shape)),
                        f"Ejecta profile for {model_name} contains NaN or Inf values.",
                    )

    def test_ejecta_profile_uses_bearings(self):
        # The realmoon rim varies with bearing, so its ejecta profile near the rim must too. Passing the bearings through to the
        # profile function previously failed silently and every point was evaluated at bearing 0.
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            morphology = Morphology.maker(morphology="realmoon", surface=surface)
            crater = morphology.Crater.maker(diameter=50e3, location=(0.0, 0.0))
            rvals = np.full(36, crater.radius)
            bvals = np.linspace(0.0, 350.0, 36)
            ejecta = morphology.ejecta_profile(crater, rvals, bvals)
            ejecta_at_zero = morphology.ejecta_profile(crater, rvals, np.zeros_like(bvals))
            self.assertFalse(np.array_equal(ejecta, ejecta_at_zero), "realmoon ejecta profile ignores bearings")

    def test_ejecta_distribution_executes(self):
        # The ray (dorays=True) path calls ejecta_profile with bearings as a keyword argument
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(morphology=model_name, surface=surface)
                crater = morphology.Crater.maker(diameter=50e3, location=(0.0, 0.0))
                rvals = np.linspace(0.0, 5.0 * crater.radius, 50)
                bvals = np.linspace(0.0, 350.0, 50)
                thickness, intensity = morphology.ejecta_distribution(crater, rvals, bvals)
                self.assertEqual(thickness.shape, rvals.shape)
                self.assertTrue(np.all(np.isfinite(thickness)) and np.all(np.isfinite(intensity)))

    def test_form_craters_matches_sequential(self):
        # Forming a batch of non-overlapping craters in parallel must give exactly the same surface as forming them one at a time
        locations = [(0.0, 0.0), (90.0, 0.0), (180.0, 0.0), (-90.0, 0.0), (0.0, 60.0), (0.0, -60.0)]
        for model_name in morphology_models:
            results = []
            for batched in (False, True):
                with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                    surface = Surface.maker(
                        simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False
                    )
                    morphology = Morphology.maker(morphology=model_name, surface=surface, rng_seed=4321)
                    craters = [morphology.Crater.maker(diameter=100e3, location=loc) for loc in locations]
                    face_sets = [set(np.asarray(morphology.overlap_function(c)[1]).tolist()) for c in craters]
                    self.assertTrue(
                        all(not (face_sets[i] & face_sets[j]) for i in range(len(face_sets)) for j in range(i)),
                        "test craters must have non-overlapping regions",
                    )
                    if batched:
                        morphology.form_craters(craters)
                    else:
                        for crater in craters:
                            morphology.form_crater(crater)
                    variables = surface.uxds.variables
                    results.append(
                        [variables[name].data.copy() for name in ("face_elevation", "node_elevation", "ejecta_thickness")]
                    )
            for sequential, batched in zip(*results, strict=True):
                np.testing.assert_array_equal(batched, sequential, err_msg=f"{model_name}: batched result differs")

    def test_form_craters_rejects_overlapping_regions(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(morphology=model_name, surface=surface)
                craters = [morphology.Crater.maker(diameter=150e3, location=(0.0, 0.0)) for _ in range(2)]
                with self.assertRaises(ValueError):
                    morphology.form_craters(craters)

    def test_crater_depth_surface(self):
        # Tests that the surface elevations are expected
        diameter_list = [200e3, 500e3, 1000e3]
        delta_vals = [0.3, 0.3, 0.2]

        surface_args = {
            "icosphere": {"gridlevel": 6},
            "arbitrary_resolution": {
                "pix": 30e3,
            },
            "hireslocal": {
                "pix": 1.0e3,
                "local_location": (0, 0),
                "local_radius": 100e3,
            },
        }

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            for name, args in surface_args.items():
                sim = Simulation(simdir=simdir, surface=name, ask_overwrite=False, rng_seed=8675309, **args)
                for diameter, delta in zip(diameter_list, delta_vals, strict=False):
                    sim.reset()
                    # verify that the surface is flat
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.min(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.min(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.max(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.max(),
                        0.0,
                        delta=1e0,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )

                    crater = sim.Crater.maker(diameter=diameter, location=(0, 0), peak_height=0.0)
                    sim.emplace(crater)

                    # Verify that the crater depth and rim heights are close to the expected values
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.min() / crater.floor_elevation,
                        1.0,
                        delta=delta,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.min() / crater.floor_elevation,
                        1.0,
                        delta=delta,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.node_elevation.max() / crater.rim_height,
                        1.0,
                        delta=2 * delta,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )
                    self.assertAlmostEqual(
                        sim.surface.face_elevation.max() / crater.rim_height,
                        1.0,
                        delta=2 * delta,
                        msg=f"Failed for {name} with diameter {diameter}",
                    )


if __name__ == "__main__":
    unittest.main()
