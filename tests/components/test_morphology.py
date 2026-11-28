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

    def test_peak_ring_layout(self):
        # Central peaks open into peak rings from 200 to 300 km, at a nearly constant fraction of the rim radius
        from cratermaker.components.morphology.basicmoon import (
            _PEAK_RING_D_REF,
            _PEAK_RING_EXPONENT,
            _PEAK_RING_RATIO,
            peak_ring_layout,
        )

        self.assertEqual(peak_ring_layout(150e3), (0.0, 0.0))
        fraction, ring_radius = peak_ring_layout(_PEAK_RING_D_REF)
        self.assertEqual(fraction, 1.0)
        self.assertAlmostEqual(ring_radius / (_PEAK_RING_D_REF / 2), _PEAK_RING_RATIO)
        # The ring radius grows monotonically and continuously with diameter
        diameters = np.geomspace(150e3, 1500e3, 400)
        radii = np.array([peak_ring_layout(d)[1] for d in diameters])
        self.assertTrue(np.all(np.diff(radii) >= 0.0))
        self.assertLess(np.max(np.diff(radii)), 0.02 * radii[-1])
        _, ring_radius = peak_ring_layout(900e3)
        self.assertAlmostEqual(ring_radius / 450e3, _PEAK_RING_RATIO * (900e3 / _PEAK_RING_D_REF) ** _PEAK_RING_EXPONENT)

    def test_ring_layout_capped_at_orientale(self):
        # Beyond Orientale, the ring ratios are not extrapolated and the scatter of the ring radii stays fixed in m
        from cratermaker.components.morphology.basicmoon import _RING_DATA_MAX_DIAMETER, peak_ring_layout, scarp_layout

        D0 = _RING_DATA_MAX_DIAMETER
        self.assertAlmostEqual(peak_ring_layout(3 * D0)[1] / (1.5 * D0), peak_ring_layout(D0)[1] / (0.5 * D0))
        for layout in (lambda d, rng: peak_ring_layout(d, rng=rng)[1], lambda d, rng: scarp_layout(d, "multiring", rng=rng)[0]):
            spread = []
            for d in (D0, 3 * D0):
                rng = np.random.default_rng(1)
                radii = np.array([layout(d, rng) for _ in range(2000)])
                spread.append(np.std(radii - np.mean(radii)))
            self.assertLess(abs(spread[1] / spread[0] - 1.0), 0.1)

    def test_ring_scarp(self):
        # Only multiring basins have a ring scarp; its height grows in from 0 as they grow from 500 to 650 km
        from cratermaker.components.morphology.basicmoon import _SCARP, scarp_layout

        self.assertEqual(scarp_layout(400e3, "complex"), (0.0, 0.0, 0.0))
        self.assertEqual(scarp_layout(500e3, "multiring")[1], 0.0)
        heights = [scarp_layout(d, "multiring")[1] / d for d in np.geomspace(500e3, 700e3, 50)]
        self.assertTrue(np.all(np.diff(heights) >= 0.0))
        self.assertAlmostEqual(heights[-1], _SCARP["height_over_D"])

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            morphology = Morphology.maker(morphology="basicmoon", surface=surface, rng_seed=7)
            crater = morphology.Crater.maker(diameter=1000e3, location=(0.0, 0.0))
            self.assertTrue(crater.has_scarp)
            self.assertGreater(crater.peak_ring_radius, 0.0)
            # The scarp is a ridge added to the profile: its crest is at scarp_radius, scarp_height above the profile without it
            plain = morphology.Crater.maker(crater=crater, scarp_radius=0.0, scarp_height=0.0, scarp_width=0.0)
            self.assertFalse(plain.has_scarp)
            r = np.linspace(0.0, 1.5 * crater.radius, 3001)
            ridge = morphology.crater_profile(crater, r) - morphology.crater_profile(plain, r)
            self.assertAlmostEqual(r[np.argmax(ridge)] / crater.scarp_radius, 1.0, places=2)
            self.assertAlmostEqual(ridge.max() / crater.scarp_height, 1.0, places=3)
            # The steep inner face is narrower than the back slope
            inner = r[(ridge > 0.5 * crater.scarp_height) & (r < crater.scarp_radius)]
            outer = r[(ridge > 0.5 * crater.scarp_height) & (r > crater.scarp_radius)]
            self.assertLess(crater.scarp_radius - inner.min(), outer.max() - crater.scarp_radius)
            self.assertEqual(ridge[r > crater.scarp_radius + 4 * crater.scarp_width].max(), 0.0)

    def test_legacy_ring_rows_are_dropped(self):
        # Files written when the rings of multiring basins were separate ring craters still load: their ring rows are skipped
        from cratermaker.components.crater import _drop_legacy_rings

        self.assertIsNone(_drop_legacy_rings({"id": 2, "morphology_type": "ring", "parent": 1, "isring": True}))
        self.assertEqual(_drop_legacy_rings({"id": 1, "morphology_type": "multiring", "isring": False}),
                         {"id": 1, "morphology_type": "multiring"})

    def test_morphology_type_classified_by_size(self):
        # Without an explicit morphology_type, craters are classified from their size
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            for model_name in morphology_models:
                morphology = Morphology.maker(morphology=model_name, surface=surface, rng_seed=7)
                # Only multiring basins have a ring scarp (the Outer Rook ring of Orientale); the terraces of complex craters
                # are variation about the axisymmetric profile
                expected = {1e3: ("simple", False), 60e3: ("complex", False), 600e3: ("multiring", True)}
                for diameter, (morphology_type, has_scarp) in expected.items():
                    crater = morphology.Crater.maker(diameter=diameter, location=(0.0, 0.0))
                    self.assertEqual(crater.morphology_type, morphology_type, f"{model_name} D={diameter}")
                    self.assertEqual(crater.has_scarp, has_scarp, f"{model_name} D={diameter}")
                explicit = morphology.Crater.maker(diameter=15e3, location=(0.0, 0.0), morphology_type="complex")
                self.assertEqual(explicit.morphology_type, "complex")

    def test_realmoon_ejecta_never_negative(self):
        # The field is shared between the uplift and the ejecta in proportion to their parts of the profile's height, so the
        # ejecta part is never negative (and never clipped when the crater forms), and the crater and ejecta parts sum to the total
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            morphology = Morphology.maker(morphology="realmoon", surface=surface, rng_seed=11, hfield_pix=1e3)
            for diameter in (60e3, 250e3):
                crater = morphology.Crater.maker(diameter=diameter, location=(0.0, 0.0))
                r = np.repeat(np.linspace(0.0, 2.0 * crater.radius, 400), 72)
                b = np.tile(np.arange(72) * 5.0, 400)
                total = morphology.crater_profile(crater, r, bearings=b, include_crater=True, include_ejecta=True)
                crater_part = morphology.crater_profile(crater, r, bearings=b, include_crater=True, include_ejecta=False)
                ejecta_part = morphology.crater_profile(crater, r, bearings=b, include_crater=False, include_ejecta=True)
                self.assertGreaterEqual(ejecta_part.min(), 0.0)
                # The ejecta takes at most the whole field, so it is never thicker than the crater is high
                self.assertLessEqual(ejecta_part.max(), total.max())
                np.testing.assert_allclose(crater_part + ejecta_part, total, atol=1e-6)

    def test_realmoon_wall_levels(self):
        # Complex craters get stacked wall levels whose outermost one is the rim outline; simple craters get none
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            morphology = Morphology.maker(morphology="realmoon", surface=surface, rng_seed=11, hfield_pix=500.0)
            crater = morphology.Crater.maker(diameter=100e3, location=(0.0, 0.0))
            self.assertEqual(crater.morphology_type, "complex")
            nominal, levels, params = crater.wall
            self.assertGreaterEqual(len(nominal), 2)
            self.assertEqual(nominal[0], 1.0)
            self.assertTrue(np.all(np.diff(nominal) < 0))
            self.assertGreaterEqual(nominal[-1], crater.floor_radius / crater.radius)
            rim = crater.hfield[0]
            self.assertEqual(levels.shape, (3, len(nominal), len(rim)))
            np.testing.assert_allclose(rim, levels[0].max(axis=0))
            self.assertTrue(np.all((levels[1] >= 0) & (levels[1] <= 1)))
            r = np.linspace(0, 2 * crater.radius, 200)
            self.assertTrue(np.all(np.isfinite(morphology.crater_profile(crater, r, bearings=np.full(200, 30.0)))))
            simple = morphology.Crater.maker(diameter=2e3, location=(0.0, 0.0))
            self.assertEqual(len(simple.wall[0]), 0)

    def test_realmoon_copy_keeps_field(self):
        # Copying a crater (with the constructor or with maker) keeps the field it has already drawn, without drawing any random
        # numbers; a copy with a different size is a different crater, which draws its own
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(simdir=simdir, target=self.target, reset=True, gridlevel=self.gridlevel, ask_overwrite=False)
            morphology = Morphology.maker(morphology="realmoon", surface=surface, rng_seed=7)
            crater = morphology.Crater.maker(diameter=600e3, location=(10.0, 5.0))
            rim, field, _ = crater.hfield
            rng_state = morphology.rng.bit_generator.state
            copies = [type(crater)(crater=crater, morphology=morphology), morphology.Crater.maker(crater=crater)]
            for c in copies:
                self.assertEqual(c.id, crater.id)
                self.assertEqual(c.morphology_type, crater.morphology_type)
                self.assertEqual(c.rim_height, crater.rim_height)
                np.testing.assert_array_equal(c.hfield[0], rim)
                np.testing.assert_array_equal(c.hfield[1], field)
            self.assertEqual(morphology.rng.bit_generator.state, rng_state, "copying drew random numbers")
            resized = morphology.Crater.maker(crater=crater, diameter=5e3)
            self.assertEqual(resized.morphology_type, "simple")
            self.assertIsNone(getattr(resized._var, "_hfield", None))

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
                    # The rim height can only be measured where the surface contains the rim and its mesh resolves it (basins
                    # have narrow rims)
                    if crater.rim_width >= 2 * sim.surface.pix and getattr(sim.surface, "local_radius", np.inf) > crater.radius:
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


class TestRealMoon(unittest.TestCase):
    """The realistic lunar morphology: the basic crater warped by a rim outline, plus a field of variation (the H-field)."""

    def test_ejecta_balances_crater_without_profile_peak(self):
        """
        The ejecta of a crater balances its crater part. A realmoon peak takes its volume from the floor, so the balance is with
        the profile without the basic Gaussian peak; a basicmoon peak is part of the profile.
        """
        from types import SimpleNamespace

        from cratermaker.components.morphology.basicmoon import _PROFILE_PARAMETERS

        for name in ("basicmoon", "realmoon"):
            morphology = Morphology.maker(name, rng_seed=11)
            for diameter in (12e3, 40e3, 90e3, 200e3):
                crater = morphology.Crater.maker(final_diameter=diameter, morphology=morphology, monte_carlo=False)
                params = SimpleNamespace(**{k: getattr(crater, k) for k in _PROFILE_PARAMETERS})
                if name == "realmoon":
                    params.peak_height = 0.0
                excavated = morphology.estimate_volume(params, include_crater=True, include_ejecta=False)
                deposit = crater.frac_ejrim * 2 * np.pi * crater.rim_height * crater.radius**2 / (-crater.ejprofile - 2)
                with self.subTest(morphology=name, diameter=diameter):
                    self.assertLess(abs(deposit + excavated) / -excavated, 1e-6)

    def test_volume_conserved_and_field_kept_on_copy(self):
        import cratermaker

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(
                simdir=simdir,
                surface="hireslocal",
                morphology="realmoon",
                local_location=(0, 0),
                local_radius=40e3,
                pix=400.0,
                ask_overwrite=False,
                reset=True,
                rng_seed=5,
                save_observed=False,
                save_emplaced=False,
            )
            area = np.asarray(sim.surface.face_area)
            before = np.array(sim.surface.face_elevation)
            crater = sim.Crater.maker(diameter=25e3, location=(0, 0))
            rim, field, rho_max = crater.hfield
            copy = type(crater)(crater=crater)
            np.testing.assert_array_equal(copy.hfield[1], field)
            self.assertEqual(rho_max, 2.0)
            self.assertTrue(np.all(rim > 0))
            self.assertGreater(np.std(field), 0.0)
            sim.emplace(crater)
            change = np.asarray(sim.surface.face_elevation) - before
            self.assertGreater(np.sum(np.abs(change) * area), 0.0)
            self.assertLess(abs(np.sum(change * area)) / np.sum(np.abs(change) * area), 1e-10)

    def test_field_resolution_follows_the_mesh(self):
        from cratermaker.components.morphology import Morphology

        fine = Morphology.maker("realmoon", hfield_pix=20.0, rng_seed=1)
        coarse = Morphology.maker("realmoon", hfield_pix=2000.0, rng_seed=1)
        f = fine.Crater.maker(diameter=20e3, location=(0, 0)).hfield[1]
        c = coarse.Crater.maker(diameter=20e3, location=(0, 0)).hfield[1]
        self.assertGreater(f.shape[0], c.shape[0])
        self.assertGreater(f.shape[1], c.shape[1])

    def test_peak_center(self):
        """A realmoon peak's center is offset at random unless given; a given offset is kept and the peak is built around it."""
        import cratermaker

        D = 60e3
        changes = {}
        for offset in (True, False):
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                sim = cratermaker.Simulation(
                    simdir=simdir,
                    surface="hireslocal",
                    morphology="realmoon",
                    local_location=(0, 0),
                    local_radius=1.5 * D,
                    pix=D / 200,
                    ask_overwrite=False,
                    reset=True,
                    rng_seed=3,
                    save_observed=False,
                    save_emplaced=False,
                )
                random = sim.Crater.maker(diameter=D, location=(0, 0), morphology_type="complex")
                self.assertGreater(random.peak_center_distance, 0.0)
                distance = 0.3 * random.radius
                crater = sim.Crater.maker(
                    diameter=D, location=(0, 0), morphology_type="complex", peak_center_distance=distance, peak_center_bearing=90.0
                )
                self.assertEqual(crater.peak_center_distance, distance)
                self.assertEqual(crater.peak_center_bearing, 90.0)
                copy = sim.Crater.maker(crater=crater)
                self.assertEqual((copy.peak_center_distance, copy.peak_center_bearing), (distance, 90.0))
                if not offset:
                    crater = sim.Crater.maker(crater=crater, peak_height=0.0)
                region = sim.surface.extract_region((0, 0), 0.75 * D)
                sim.emplace(crater)
                b = np.radians(np.asarray(region.face_bearing))
                r = np.asarray(region.face_distance)
                changes[offset] = (r * np.sin(b), r * np.cos(b), np.asarray(region.face_area), np.array(region.face_elevation))
        x, y, area, with_peak = changes[True]
        peak = np.clip(with_peak - changes[False][3], 0.0, None) * area
        cx, cy = np.sum(peak * x) / np.sum(peak), np.sum(peak * y) / np.sum(peak)
        # The peak is built around the given center, 0.3 R to the east
        self.assertGreater(cx, 0.5 * distance)
        self.assertLess(cx, 1.5 * distance)
        self.assertLess(abs(cy), 0.5 * distance)

    def _peak_sim(self, simdir, D, rng_seed):
        import cratermaker

        return cratermaker.Simulation(
            simdir=simdir,
            surface="hireslocal",
            morphology="realmoon",
            local_location=(0, 0),
            local_radius=1.5 * D,
            pix=D / 200,
            ask_overwrite=False,
            reset=True,
            rng_seed=rng_seed,
            save_observed=False,
            save_emplaced=False,
        )

    def test_peak_elongation(self):
        """A given elongation stretches the peak along the given bearing."""
        D = 60e3
        fixed = dict(
            diameter=D, location=(0, 0), morphology_type="complex", rim_height=1000.0, floor_elevation=-3000.0, floor_radius=18e3,
            rim_width=3e3, wall_curvature=0.2, floor_blend=0.5, peak_width=3240.0, peak_center_distance=0.0,
            peak_center_bearing=0.0, peak_elongation=3.0, peak_seed=11, hfield_seed=12, ray_seed=13,
        )
        moments = {}
        for bearing in (0.0, 90.0):
            change = {}
            for height in (1000.0, 0.0):
                with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                    sim = self._peak_sim(simdir, D, 5)
                    crater = sim.Crater.maker(peak_height=height, peak_elongation_bearing=bearing, **fixed)
                    self.assertEqual((crater.peak_elongation, crater.peak_elongation_bearing), (3.0, bearing))
                    region = sim.surface.extract_region((0, 0), 0.75 * D)
                    before = np.array(region.face_elevation)
                    sim.emplace(crater)
                    change[height] = np.array(region.face_elevation) - before
            b = np.radians(np.asarray(region.face_bearing))
            r = np.asarray(region.face_distance)
            x, y = r * np.sin(b), r * np.cos(b)
            # The peak, on the floor (removing it also changes the volume balance of the profile, and so the ejecta)
            peak = np.where(r < 18e3, np.clip(change[1000.0] - change[0.0], 0.0, None), 0.0) * np.asarray(region.face_area)
            moments[bearing] = (np.sum(peak * x**2), np.sum(peak * y**2))
        # Long axis north-south: more spread in y than in x; east-west: the reverse
        self.assertGreater(moments[0.0][1], 1.5 * moments[0.0][0])
        self.assertGreater(moments[90.0][0], 1.5 * moments[90.0][1])

    def test_crater_reproducible_from_parameters(self):
        """A realmoon crater given all its parameters and seeds is the same in simulations with different random seeds."""
        from cratermaker.components.morphology.basicmoon import _PROFILE_PARAMETERS

        D = 60e3
        names = [n for n in _PROFILE_PARAMETERS if n not in ("diameter", "radius")] + ["hfield_seed", "ray_seed"]
        changes = []
        given = None
        for rng_seed in (21, 22):
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                sim = self._peak_sim(simdir, D, rng_seed)
                if given is None:
                    first = sim.Crater.maker(diameter=D, location=(0, 0), morphology_type="complex", ray_seed=31)
                    given = {n: getattr(first, n) for n in names}
                crater = sim.Crater.maker(diameter=D, location=(0, 0), morphology_type="complex", **given)
                region = sim.surface.extract_region((0, 0), 0.75 * D)
                before = np.array(region.face_elevation)
                sim.emplace(crater)
                changes.append(np.array(region.face_elevation) - before)
        np.testing.assert_allclose(changes[0], changes[1], atol=1e-6)

    def test_central_peak(self):
        """A complex crater gets a central peak of its height that changes nothing beyond the floor and conserves volume."""
        import cratermaker

        D = 60e3
        surfaces = {}
        for peak in (True, False):
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                sim = cratermaker.Simulation(
                    simdir=simdir,
                    surface="hireslocal",
                    morphology="realmoon",
                    local_location=(0, 0),
                    local_radius=1.5 * D,
                    pix=D / 200,
                    ask_overwrite=False,
                    reset=True,
                    rng_seed=7,
                    save_observed=False,
                    save_emplaced=False,
                )
                crater = sim.Crater.maker(diameter=D, location=(0, 0), morphology_type="complex")
                self.assertGreater(crater.peak_height, 0.0)
                from cratermaker.components.morphology.realmoon import _PEAK_WIDTH
                self.assertAlmostEqual(crater.peak_width, _PEAK_WIDTH * crater.radius)
                if not peak:
                    crater = sim.Crater.maker(crater=crater, peak_height=0.0)
                    self.assertEqual(crater.peak_height, 0.0)
                region = sim.surface.extract_region((0, 0), 0.75 * D)
                sim.emplace(crater)
                surfaces[peak] = (
                    np.asarray(region.face_distance),
                    np.asarray(region.face_area),
                    np.array(region.face_elevation),
                    crater,
                )
        r, area, with_peak, crater = surfaces[True]
        without = surfaces[False][2]
        change = with_peak - without
        # The relaxed peak reaches its height above the floor around it, within the tolerance of the amplitude iteration. The
        # floor is lowered by the layer that holds the peak's volume, thickest at the center, so the height is measured from the
        # lowest point of the floor near the peak (where there are no columns, the change is that of the layer alone)
        near = r < 2.5 * crater.peak_width
        self.assertAlmostEqual((change.max() - change[near].min()) / crater.peak_height, 1.0, delta=0.05)
        # Beyond the floor only the faint tail of the (randomly offset and elongated) envelope remains, and the peak's volume is
        # taken from the floor
        outside = np.abs(change[r > crater.floor_radius * 1.05])
        print(f"largest change beyond the floor: {outside.max():.3g} m of a {crater.peak_height:.0f} m peak")
        self.assertLess(outside.max(), 1e-2 * crater.peak_height)
        self.assertLess(abs(np.sum(change * area)) / np.sum(np.abs(change) * area), 1e-8)


if __name__ == "__main__":
    unittest.main()
