import tempfile
import unittest

# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

import cratermaker
from cratermaker import Crater, Target

warnings.filterwarnings("ignore", category=FutureWarning, module="xarray")
warnings.filterwarnings("ignore", category=FutureWarning, module="uxarray")


def warning_with_breakpoint(message, category, filename, lineno, file=None, line=None):
    print(f"{category} Warning in line {lineno} of {filename} : {message}")
    return


warnings.simplefilter("always")  # Always trigger the warnings
warnings.showwarning = warning_with_breakpoint


class TestSimulation(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        target = Target(name="Moon")
        self.pix = target.radius / 10.0
        self.gridlevel = 5

    def test_simulation_defaults(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False)
            self.assertEqual(sim.target.name, "Moon")

    def test_simulation_save(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Test basic save operation
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False)
            sim.save()

            # Test that variables are saved correctly
            sim.surface.update_elevation(1.0)
            np.testing.assert_array_equal(
                sim.surface.uxds["node_elevation"].values,
                np.ones(sim.surface.uxds.uxgrid.n_node),
            )
            np.testing.assert_array_equal(
                sim.surface.uxds["face_elevation"].values,
                np.ones(sim.surface.uxds.uxgrid.n_face),
            )

            sim.save()
            filename = (
                Path(sim.surface.output_dir)
                / f"{sim.surface._output_file_prefix}{sim.interval:06d}.{sim.surface._output_file_extension}"
            )
            self.assertTrue(filename.exists())
            with xr.open_dataset(filename, engine="h5netcdf") as ds:
                ds = ds.isel(interval=-1)
                np.testing.assert_array_equal(ds["node_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_node))
                np.testing.assert_array_equal(ds["face_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_face))

        return

    def test_emplace(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            cdiam = 2 * self.pix
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, save_actions=None)
            crater = cratermaker.Crater.maker(diameter=cdiam)
            sim.emplace(crater)
            pdiam = crater.projectile_diameter

            sim.emplace(diameter=cdiam)
            sim.emplace(projectile_diameter=pdiam)
        return

    def test_populate(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, save_actions=None)
            # Test that populate will work even if no craters are returned
            sim.populate(age=1e-6)
            sim.populate(age=3.8e3)
        return

    def test_run(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            common_args = {
                "simdir": simdir,
                "gridlevel": self.gridlevel,
                "reset": True,
                "ask_overwrite": False,
                "save_actions": None,
            }
            sim = cratermaker.Simulation(do_subpixel_degradation=True, **common_args)
            sim.run(age=1000)

            sim = cratermaker.Simulation(**common_args)
            sim.run(age=1000, time_interval=100)

            # Test that the simulation doesn't fail when the age doesn't divide evenly by the time_interval
            sim = cratermaker.Simulation(**common_args)
            sim.run(age=1010, time_interval=100)

            # Test runs with a single time interval between time_start and time_end
            sim = cratermaker.Simulation(**common_args)
            sim.run(time_start=1010, time_end=1000, time_interval=10)

            # Tests if a run can continue where a previous run left off
            sim.run(time_start=1000, time_end=900, time_interval=10)

            # Tests that runs can continue without passing time_start explicitly
            sim.run(time_end=880, time_interval=10)
        return

    def test_invalid_run_args(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, save_actions=None)

            # Test case: Neither the age nor the diameter_number argument is provided
            with self.assertRaises(ValueError):
                sim.run()

            # Test case: Both the age and diameter_number arguments are provided
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, diameter_number=(300e3, 80))

            # Test case: Both the time_end and diameter_number_end arguments are provided
            with self.assertRaises(ValueError):
                sim.run(time_end=3.0e3, diameter_number_end=(300e3, 80))

            # Test case: The age argument is provided but is not a scalar
            with self.assertRaises(ValueError):
                sim.run(age=[3.8e3])

            # Test case: The time_end argument is provided but is not a scalar
            with self.assertRaises(ValueError):
                sim.run(time_end=[3.0e3])

            # Test case: The time_interval is provided but is not a positive scalar
            with self.assertRaises(ValueError):
                sim.run(time_start=3.8e3, time_interval=-100.0)

            # Test case: The time_interval provided is negative, or is greater than age - time_end
            with self.assertRaises(ValueError):
                sim.run(time_start=3.8e3, time_end=3.0e3, time_interval=1000.0)

            # Test case: The diameter_number argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number=(300e3, -80))

            # Test case: The diameter_number_end argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number_end=(300e3, -80))

            # Test case: The diameter_number_interval argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number_interval=(300e3, -80))

            # Test case: The ninterval is provided but is not an integer or is less than 1
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, ninterval=0)

            # Test case: The ninterval is provided and  time_interval is also provided
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, ninterval=100, time_interval=100.0)

        return

    def test_simulation_to_config(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # First simulation: no target passed, should default to "Moon"
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel, reset=True, ask_overwrite=False)
            self.assertIsInstance(sim.target, Target)
            self.assertEqual(sim.target.name, "Moon")
            del sim

            # Second simulation: override target with "Mars"
            sim = cratermaker.Simulation(simdir=simdir, target="Mars", reset=False, ask_overwrite=False)
            self.assertEqual(sim.target.name, "Mars")
            self.assertEqual(sim.production.version, "Mars")
            del sim

            # Third simulation: no target passed, should read "Mars" from config
            sim = cratermaker.Simulation(simdir=simdir, reset=False, ask_overwrite=False)
            self.assertEqual(sim.target.name, "Mars")
            del sim
        return

    def test_save_actions(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            common_args = {"simdir": simdir, "gridlevel": self.gridlevel, "reset": True, "ask_overwrite": False}

            # Test default behavior
            sim = cratermaker.Simulation(**common_args)
            default_sim_save_actions = sim.save_actions.copy()
            self.assertEqual(len(default_sim_save_actions), 1)
            sim = cratermaker.Simulation(save_actions="default", **common_args)
            self.assertEqual(sim.save_actions, default_sim_save_actions)
            sim = cratermaker.Simulation(save_actions=None, **common_args)
            self.assertEqual(len(sim.save_actions), 0)

            new_sim_save_action = {
                "plot": {
                    "include_counting": True,
                    "show": False,
                    "save": True,
                }
            }
            new_counting_save_action = {"export": {"driver": "SCC", "crater_type": "both"}}
            ninterval = 4

            old_sim_action_len = len(sim.save_actions)
            old_count_action_len = len(sim.counting.save_actions)
            sim.add_save_action(new_sim_save_action)
            sim.counting.add_save_action(new_counting_save_action)
            self.assertEqual(len(sim.save_actions), old_sim_action_len + 1)
            self.assertEqual(len(sim.counting.save_actions), old_count_action_len + 1)
            sim.run(age=4000, ninterval=ninterval)
            img_out = list(sim.counting.plot_dir.glob(f"*.{sim.counting._output_image_file_extension}"))
            self.assertEqual(len(img_out), ninterval + 1)
            scc_out = list(sim.counting.export_dir.glob("*.scc"))
            self.assertEqual(len(scc_out), 2 * (ninterval + 1))
        return

    def test_volume_conservation(self):
        # Craters and their ejecta, and the degradation they cause, only move material around, so a simulation should not change
        # the total volume of the surface beyond rounding error. Craters that the mesh resolves so poorly that their sampled shape
        # removes no material are skipped, as their ejecta could not balance them. Multiring basins are not held to volume
        # conservation (their final surface keeps little of their excavation), so craters are limited to sizes below the
        # multiring transition. (Modeled on example 1.2.)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(
                target="Moon",
                morphology="basicmoon",
                largest_crater=400e3,
                gridlevel=self.gridlevel,
                save_observed=False,
                save_emplaced=False,
                simdir=simdir,
                ask_overwrite=False,
                reset=True,
                rng_seed=3234653,
            )
            area = np.asarray(sim.surface.face_area)
            initial_elevation = np.array(sim.surface.face_elevation, dtype=float)
            sim.run(age=4310)
            change = np.asarray(sim.surface.face_elevation) - initial_elevation
            net = np.sum(change * area)
            gross = np.sum(np.abs(change) * area)
            self.assertGreater(gross, 0.0)
            self.assertLess(abs(net) / gross, 1e-12, f"net volume change {net:.3e} m3 out of {gross:.3e} m3 moved")

    def test_count_superdomain(self):
        # On a hireslocal surface, a crater centered outside the local region is only tagged for counting when count_superdomain
        # is set, either for the simulation (passed to the Counting model) or for one emplace call, and the setting is kept when
        # the simulation is resumed
        def offside_crater(sim, **kwargs):
            local_radius = sim.surface.local.radius
            sim.emplace(radius=800e3, relative_location={"distance": local_radius + 200e3, "bearing": 90.0}, **kwargs)
            return len(sim.counting.observed)

        settings = {"pix": 10e3, "local_radius": 250e3, "ask_overwrite": False, "rng_seed": 234471661}
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(surface="hireslocal", simdir=simdir, reset=True, **settings)
            self.assertFalse(sim.counting.count_superdomain)
            self.assertEqual(offside_crater(sim), 0)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(surface="hireslocal", simdir=simdir, reset=True, **settings)
            self.assertEqual(offside_crater(sim, count_superdomain=True), 1)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(surface="hireslocal", simdir=simdir, reset=True, count_superdomain=True, **settings)
            self.assertTrue(sim.counting.count_superdomain)
            self.assertEqual(offside_crater(sim), 1)
            sim.save()
            resumed = cratermaker.Simulation(simdir=simdir, reset=False, ask_overwrite=False)
            self.assertTrue(resumed.counting.count_superdomain)

    def test_resume_reproduces_run(self):
        # A run resumed from its saved state continues exactly as an uninterrupted run would: the state of the random number
        # generator saved with the simulation is restored, not replaced by one freshly seeded (and advanced by the setup)
        kw = dict(morphology="basicmoon", rng_seed=230982307, gridlevel=4, ask_overwrite=False, save_actions=[], do_counting=True)

        def interval_craters(simdir, interval):
            with xr.open_dataset(Path(simdir) / "counting" / f"emplaced_craters{interval:06d}.nc") as ds:
                return np.sort(ds["diameter"].values.ravel())

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as continuous, tempfile.TemporaryDirectory(
            ignore_cleanup_errors=True
        ) as resumed:
            sim = cratermaker.Simulation(simdir=continuous, reset=True, **kw)
            sim.run(time_start=4310, time_end=4290, time_interval=10)
            final_state = sim.rng.bit_generator.state

            sim = cratermaker.Simulation(simdir=resumed, reset=True, **kw)
            sim.run(time_start=4310, time_end=4300, time_interval=10)
            saved_state = sim.rng.bit_generator.state
            del sim
            sim = cratermaker.Simulation(simdir=resumed, reset=False, **kw)
            self.assertEqual(sim.rng.bit_generator.state, saved_state)
            sim.run(time_start=4310, time_end=4290, time_interval=10)
            self.assertEqual(sim.rng.bit_generator.state, final_state)
            np.testing.assert_array_equal(interval_craters(continuous, 2), interval_craters(resumed, 2))

    def test_crater_size_limits(self):
        # The smallest and largest crater sizes given to the constructor are applied, limited to the surface and the target, and
        # kept when the simulation is resumed (unless new ones are given)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(
                target="Moon",
                gridlevel=self.gridlevel,
                simdir=simdir,
                ask_overwrite=False,
                reset=True,
                smallest_crater=100e3,
                largest_crater=400e3,
            )
            self.assertEqual(sim.smallest_crater, 100e3)
            self.assertEqual(sim.largest_crater, 400e3)
            sim.save()
            resumed = cratermaker.Simulation(simdir=simdir, ask_overwrite=False, reset=False)
            self.assertEqual(resumed.smallest_crater, 100e3)
            self.assertEqual(resumed.largest_crater, 400e3)
            overridden = cratermaker.Simulation(simdir=simdir, ask_overwrite=False, reset=False, largest_crater=300e3)
            self.assertEqual(overridden.largest_crater, 300e3)
            capped = cratermaker.Simulation(
                target="Moon", gridlevel=self.gridlevel, simdir=simdir, ask_overwrite=False, reset=True, largest_crater=1e9
            )
            self.assertEqual(capped.largest_crater, capped.target.diameter)

    def test_quasimc_reproducible(self):
        # The quasi-Monte Carlo craters must be created with the simulation's random number generator, so that two simulations
        # with the same seed produce the same craters (including the randomly drawn projectile properties).
        attrs = ("id", "time", "diameter", "transient_diameter", "projectile_velocity", "projectile_angle", "orientation")
        results = []
        for _ in range(2):
            with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
                quasimc_file = Path(simdir) / "qmc.csv"
                with Path.open(quasimc_file, "w") as f:
                    f.write(
                        "name,latitude,longitude,diameter,production_time,production_time_stdev,production_D,production_N,production_N_stdev,production_sequence\n"
                    )
                    f.write("South Pole-Aitken,-53,191,2400000,,,20,999,,0\n")
                    f.write("Nectaris,-15.6,35.1,885000,,,20,172,20,40\n")
                    f.write("Imbrium,37,341.5,1321000,3922,12,,,,100\n")
                sim = cratermaker.Simulation(
                    simdir=simdir,
                    gridlevel=self.gridlevel,
                    reset=True,
                    ask_overwrite=False,
                    quasimc_file=quasimc_file,
                    rng_seed=8675309,
                )
                from_file = [tuple(getattr(c, a) for a in attrs) for c in sim.quasimc_craters]
                merged = sim.populate(time_start=4500, time_end=3800, emplace=False)
                named = [tuple(getattr(c, a) for a in attrs) for c in merged if c.name is not None]
                results.append((from_file, named))
        self.assertEqual(len(results[0][1]), 3)
        self.assertEqual(results[0], results[1])

    def test_quasimc(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            quasimc_file = Path(simdir) / "qmc.csv"
            with Path.open(quasimc_file, "w") as f:
                f.write(
                    "name,latitude,longitude,diameter,production_time,production_time_stdev,production_D,production_N,production_N_stdev,production_sequence\n"
                )
                f.write("South Pole-Aitken,-53,191,2400000,,,20,999,,0\n")
                f.write("Vaporum,14.2,3.1,410000,,,,,,5\n")
                f.write("Fecunditatis,-4.6,52,690000,,,90,10,4,10\n")
                f.write("Nectaris,-15.6,35.1,885000,,,20,172,20,40\n")
                f.write("Imbrium,37,341.5,1321000,3922,12,,,,100\n")

            sim = cratermaker.Simulation(
                simdir=simdir, gridlevel=self.gridlevel, reset=True, ask_overwrite=False, quasimc_file=quasimc_file
            )
            self.assertEqual(len(sim.quasimc_craters), 5)
            self.assertTrue(all(isinstance(crater, Crater) for crater in sim.quasimc_craters))
            self.assertTrue(
                all(
                    crater.name in ["South Pole-Aitken", "Vaporum", "Fecunditatis", "Nectaris", "Imbrium"]
                    for crater in sim.quasimc_craters
                )
            )
            self.assertTrue(all(crater.time is not None for crater in sim.quasimc_craters))
            min_diameter = min([crater.diameter for crater in sim.quasimc_craters])
            self.assertEqual(sim.largest_crater, min_diameter)

            sim.quasimc_craters.append(
                sim.Crater.maker(
                    name="Copernicus",
                    location=(339.9214, 9.6209),
                    diameter=96070,
                    production_time=(800.0, 15.0),
                    production_sequence=200,
                )
            )

            self.assertEqual(len(sim.quasimc_craters), 6)
            self.assertTrue(all(isinstance(crater, Crater) for crater in sim.quasimc_craters))
            self.assertTrue(
                all(
                    crater.name in ["South Pole-Aitken", "Vaporum", "Fecunditatis", "Nectaris", "Imbrium", "Copernicus"]
                    for crater in sim.quasimc_craters
                )
            )
            self.assertTrue(all(crater.time is not None for crater in sim.quasimc_craters))
            min_diameter = min([crater.diameter for crater in sim.quasimc_craters])
            self.assertEqual(sim.largest_crater, min_diameter)

            with Path.open(quasimc_file, "w") as f:
                f.write(
                    "name,latitude,longitude,diameter,production_time,production_time_stdev,production_D,production_N,production_N_stdev,production_sequence\n"
                )
                f.write("South Pole-Aitken,-53,191,2400000,,,20,999,,0\n")
                f.write("Vaporum,14.2,3.1,410000,,,,,,5\n")
                f.write("Fecunditatis,-4.6,52,690000,,,90,10,4,10\n")
                f.write("Nectaris,-15.6,35.1,885000,,,20,172,20,40\n")
                f.write("Imbrium,37,341.5,1321000,3922,12,,,,100\n")
                f.write("Tycho,-43.2958,348.7847,85294,109,4,,,,300\n")

            sim.quasimc_file = quasimc_file
            self.assertEqual(len(sim.quasimc_craters), 6)
            self.assertTrue(
                all(
                    crater.name in ["South Pole-Aitken", "Vaporum", "Fecunditatis", "Nectaris", "Imbrium", "Tycho"]
                    for crater in sim.quasimc_craters
                )
            )
            self.assertTrue(all(crater.time is not None for crater in sim.quasimc_craters))
            min_diameter = min([crater.diameter for crater in sim.quasimc_craters])
            self.assertEqual(sim.largest_crater, min_diameter)

            sim.largest_crater = 200e3
            self.assertEqual(sim.largest_crater, 200e3)
            self.assertEqual(sim.production.diameter_range[1], 200e3)

        return

    def test_positive_only_ejecta(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(
                simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, save_actions=None, save_observed=False
            )
            sim.run(age=4310)
            self.assertFalse(any(sim.surface.ejecta_thickness < 0))

        return

    def test_repeatable_crater_copy(self):
        """
        Tests two things: Copied+modified craters retain their old non-modified values and repeated craters called with the same simulation rng_seed are the same.
        """
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            crater_vals = []
            diffs = ["location", "projectile_velocity", "orientation", "id"]
            sames = ["diameter", "time"]
            ignores = ["ejecta_rmax", "emplaceable"]
            for sim_num in range(2):
                sim = cratermaker.Simulation(
                    simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, reset=True, rng_seed=2345253
                )

                for iteration in range(2):
                    crater = sim.Crater.maker(diameter=5000.0, time=1000.0)
                    new_crater = sim.Crater.maker(crater=crater, time=2000.0)
                    crater_dict = crater.as_dict()
                    crater_vals.append(crater_dict)
                    new_crater_dict = new_crater.as_dict()
                    test_dict = {
                        k: crater_dict[k] != new_crater_dict[k]
                        if k == "time" or k == "id"
                        else crater_dict[k] == new_crater_dict[k]
                        for k in crater_dict
                        if k not in ignores
                    }
                    mismatch = {k: v for k, v in test_dict.items() if not v}
                    self.assertTrue(
                        all(test_dict.values()),
                        msg=f"Failed for sim_num {sim_num} iteration {iteration}. Mismatched keys: {mismatch.keys()}",
                    )

                self.assertTrue(
                    all(crater_vals[-2][k] == crater_vals[-1][k] for k in sames),
                    msg=f"Failed same values check on sim_num {sim_num}",
                )
                self.assertTrue(
                    all(crater_vals[-2][k] != crater_vals[-1][k] for k in diffs),
                    msg=f"Failed different values check on sim_num {sim_num}",
                )

            combo = diffs + sames
            self.assertTrue(
                all(crater_vals[0][k] == crater_vals[2][k] for k in combo),
                msg="Failed repeatability of iteration 0 between the two sims",
            )
            self.assertTrue(
                all(crater_vals[1][k] == crater_vals[3][k] for k in combo),
                msg="Failed repeatability of iteration 1 between the two sims",
            )


if __name__ == "__main__":
    unittest.main()
