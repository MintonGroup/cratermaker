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
            with xr.open_dataset(filename) as ds:
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
                simdir=simdir, gridlevel=self.gridlevel, ask_overwrite=False, save_actions=None, do_counting=False
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
                    self.assertTrue(
                        all(
                            crater_dict[k] != new_crater_dict[k]
                            if k == "time" or k == "id"
                            else crater_dict[k] == new_crater_dict[k]
                            for k in crater_dict
                        ),
                        msg=f"Failed for sim_num {sim_num} iteration {iteration}",
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
