import tempfile
import unittest

# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

import cratermaker
from cratermaker import Target
from cratermaker.constants import _COMBINED_DATA_FILE_NAME, _EXPORT_DIR

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
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            self.assertEqual(sim.target.name, "Moon")

    def test_simulation_save(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Test basic save operation
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
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

            filename = Path(sim.data_dir) / _COMBINED_DATA_FILE_NAME.replace(
                ".nc", f"{sim.interval_number:06d}.nc"
            )
            self.assertTrue(filename.exists())
            with xr.open_dataset(filename) as ds:
                ds = ds.isel(time=-1)
                np.testing.assert_array_equal(
                    ds["node_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_node)
                )
                np.testing.assert_array_equal(
                    ds["face_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_face)
                )

            # Test saving combined data
            sim.save(combine_data_files=True)
            filename = _COMBINED_DATA_FILE_NAME

            filename = Path(sim.data_dir) / _COMBINED_DATA_FILE_NAME
            self.assertTrue(filename.exists())
            with xr.open_dataset(filename) as ds:
                ds = ds.isel(time=-1)
                np.testing.assert_array_equal(
                    ds["node_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_node)
                )
                np.testing.assert_array_equal(
                    ds["face_elevation"].values, np.ones(sim.surface.uxds.uxgrid.n_face)
                )

        return

    def test_simulation_export_vtk(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            # Test with default parameters
            default_out_dir = Path(sim.simdir) / _EXPORT_DIR
            expected_files = ["surface000000.vtp"]
            sim.export("vtp")
            self.assertTrue(Path(default_out_dir).is_dir())
            for f in expected_files:
                self.assertTrue((Path(default_out_dir / f).exists()))

    def test_emplace(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            cdiam = 2 * self.pix
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            crater = cratermaker.Crater.maker(final_diameter=cdiam)
            sim.emplace(crater)
            pdiam = crater.projectile_diameter

            sim.emplace(final_diameter=cdiam)
            sim.emplace(projectile_diameter=pdiam)
        return

    def test_populate(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            # Test that populate will work even if no craters are returned
            sim.populate(age=1e-6)
            sim.populate(age=3.8e3)
        return

    def test_run(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            sim.run(age=1000)

            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            sim.run(age=1000, age_interval=100)

            # Test that the simulation doesn't fail when the age doesn't divide evenly by the age_interval
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            sim.run(age=1010, age_interval=100)

    def test_invalid_run_args(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)

            # Test case: Neither the age nor the diameter_number argument is provided
            with self.assertRaises(ValueError):
                sim.run()

            # Test case: Both the age and diameter_number arguments are provided
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, diameter_number=(300e3, 80))

            # Test case: Both the age_end and diameter_number_end arguments are provided
            with self.assertRaises(ValueError):
                sim.run(age_end=3.0e3, diameter_number_end=(300e3, 80))

            # Test case: The age argument is provided but is not a scalar
            with self.assertRaises(ValueError):
                sim.run(age=[3.8e3])

            # Test case: The age_end argument is provided but is not a scalar
            with self.assertRaises(ValueError):
                sim.run(age_end=[3.0e3])

            # Test case: The age_interval is provided but is not a positive scalar
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, age_interval=-100.0)

            # Test case: The age_interval provided is negative, or is greater than age - age_end
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, age_end=3.0e3, age_interval=1000.0)

            # Test case: The diameter_number argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number=(300e3, -80))

            # Test case: The diameter_number_end argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number_end=(300e3, -80))

            # Test case: The diameter_number_interval argument is not a pair of values, or any of them are less than 0
            with self.assertRaises(ValueError):
                sim.run(diameter_number_interval=(300e3, -80))

            # Test case: The age_interval and diameter_number_interval arguments are both provided
            with self.assertRaises(ValueError):
                sim.run(
                    age=3.8e3, age_interval=100.0, diameter_number_interval=(300e3, 80)
                )

            # Test case: The diameter_number_interval provided is negative, or is greater than diameter_number - diameter_number_end
            with self.assertRaises(ValueError):
                sim.run(
                    diameter_number=(300e3, 80),
                    diameter_number_end=(300e3, 30),
                    diameter_number_interval=(300e3, 100),
                )

            # Test case: The ninterval is provided but is not an integer or is less than 1
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, ninterval=0)

            # Test case: The ninterval is provided and either age_interval or diameter_number_interval is also provided
            with self.assertRaises(ValueError):
                sim.run(age=3.8e3, ninterval=100, age_interval=100.0)

        return

    def test_simulation_to_config(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # First simulation: no target passed, should default to "Moon"
            sim = cratermaker.Simulation(simdir=simdir, gridlevel=self.gridlevel)
            self.assertIsInstance(sim.target, Target)
            self.assertEqual(sim.target.name, "Moon")
            del sim

            # Second simulation: override target with "Mars"
            sim = cratermaker.Simulation(simdir=simdir, target="Mars", resume_old=True)
            self.assertEqual(sim.target.name, "Mars")
            del sim

            # Third simulation: no target passed, should read "Mars" from config
            sim = cratermaker.Simulation(simdir=simdir, resume_old=True)
            self.assertEqual(sim.target.name, "Mars")
            del sim
        return


if __name__ == "__main__":
    unittest.main()
