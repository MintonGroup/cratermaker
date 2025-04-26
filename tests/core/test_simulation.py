import unittest
import cratermaker 
from cratermaker import Target
import tempfile
import os
from pathlib import Path
import numpy as np
import xarray as xr
from cratermaker.constants import _COMBINED_DATA_FILE_NAME, _EXPORT_DIR
# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
import warnings
warnings.filterwarnings("ignore",category=FutureWarning,module="xarray")
warnings.filterwarnings("ignore",category=FutureWarning,module="uxarray")

def warning_with_breakpoint(message, category, filename, lineno, file=None, line=None):
    print(f"{category} Warning in line {lineno} of {filename} : {message}")
    return
warnings.simplefilter("always")  # Always trigger the warnings
warnings.showwarning = warning_with_breakpoint

class TestSimulation(unittest.TestCase):
    
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.simdir = self.temp_dir.name
        print(f"Temporary directory created: {self.simdir}")
        target = Target(name="Moon")
        self.pix = target.radius / 10.0
        self.cwd = Path.cwd()
        self.gridlevel = 5
        os.chdir(self.temp_dir.name)
        
    def tearDown(self):
        # Clean up temporary directory
        os.chdir(self.cwd)  # Change back to the original working directory
        self.temp_dir.cleanup() 
        return           

    def test_simulation_defaults(self):
        sim = cratermaker.Simulation(simdir=self.simdir, gridlevel=self.gridlevel)
        self.assertEqual(sim.target.name, "Moon")
        
    def test_simulation_save(self):
        # Test basic save operation
        sim = cratermaker.Simulation(simdir=self.simdir,gridlevel=self.gridlevel)
        sim.save()
       
        # Test that variables are saved correctly
        sim.surf.set_elevation(1.0)
        np.testing.assert_array_equal(sim.surf.data["node_elevation"].values, np.ones(sim.surf.data.uxgrid.n_node)) 
        np.testing.assert_array_equal(sim.surf.data["face_elevation"].values, np.ones(sim.surf.data.uxgrid.n_face)) 
        
        sim.save()
        
        filename = os.path.join(sim.data_dir,_COMBINED_DATA_FILE_NAME.replace(".nc", f"{sim.interval_number:06d}.nc"))
        self.assertTrue(os.path.exists(filename))
        with xr.open_dataset(filename) as ds:
            ds = ds.isel(time=-1)
            np.testing.assert_array_equal(ds["node_elevation"].values, np.ones(sim.surf.data.uxgrid.n_node))
            np.testing.assert_array_equal(ds["face_elevation"].values, np.ones(sim.surf.data.uxgrid.n_face))
    
        # Test saving combined data
        sim.save(combine_data_files=True)
        filename = _COMBINED_DATA_FILE_NAME
        
        filename = os.path.join(sim.data_dir,_COMBINED_DATA_FILE_NAME)
        self.assertTrue(os.path.exists(filename))
        with xr.open_dataset(filename) as ds:
            ds = ds.isel(time=-1)
            np.testing.assert_array_equal(ds["node_elevation"].values, np.ones(sim.surf.data.uxgrid.n_node))
            np.testing.assert_array_equal(ds["face_elevation"].values, np.ones(sim.surf.data.uxgrid.n_face))
    
        return 
        
    def test_simulation_export_vtk(self):
      
        sim = cratermaker.Simulation(simdir=self.simdir,gridlevel=self.gridlevel)
        # Test with default parameters
        default_out_dir = os.path.join(sim.simdir, _EXPORT_DIR)
        expected_files = ["surf000000.vtp"]
        sim.export_vtk()
        self.assertTrue(os.path.isdir(default_out_dir))
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(default_out_dir, f)))
        
    def test_emplace_crater(self):
        cdiam = 2*self.pix
        sim = cratermaker.Simulation(simdir=self.simdir,gridlevel=self.gridlevel)
        sim.emplace_crater(final_diameter=cdiam)
        pdiam = sim.crater.projectile_diameter
        
        sim.emplace_crater(final_diameter=cdiam)
        sim.emplace_crater(projectile_diameter=pdiam)
        return
    
    def test_populate(self):
        sim = cratermaker.Simulation(simdir=self.simdir,gridlevel=self.gridlevel)
        # Test that populate will work even if no craters are returned
        sim.populate(age=1e-6)
        return
    
    def test_invalid_run_args(self):
        sim = cratermaker.Simulation(simdir=self.simdir,gridlevel=self.gridlevel)

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
            sim.run(age=3.8e3, age_interval=100.0, diameter_number_interval=(300e3, 80))

        # Test case: The diameter_number_interval provided is negative, or is greater than diameter_number - diameter_number_end
        with self.assertRaises(ValueError):
            sim.run(diameter_number=(300e3, 80), diameter_number_end=(300e3, 30), diameter_number_interval=(300e3, 100))

        # Test case: The ninterval is provided but is not an integer or is less than 1
        with self.assertRaises(ValueError):
            sim.run(age=3.8e3, ninterval=0)

        # Test case: The ninterval is provided and either age_interval or diameter_number_interval is also provided
        with self.assertRaises(ValueError):
            sim.run(age=3.8e3, ninterval=100, age_interval=100.0)

        return

if __name__ == '__main__':
    unittest.main()
