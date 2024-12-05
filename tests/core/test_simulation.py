import unittest
import cratermaker 
from cratermaker import Target
import tempfile
import os
import numpy as np
import xarray as xr
from cratermaker.core.surface import _COMBINED_DATA_FILE_NAME
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
        self.target = Target(name="Moon") 
        self.pix = self.target.radius / 10.0
        os.chdir(self.temp_dir.name) 
        
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup() 
        return           

    def test_simulation_defaults(self):
        sim = cratermaker.Simulation(pix=self.pix)
        self.assertEqual(sim.target.name, "Moon")
        
    def test_simulation_save(self):
        # Test basic save operation
        sim = cratermaker.Simulation(pix=self.pix, target=self.target)
        sim.save()
       
        # Test that variables are saved correctly
        sim.surf.set_elevation(1.0)
        np.testing.assert_array_equal(sim.surf["node_elevation"].values, np.ones(sim.surf.uxgrid.n_node)) 
        np.testing.assert_array_equal(sim.surf["face_elevation"].values, np.ones(sim.surf.uxgrid.n_face)) 
        
        sim.save()
        
        filename = os.path.join(sim.data_dir,_COMBINED_DATA_FILE_NAME.replace(".nc", f"{sim.interval_number:06d}.nc"))
        self.assertTrue(os.path.exists(filename))
        with xr.open_dataset(filename) as ds:
            ds = ds.isel(Time=-1)
            np.testing.assert_array_equal(ds["node_elevation"].values, np.ones(sim.surf.uxgrid.n_node))
            np.testing.assert_array_equal(ds["face_elevation"].values, np.ones(sim.surf.uxgrid.n_face))
    
        # Test saving combined data
        sim.save(combine_data_files=True)
        filename = _COMBINED_DATA_FILE_NAME
        
        filename = os.path.join(sim.data_dir,_COMBINED_DATA_FILE_NAME)
        self.assertTrue(os.path.exists(filename))
        with xr.open_dataset(filename) as ds:
            ds = ds.isel(Time=-1)
            np.testing.assert_array_equal(ds["node_elevation"].values, np.ones(sim.surf.uxgrid.n_node))
            np.testing.assert_array_equal(ds["face_elevation"].values, np.ones(sim.surf.uxgrid.n_face))
    
        return 
        
    def test_simulation_export_vtk(self):
      
        sim = cratermaker.Simulation(pix=self.pix, target=self.target) 
        # Test with default parameters
        default_out_dir = os.path.join(sim.simdir, "vtk_files")
        expected_files = ["staticFieldsOnCells.vtp","staticFieldsOnVertices.vtp","timeDependentFieldsOnCells.pvd","timeDependentFieldsOnVertices.pvd"]
        sim.export_vtk()
        self.assertTrue(os.path.isdir(default_out_dir))
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(default_out_dir, f)))
            
        # Test with custom output directory
        custom_out_dir = os.path.join(sim.simdir, "custom_vtk_files")
        sim.export_vtk(out_dir=custom_out_dir)
        self.assertTrue(os.path.isdir(custom_out_dir))
        for f in expected_files:
            self.assertTrue(os.path.exists(os.path.join(custom_out_dir, f)))        
        
    def test_emplace_crater(self):
        cdiam = 2*self.pix
        sim = cratermaker.Simulation(pix=self.pix)
        sim.emplace_crater(diameter=cdiam)
        pdiam = sim.projectile.diameter
        
        sim.emplace_crater(diameter=cdiam)
        sim.emplace_crater(diameter=pdiam, from_projectile=True)
        return
    
    def test_populate(self):
        sim = cratermaker.Simulation(pix=self.pix)
        # Test that populate will work even if no craters are returned
        sim.populate(age=1e-6)
        return
    
    def test_invalid_run_args(self):
        sim = cratermaker.Simulation(pix=self.pix)

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
