import unittest
from unittest.mock import patch, MagicMock
import cratermaker 
from cratermaker import Target
import tempfile
import os
import numpy as np
import xarray as xr
from cratermaker.core.surface import  _COMBINED_DATA_FILE_NAME

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
        self.assertEqual(sim.target.material.name, "Soft Rock")
        
    def test_simulation_save(self):
        # Test basic save operation
        sim = cratermaker.Simulation(pix=self.pix, target=self.target)
        sim.save()
       
        # Test that variables are saved correctly
        sim.surf.set_elevation(1.0)
        np.testing.assert_array_equal(sim.surf["elevation"].values, np.ones(sim.surf.uxgrid.n_node)) 
        sim.save()
        
        with xr.open_dataset(sim.surf.elevation_file) as ds:
            np.testing.assert_array_equal(ds["elevation"].values, np.ones(sim.surf.uxgrid.n_node))
        
        # Test saving combined data
        sim.save(combine_data_files=True)
        combined_file = os.path.join(sim.surf.data_dir, _COMBINED_DATA_FILE_NAME)
        self.assertTrue(os.path.exists(combined_file))
        
        with xr.open_dataset(combined_file) as ds:
            np.testing.assert_array_equal(ds["elevation"].values, np.ones(sim.surf.uxgrid.n_node))
       
        return 
        
    def test_simulation_export_vtk(self):
        sim = cratermaker.Simulation(pix=self.pix, target=self.target) 
        
        # Test with default parameters
        default_out_dir = os.path.join(sim.simdir, "vtk_files")
        expected_files = ["staticFieldsOnCells.vtp", "staticFieldsOnVertices.vtp"]
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
            

if __name__ == '__main__':
    unittest.main()
