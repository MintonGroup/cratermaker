import unittest
from unittest.mock import patch, MagicMock
import cratermaker 
from cratermaker import Target
import tempfile
import os
import numpy as np
import xarray as xr
from cratermaker.core.surface import _DATA_DIR, _GRID_FILE_NAME, _GRID_TEMP_DIR

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
        np.testing.assert_array_equal(sim.surf['elevation_face'].values, np.ones(sim.surf.uxgrid.n_face)) 
        np.testing.assert_array_equal(sim.surf['elevation_node'].values, np.ones(sim.surf.uxgrid.n_node)) 
        sim.save()
        
        with xr.open_dataset(os.path.join(sim.surf.data_dir, "elevation_node.nc")) as ds:
            np.testing.assert_array_equal(ds['elevation_node'].values, np.ones(sim.surf.uxgrid.n_node))
        
        with xr.open_dataset(os.path.join(sim.surf.data_dir, "elevation_face.nc")) as ds:
            np.testing.assert_array_equal(ds['elevation_face'].values, np.ones(sim.surf.uxgrid.n_face))
        
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
