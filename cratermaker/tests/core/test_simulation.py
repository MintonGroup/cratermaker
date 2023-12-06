import unittest
from unittest.mock import patch, MagicMock
import cratermaker 
from cratermaker import Target
import tempfile
import os
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
        sim = cratermaker.Simulation(pix=self.pix, target=self.target)
        sim.save()

if __name__ == '__main__':
    unittest.main()
