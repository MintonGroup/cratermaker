import unittest
import cratermaker 
import tempfile
import os
import warnings

class TestRealistic(unittest.TestCase):

    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.target = cratermaker.Target(name="Moon") 
        self.gridlevel = 4
        os.chdir(self.temp_dir.name) 
        
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup() 
        return           

    def test_realistic(self):
        sim = cratermaker.Simulation(gridlevel=self.gridlevel)
        sim.apply_noise(model="ridged")

if __name__ == '__main__':
    unittest.main()
