import unittest
import cratermaker 
import tempfile
import os
import warnings
try:
    import mpas_tools
    MPAS_TOOLS_AVAILABLE = True
except ModuleNotFoundError:
    MPAS_TOOLS_AVAILABLE = False

class TestRealistic(unittest.TestCase):

    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.target = cratermaker.Target(name="Moon") 
        self.pix = self.target.radius / 10.0
        os.chdir(self.temp_dir.name) 
        
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup() 
        return           

    def test_realistic(self):
        if MPAS_TOOLS_AVAILABLE:
            sim = cratermaker.Simulation(pix=self.pix)
            sim.apply_noise(model="ridged")
        else:
            warnings.warn("mpas_tools not available. Skipping test_realistic.")
        

if __name__ == '__main__':
    unittest.main()
