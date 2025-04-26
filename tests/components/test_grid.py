import unittest
import os
from pathlib import Path
import tempfile
from cratermaker import Target, get_grid_type, available_grid_types
from cratermaker.constants import  _GRID_FILE_NAME

gridtypes = available_grid_types()

class TestGrid(unittest.TestCase):
    """
    A collection of unit tests for the GridMaker classes in the cratermaker project.

    Parameters
    ----------
    temp_dir : TemporaryDirectory
        A temporary directory for testing file generation and I/O.
    grid_file : str
        Path to the temporary grid file.
    target : Target
        Target object representing a celestial body.
    pix : float
        Pixel size or resolution of the grid.

    """    
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.grid_file = Path(self.temp_dir.name) / _GRID_FILE_NAME
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4
        self.cwd = Path.cwd()
        os.chdir(self.temp_dir.name)
        
        return
    
    def tearDown(self):
        # Clean up temporary directory
        os.chdir(self.cwd)
        self.temp_dir.cleanup() 
        return

    def test_generate_grid(self):
        gridargs = {
            "icosphere": {
                "level" :self.gridlevel, 
                "radius": self.target.radius
                },
            "arbitrary_resolution": {
                "pix": self.pix, 
                "radius": self.target.radius
                },
            "hireslocal": {
                "pix": self.pix, 
                "radius": self.target.radius, 
                "local_location": (0, 0), 
                "local_radius": 100e3, 
                "superdomain_scale_factor": 10
                }
            }

        for name, args in gridargs.items():
            grid = get_grid_type(name)(**args)
            grid.generate_grid()
            self.assertTrue(os.path.exists(grid.file))
        
        return

if __name__ == '__main__':
    unittest.main()
