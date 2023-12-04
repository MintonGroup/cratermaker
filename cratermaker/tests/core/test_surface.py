import unittest
import os
import shutil
import numpy as np
from cratermaker import Target
from cratermaker import Surface, initialize_surface, generate_grid, generate_data_file
from cratermaker.core.surface import _DATA_DIR
# Import other necessary modules and functions

class TestSurface(unittest.TestCase):
    
    def setUp(self):
        # Initialize a target and surface for testing
        self.target = Target(name="Moon")
        self.pix = np.sqrt(4 * np.pi * self.target.radius**2) * 1e-2 # Make a smaller mesh for testing  
        self.surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
    
    def test_initialize_surface(self):
        # Delete any existing surface data file and re-generate it
        directory = os.path.join(os.getcwd(), _DATA_DIR)
        if os.path.exists(directory) and os.path.isdir(directory):
            try:
                shutil.rmtree(directory)
            except Exception as error:
                print(f"Error: {directory} : {error}")

        # Initializing it first should run the jigsaw mesh generator
        self.surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
        
        # Try initializing the surface again with the same parameters. This should find the existing grid file and load it 
        self.surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
        return

    def test_set_elevation(self):
        # Test set_elevation with valid and invalid inputs
        pass

    def test_calculate_haversine_distance(self):
        # Test calculate_haversine_distance with example values
        pass

    def test_get_cell_distance(self):
        # Test get_cell_distance functionality
        pass

    def test_calculate_initial_bearing(self):
        # Test calculate_initial_bearing with specific cases
        pass

    # Continue with other test methods for each functionality



    def test_generate_grid(self):
        # Test generate_grid with different parameters
        pass

    # Add more tests as needed

if __name__ == '__main__':
    unittest.main()