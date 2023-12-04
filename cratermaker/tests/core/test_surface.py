import unittest
import os
import shutil
import numpy as np
import tempfile
from cratermaker import Target
from cratermaker import Surface, initialize_surface, generate_grid, generate_data
from cratermaker.core.surface import _DATA_DIR, _GRID_FILE_NAME, _GRID_TEMP_DIR
# Import other necessary modules and functions

class TestSurface(unittest.TestCase):
    
    def setUp(self):
        # Initialize a target and surface for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.grid_file = os.path.join(self.temp_dir.name, _GRID_FILE_NAME)
        self.grid_temp_dir = os.path.join(self.temp_dir.name, _GRID_TEMP_DIR)
        os.mkdir(self.grid_temp_dir)
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        os.chdir(self.temp_dir.name)
        
        return
    
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup() 
        return

    def test_generate_grid(self):
        # Generate grid
        generate_grid(target=self.target, 
                      pix=self.pix, 
                      grid_file=self.grid_file,
                      grid_temp_dir=self.grid_temp_dir)

        # Check if grid file is created
        self.assertTrue(os.path.exists(self.grid_file))
        return

    def test_initialize_surface(self):
        # Delete any existing surface data file and re-generate it
        directory = os.path.join(os.getcwd(), _DATA_DIR)
        if os.path.exists(directory) and os.path.isdir(directory):
            try:
                shutil.rmtree(directory)
            except Exception as error:
                print(f"Error: {directory} : {error}")

        # Initializing it first should run the jigsaw mesh generator
        surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
        
        # Try initializing the surface again with the same parameters. This should find the existing grid file and load it 
        surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=False)
        
        # Test regridding if the parameters change
        n_face_orig = surf.uxgrid.n_face
        
        surf = initialize_surface(pix=2*self.pix, target=self.target, reset_surface=False)
        self.assertGreater(n_face_orig, surf.uxgrid.n_face)
       
        # Test different target values
        surf = initialize_surface(pix=self.pix, target=Target(name="Mars"), reset_surface=False)
        surf = initialize_surface(pix=self.pix, target="Mercury", reset_surface=False)
        
        # Test bad values
        with self.assertRaises(TypeError):
            surf = initialize_surface(pix=self.pix, target=1, reset_surface=False)
        with self.assertRaises(ValueError):
            surf = initialize_surface(pix=self.pix, target="Arrakis", reset_surface=False)
            surf = initialize_surface(pix=self.pix, target=Target(name="Salusa Secundus"), reset_surface=False)
        return

    def test_set_elevation(self):
        surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
        # Test with valid elevation data
        new_elev = np.random.rand(surf.uxgrid.n_face)  # Generate random elevation data
        surf.set_elevation(new_elev)

        # Check if the elevation data is correctly set
        np.testing.assert_array_equal(surf['elevation'].values, new_elev)

        # Test with invalid elevation data (wrong size)
        new_elev = np.random.rand(surf.uxgrid.n_face + 1)  # Incorrect size

        # Expect ValueError for incorrect size
        with self.assertRaises(ValueError):
            surf.set_elevation(new_elev)

        # Test setting elevation to None (should set to zero)
        surf.set_elevation(None)

        # Check if the elevation data is set to zero
        np.testing.assert_array_equal(surf['elevation'].values, np.zeros(surf.uxgrid.n_face))
        
        return
    
    def test_calculate_haversine_distance(self):
        # Example coordinates (lat/lon in radians)
        lon1, lat1 = np.radians(0), np.radians(0)  # Equator, prime meridian
        lon2, lat2 = np.radians(90), np.radians(0)  # 90 degrees East, equator
        radius = 6371000  # Earth's radius in meters

        # Known distance should be 1/4 the circumference of the Earth
        expected_distance = np.pi * radius / 2
        calculated_distance = Surface.calculate_haversine_distance(lon1, lat1, lon2, lat2, radius)

        # Compare the expected and calculated distances
        self.assertAlmostEqual(calculated_distance, expected_distance, places=1)

    def test_get_cell_distance(self):
        surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=True)
        location = (0, 0)  # Equator, prime meridian in radians
        target_radius = self.target.radius

        # Convert important points to radians
        north_pole = (0, np.radians(90))
        south_pole = (0, np.radians(-90))
        antipode = (np.radians(180), 0)

        # Function to find nearest cell index
        def find_nearest_index(point):
            distances = surf.calculate_haversine_distance(surf.uxgrid.face_lon, surf.uxgrid.face_lat, *point, target_radius)
            return np.argmin(distances.data)

        north_pole_idx = find_nearest_index(north_pole)
        south_pole_idx = find_nearest_index(south_pole)
        antipode_idx = find_nearest_index(antipode)

        # Test distances
        distances = surf.get_cell_distance(location, target_radius)
        self.assertAlmostEqual(distances[north_pole_idx], target_radius * np.pi / 2, delta=2*self.pix)
        self.assertAlmostEqual(distances[south_pole_idx], target_radius * np.pi / 2, delta=2*self.pix)
        self.assertAlmostEqual(distances[antipode_idx], target_radius * np.pi, delta=2*self.pix)

    def test_calculate_initial_bearing(self):
        # Example coordinates (lat/lon in radians)
        lon1, lat1 = np.radians(0), np.radians(0)  # Equator, prime meridian
        lon2, lat2 = np.radians(90), np.radians(0)  # 90 degrees East, equator

        # The bearing from (0, 0) to (90, 0) should be 90 degrees in radians
        expected_bearing = np.radians(90)
        calculated_bearing = Surface.calculate_initial_bearing(lon1, lat1, lon2, lat2)

        # Compare the expected and calculated bearings
        self.assertAlmostEqual(calculated_bearing, expected_bearing, places=1)


if __name__ == '__main__':
    unittest.main()