import unittest
import os
import shutil
import numpy as np
import tempfile
from cratermaker import Target
from cratermaker import Surface
from cratermaker.core.surface import UniformGrid, HiResLocalGrid, _DATA_DIR, _GRID_FILE_NAME, _GRID_TEMP_DIR
from cratermaker.utils.montecarlo import get_random_location
from cratermaker.utils.general_utils import normalize_coords


class TestSurface(unittest.TestCase):
    """
    A collection of unit tests for the Surface class in the cratermaker project.

    Parameters
    ----------
    temp_dir : TemporaryDirectory
        A temporary directory for testing file generation and I/O.
    grid_file : str
        Path to the temporary grid file.
    grid_temp_dir : str
        Path to the temporary directory for grid generation.
    target : Target
        Target object representing a celestial body.
    pix : float
        Pixel size or resolution of the grid.

    """    
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
        grid_strategy = UniformGrid(pix=self.pix, radius=self.target.radius)
        grid_strategy.generate_grid(grid_file=self.grid_file, grid_temp_dir=self.grid_temp_dir)
        self.assertTrue(os.path.exists(self.grid_file))
        
        grid_strategy = HiResLocalGrid(pix=self.pix, radius=self.target.radius, local_location=(0, 0), local_radius=100e3, superdomain_scale_factor=10)
        grid_strategy.generate_grid(grid_file=self.grid_file, grid_temp_dir=self.grid_temp_dir)
        self.assertTrue(os.path.exists(self.grid_file))
        return

    def test_initialize_surface(self):
        directory = os.path.join(os.getcwd(), _DATA_DIR)
        if os.path.exists(directory) and os.path.isdir(directory):
            try:
                shutil.rmtree(directory)
            except Exception as error:
                print(f"Error: {directory} : {error}")

        # Initializing it first should run the jigsaw mesh generator
        surf = Surface.initialize(pix=self.pix, target=self.target, reset_surface=True)
        
        # Try initializing the surface again with the same parameters. This should find the existing grid file and load it 
        surf = Surface.initialize(pix=self.pix, target=self.target, reset_surface=False)
        
        # Test regridding if the parameters change
        n_face_orig = surf.uxgrid.n_face
    
        surf = Surface.initialize(pix=2*self.pix, target=self.target, reset_surface=False)
        self.assertGreater(n_face_orig, surf.uxgrid.n_face)
    
        # Test different target values
        surf = Surface.initialize(pix=self.pix, target=Target(name="Mars"), reset_surface=False)
        surf = Surface.initialize(pix=self.pix, target="Mercury", reset_surface=False)
        
        # Test bad values
        with self.assertRaises(TypeError):
            surf = Surface.initialize(pix=self.pix, target=1, reset_surface=False)
        with self.assertRaises(ValueError):
            surf = Surface.initialize(pix=self.pix, target="Arrakis", reset_surface=False)
        with self.assertRaises(ValueError):
            surf = Surface.initialize(pix=self.pix, target=Target(name="Salusa Secundus"), reset_surface=False)
        return


    def test_set_elevation(self):
        surf = Surface.initialize(pix=self.pix, target=self.target, reset_surface=True)
        # Test with valid elevation data
        new_elev = np.random.rand(surf.uxgrid.n_node)  # Generate random elevation data
        surf.set_elevation(new_elev)

        # Test with invalid elevation data (wrong size)
        new_elev = np.random.rand(surf.uxgrid.n_node + 1)  # Incorrect size

        # Expect ValueError for incorrect size
        with self.assertRaises(ValueError):
            surf.set_elevation(new_elev)

        # Test setting elevation to None (should set to zero)
        surf.set_elevation(None)

        # Check if the elevation data is set to zero
        np.testing.assert_array_equal(surf['node_elevation'].values, np.zeros(surf.uxgrid.n_node))
        np.testing.assert_array_equal(surf['face_elevation'].values, np.zeros(surf.uxgrid.n_face))
        
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


    def test_get_face_distance(self):
        surf = Surface.initialize(pix=self.pix, target=self.target, reset_surface=True)
        
        location = get_random_location()
        lon = location[0]
        lat = location[1]

        # Convert important points to radians
        north = (lon, lat + 90)
        south = (lon, lat - 90)
        antipode = (lon-180, -lat)
        
        north = normalize_coords(north)
        south = normalize_coords(south)
        antipode = normalize_coords(antipode)
        
        _,north_idx = surf.find_nearest_index(north)
        _,south_idx = surf.find_nearest_index(south)
        _,antipode_idx = surf.find_nearest_index(antipode)
        
        north_distance = surf.target.radius * np.pi / 2 
        south_distance = surf.target.radius * np.pi / 2
        antipode_distance = surf.target.radius * np.pi
        delta = 2*self.pix

        # Test distances
        _,distances = surf.get_distance(location)
        self.assertAlmostEqual(distances[north_idx], north_distance, delta=delta, msg=f"North face distance ratio: {distances[north_idx].item()/north_distance}")
        self.assertAlmostEqual(distances[south_idx], south_distance, delta=delta, msg=f"South face distance ratio: {distances[south_idx].item()/south_distance}")
        self.assertAlmostEqual(distances[antipode_idx], antipode_distance, delta=delta, msg=f"Antipode face distance ratio: {distances[antipode_idx].item()/antipode_distance}")
        
        
    def test_get_node_distance(self):
        surf = Surface.initialize(pix=self.pix, target=self.target, reset_surface=True)
        
        location = get_random_location()
        lon = location[0]
        lat = location[1]

        # Convert important points to radians
        north = (lon, lat + 90)
        south = (lon, lat - 90)
        antipode = (lon+180, -lat)
        
        north = normalize_coords(north)
        south = normalize_coords(south)
        antipode = normalize_coords(antipode)
        
        north_idx,_ = surf.find_nearest_index(north)
        south_idx,_ = surf.find_nearest_index(south)
        antipode_idx,_ = surf.find_nearest_index(antipode)

        north_distance = surf.target.radius * np.pi / 2 
        south_distance = surf.target.radius * np.pi / 2
        antipode_distance = surf.target.radius * np.pi
        delta = 2*self.pix

        # Test distances
        node_distances,_ = surf.get_distance(location)
        self.assertAlmostEqual(node_distances[north_idx], north_distance, delta=delta, msg=f"North node distance ratio: {node_distances[north_idx].item()/north_distance}")
        self.assertAlmostEqual(node_distances[south_idx], south_distance, delta=delta, msg=f"South node distance ratio: {node_distances[south_idx].item()/south_distance}")
        self.assertAlmostEqual(node_distances[antipode_idx], antipode_distance, delta=delta, msg=f"Antipode node distance ratio: {node_distances[antipode_idx].item()/antipode_distance}")
        
        
    def test_calculate_initial_bearing(self):
        # Example coordinates (lat/lon in radians)
        lon1, lat1 = np.radians(0), np.radians(0)  # Equator, prime meridian
        lon2, lat2 = np.radians(90), np.radians(0)  # 90 degrees East, equator

        # The bearing from (0, 0) to (90, 0) should be 90 degrees in radians
        expected_bearing = np.radians(90)
        calculated_bearing = Surface.calculate_initial_bearing(lon1, lat1, lon2, lat2)

        # Compare the expected and calculated bearings
        self.assertAlmostEqual(calculated_bearing, expected_bearing, places=1)


    #TODO: Fix this test case since the method waws
    # def test_average_surface(self):
    #     """
    #     Test the average_surface function. This will generate a random noisy surface, compute an average over a small region, and then test that the surface elevations relative to the average are close to zero.
    #     """
    #     sim = Simulation(pix=self.pix)
    #     sim.apply_noise(model="ridged")
    #     location = get_random_location()
    #     region_radius = 100e3

    #     sim.surf.get_reference_surface(location, region_radius)
        
    #     # Find cells within the crater radius
    #     if 'face_crater_distance' in sim.surf:
    #         cells_within_radius = sim.surf['face_crater_distance'] <= region_radius
    #     else:
    #         _, cells_within_radius = sim.surf.get_distance(location)
    #         cells_within_radius = cells_within_radius <= region_radius
        
    #     vert_vars = ['face_x', 'face_y', 'face_z'] 
    #     region_mesh = sim.surf.uxgrid._ds[vert_vars].where(cells_within_radius, drop=True)
    #     region_elevation = sim.surf['face_elevation'].where(cells_within_radius, drop=True) 
    #     region_surf = elevation_to_cartesian(region_mesh, region_elevation)
        
    #     cap_mult = np.linalg.norm(sim.surf.reference_surface_vector) / sim.target.radius
        
    #     region_mesh *= cap_mult 
    #     region_mesh['face_x'] += sim.surf.reference_surface_center[0] 
    #     region_mesh['face_y'] += sim.surf.reference_surface_center[1] 
    #     region_mesh['face_z'] += sim.surf.reference_surface_center[2] 

    #     # Fetch x, y, z values of the mesh within the region
    #     region_delta = region_surf - region_mesh

    #     # Fetch the areas of the cells within the region
    #     cell_areas = sim.surf['face_areas'].where(cells_within_radius, drop=True)

    #     # Calculate the weighted average of the x, y, and z coordinates to get the average surface vector
    #     weighted_x = (region_delta['face_x'] * cell_areas).sum() / cell_areas.sum()
    #     weighted_y = (region_delta['face_y'] * cell_areas).sum() / cell_areas.sum()
    #     weighted_z = (region_delta['face_z'] * cell_areas).sum() / cell_areas.sum()
    #     average_center = np.array([weighted_x.item(), weighted_y.item(), weighted_z.item()])
    #     average_center = np.linalg.norm(average_center)
    #     self.assertAlmostEqual(average_center / sim.target.radius, 0.0, delta=1e-5)
        

if __name__ == '__main__':
     
    unittest.main()