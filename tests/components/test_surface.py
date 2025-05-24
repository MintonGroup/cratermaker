import tempfile
import unittest
from pathlib import Path

import numpy as np

from cratermaker import Surface, Target
from cratermaker.utils.general_utils import normalize_coords
from cratermaker.utils.montecarlo_utils import get_random_location

surfacetypes = Surface.available()


class TestSurface(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.target = Target(name="Moon")
        self.pix = self.target.radius / 10.0
        self.gridlevel = 4

        return

    def test_initialize_surface(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Initializing it first should run the mesh generator
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )
            del surface

            # Try initializing the surface again with the same parameters. This should find the existing grid file and load it
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=False
            )

            # Test regridding if the parameters change
            n_face_orig = surface.uxgrid.n_face
            del surface

            surface = Surface.maker(
                simdir=simdir,
                gridlevel=self.gridlevel - 1,
                target=self.target,
                reset=False,
            )
            self.assertGreater(n_face_orig, surface.uxds.uxgrid.n_face)
            del surface

            # Test different target values
            surface = Surface.maker(
                simdir=simdir,
                gridlevel=self.gridlevel,
                target=Target(name="Mars"),
                reset=False,
            )
            del surface

            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target="Mercury", reset=False
            )
            del surface

            # Test bad values
            with self.assertRaises(TypeError):
                Surface.maker(
                    simdir=simdir, gridlevel=self.gridlevel, target=1, reset=False
                )
            with self.assertRaises(ValueError):
                Surface.maker(
                    simdir=simdir,
                    gridlevel=self.gridlevel,
                    target="Arrakis",
                    reset=False,
                )
            with self.assertRaises(ValueError):
                Surface.maker(
                    simdir=simdir,
                    gridlevel=self.gridlevel,
                    target=Target(name="Salusa Secundus"),
                    reset=False,
                )
        return

    def test_update_elevation(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )
            # Test with valid elevation data
            new_elev = np.random.rand(
                surface.uxds.uxgrid.n_node
            )  # Generate random elevation data
            surface.update_elevation(new_elev)

            # Test with invalid elevation data (wrong size)
            new_elev = np.random.rand(surface.uxds.uxgrid.n_node + 1)  # Incorrect size

            # Expect ValueError for incorrect size
            with self.assertRaises(ValueError):
                surface.update_elevation(new_elev)

            surface.update_elevation(0.0, overwrite=True)

            # Check if the elevation data is set to zero
            np.testing.assert_array_equal(
                surface.uxds["node_elevation"].values,
                np.zeros(surface.uxds.uxgrid.n_node),
            )
            np.testing.assert_array_equal(
                surface.uxds["face_elevation"].values,
                np.zeros(surface.uxds.uxgrid.n_face),
            )

        return

    def test_calculate_haversine_distance(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Example coordinates (lat/lon in radians)
            lon1, lat1 = np.radians(0), np.radians(0)  # Equator, prime meridian
            lon2, lat2 = np.radians(90), np.radians(0)  # 90 degrees East, equator
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target="Earth"
            )

            # Known distance should be 1/4 the circumference of the Earth
            expected_distance = np.pi * surface.radius / 2
            calculated_distance = surface.calculate_haversine_distance(
                lon1, lat1, lon2, lat2
            )[0]

            # Compare the expected and calculated distances
            self.assertAlmostEqual(calculated_distance, expected_distance, places=1)

    def test_get_face_distance(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )

            location = get_random_location()[0]
            lon = location[0]
            lat = location[1]

            # Convert important points to radians
            north = (lon, lat + 90)
            south = (lon, lat - 90)
            antipode = (lon - 180, -lat)

            north = normalize_coords(north)
            south = normalize_coords(south)
            antipode = normalize_coords(antipode)

            north_idx, _ = surface.find_nearest_index(north)
            south_idx, _ = surface.find_nearest_index(south)
            antipode_idx, _ = surface.find_nearest_index(antipode)

            north_distance = surface.target.radius * np.pi / 2
            south_distance = surface.target.radius * np.pi / 2
            antipode_distance = surface.target.radius * np.pi
            delta = 2 * self.pix

            # Test distances
            distances, _ = surface.get_distance(location)
            self.assertAlmostEqual(
                distances[north_idx],
                north_distance,
                delta=delta,
                msg=f"North face distance ratio: {distances[north_idx].item() / north_distance}",
            )
            self.assertAlmostEqual(
                distances[south_idx],
                south_distance,
                delta=delta,
                msg=f"South face distance ratio: {distances[south_idx].item() / south_distance}",
            )
            self.assertAlmostEqual(
                distances[antipode_idx],
                antipode_distance,
                delta=delta,
                msg=f"Antipode face distance ratio: {distances[antipode_idx].item() / antipode_distance}",
            )

    def test_get_node_distance(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )

            location = get_random_location()[0]
            lon = location[0]
            lat = location[1]

            # Convert important points to radians
            north = (lon, lat + 90)
            south = (lon, lat - 90)
            antipode = (lon + 180, -lat)

            north = normalize_coords(north)
            south = normalize_coords(south)
            antipode = normalize_coords(antipode)

            _, north_idx = surface.find_nearest_index(north)
            _, south_idx = surface.find_nearest_index(south)
            _, antipode_idx = surface.find_nearest_index(antipode)

            north_distance = surface.target.radius * np.pi / 2
            south_distance = surface.target.radius * np.pi / 2
            antipode_distance = surface.target.radius * np.pi
            delta = 2 * self.pix

            # Test distances
            _, node_distances = surface.get_distance(location)
            self.assertAlmostEqual(
                node_distances[north_idx],
                north_distance,
                delta=delta,
                msg=f"North node distance ratio: {node_distances[north_idx].item() / north_distance}",
            )
            self.assertAlmostEqual(
                node_distances[south_idx],
                south_distance,
                delta=delta,
                msg=f"South node distance ratio: {node_distances[south_idx].item() / south_distance}",
            )
            self.assertAlmostEqual(
                node_distances[antipode_idx],
                antipode_distance,
                delta=delta,
                msg=f"Antipode node distance ratio: {node_distances[antipode_idx].item() / antipode_distance}",
            )

    def test_calculate_initial_bearing(self):
        # Example coordinates (lat/lon in radians)
        lon1, lat1 = np.radians(0), np.radians(0)  # Equator, prime meridian
        lon2, lat2 = np.radians(90), np.radians(0)  # 90 degrees East, equator

        # The bearing from (0, 0) to (90, 0) should be 90 degrees in radians
        expected_bearing = np.radians(90)
        calculated_bearing = Surface.calculate_initial_bearing(lon1, lat1, lon2, lat2)

        # Compare the expected and calculated bearings
        self.assertAlmostEqual(calculated_bearing, expected_bearing, places=1)

    # def test_get_random_on_face(self):
    #     # Tests that the random location is within the face we expect
    #     surface = Surface.maker(gridlevel=self.gridlevel*2, target=self.target, reset=True)
    #     n_per_face = 10
    #     for i in surface.n_face:
    #         original_face_index = i.values.item()
    #         for _ in range(n_per_face):
    #             location = surface.get_random_location_on_face(original_face_index)
    #             new_face_index, _ = surface.find_nearest_index(location)
    #             self.assertEqual(original_face_index, new_face_index)

    def test_face_surface_values(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            # Tests that the face_surface generates the correct values
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )
            total_area_1 = surface.uxgrid.calculate_total_face_area()
            total_area_2 = surface.face_areas.sum().item()
            ratio = np.sqrt(total_area_2 / total_area_1) / self.target.radius
            self.assertAlmostEqual(ratio, 1.0, places=2)
            self.assertAlmostEqual(
                total_area_2 / (4 * np.pi * self.target.radius**2), 1.0, places=2
            )

    def test_generate_grid(self):
        gridargs = {
            "icosphere": {
                "level": self.gridlevel,
            },
            "arbitrary_resolution": {
                "pix": self.pix,
            },
            "hireslocal": {
                "pix": self.pix,
                "local_location": (0, 0),
                "local_radius": 100e3,
                "superdomain_scale_factor": 10,
            },
        }

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            for name, args in gridargs.items():
                surface = Surface.maker(name, simdir=simdir, **args)
                self.assertTrue(Path(surface.grid_file).exists())

        return

    def test_add_data(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as simdir:
            test_value = 76.0
            surface = Surface.maker(
                simdir=simdir, gridlevel=self.gridlevel, target=self.target, reset=True
            )
            region_view = surface.extract_region(location=(0, 0), region_radius=100e3)

            for obj, uxds in zip(
                [surface, region_view], [surface.uxds, region_view.surface.uxds]
            ):
                surface.reset()
                obj.add_data(
                    name="test_data",
                    long_name="data for testing",
                    units="furlongs per fortnight",
                    data=test_value,
                )
                n_face = obj.n_face
                n_node = obj.n_node
                face_indices = obj.face_indices
                node_indices = obj.node_indices
                # Test that the attributes were added correctly
                self.assertIn("test_data", uxds)
                self.assertEqual(uxds["test_data"].long_name, "data for testing")
                self.assertEqual(uxds["test_data"].units, "furlongs per fortnight")
                np.testing.assert_array_equal(
                    uxds["test_data"].sel(n_face=face_indices).values,
                    np.full(n_face, test_value),
                )

                obj.add_data(
                    name="test_data",
                    data=-test_value,
                    long_name="this should be ignored",
                    units="slug angstrom^2 / fortnight^2",
                    overwrite=True,
                )
                # test that the attributes were ignored but the data was added
                self.assertIn("test_data", uxds)
                self.assertEqual(uxds["test_data"].long_name, "data for testing")
                self.assertEqual(uxds["test_data"].units, "furlongs per fortnight")
                np.testing.assert_array_equal(
                    uxds["test_data"].sel(n_face=face_indices).values,
                    np.full(n_face, -test_value),
                )

                # Add scalar node data
                obj.add_data(name="scalar_node", data=test_value, isfacedata=False)
                self.assertIn("scalar_node", uxds)
                np.testing.assert_array_equal(
                    uxds["scalar_node"].sel(n_node=node_indices).values,
                    np.full(n_node, test_value),
                )

                # Add array face data
                data_array = np.arange(n_face)
                obj.add_data(name="array_face", data=data_array)
                np.testing.assert_array_equal(
                    uxds["array_face"].sel(n_face=face_indices).values, data_array
                )

                # Overwrite behavior
                obj.add_data(name="array_face", data=np.ones(n_face), overwrite=False)
                np.testing.assert_array_equal(
                    uxds["array_face"].sel(n_face=face_indices).values,
                    data_array + 1,
                )

                obj.add_data(name="array_face", data=np.zeros(n_face), overwrite=True)
                np.testing.assert_array_equal(
                    uxds["array_face"].sel(n_face=face_indices).values,
                    np.zeros(n_face),
                )


if __name__ == "__main__":
    unittest.main()
