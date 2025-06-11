from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
import threading
import warnings
from abc import abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import uxarray as uxr
import xarray as xr
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import OptimizeWarning, curve_fit
from uxarray import INT_FILL_VALUE, UxDataArray, UxDataset

from cratermaker._cratermaker import surface_functions
from cratermaker.constants import (
    _COMBINED_DATA_FILE_NAME,
    _GRID_FILE_NAME,
    _SMALLFAC,
    _SURFACE_DIR,
    _VSMALL,
    FloatLike,
)
from cratermaker.utils import export
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import (
    format_large_units,
    validate_and_normalize_location,
)
from cratermaker.utils.montecarlo_utils import get_random_location_on_face

if TYPE_CHECKING:
    from uxarray.grid.neighbors import SpatialHash

    from cratermaker.components.target import Target

surface_lock = threading.Lock()


class Surface(ComponentBase):
    _registry: dict[str, type[Surface]] = {}

    """
    Used for handling surface-related data and operations in the cratermaker project.

    It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations.
    The Surface class extends UxDataset for the cratermaker project.

    Parameters
    ----------
    target : Target, optional
        The target body or name of a known target body for the impact simulation.
    reset : bool, optional
        Flag to indicate whether to reset the surface. Default is the value of `regrid`
    regrid : bool, optional
        Flag to indicate whether to regrid the surface. Default is False.
    simdir : str | Path
        The main project simulation directory. Defaults to the current working directory if None.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        target: Target | str | None = None,
        simdir: str | Path | None = None,
        **kwargs,
    ):
        from cratermaker.components.target import Target

        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_uxds", None)
        object.__setattr__(self, "_pix", None)
        object.__setattr__(self, "_pix_mean", None)
        object.__setattr__(self, "_pix_std", None)
        object.__setattr__(self, "_pix_min", None)
        object.__setattr__(self, "_pix_max", None)
        object.__setattr__(self, "_area", None)
        object.__setattr__(self, "_edge_tree", None)
        object.__setattr__(self, "_edge_face_distance", None)
        object.__setattr__(self, "_edge_lengths", None)
        object.__setattr__(self, "_edge_indices", None)
        object.__setattr__(self, "_node_tree", None)
        object.__setattr__(self, "_face_tree", None)
        object.__setattr__(self, "_face_area", None)
        object.__setattr__(self, "_face_size", None)
        object.__setattr__(self, "_face_bin_indices", None)
        object.__setattr__(self, "_face_bin_argmin", None)
        object.__setattr__(self, "_face_bin_argmax", None)
        object.__setattr__(self, "_face_bin_area", None)
        object.__setattr__(self, "_face_x", None)
        object.__setattr__(self, "_face_y", None)
        object.__setattr__(self, "_face_z", None)
        object.__setattr__(self, "_node_x", None)
        object.__setattr__(self, "_node_y", None)
        object.__setattr__(self, "_node_z", None)
        object.__setattr__(self, "_smallest_length", None)

        super().__init__(simdir=simdir, **kwargs)

        self._data_variable_init = {
            "node_elevation": {
                "units": "m",
                "long_name": "elevation of nodes",
                "initial_value": 0.0,
                "isfacedata": False,
            },
            "face_elevation": {
                "units": "m",
                "long_name": "elevation of faces",
                "initial_value": 0.0,
                "isfacedata": True,
            },
        }

        self.target = Target.maker(target, **kwargs)

        return

    def __str__(self) -> str:
        base = super().__str__()
        return (
            f"{base}\nTarget: {self.target.name}\nGrid File: {self.grid_file}\n"
            f"Number of faces: {self.n_face}\n"
            f"Number of nodes: {self.n_node}"
        )

    def __del__(self):
        try:
            if hasattr(self, "_uxds") and hasattr(self._uxds, "close"):
                if hasattr(self._uxds, "uxgrid") and hasattr(self._uxds.uxgrid, "_ds"):
                    self._uxds.uxgrid._ds.close()
                self._uxds.close()
        except Exception as e:
            warnings.warn(f"An error occurred while closing the dataset: {e}", RuntimeWarning, stacklevel=2)
        try:
            if hasattr(self, "_grid") and hasattr(self._grid, "uxgrid") and hasattr(self._grid.uxgrid, "_ds"):
                self._grid.uxgrid._ds.close()
        except Exception as e:
            warnings.warn(f"An error occurred while closing the dataset: {e}", RuntimeWarning, stacklevel=2)

    @classmethod
    def maker(
        cls: Surface,
        surface: str | Surface | None = None,
        target: Target | str | None = None,
        reset: bool = False,
        regrid: bool = False,
        simdir: str | Path | None = None,
        **kwargs,
    ) -> Surface:
        """
        Factory method to create a Surface instance from a grid file.

        Parameters
        ----------
        surface : str or Surface, optional
            The name of the type of grid used for the surface. Default is "icosphere".
        target : Target, optional
            The target body or name of a known target body for the impact simulation.
        reset : bool, optional
            Flag to indicate whether to reset the surface. Default is True.
        regrid : bool, optional
            Flag to indicate whether to regrid the surface. Default is False.
        simdir : str | Path
            The main project simulation directory. Defaults to the current working directory if None.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        Surface
            An initialized Surface object.
        """
        if surface is None:
            surface = "icosphere"

        surface = super().maker(
            component=surface,
            target=target,
            reset=reset,
            regrid=regrid,
            simdir=simdir,
            **kwargs,
        )
        return surface

    def reset(self, **kwargs: Any) -> None:
        """
        Reset the surface to its initial state.
        """
        # Remove all old data from the dataset
        varlist = list(self.uxds.data_vars)
        for name in varlist:
            if name not in self._data_variable_init:
                del self.uxds[name]
        for name, entry in self._data_variable_init.items():
            self._add_new_data(
                name=name,
                data=entry["initial_value"],
                long_name=entry["long_name"],
                units=entry["units"],
                isfacedata=entry["isfacedata"],
                save_to_file=True,
            )

        return

    def extract_region(self, location: tuple[FloatLike, FloatLike], region_radius: FloatLike):
        """
        Extract a regional grid based on a given location and radius.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.
        region_radius : float
            The radius of the region to extract in meters.

        Returns
        -------
        LocalSurface
            A LocalSurface object containing a view of the regional grid.

        """
        region_angle = np.rad2deg(region_radius / self.radius)
        if len(location) == 1:
            location = location.item()
        location = validate_and_normalize_location(location)
        coords = np.asarray(location)

        if region_angle < 180.0:
            face_indices = self.face_tree.query_radius(coords, region_angle)
            if len(face_indices) == 0:
                return None

            # First select edges and nodes that are attached to these faces
            edge_indices = np.unique(self.face_edge_connectivity[face_indices].ravel())
            edge_indices = edge_indices[edge_indices != INT_FILL_VALUE]

            node_indices = np.unique(self.face_node_connectivity[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]

            # Now add in all faces that are connected to anything inside the region, so that the outermost border of the local region has a buffer of faces
            # These are needed for diffusion calculations
            neighbor_faces = self.face_face_connectivity[face_indices]
            neighbor_faces = neighbor_faces[neighbor_faces != INT_FILL_VALUE]
            node_faces = self.node_face_connectivity[node_indices]
            node_faces = node_faces[node_faces != INT_FILL_VALUE]
            edge_faces = self.edge_face_connectivity[edge_indices]
            edge_faces = edge_faces[edge_faces != INT_FILL_VALUE]
            face_indices = np.unique(np.concatenate((face_indices, neighbor_faces, node_faces, edge_faces)))
        else:  # This is the entire surface
            face_indices = slice(None)
            edge_indices = slice(None)
            node_indices = slice(None)

        return LocalSurface(
            surface=self,
            face_indices=face_indices,
            node_indices=node_indices,
            edge_indices=edge_indices,
            location=location,
            region_radius=region_radius,
        )

    def add_data(
        self,
        name: str,
        data: FloatLike | NDArray,
        long_name: str | None = None,
        units: str | None = None,
        isfacedata: bool = True,
        overwrite: bool = False,
    ) -> None:
        """
        Adds new data to the surface.

        If the data variable already exists, it will be overwritten if `overwrite` is set to True.

        Parameters
        ----------
        name : str
            Name of the data variable. This will also be used as the data file name.
        data : scalar or array-like
            Data file to be saved. If data is a scalar, then the data file will be filled with that value. If data is an array, then the data file will be filled with the array values. The data array must have the same size as the number of faces or nodes in the grid.
        long_name : str, optional
            Long name of the data variable that will be saved as an attribute if this is new data. If the data already exists on the surface, this will be ignored.
        units : str, optional
            Units of the data variable that will be saved as an attribute if this is new data. If the data already exists on the surface, this will be ignored.
        isfacedata : bool, optional, default True
            Flag to indicate whether the data is face data or node data. This is only needed if `data` is a scalar, otherwise it is ignored
        overwrite : bool, optional, default False
            By default, new data is added to the old data. This flag indicates that the data should be overwritten, replacing any old data with the new data.


        Returns
        -------
        None
        """
        return self._full().add_data(
            name=name,
            data=data,
            long_name=long_name,
            units=units,
            isfacedata=isfacedata,
            overwrite=overwrite,
        )

    def update_elevation(
        self,
        new_elevation: ArrayLike | FloatLike,
        overwrite: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Update the elevation data for the target's surface mesh. This method will determine whether to update the node or face data (or both) depending on the size of the input data. If a scalar is passed, both node and face elevation will be updated to that value. By default, the elevation data will be added to the existing data. If `overwrite` is set to True, the existing data will be replaced with the new data.

        Parameters
        ----------
        new_elevation : ArrayLike | FloatLike
            Elevation to be added (or replaced, if overwrite is True). This can be a scalar, an array with the same size as the number of faces, an array with the same size as the number of nodes, or an array with the same size as the number of faces + nodes.
        overwrite : bool, optional
            If True, the existing data will be replaced with the new data. Default is False.
        **kwargs : Any
            Additional keyword arguments.

        Notes
        -----
        When passing combined data, the first part of the array will be used for face elevation and the second part for node elevation.
        """
        return self._full().update_elevation(new_elevation=new_elevation, overwrite=overwrite, **kwargs)

    def apply_diffusion(self, kdiff: FloatLike | NDArray) -> None:
        """
        Apply diffusion to the surface.

        Parameters
        ----------
        kdiff : float or array-like
            The degradation state of the surface, which is the product of diffusivity and time. It can be a scalar or an array of the same size as the number of faces in the grid.
            If it is a scalar, the same value is applied to all faces. If it is an array, it must have the same size as the number of faces in the grid.
            The value of kdiff must be greater than 0.0.

        """
        return self._full().apply_diffusion(kdiff)

    def slope_collapse(self, critical_slope_angle: FloatLike = 35.0) -> None:
        """
        Collapse all slopes larger than the critical slope angle.

        Parameters
        ----------
        critical_slope_angle : float
            The critical slope angle (angle of repose) in degrees.
        """
        return self._full().slope_collapse(critical_slope_angle)

    def compute_slope(self) -> NDArray[np.float64]:
        """
        Compute the slope of the surface.

        Returns
        -------
        NDArray[np.float64]
            The slope of all faces in degrees.
        """
        return self._full().compute_slope()

    def apply_noise(
        self,
        model: str = "turbulence",
        noise_width: FloatLike = 1000e3,
        noise_height: FloatLike = 1e3,
        **kwargs: Any,
    ) -> None:
        """
        Apply noise to the surface.

        Parameters
        ----------
        model : str
            The noise model to use. Options are "turbulence"
        noise_width : float
            The width of the noise in meters.
        noise_height : float
            The height of the noise in meters.
        kwargs : Any
            Additional arguments to pass to the noise model.
        """
        return self._full().apply_noise(model=model, noise_width=noise_width, noise_height=noise_height, **kwargs)

    def calculate_face_and_node_distances(self, location: tuple[float, float]) -> tuple[NDArray, NDArray]:
        """
        Computes the distances from a given location to all faces and nodes.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        NDArray
            Array of distances for each face in meters.
        NDArray
            Array of distances for each node in meters.
        """
        return self._full().calculate_face_and_node_distances(location)

    def calculate_face_and_node_bearings(self, location: tuple[float, float]) -> tuple[NDArray, NDArray]:
        """
        Computes the initial bearing from a given location to all faces and nodes.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        NDArray
            Array of initial bearings for each face in radians.
        NDArray
            Array of initial bearings for each node in radians.

        Notes
        -----
        This is intended to be used as a helper to calculate_face_and_node_bearings.
        """
        return self._full().calculate_face_and_node_bearings(location)

    def find_nearest_node(self, location):
        """
        Find the index of the nearest node to a given point.

        This method calculates the Haversine distance from the given point to each node in the grid, and returns the index of the node with the minimum distance.

        Parameters
        ----------
        location : tuple
            A tuple containing two elements: (longitude, latitude) in degrees.

        Returns
        -------
        int
            The index of the nearest node in the grid to the given point.

        Notes
        -----
        The method uses the ball tree query method that is included in the UxArray.Grid class.
        """
        location = validate_and_normalize_location(location)
        if len(location) == 1:
            location = location.item()
        coords = np.asarray(location)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", Warning)
            node_ind = self.node_tree.query(coords=coords, k=1, return_distance=False)

        return node_ind.item()

    def find_nearest_face(self, location):
        """
        Find the index of the nearest face to a given point.

        This method calculates the Haversine distance from the given point to each face in the grid, and returns the index of the face with the minimum distance.

        Parameters
        ----------
        location : tuple
            A tuple containing two elements: (longitude, latitude) in degrees.

        Returns
        -------
        int
            The index of the nearest face in the grid to the given point.

        Notes
        -----
        The method uses the ball tree query method that is included in the UxArray.Grid class.  It is only approximate, as it looks for whichever face center is closest to the location. This means that it can return a neighboring face, rather than the face that contains the point.
        """
        location = validate_and_normalize_location(location)
        if len(location) == 1:
            location = location.item()
        coords = np.asarray(location)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            face_ind = self.face_tree.query(coords=coords, k=1, return_distance=False)

        return face_ind.item()

    def interpolate_node_elevation_from_faces(self) -> None:
        """
        Update node elevations by area-weighted averaging of adjacent face elevations.

        For each node, the elevation is computed as the area-weighted average of the elevations
        of the surrounding faces.

        Returns
        -------
        None
        """
        face_elevation = self.face_elevation
        face_area = self.face_area
        node_face_conn = self.node_face_connectivity

        node_elevation = np.zeros(self.n_node, dtype=np.float64)

        for node_id in range(self.n_node):
            connected_faces = node_face_conn[node_id]
            valid = connected_faces != INT_FILL_VALUE
            faces = connected_faces[valid]

            if faces.size == 0:
                continue

            areas = face_area[faces]
            elevations = face_elevation[faces]

            total_area = np.sum(areas)
            if total_area > 0:
                node_elevation[node_id] = np.sum(elevations * areas) / total_area

        self.node_elevation = node_elevation

    def get_random_location_on_face(self, face_index: int, **kwargs) -> float | tuple[float, float] | ArrayLike:
        """
        Generate a random coordinate within a given face of a the mesh.

        Parameters
        ----------
        grid : uxarray.Grid
            The grid object containing the mesh information.
        face_index : int | NDArray[np.int64]
            The index or array of indices of the face within the grid to obtain the random sample.
        size : int or tuple of ints, optional
            The number of samples to generate. If size is None (the default), a single tuple is returned. If size is greater than 1,
            then an array of tuples is returned. The default is 1.

        Returns
        -------
        (lon,lat) or ndarray[(lon,lat)] of given size
            A pair or array of pairs of longitude and latitude values in degrees.

        Notes
        -----
        This method is a wrapper for :func:`cratermaker.utils.montecarlo_utils.get_random_location_on_face`.
        """
        return get_random_location_on_face(self.uxgrid, face_index, rng=self.rng, **kwargs)

    def elevation_to_cartesian(self, element="face") -> NDArray[np.float64]:
        """
        Convert either the face or node elevations to Cartesian coordinates.

        Parameters
        ----------
        element : str, optional
            The type of element to convert. Can be "face" or "node". Default is "face".

        Returns
        -------
        NDArray[np.float64, 3]
            The Cartesian coordinates of the face elevations.
        """
        return self._full().elevation_to_cartesian(element=element)

    def save(
        self,
        combine_data_files: bool = False,
        interval_number: int = 0,
        time_variables: dict | None = None,
        **kwargs,
    ) -> None:
        """
        Save the surface data to the specified directory. Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the time dimension. If 'interval_number' is included as a key in `time_variables`, then this will be appended to the data file name.

        Parameters
        ----------
        combine_data_files : bool, optional
            If True, combine all data variables into a single NetCDF file, otherwise each variable will be saved to its own NetCDF file. Default is False.
        interval_number : int, optional
            Interval number to append to the data file name. Default is 0.
        time_variables : dict, optional
            Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the time dimension. Default is None.
        """
        do_not_save = ["face_area"]

        self.data_dir.mkdir(parents=True, exist_ok=True)

        if time_variables is None:
            time_variables = {"elapsed_time": float(interval_number)}
        else:
            if not isinstance(time_variables, dict):
                raise TypeError("time_variables must be a dictionary")

        self.uxds.close()

        ds = self.uxds.expand_dims(dim="time").assign_coords({"time": [interval_number]})
        for k, v in time_variables.items():
            ds[k] = xr.DataArray(data=[v], name=k, dims=["time"], coords={"time": [interval_number]})

        drop_vars = [k for k in ds.data_vars if k in do_not_save]
        if len(drop_vars) > 0:
            ds = ds.drop_vars(drop_vars)

        self._save_data(ds, interval_number, combine_data_files)

        save_geometry = interval_number == 0
        self.export(
            format="vtp", interval_number=interval_number, time_variables=time_variables, save_geometry=save_geometry, **kwargs
        )

        return

    def export(self, format="vtp", **kwargs) -> None:
        """
        Export the surface mesh to a file in the specified format. Currently only VTK is supported.
        """
        if format == "vtp" or format == "vtk":
            export.to_vtk(self, **kwargs)
        else:
            raise ValueError(f"Unsupported export format: {format}")

    @staticmethod
    def _sphere_function(coords, x_c, y_c, z_c, r):
        """
        Compute the sphere function.

        Parameters
        ----------
        coords : ndarray
            Array of x, y, and z coordinates.
        x_c, y_c, z_c : float
            Center coordinates of the sphere.
        r : float
            Radius of the sphere.

        Returns
        -------
        ndarray
            The values of the sphere function at the given coordinates.
        """
        x, y, z = coords.T
        return (x - x_c) ** 2 + (y - y_c) ** 2 + (z - z_c) ** 2 - r**2

    def _calculate_distance(
        self,
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: ArrayLike,
        lat2: ArrayLike,
    ) -> NDArray[np.float64]:
        """
        Calculate the great circle distance between one point and one or more other points in meters.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike or ArrayLike
            Longitude of the second point or array of points in radians.
        lat2 : FloatLike or ArrayLike
            Latitude of the second point or array of points in radians.

        Returns
        -------
        NDArray
            Great circle distance between the two points in meters.

        Notes
        -----
        This is a wrapper for a compiled Rust function and is intended to be used as a helper to calculate_face_and_node_distances.
        """
        # Validate that lon1 and lat1 are single points
        if not np.isscalar(lon1):
            if lon1.size != 1:
                raise ValueError("lon1 must be a single point")
            lon1 = lon1.item()
        if not np.isscalar(lat1):
            if lat1.size != 1:
                raise ValueError("lat1 must be a single point")
            lat1 = lat1.item()
        if np.isscalar(lon2):
            lon2 = np.array([lon2])
        if np.isscalar(lat2):
            lat2 = np.array([lat2])

        return surface_functions.calculate_distance(lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2, radius=self.radius)

    @staticmethod
    def _calculate_bearing(
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: FloatLike | ArrayLike,
        lat2: FloatLike | ArrayLike,
    ) -> NDArray[np.float64]:
        """
        Calculate the initial bearing from one point to one or more other points in radians.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike or ArrayLike
            Longitude of the second point or array of points in radians.
        lat2 : FloatLike or ArrayLike
            Latitude of the second point or array of points in radians.

        Returns
        -------
        NDArray
            Initial bearing from the first point to the second point or points in radians.
        """
        return LocalSurface._calculate_bearing(lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2)

    @staticmethod
    def _compute_elevation_to_cartesian(position: NDArray, elevation: NDArray) -> NDArray:
        """
        Convert elevation values to Cartesian coordinates.

        Parameters
        ----------
        position : NDArray
            The position of the points in Cartesian coordinates.
        elevation : NDArray
            The elevation values to convert.

        Returns
        -------
        NDArray
            The Cartesian coordinates of the points with the given elevation.
        """
        runit = position / np.linalg.norm(position, axis=1, keepdims=True)

        return position + elevation[:, np.newaxis] * runit

    def _load_from_files(self, reset: bool = False, **kwargs: Any) -> None:
        """
        Load the grid and data files into the surface object.

        This function loads the grid file and data files from the specified directory. If the grid file does not exist, it will attempt to create a new grid.
        If the data files do not exist, it will create an empty dataset. If reset is True, it will delete all data files except the grid file.

        Parameters
        ----------
        reset : bool, optional
            Flag to indicate whether to reset the surface. Default is False.
        """
        # Get the names of all data files in the data directory that are not the grid file
        regrid = self._regrid_if_needed(**kwargs)
        reset = reset or regrid

        data_file_list = list(self.data_dir.glob("*.nc"))
        if self.grid_file in data_file_list:
            data_file_list.remove(self.grid_file)

        # if data_file_list is empty, set reset to True
        reset = reset or not data_file_list

        # If reset is True, delete all data files except the grid file
        if reset:
            for f in data_file_list:
                f.unlink()
            data_file_list = []

        try:
            with xr.open_dataset(self.grid_file) as uxgrid:
                uxgrid.load()
                if reset:  # Create an empty dataset
                    self._uxds = uxr.UxDataset()
                else:  # Read data from from existing datafiles
                    with uxr.open_mfdataset(uxgrid, data_file_list, use_dual=False) as ds:
                        self._uxds = ds.isel(time=-1).load()
                self._uxds.uxgrid = uxr.Grid.from_dataset(uxgrid)
                self._uxgrid = uxgrid
        except Exception as e:
            raise RuntimeError("Error loading grid and data files") from e

        if reset:
            self.reset(**kwargs)

        return

    def _generate_grid(self, **kwargs: Any) -> None:
        """
        Generate a tessellated mesh of a sphere of based on the particular Surface component that is being used.
        """
        self.grid_file.unlink(missing_ok=True)

        points = self._generate_face_distribution(**kwargs)

        threshold = min(10 ** np.floor(np.log10(self.pix / self.radius)), 1e-6)
        uxgrid = uxr.Grid.from_points(points, method="spherical_voronoi", threshold=threshold)
        uxgrid.attrs["_id"] = self._id
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            uxgrid.to_xarray().to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())

        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name, self.grid_file)
        self._uxgrid = uxgrid

        regrid = self._regrid_if_needed(**kwargs)
        if regrid:
            raise ValueError("Grid file does not match the expected parameters.")
        self._compute_face_size(uxgrid)

        return

    def _regrid_if_needed(self, regrid: bool = False, **kwargs: Any) -> bool:
        """
        Check if the existing grid matches the desired parameters determine if regridding is necessary.

        This function checks if a grid file exists and matches the specified parameters based on a unique hash generated from these
        parameters. If the grid does not exist or does not match the parameters it generates a new grid and returns True.

        Parameters
        ----------
        regrid: bool, optional
            Flag to force regridding even if the grid file exists. Default is False.

        Returns
        -------
        bool
            A boolean indicating whether the grid should be regenerated.
        """
        # Find out if the file exists, if it does't we'll need to make a new grid
        self.data_dir.mkdir(parents=True, exist_ok=True)
        regrid = regrid or not Path(self.grid_file).exists()

        if not regrid:
            try:
                with xr.open_dataset(self.grid_file) as ds:
                    ds.load()
                    uxgrid = uxr.Grid.from_dataset(ds)
                    old_id = uxgrid.attrs.get("_id")
                    regrid = old_id != self._id
            except Exception:
                # Failed to open an old file for whatever reason, so we'll need to regrid
                regrid = True

        if regrid:
            print("Creating a new grid")
            try:
                self._generate_grid(**kwargs)
            except Exception as e:
                raise RuntimeError("Failed to create a new grid") from e

        return regrid

    def _full(self):
        return LocalSurface(self, face_indices=slice(None), node_indices=slice(None), edge_indices=slice(None))

    def _save_data(self, ds: xr.Dataset | xr.DataArray, interval_number: int = 0, combine_data_files: bool = False) -> None:
        """
        Save the data to the specified directory.

        If `combine_data_files` is True, then all data variables are saved to a single NetCDF
        file. If False, then only the data variables for the current interval are saved to a NetCDF file with the interval number
        appended.

        Parameters
        ----------
        ds : xr.Dataset or xr.DataArray
            The data to be saved.
        interval_number : int, Default is 0.
            Interval number to append to the data file name. Default is 0.
        combine_data_files : bool, Default is False.
            If True, combine all data variables into a single NetCDF file, otherwise each variable will be saved to its own NetCDF file.

        Notes
        -----
        This function first saves to a temporary file and then moves that file to the final destination. This is done to avoid file
        locking issues with NetCDF files.
        """
        if isinstance(ds, xr.DataArray):
            ds = ds.to_dataset()

        if "time" not in ds.dims:
            ds = ds.expand_dims(["time"])
        if "time" not in ds.coords:
            ds = ds.assign_coords({"time": [interval_number]})

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as temp_dir:
            if combine_data_files:
                filename = _COMBINED_DATA_FILE_NAME
            else:
                filename = _COMBINED_DATA_FILE_NAME.replace(".nc", f"{interval_number:06d}.nc")

            data_file = self.data_dir / filename
            if data_file.exists():
                with xr.open_mfdataset(data_file) as ds_file:
                    ds_file = ds.merge(ds_file, compat="override")
                    ds_file.load()
            else:
                ds_file = ds

            temp_file = Path(temp_dir) / filename

            comp = {"zlib": True, "complevel": 9}
            encoding = dict.fromkeys(ds_file.data_vars, comp)
            ds_file.to_netcdf(temp_file, encoding=encoding)
            ds_file.close()
            shutil.move(temp_file, data_file)

        return

    def _add_new_data(
        self,
        name: str,
        long_name: str | None = None,
        units: str | None = None,
        data: FloatLike | NDArray | None = None,
        isfacedata: bool = True,
        save_to_file: bool = False,
        interval_number: int = 0,
        combine_data_files: bool = False,
    ) -> None:
        """
        Generate either a node or face data variable and optionally save it to a file. If the data variable already exists, it will be overwritten.

        Parameters
        ----------
        name : str
            Name of the data variable. This will also be used as the data file name.
        long_name : str, optional
            Long name of the data variable that will be saved as an attribute.
        units : str, optional
            Units of the data variable that will be saved as an attribute.
        data : scalar or array-like
            Data file to be saved. If data is a scalar, then the data file will be filled with that value. If data is an array, then the data file will be filled with the array values. The data array must have the same size as the number of faces or nodes in the grid.
        isfacedata : bool, optional
            Flag to indicate whether the data is face data or node data. Default is True.
        save_to_file: bool, optional
            Specify whether the data should be saved to a file. Default is False.
        interval_number : int, optional, default 0
            The interval number to use when saving the data to the data file.
        combine_data_files : bool, optional
            If True, combine the current data with the existing data for previous intervals in the data file. Default is False.

        Returns
        -------
        None
        """
        with xr.open_dataset(self.grid_file) as ds:
            ds.load()
            uxgrid = uxr.Grid.from_dataset(ds)
        if long_name is None and units is None:
            attrs = None
        else:
            attrs = {}
            if long_name:
                attrs["long_name"] = long_name
            if units:
                attrs["units"] = units

        if isfacedata:
            dims = ["n_face"]
            size = uxgrid.n_face
        else:
            dims = ["n_node"]
            size = uxgrid.n_node

        if data is None:
            data = np.zeros(size, dtype=np.float64)
        elif np.isscalar(data):
            data = np.full(size, data)
        else:
            if data.size != size:
                raise ValueError("data must have the same size as the number of faces or nodes in the grid")
        uxda = UxDataArray(
            data=data,
            dims=dims,
            attrs=attrs,
            name=name,
            uxgrid=uxgrid,
        )

        self._uxds[name] = uxda

        if save_to_file:
            self._save_data(uxda, interval_number, combine_data_files)
        return

    @abstractmethod
    def _generate_face_distribution(self, **kwargs: Any) -> tuple[NDArray, NDArray, NDArray]: ...

    def _compute_face_size(self, uxgrid: UxDataset | None = None) -> None:
        """
        Compute the effective pixel size of the mesh based on the face areas.

        Parameters
        ----------
        uxgrid : UxDataset
            The grid object containing the mesh information. If not set, it will be retrieved self.uxds.uxgrid

        """
        if uxgrid is None:
            if self.uxgrid is None:
                return
            else:
                uxgrid = self.uxgrid
        self._face_area = uxgrid.face_areas.values * self.radius**2
        self._face_size = np.sqrt(self._face_area)
        self._pix_mean = float(self._face_size.mean().item())
        self._pix_std = float(self._face_size.std().item())
        self._pix_min = float(self._face_size.min().item())
        self._pix_max = float(self._face_size.max().item())
        return

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self._radius]

    @property
    def _id(self):
        """
        The hash id of the grid. This is used for determining if the grid needs to be regridded.
        """
        combined = ":".join(str(v) for v in self._hashvars)
        hash_object = hashlib.sha256(combined.encode())
        return hash_object.hexdigest()

    @property
    def uxds(self) -> UxDataset:
        """
        The data associated with the surface. This is an instance of UxDataset.
        """
        return self._uxds

    @property
    def uxgrid(self):
        """
        The grid object.
        """
        if self.uxds is not None:
            return self.uxds.uxgrid

    @property
    def data_dir(self):
        """
        Directory for data files.
        """
        return self.simdir / _SURFACE_DIR

    @property
    def gridtype(self):
        """
        The name of the grid type.
        """
        return self._component_name

    @property
    def grid_file(self):
        """
        Path to the grid file.
        """
        return self.simdir / _SURFACE_DIR / _GRID_FILE_NAME

    @property
    def target(self):
        """
        The target body for the impact simulation. Set during initialization.
        """
        return self._target

    @target.setter
    def target(self, value):
        from cratermaker.components.target import Target

        self._target = Target.maker(value)
        return

    @property
    def face_elevation(self) -> NDArray[np.float64]:
        """
        The elevation of the faces.
        """
        return self.uxds["face_elevation"].values

    @face_elevation.setter
    def face_elevation(self, value: NDArray) -> None:
        """
        Set the elevation of the faces.

        Parameters
        ----------
        value : NDArray
            The elevation values to set for the faces.
        """
        value = np.asarray(value, dtype=np.float64)
        if value.size != self.n_face:
            raise ValueError(f"Value must have size {self.n_face}, got {value.size} instead.")
        self.uxds["face_elevation"][:] = value
        return

    @property
    def node_elevation(self) -> NDArray[np.float64]:
        """
        The elevation of the nodes.
        """
        return self.uxds["node_elevation"].values

    @node_elevation.setter
    def node_elevation(self, value: NDArray) -> None:
        """
        Set the elevation of the nodes.

        Parameters
        ----------
        value : NDArray
            The elevation values to set for the nodes.
        """
        value = np.asarray(value, dtype=np.float64)
        if value.size != self.n_node:
            raise ValueError(f"Value must have size {self.n_node}, got {value.size} instead.")
        self.uxds["node_elevation"][:] = value
        return

    @property
    def pix_mean(self) -> float:
        """
        The mean pixel size of the mesh.
        """
        if self._pix_mean is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_mean

    @property
    def pix_std(self) -> float:
        """
        The standard deviation of the pixel size of the mesh.
        """
        if self._pix_std is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_std

    @property
    def pix_min(self) -> float:
        """
        The minimum pixel size of the mesh.
        """
        if self._pix_min is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_min

    @property
    def pix_max(self) -> float:
        """
        The maximum pixel size of the mesh.
        """
        if self._pix_max is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_max

    @property
    def radius(self):
        """
        Radius of the target body.
        """
        return self.target.radius

    @property
    def area(self) -> float:
        """
        Total surface area of the target body.
        """
        if self._area is None:
            self._area = float(self.face_area.sum())
        return self._area

    @property
    def face_area(self) -> NDArray[np.float64]:
        """
        The areas of each face.

        Notes
        -----
        Unlike uxarray.Grid.face_area, this is in meters squared rather than normalized to a unit sphere.

        """
        if self._face_area is None:
            self._compute_face_size()
        return self._face_area

    @property
    def face_size(self) -> NDArray[np.float64]:
        """
        The effective size of each face in meters.

        This is simply the square root of the face area, but is useful for certain comparisons and is equivalent to the `pix` variable from CTEM
        """
        if self._face_size is None:
            self._compute_face_size()
        return self._face_size

    @property
    def smallest_length(self) -> float:
        """
        The smallest length of the mesh.
        """
        if self._smallest_length is None:
            self._smallest_length = float(np.min(self.face_size) * _SMALLFAC)
        return self._smallest_length

    @property
    def face_lat(self) -> NDArray[np.float64]:
        """
        Latitude of the center of each face in degrees.
        """
        return self.uxgrid.face_lat.values

    @property
    def face_lon(self) -> NDArray[np.float64]:
        """
        Longitude of the center of each face in degrees.
        """
        return self.uxgrid.face_lon.values

    @property
    def face_x(self) -> NDArray[np.float64]:
        """
        Cartesian x location of the center of each face in meters.
        """
        if self._face_x is None:
            self._face_x = self.uxgrid.face_x.values * self.radius
        return self._face_x

    @property
    def face_y(self) -> NDArray[np.float64]:
        """
        Cartesian y location of the center of each face in meters.
        """
        if self._face_y is None:
            self._face_y = self.uxgrid.face_y.values * self.radius
        return self._face_y

    @property
    def face_z(self) -> NDArray[np.float64]:
        """
        Cartesian z location of the center of each face in meters.
        """
        if self._face_z is None:
            self._face_z = self.uxgrid.face_z.values * self.radius
        return self._face_z

    @property
    def node_lat(self) -> NDArray[np.float64]:
        """
        Latitude of each node in degrees.
        """
        return self.uxgrid.node_lat.values

    @property
    def node_lon(self) -> NDArray[np.float64]:
        """
        Longitude of each node in degrees.
        """
        return self.uxgrid.node_lon.values

    @property
    def node_x(self) -> NDArray[np.float64]:
        """
        Cartesian x location of each node in meters.
        """
        if self._node_x is None:
            self._node_x = self.uxgrid.node_x.values * self.radius
        return self._node_x

    @property
    def node_y(self) -> NDArray[np.float64]:
        """
        Cartesian y location of each node in meters.
        """
        if self._node_y is None:
            self._node_y = self.uxgrid.node_y.values * self.radius
        return self._node_y

    @property
    def node_z(self) -> NDArray[np.float64]:
        """
        Cartesian z location of each node in meters.
        """
        if self._node_z is None:
            self._node_z = self.uxgrid.node_z.values * self.radius
        return self._node_z

    @property
    def face_indices(self) -> NDArray[np.int64]:
        """
        The indices of the faces of the surface.
        """
        return self.uxds.n_face.values

    @property
    def node_indices(self) -> NDArray[np.int64]:
        """
        The indices of the nodes of the surface.
        """
        return self.uxds.n_node.values

    @property
    def edge_indices(self) -> NDArray[np.int64]:
        """
        The indices of the edges of the surface.
        """
        if self._edge_indices is None:
            self._edge_indices = np.arange(self.uxgrid.n_edge, dtype=np.int64)
        return self._edge_indices

    def _compute_face_bins(self) -> None:
        """
        Compute the face bins based on the face areas. This is used to bin faces by their area for crater generation.
        """
        min_area = self.face_area.min()
        max_area = self.face_area.max()
        max_bin_index = np.ceil(np.log2(max_area / min_area)).astype(int)
        bins = [[] for _ in range(max_bin_index)]

        for face_index, area in enumerate(self.face_area):
            bin_index = np.floor(np.log2(area / min_area)).astype(int)
            bins[bin_index].append(face_index)

        self._face_bin_indices = [np.array(bins[i]) for i in range(max_bin_index) if len(bins[i]) > 0]

        self._face_bin_area = [np.sum(self.face_area[face_indices]) for face_indices in self.face_bin_indices]

        self._face_bin_argmin = [
            int(face_indices[np.argmin(self.face_area[face_indices])]) for face_indices in self._face_bin_indices
        ]

        self._face_bin_argmax = [
            int(face_indices[np.argmax(self.face_area[face_indices])]) for face_indices in self._face_bin_indices
        ]
        return

    @property
    def face_bin_indices(self) -> list[NDArray]:
        """
        Faces binned based on their area.

        All faces within a factor of 2 in area are in the same bin. This property returns a list of face indices lists for each bin.
        The keys are the bin indices, and the values are lists of face indices of faces within that bin.

        This is used when generating craters on surfaces with varying face sizes, so that the smallest crater is sized for the smallest face of a particular bin, rather than for the entire surface.
        """
        if self._face_bin_indices is None:
            self._compute_face_bins()

        return self._face_bin_indices

    @property
    def face_bin_area(self) -> list[float]:
        """
        The total area of all faces in each bin.
        """
        if self._face_bin_area is None:
            self._compute_face_bins()

        return self._face_bin_area

    @property
    def face_bin_argmin(self) -> list[int]:
        """
        The index of the smallest face in each bin.
        """
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return self._face_bin_argmin

    @property
    def face_bin_argmax(self) -> list[int]:
        """
        The index of the largest face in each bin.
        """
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return self._face_bin_argmax

    @property
    def face_bin_min_areas(self) -> list[float]:
        """
        The area of the smallest face in each bin.
        """
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return [float(self.face_area[face_index]) for face_index in self.face_bin_argmin]

    @property
    def face_bin_max_areas(self) -> list[float]:
        """
        The area of the largest face in each bin.
        """
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return [float(self.face_area[face_index]) for face_index in self.face_bin_argmax]

    @property
    def face_bin_min_sizes(self) -> list[float]:
        """
        The effective size of the smallest face in each bin.
        """
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return [float(self.face_size[face_index]) for face_index in self.face_bin_argmin]

    @property
    def face_bin_max_sizes(self) -> list[float]:
        """
        The effective size of the largest face in each bin.
        """
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return [float(self.face_size[face_index]) for face_index in self.face_bin_argmax]

    @property
    def n_face(self) -> int:
        """
        Total number of faces.
        """
        return int(self.uxgrid.n_face)

    @property
    def n_node(self) -> int:
        """
        Total number of nodes.
        """
        return int(self.uxgrid.n_node)

    @property
    def n_edge(self) -> int:
        """
        Total number of edges.
        """
        return int(self.uxgrid.n_edge)

    @property
    def n_nodes_per_face(self) -> NDArray[np.int64]:
        """
        The number of nodes that make up each face.

        Dimensions: `(n_node, )`
        """
        return self.uxgrid.n_nodes_per_face.values

    @property
    def n_max_face_faces(self) -> int:
        """
        The maximum number of faces that surround a face.
        """
        return int(self.uxgrid.n_max_face_faces)

    @property
    def edge_face_distance(self) -> NDArray[np.float64]:
        """
        Distances between the centers of the faces that saddle each edge in meters.
        """
        if self._edge_face_distance is None:
            self._edge_face_distance = surface_functions.compute_edge_distances(
                self.edge_face_connectivity, self.face_lon, self.face_lat, self.radius
            )
        return self._edge_face_distance

    @property
    def edge_length(self) -> NDArray[np.float64]:
        """
        Lengths of each edge in meters.

        Notes
        -----
        This is computed as the distance between the nodes that define the edge
        """
        if self._edge_lengths is None:
            self._edge_lengths = surface_functions.compute_edge_distances(
                self.edge_node_connectivity, self.node_lon, self.node_lat, self.radius
            )
        return self._edge_lengths

    @property
    def face_node_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the nodes that make up each face.

        Dimensions: `(n_face, n_max_face_nodes)`

        Nodes are in counter-clockwise order.
        """
        return self.uxgrid.face_node_connectivity.values

    @property
    def face_face_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the faces that surround each face.

        Dimensions: `(n_face, n_max_face_faces)`
        """
        return self.uxgrid.face_face_connectivity.values

    @property
    def face_edge_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the edges that surround each face.

        Dimensions (n_face, n_max_face_edges)

        """
        return self.uxgrid.face_edge_connectivity.values

    @property
    def node_face_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the faces that surround each node.

        Dimensions: `(n_node, n_max_node_faces)`
        """
        return self.uxgrid.node_face_connectivity.values

    @property
    def edge_face_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the faces that saddle each edge.

        Dimensions (n_edge, 2)

        """
        return self.uxgrid.edge_face_connectivity.values

    @property
    def edge_node_connectivity(self) -> NDArray[np.int64]:
        """
        Indices of the nodes that define each edge.

        Dimensions (n_edge, 2)

        """
        return self.uxgrid.edge_node_connectivity.values

    @property
    def face_tree(self):
        if self._face_tree is None:
            self._face_tree = self.uxgrid.get_ball_tree(
                "face centers",
                distance_metric="haversine",
                coordinate_system="spherical",
                reconstruct=True,
            )

        return self._face_tree

    @property
    def node_tree(self):
        if self._node_tree is None:
            self._node_tree = self.uxgrid.get_ball_tree(
                "nodes",
                distance_metric="haversine",
                coordinate_system="spherical",
                reconstruct=True,
            )

        return self._node_tree

    @property
    def edge_tree(self):
        if self._edge_tree is None:
            self._edge_tree = self.uxgrid.get_ball_tree(
                "edge centers",
                distance_metric="haversine",
                coordinate_system="spherical",
                reconstruct=True,
            )

        return self._edge_tree


class LocalSurface:
    """
    Generates a regional view of a subset of the surface mesh without making copies of any of the data.

    Parameters
    ----------
    surface : Surface
        The surface object that contains the mesh data.
    face_indices : NDArray | slice
        The indices of the faces to include in the view.
    node_indices : NDArray | slice | None, optional
        The indices of the nodes to include in the view. If None, all nodes connected to the faces will be extracted when required
    edge_indices : NDArray | slice | None, optional
        The indices of the edges to include in the view. If None, all edges connected to the faces will be extracted when required.
    location : tuple[float, float] | None, optional
        The location of the center of the view in degrees. This is intended to be passed via the extract_region method of Surface.
    region_radius : FloatLike | None, optional
        The radius of the region to include in the view in meters. This is intended to be passed via the extract_region method of Surface.
    """

    def __init__(
        self,
        surface: Surface,
        face_indices: NDArray | slice,
        node_indices: NDArray | slice | None = None,
        edge_indices: NDArray | slice | None = None,
        location: tuple[float, float] | None = None,
        region_radius: FloatLike | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_surface", Surface.maker(surface))
        object.__setattr__(self, "_face_indices", face_indices)
        object.__setattr__(self, "_node_indices", node_indices)
        object.__setattr__(self, "_edge_indices", edge_indices)
        object.__setattr__(self, "_region_radius", region_radius)
        object.__setattr__(self, "_location", None)
        object.__setattr__(self, "_area", None)
        object.__setattr__(self, "_n_edge", None)
        object.__setattr__(self, "_n_face", None)
        object.__setattr__(self, "_n_node", None)
        object.__setattr__(self, "_face_distance", None)
        object.__setattr__(self, "_node_distance", None)
        object.__setattr__(self, "_face_bearing", None)
        object.__setattr__(self, "_node_bearing", None)
        object.__setattr__(self, "_edge_face_connectivity", None)
        object.__setattr__(self, "_face_edge_connectivity", None)
        object.__setattr__(self, "_face_node_connectivity", None)
        object.__setattr__(self, "_face_face_connectivity", None)
        object.__setattr__(self, "_node_face_connectivity", None)

        if location is not None:
            self._location = validate_and_normalize_location(location)
            self._face_distance, self._node_distance = self.calculate_face_and_node_distances()
            self._face_bearing, self._node_bearing = self.calculate_face_and_node_bearings()
        return

    def __str__(self) -> str:
        """
        String representation of the LocalSurface object.
        """
        base = "<LocalSurface>"
        if self.location:
            base += f"\nLocation: {self.location[0]:.2f}°, {self.location[1]:.2f}°"

        if self.region_radius:
            base += f"\nRegion Radius: {format_large_units(self.region_radius, quantity='length')}"

        return f"{base}\nNumber of faces: {self.n_face}\nNumber of nodes: {self.n_node}"

    def add_data(
        self,
        name: str,
        data: FloatLike | NDArray,
        long_name: str | None = None,
        units: str | None = None,
        isfacedata: bool = True,
        overwrite: bool = False,
    ) -> None:
        """
        Adds new data to the surface.

        Parameters
        ----------
        name : str
            Name of the data variable. This will also be used as the data file name.
        data : scalar or array-like
            Data file to be saved. If data is a scalar, then the data file will be filled with that value. If data is an array, then the data file will be filled with the array values. The data array must have the same size as the number of faces or nodes in the grid.
        long_name : str, optional
            Long name of the data variable that will be saved as an attribute if this is new data. If the data already exists on the surface, this will be ignored.
        units : str, optional
            Units of the data variable that will be saved as an attribute if this is new data. If the data already exists on the surface, this will be ignored.
        isfacedata : bool, optional, default True
            Flag to indicate whether the data is face data or node data. This is only needed if `data` is a scalar, otherwise it is ignored
        overwrite : bool, optional, default False
            By default, new data is added to the old data. This flag indicates that the data should be overwritten, replacing any old data with the new data.

        Returns
        -------
        None
        """
        # Check if the data is a scalar or an array
        if np.isscalar(data):
            n = self.n_face if isfacedata else self.n_node
            data = np.full(n, data)
        elif isinstance(data, list):
            data = np.array(data)
        else:
            data = np.asarray(data)
        if data.size == self.n_face:
            isfacedata = True
            indices = self.face_indices
            dim_name = "n_face"
        elif data.size == self.n_node:
            isfacedata = False
            indices = self.node_indices
            dim_name = "n_node"
        else:
            raise ValueError("data must be a scalar or an array with the same size as the number of faces or nodes in the grid")

        if name not in self.surface.uxds.data_vars:
            self.surface._add_new_data(name, data=0.0, long_name=long_name, units=units, isfacedata=isfacedata)

        # This prevents concurrent writes to the same data variable when used in
        with surface_lock:
            if overwrite:
                self.surface.uxds[name].loc[{dim_name: indices}] = data
            else:
                self.surface.uxds[name].loc[{dim_name: indices}] += data

        return

    def update_elevation(
        self,
        new_elevation: ArrayLike | FloatLike,
        overwrite: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Update the elevation data for the target's surface mesh. This method will determine whether to update the node or face data (or both) depending on the size of the input data. If a scalar is passed, both node and face elevation will be updated to that value. By default, the elevation data will be added to the existing data. If `overwrite` is set to True, the existing data will be replaced with the new data.

        Parameters
        ----------
        new_elevation : ArrayLike | FloatLike
            Elevation to be added (or replaced, if overwrite is True). This can be a scalar, an array with the same size as the number of faces, an array with the same size as the number of nodes, or an array with the same size as the number of faces + nodes.
        overwrite : bool, optional
            If True, the existing data will be replaced with the new data. Default is False.
        **kwargs : Any
            Additional keyword arguments.

        Notes
        -----
        When passing combined data, the first part of the array will be used for face elevation and the second part for node elevation.
        """

        def raise_invalid_elevation_error():
            raise ValueError(
                "new_elev must be None, a scalar, or an array with the same size as the number of nodes, faces, or nodes+faces"
            )

        try:
            new_elevation = np.asarray(new_elevation)

            if np.asarray(new_elevation).size == 1:
                update_node = True
                update_face = True
                new_face_elev = np.full(self.n_face, new_elevation.item())
                new_node_elev = np.full(self.n_node, new_elevation.item())
            else:
                if new_elevation.size == self.n_face:
                    update_face = True
                    update_node = False
                    new_face_elev = new_elevation
                elif new_elevation.size == self.n_node:
                    update_face = False
                    update_node = True
                    new_node_elev = new_elevation
                elif new_elevation.size == self.n_face + self.n_node:
                    update_face = True
                    update_node = True
                    new_face_elev = new_elevation[: self.n_face]
                    new_node_elev = new_elevation[self.n_face :]
                else:
                    raise_invalid_elevation_error()
        except Exception as e:
            raise ValueError("new_elev must be None, a scalar, or an array") from e

        if update_face:
            self.add_data(name="face_elevation", data=new_face_elev, overwrite=overwrite)
        if update_node:
            self.add_data(name="node_elevation", data=new_node_elev, overwrite=overwrite)

        return

    def apply_diffusion(self, kdiff: FloatLike | NDArray) -> None:
        """
        Apply diffusion to the surface.

        Parameters
        ----------
        kdiff : float or array-like
            The degradation state of the surface, which is the product of diffusivity and time. It can be a scalar or an array of the same size as the number of faces in the grid.
            If it is a scalar, the same value is applied to all faces. If it is an array, it must have the same size as the number of faces in the grid.
            The value of kdiff must be greater than 0.0.

        """
        if np.isscalar(kdiff):
            kdiff = np.full(self.n_face, kdiff)
        elif kdiff.size != self.n_face:
            raise ValueError("kdiff must be a scalar or an array with the same size as the number of faces in the grid")
        if np.any(kdiff < 0.0):
            raise ValueError("kdiff must be greater than 0.0")
        kdiff = np.asarray(kdiff, dtype=np.float64)
        kdiffmax = np.max(kdiff)

        if abs(kdiffmax) < _VSMALL:
            return

        delta_face_elevation = surface_functions.apply_diffusion(
            face_kappa=kdiff,
            face_elevation=self.face_elevation,
            face_area=self.face_area,
            edge_face_connectivity=self.edge_face_connectivity,
            edge_face_distance=self.edge_face_distance,
            edge_length=self.edge_length,
        )
        self.update_elevation(delta_face_elevation)
        self.add_data("ejecta_thickness", delta_face_elevation)
        self.interpolate_node_elevation_from_faces()
        return

    def slope_collapse(self, critical_slope_angle: FloatLike = 35.0) -> NDArray:
        """
        Collapse all slopes larger than the critical slope angle.

        Parameters
        ----------
        critical_slope_angle : float
            The critical slope angle (angle of repose) in degrees.
        """
        try:
            critical_slope = np.tan(np.deg2rad(critical_slope_angle))
        except ValueError as e:
            raise ValueError("critical_slope_angle must be between 0 and 90 degrees") from e

        delta_face_elevation = surface_functions.slope_collapse(
            critical_slope=critical_slope,
            face_elevation=self.face_elevation,
            face_area=self.face_area,
            edge_face_connectivity=self.edge_face_connectivity,
            face_edge_connectivity=self.face_edge_connectivity,
            edge_face_distance=self.edge_face_distance,
            edge_length=self.edge_length,
        )
        self.update_elevation(delta_face_elevation)
        self.add_data("ejecta_thickness", delta_face_elevation)
        self.interpolate_node_elevation_from_faces()

    def compute_slope(self) -> NDArray[np.float64]:
        """
        Compute the slope of the surface.

        Returns
        -------
        NDArray[np.float64]
            The slope of all faces in degrees.
        """
        slope = surface_functions.compute_slope(
            face_elevation=self.face_elevation,
            edge_face_connectivity=self.edge_face_connectivity,
            face_edge_connectivity=self.face_edge_connectivity,
            edge_face_distance=self.edge_face_distance,
            edge_length=self.edge_length,
        )

        return np.rad2deg(np.arctan(slope))

    def apply_noise(
        self,
        model: str = "turbulence",
        noise_width: FloatLike = 1000e3,
        noise_height: FloatLike = 1e3,
        **kwargs: Any,
    ) -> None:
        """
        Apply noise to the node elevations of the surface view.

        Parameters
        ----------
        noise_width : float
            The spatial wavelength of the noise.
        noise_height : float
            The amplitude of the noise.
        """
        num_octaves = kwargs.pop("num_octaves", 12)
        anchor = kwargs.pop(
            "anchor",
            self.surface.rng.uniform(0, 2 * np.pi, size=(num_octaves, 3)),
        )
        x = np.concatenate([self.face_x, self.node_x]) / self.surface.radius
        y = np.concatenate([self.face_y, self.node_y]) / self.surface.radius
        z = np.concatenate([self.face_z, self.node_z]) / self.surface.radius

        if model == "turbulence":
            noise = surface_functions.turbulence_noise(
                x=x,
                y=y,
                z=z,
                noise_height=noise_height / self.surface.radius,
                noise_width=noise_width / self.surface.radius,
                freq=2.0,
                pers=0.5,
                anchor=anchor,
                seed=self.surface.rng.integers(0, 2**32 - 1),
            )
        else:
            raise ValueError(f"Unknown noise model: {model}")
        # Compute the weighted mean to ensure volume conservation
        mean = np.sum(noise[: self.n_face] * self.face_area) / self.area
        noise -= mean
        self.update_elevation(noise)
        return

    def calculate_face_and_node_distances(
        self, location: tuple[float, float] | None = None
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Computes the distances between nodes and faces and a given location.

        Parameters
        ----------
        location : tuple[float, float], option
            tuple containing the longitude and latitude of the location in degrees. If None, the location of the view center is used if it is set.

        Returns
        -------
        NDArray
            Array of face distances in meters.
        NDArray
            Array of node distances in meters.
        """
        if location is None:
            if self.location is None:
                raise ValueError("location must be set.")
            location = self.location

        if len(location) == 1:
            location = location.item()
        if len(location) != 2:
            raise ValueError("location must be a single pair of (longitude, latitude).")
        location = validate_and_normalize_location(location)
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.node_lon)
        node_lat2 = np.deg2rad(self.node_lat)
        face_lon2 = np.deg2rad(self.face_lon)
        face_lat2 = np.deg2rad(self.face_lat)
        return self._calculate_distance(lon1, lat1, face_lon2, face_lat2), self._calculate_distance(
            lon1, lat1, node_lon2, node_lat2
        )

    def calculate_face_and_node_bearings(self, location: tuple[float, float] | None = None) -> tuple[NDArray, NDArray]:
        """
        Computes the initial bearing between nodes and faces and a given location.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        NDArray
            Array of initial bearings for each face in radians.
        NDArray
            Array of initial bearings for each node in radians.
        """
        if location is None:
            if self.location is None:
                raise ValueError("location must be set.")
            location = self.location

        if len(location) == 1:
            location = location.item()
        if len(location) != 2:
            raise ValueError("location must be a single pair of (longitude, latitude).")
        location = validate_and_normalize_location(location)
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.node_lon)
        node_lat2 = np.deg2rad(self.node_lat)
        face_lon2 = np.deg2rad(self.face_lon)
        face_lat2 = np.deg2rad(self.face_lat)
        return (
            surface_functions.calculate_bearing(lon1, lat1, face_lon2, face_lat2),
            surface_functions.calculate_bearing(lon1, lat1, node_lon2, node_lat2),
        )

    def interpolate_node_elevation_from_faces(self) -> None:
        """
        Update node elevations by area-weighted averaging of adjacent face elevations.

        For each node, the elevation is computed as the area-weighted average of the elevations
        of the surrounding faces.

        Returns
        -------
        None
        """
        node_elevation = surface_functions.interpolate_node_elevation_from_faces(
            face_area=self.face_area,
            face_elevation=self.face_elevation,
            node_face_connectivity=self.node_face_connectivity,
        )
        self.update_elevation(node_elevation, overwrite=True)
        return

    def get_reference_surface(self, reference_radius: float, **kwargs: Any) -> NDArray[np.float64]:
        """
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.

        Parameters
        ----------
        reference_radius : float
            The radius of the reference region to compute the average over in meters.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        NDArray[np.float64]
            The face and node elevation points of the reference sphere, or the original elevation points is th reference region is too small
        """

        def _find_reference_elevations(region_coords, region_elevation):
            # Perform the curve fitting to get the best fitting spherical cap for the reference surface

            region_surf = self.surface._compute_elevation_to_cartesian(region_coords, region_elevation)

            # Initial guess for the sphere center and radius
            guess_radius = 1.0 + region_elevation.mean()
            initial_guess = [0, 0, 0, guess_radius]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                try:
                    bounds = ([-1.0, -1.0, -1.0, 0.5], [1.0, 1.0, 1.0, 2.0])
                    fit_result, _ = curve_fit(
                        self.surface._sphere_function,
                        region_surf,
                        np.zeros_like(region_elevation),
                        p0=initial_guess,
                        bounds=bounds,
                    )
                except Exception:
                    fit_result = initial_guess

            reference_sphere_center = fit_result[:3]
            reference_sphere_radius = fit_result[3]

            # Find the point along the original vector that intersects the sphere
            f_vec = region_coords / self.surface.radius
            a = f_vec[:, 0] ** 2 + f_vec[:, 1] ** 2 + f_vec[:, 2] ** 2
            b = -2 * (
                f_vec[:, 0] * reference_sphere_center[0]
                + f_vec[:, 1] * reference_sphere_center[1]
                + f_vec[:, 2] * reference_sphere_center[2]
            )
            c = np.dot(reference_sphere_center, reference_sphere_center) - reference_sphere_radius**2
            sqrt_term = b**2 - 4 * a * c
            valid = ~np.isnan(a) & (sqrt_term >= 0.0)

            # Initialize t with default value
            t = np.full_like(a, 1.0)

            # Calculate square root only for valid terms
            sqrt_valid_term = np.sqrt(np.where(valid, sqrt_term, 0))

            # Apply the formula only where valid
            t = np.where(valid, (-b + sqrt_valid_term) / (2 * a), t)
            if np.any(t[valid] < 0):
                t = np.where(valid & (t < 0), (-b - sqrt_valid_term) / (2 * a), t)

            elevations = self.surface.radius * (t * np.linalg.norm(f_vec, axis=1) - 1)
            return elevations

        # Find cells within the crater radius
        faces_within_region = self.face_distance <= reference_radius
        nodes_within_region = self.node_distance <= reference_radius
        points_within_region = np.concatenate([faces_within_region, nodes_within_region])

        elevation = np.concatenate([self.face_elevation, self.node_elevation])

        n_reference = np.sum(faces_within_region) + np.sum(nodes_within_region)

        # if there are not enough faces or nodes within the radius, return the original elevations and use that as the reference
        if n_reference < 5 or not nodes_within_region.any():
            return elevation

        combo_grid = np.column_stack(
            (
                np.concatenate([self.face_x, self.node_x]),
                np.concatenate([self.face_y, self.node_y]),
                np.concatenate([self.face_z, self.node_z]),
            )
        )
        # find_refence_elevations expects coordinates to be on the unit sphere, so we need to normalize
        region_coords = combo_grid[points_within_region] / self.surface.radius
        region_elevation = elevation[points_within_region] / self.surface.radius

        reference_elevation = elevation
        reference_elevation[points_within_region] = _find_reference_elevations(region_coords, region_elevation)

        return reference_elevation

    def extract_subregion(self, subregion_radius: FloatLike):
        """
        Extract a subset of the LocalSurface region with a smaller radius than the original region.

        Parameters
        ----------
        subregion_radius : float
            The radius of the subregion to extract in meters.

        Returns
        -------
        LocalSurface
            A LocalSurface object containing a view of the regional grid.

        """
        if subregion_radius > self.region_radius:
            raise ValueError("subregion_radius must be smaller than the original region radius")

        if subregion_radius < self.region_radius:
            if isinstance(self.face_indices, slice):
                face_indices = np.arange(self.n_face)[self.face_indices][self.face_distance <= subregion_radius]
            else:
                face_indices = self.face_indices[self.face_distance <= subregion_radius]
            if face_indices.size == 0:
                return None

            # First select edges and nodes that are attached to these faces
            edge_indices = np.unique(self.surface.face_edge_connectivity[face_indices].ravel())
            edge_indices = edge_indices[edge_indices != INT_FILL_VALUE]

            node_indices = np.unique(self.surface.face_node_connectivity[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]

            # Now add in all faces that are connected to anything inside the region, so that the outermost border of the local region has a buffer of faces
            # These are needed for diffusion calculations
            neighbor_faces = self.surface.face_face_connectivity[face_indices]
            neighbor_faces = neighbor_faces[neighbor_faces != INT_FILL_VALUE]
            node_faces = self.surface.node_face_connectivity[node_indices]
            node_faces = node_faces[node_faces != INT_FILL_VALUE]
            edge_faces = self.surface.edge_face_connectivity[edge_indices]
            edge_faces = edge_faces[edge_faces != INT_FILL_VALUE]
            face_indices = np.unique(np.concatenate((face_indices, neighbor_faces, node_faces, edge_faces)))
        else:
            face_indices = self.face_indices
            edge_indices = self.edge_indices
            node_indices = self.node_indices

        return LocalSurface(
            surface=self.surface,
            face_indices=face_indices,
            node_indices=node_indices,
            edge_indices=edge_indices,
            location=self.location,
            region_radius=subregion_radius,
        )

    def compute_volume(self, elevation: NDArray) -> NDArray:
        """
        Compute the volume of an array of elevation points.

        Parameters
        ----------
        elevation : NDArray
            The elevation values to compute. This should be the same size as the number of faces in the grid.

        Returns
        -------
        float
            The volume of the elevation points
        """
        if elevation.size != self.n_face:
            raise ValueError("elevation must be an array with the same size as the number of faces in the grid")
        return np.sum(elevation * self.face_area)

    def elevation_to_cartesian(self, element="face") -> NDArray[np.float64]:
        """
        Convert either the face or node elevations to Cartesian coordinates.

        Parameters
        ----------
        element : str, optional
            The type of element to convert. Can be "face" or "node". Default is "face".

        Returns
        -------
        NDArray[np.float64, 3]
            The Cartesian coordinates of the face elevations.
        """
        if element not in ("face", "node"):
            raise ValueError("element must be either 'face' or 'node'")
        if element == "face":
            position = np.column_stack((self.face_x, self.face_y, self.face_z))
            elevation = self.face_elevation
        elif element == "node":
            position = np.column_stack((self.node_x, self.node_y, self.node_z))
            elevation = self.node_elevation

        return self.surface._compute_elevation_to_cartesian(position, elevation)

    def _calculate_distance(
        self,
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: FloatLike | ArrayLike,
        lat2: FloatLike | ArrayLike,
    ) -> NDArray[np.float64]:
        """
        Calculate the great circle distance between one point and one or more other points in meters.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike or ArrayLike
            Longitude of the second point or array of points in radians.
        lat2 : FloatLike or ArrayLike
            Latitude of the second point or array of points in radians.

        Returns
        -------
        NDArray
            Great circle distance between the two points in meters.

        Notes
        -----
        This is a wrapper for a compiled Rust function and is intended to be used as a helper to calculate_face_and_node_distances.
        """
        return self.surface._calculate_distance(lon1, lat1, lon2, lat2)

    @staticmethod
    def _calculate_bearing(
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: FloatLike | ArrayLike,
        lat2: FloatLike | ArrayLike,
    ) -> NDArray[np.float64]:
        """
        Calculate the initial bearing from one point to one or more other points in radians.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike or ArrayLike
            Longitude of the second point or array of points in radians.
        lat2 : FloatLike or ArrayLike
            Latitude of the second point or array of points in radians.

        Returns
        -------
        NDArray
            Initial bearing from the first point to the second point or points in radians.

        Notes
        -----
        This is intended to be used as a helper to calculate_face_and_node_bearings.
        """
        # Calculate differences in coordinates
        dlon = np.mod(lon2 - lon1 + np.pi, 2 * np.pi) - np.pi

        # Haversine formula calculations
        x = np.sin(dlon) * np.cos(lat2)
        y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
        initial_bearing = np.arctan2(x, y)

        # Normalize bearing to 0 to 2*pi
        initial_bearing = (initial_bearing + 2 * np.pi) % (2 * np.pi)

        return initial_bearing

    @staticmethod
    def _remap_connectivity_to_local(
        connectivity_array: np.ndarray,
        row_indices: np.ndarray,
        value_indices: np.ndarray,
        total_value_size: int,
        fill_value: int = INT_FILL_VALUE,
    ) -> np.ndarray:
        """
        Remap a 2D connectivity array where the rows are indexed by one entity (e.g., edges) and the values are global references to another entity (e.g., faces).

        Only rows in `row_indices` are retained, and their values are remapped from global to local indices using `value_indices`.

        Parameters
        ----------
        connectivity_array : np.ndarray
            The global 2D connectivity array of shape (N, K).
        row_indices : np.ndarray
            Global indices indicating which rows to retain.
        value_indices : np.ndarray
            Global indices of values to remap (e.g., face indices).
        total_value_size : int
            Total number of global items in the value domain (e.g., total number of faces).
        fill_value : int, optional
            Fill value indicating invalid entries (default is INT_FILL_VALUE).

        Returns
        -------
        np.ndarray
            A new 2D connectivity array of shape (len(row_indices), K) with local indices or fill_value where not applicable.
        """
        if isinstance(row_indices, slice) and row_indices == slice(None):  # This is a global slice. Just return the original array
            return connectivity_array
        if isinstance(value_indices, slice):
            value_indices = np.arange(total_value_size)[value_indices]

        # Build global-to-local index mapping
        mapping = np.full(total_value_size, fill_value, dtype=np.int64)
        mapping[value_indices] = np.arange(value_indices.size)

        subset = connectivity_array[row_indices]
        remapped = np.full_like(subset, fill_value)
        valid = subset != fill_value
        remapped[valid] = mapping[subset[valid]]

        return remapped

    @property
    def surface(self) -> Surface:
        """
        The surface object that contains the mesh data.
        """
        return self._surface

    @property
    def n_edge(self) -> int:
        """
        The number of edges in the view.
        """
        if self._n_edge is None:
            if isinstance(self._edge_indices, slice):
                self._n_edge = int(self._surface._uxds.uxgrid.n_edge[self._edge_indices].size)
            else:
                self._n_edge = self._edge_indices.size
        return self._n_edge

    @property
    def n_face(self) -> int:
        """
        The number of faces in the view.
        """
        if self._n_face is None:
            if isinstance(self.face_indices, slice):
                self._n_face = int(self.surface.face_elevation[self.face_indices].size)
            elif isinstance(self.face_indices, np.ndarray):
                self._n_face = int(self.face_indices.size)

        return self._n_face

    @property
    def n_node(self) -> int:
        """
        The number of nodes in the view.
        """
        if self._n_node is None:
            if isinstance(self.node_indices, slice):
                self._n_node = int(self.surface.node_elevation[self.node_indices].size)
            elif isinstance(self.node_indices, np.ndarray):
                self._n_node = int(self.node_indices.size)

        return self._n_node

    @property
    def face_elevation(self) -> NDArray:
        """
        The elevation of the faces.
        """
        return self.surface.face_elevation[self.face_indices]

    @face_elevation.setter
    def face_elevation(self, value: NDArray) -> None:
        """
        Set the elevation of the faces.

        Parameters
        ----------
        value : NDArray
            The elevation values to set for the faces.
        """
        if value.size != self.n_face:
            raise ValueError(f"Value must have size {self.n_face}, got {value.size} instead.")
        self.surface.face_elevation[self.face_indices] = value
        return

    @property
    def node_elevation(self) -> NDArray:
        """
        The elevation of the nodes.
        """
        return self.surface.node_elevation[self.node_indices]

    @node_elevation.setter
    def node_elevation(self, value: NDArray) -> None:
        """
        Set the elevation of the nodes.

        Parameters
        ----------
        value : NDArray
            The elevation values to set for the nodes.
        """
        if value.size != self.n_node:
            raise ValueError(f"Value must have size {self.n_node}, got {value.size} instead.")
        self.surface.node_elevation[self.node_indices] = value
        return

    @property
    def location(self) -> tuple[float, float]:
        """
        The location of the center of the view.
        """
        return self._location

    @property
    def region_radius(self) -> FloatLike:
        """
        The radius of the region to include in the view in meters.
        """
        return self._region_radius

    @property
    def face_bearing(self) -> NDArray:
        """
        The initial bearing from the location to the faces.
        """
        return self._face_bearing

    @property
    def node_bearing(self) -> NDArray:
        """
        The initial bearing from the location to the nodes.
        """
        return self._node_bearing

    @property
    def face_distance(self) -> NDArray:
        """
        The distance from the location to the faces.
        """
        return self._face_distance

    @property
    def node_distance(self) -> NDArray:
        """
        The distance from the location to the nodes.
        """
        return self._node_distance

    @property
    def face_area(self) -> NDArray:
        """
        The areas of the faces.
        """
        return self.surface.face_area[self.face_indices]

    @property
    def face_lat(self) -> NDArray:
        """
        Latitude of the center of the faces in degrees.
        """
        return self.surface.face_lat[self.face_indices]

    @property
    def face_lon(self) -> NDArray:
        """
        Longitude of the center of the faces in degrees.
        """
        return self.surface.face_lon[self.face_indices]

    @property
    def face_x(self) -> NDArray:
        """
        Cartesian x location of the center of the faces in meters.
        """
        return self.surface.face_x[self.face_indices]

    @property
    def face_y(self) -> NDArray:
        """
        Cartesian y location of the center of the faces in meters.
        """
        return self.surface.face_y[self.face_indices]

    @property
    def face_z(self) -> NDArray:
        """
        Cartesian z location of the center of the faces in meters.
        """
        return self.surface.face_z[self.face_indices]

    @property
    def edge_face_distance(self) -> NDArray:
        """
        Distances between the edges and the faces.

        Dimensions: `(n_edge)`
        """
        return self.surface.edge_face_distance[self.edge_indices]

    @property
    def edge_length(self) -> NDArray:
        """
        Lengths of the edges in meters.

        Dimensions: `(n_edge)`
        """
        return self.surface.edge_length[self.edge_indices]

    @property
    def node_lat(self) -> NDArray:
        """
        Latitude of the nodes in degrees.
        """
        return self.surface.node_lat[self.node_indices]

    @property
    def node_lon(self) -> NDArray:
        """
        Longitude of the nodes in degrees.
        """
        return self.surface.node_lon[self.node_indices]

    @property
    def node_x(self) -> NDArray:
        """
        Cartesian x location of the nodes in meters.
        """
        return self.surface.node_x[self.node_indices]

    @property
    def node_y(self) -> NDArray:
        """
        Cartesian y location of the nodes in meters.
        """
        return self.surface.node_y[self.node_indices]

    @property
    def node_z(self) -> NDArray:
        """
        Cartesian z location of the nodes in meters.
        """
        return self.surface.node_z[self.node_indices]

    @property
    def area(self) -> float:
        """
        The total area of the faces in the view.
        """
        if self._area is None:
            self._area = float(self.face_area.sum())
        return self._area

    @property
    def face_indices(self) -> NDArray:
        """
        The indices of the faces in the view.
        """
        if self._face_indices is None:
            raise ValueError("face_indices must be set to use this object.")
        return self._face_indices

    @property
    def node_indices(self) -> NDArray:
        """
        The indices of the nodes in the view.
        """
        if self._node_indices is None:
            self._node_indices = np.unique(self.surface.face_node_connectivity[self.face_indices].ravel())
            self._node_indices = self._node_indices[self._node_indices != INT_FILL_VALUE]

        return self._node_indices

    @property
    def edge_indices(self) -> NDArray:
        """
        The indices of the edges in the view.
        """
        if self._edge_indices is None:
            self._edge_indices = np.unique(self.surface.face_edge_connectivity[self.face_indices].ravel())
            self._edge_indices = self._edge_indices[self._edge_indices != INT_FILL_VALUE]
        return self._edge_indices

    @property
    def edge_face_connectivity(self) -> NDArray:
        """
        Local indices of the faces that make up the edges.

        Dimensions: `(n_edge, 2)`
        """
        if self._edge_face_connectivity is None:
            self._edge_face_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.edge_face_connectivity,
                row_indices=self.edge_indices,
                value_indices=self.face_indices,
                total_value_size=self.surface.n_face,
            )
        return self._edge_face_connectivity

    @property
    def face_edge_connectivity(self) -> NDArray:
        """
        Local indices of the edges that make up the faces.

        Dimensions: `(n_face, n_max_face_edges)`
        """
        if self._face_edge_connectivity is None:
            self._face_edge_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.face_edge_connectivity,
                row_indices=self.face_indices,
                value_indices=self.edge_indices,
                total_value_size=self.surface.n_edge,
            )
        return self._face_edge_connectivity

    @property
    def face_node_connectivity(self) -> NDArray:
        """
        Local indices of the nodes that make up the faces.

        Dimensions: `(n_face, n_max_face_nodes)`

        Nodes are in counter-clockwise order.
        """
        if self._face_node_connectivity is None:
            self._face_node_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.face_node_connectivity,
                row_indices=self.face_indices,
                value_indices=self.node_indices,
                total_value_size=self.surface.n_node,
            )
        return self._face_node_connectivity

    @property
    def face_face_connectivity(self) -> NDArray:
        """
        Indices of the faces that surround the faces.

        Dimensions: `(n_face, n_max_face_faces)`
        """
        if self._face_face_connectivity is None:
            self._face_face_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.face_face_connectivity,
                row_indices=self.face_indices,
                value_indices=self.face_indices,
                total_value_size=self.surface.n_face,
            )
        return self._face_face_connectivity

    @property
    def node_face_connectivity(self) -> NDArray:
        """
        Indices of the faces that surround the nodes.

        Dimensions: `(n_node, n_max_node_faces)`
        """
        if self._node_face_connectivity is None:
            self._node_face_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.node_face_connectivity,
                row_indices=self.node_indices,
                value_indices=self.face_indices,
                total_value_size=self.surface.n_face,
            )
        return self._node_face_connectivity

    @property
    def edge_node_connectivity(self) -> NDArray:
        """
        Local indices of the nodes that make up the edges.

        Dimensions: `(n_edge, 2)`
        """
        if self._edge_node_connectivity is None:
            self._edge_node_connectivity = self._remap_connectivity_to_local(
                connectivity_array=self.surface.edge_node_connectivity,
                row_indices=self.edge_indices,
                value_indices=self.node_indices,
                total_value_size=self.surface.n_node,
            )
        return self._edge_node_connectivity


import_components(__name__, __path__, ignore_private=True)
