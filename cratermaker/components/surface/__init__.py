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
    PairOfFloats,
)
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import validate_and_normalize_location
from cratermaker.utils.montecarlo_utils import get_random_location_on_face

if TYPE_CHECKING:
    from cratermaker.components.target import Target

surface_lock = threading.Lock()


class Surface(ComponentBase):
    _registry: dict[str, type[Surface]] = {}

    """
    This class is used for handling surface-related data and operations in the cratermaker project. It provides methods for 
    setting elevation data, calculating distances and bearings, and other surface-related computations.
    
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
        object.__setattr__(self, "_pix_mean", None)
        object.__setattr__(self, "_pix_std", None)
        object.__setattr__(self, "_pix_min", None)
        object.__setattr__(self, "_pix_max", None)
        object.__setattr__(self, "_area", None)
        object.__setattr__(self, "_node_tree", None)
        object.__setattr__(self, "_face_tree", None)
        object.__setattr__(self, "_face_areas", None)
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
        return f"{base}\nTarget: {self.target.name}\nGrid File: {self.grid_file}"

    def __del__(self):
        try:
            if hasattr(self, "_uxds") and hasattr(self._uxds, "close"):
                if hasattr(self._uxds, "uxgrid") and hasattr(self._uxds.uxgrid, "_ds"):
                    self._uxds.uxgrid._ds.close()
                self._uxds.close()
        except Exception:
            pass
        try:
            if (
                hasattr(self, "_grid")
                and hasattr(self._grid, "uxgrid")
                and hasattr(self._grid.uxgrid, "_ds")
            ):
                self._grid.uxgrid._ds.close()
        except Exception:
            pass

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

    def load_from_files(self, reset: bool = False, **kwargs: Any) -> None:
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
        regrid = self.regrid_if_needed(**kwargs)
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
                    with uxr.open_mfdataset(
                        uxgrid, data_file_list, use_dual=False
                    ) as ds:
                        self._uxds = ds.isel(time=-1).load()
                self._uxds.uxgrid = uxr.Grid.from_dataset(uxgrid)
                self._uxgrid = uxgrid
        except Exception as e:
            raise RuntimeError("Error loading grid and data files") from e

        if reset:
            self.reset(**kwargs)

        return

    def save_to_files(
        self,
        combine_data_files: bool = False,
        interval_number: int = 0,
        time_variables: dict | None = None,
        *args,
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
        do_not_save = ["face_areas"]

        self.data_dir.mkdir(parents=True, exist_ok=True)

        if time_variables is None:
            time_variables = {"elapsed_time": float(interval_number)}
        else:
            if not isinstance(time_variables, dict):
                raise TypeError("time_variables must be a dictionary")

        self.uxds.close()

        ds = self.uxds.expand_dims(dim="time").assign_coords(
            {"time": [interval_number]}
        )
        for k, v in time_variables.items():
            ds[k] = xr.DataArray(
                data=[v], name=k, dims=["time"], coords={"time": [interval_number]}
            )

        drop_vars = [k for k in ds.data_vars if k in do_not_save]
        if len(drop_vars) > 0:
            ds = ds.drop_vars(drop_vars)

        self._save_data(ds, interval_number, combine_data_files)

        return

    def generate_grid(self, **kwargs: Any) -> None:
        """
        Generate a tessellated mesh of a sphere of based on the particular Surface component that is being used.
        """
        self.grid_file.unlink(missing_ok=True)

        points = self.generate_face_distribution(**kwargs)
        uxgrid = uxr.Grid.from_points(points, method="spherical_voronoi")
        uxgrid.attrs["_id"] = self._id
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            uxgrid.to_xarray().to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())

        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name, self.grid_file)
        self._uxgrid = uxgrid

        regrid = self.regrid_if_needed(**kwargs)
        assert not regrid

        self._pix_mean, self._pix_std, self._pix_min, self._pix_max = (
            self._compute_pix_size(uxgrid)
        )

        return

    def regrid_if_needed(self, regrid: bool = False, **kwargs: Any) -> bool:
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
                self.generate_grid(**kwargs)
            except Exception as e:
                raise RuntimeError("Failed to create a new grid") from e

        return regrid

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

    def full_view(self):
        return SurfaceView(self, slice(None), slice(None))

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
        return self.full_view().apply_diffusion(kdiff)

    def slope_collapse(self, critical_slope_angle: FloatLike = 35.0) -> None:
        """
        Collapse all slopes larger than the critical slope angle.

        Parameters
        ----------
        critical_slope_angle : float
            The critical slope angle (angle of repose) in degrees.
        """
        return self.full_view().slope_collapse(critical_slope_angle)

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

        return self.full_view().apply_noise(
            model=model, noise_width=noise_width, noise_height=noise_height, **kwargs
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
        face_elevation = self.face_elevation
        face_areas = self.face_areas
        node_face_conn = self.node_face_connectivity

        node_elevation = np.zeros(self.n_node, dtype=np.float64)

        for node_id in range(self.n_node):
            connected_faces = node_face_conn[node_id]
            valid = connected_faces != INT_FILL_VALUE
            faces = connected_faces[valid]

            if faces.size == 0:
                continue

            areas = face_areas[faces]
            elevations = face_elevation[faces]

            total_area = np.sum(areas)
            if total_area > 0:
                node_elevation[node_id] = np.sum(elevations * areas) / total_area

        self.node_elevation = node_elevation

    def get_distance(self, location: tuple[float, float]) -> tuple[NDArray, NDArray]:
        """
        Computes the distances between nodes and faces and a given location.

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
        return self.full_view().get_distance(location)

    def get_initial_bearing(
        self, location: tuple[float, float]
    ) -> tuple[NDArray, NDArray]:
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
        return self.full_view().get_initial_bearing(location)

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
                raise ValueError(
                    "data must have the same size as the number of faces or nodes in the grid"
                )
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
        Adds new data

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
        return self.full_view().add_data(
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

        return self.full_view().update_elevation(
            new_elevation=new_elevation, overwrite=overwrite, **kwargs
        )

    def calculate_haversine_distance(
        self,
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: FloatLike,
        lat2: FloatLike,
    ) -> float:
        """
        Calculate the great circle distance between two points on a sphere.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike
            Longitude of the second point in radians.
        lat2 : FloatLike
            Latitude of the second point in radians.

        Returns
        -------
        float
            Great circle distance between the two points in meters.
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

        return surface_functions.calculate_haversine_distance(
            lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2, radius=self.radius
        )

    @staticmethod
    def calculate_initial_bearing(
        lon1: FloatLike, lat1: FloatLike, lon2: FloatLike, lat2: FloatLike
    ) -> float:
        """
        Calculate the initial bearing from one point to another on the surface of a sphere.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike
            Longitude of the second point in radians.
        lat2 : FloatLike
            Latitude of the second point in radians.

        Returns
        -------
        float
            Initial bearing from the first point to the second point in radians.
        """

        return SurfaceView.calculate_initial_bearing(
            lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2
        )

    def find_nearest_index(self, location):
        """
        Find the index of the nearest node and face to a given point.

        This method calculates the Haversine distance from the given point to each face in the grid,
        and returns the index of the face with the minimum distance.

        Parameters
        ----------
        location : tuple
            A tuple containing two elements: (longitude, latitude) in degrees.

        Returns
        -------
        int
            The index of the nearest face in the grid to the given point.
        int
            The index of the nearest node in the grid to the given point.

        Notes
        -----
        The method uses the ball tree query method that is included in the UxArray.Grid class.
        """

        if len(location) == 1:
            location = location.item()
        coords = np.asarray(location)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", Warning)
            node_ind = self.node_tree.query(coords=coords, k=1, return_distance=False)
            face_ind = self.face_tree.query(coords=coords, k=1, return_distance=False)

        return face_ind.item(), node_ind.item()

    @staticmethod
    def elevation_to_cartesian(position: NDArray, elevation: NDArray) -> NDArray:
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
        return SurfaceView.elevation_to_cartesian(position, elevation)

    def extract_region(
        self, location: tuple[FloatLike, FloatLike], region_radius: FloatLike
    ):
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
        SurfaceView
            A SurfaceView object containing a view of the regional grid.

        """

        region_angle = np.rad2deg(region_radius / self.radius)
        if len(location) == 1:
            location = location.item()
        coords = np.asarray(location)

        face_indices = self.face_tree.query_radius(coords, region_angle)
        if len(face_indices) == 0:
            return None

        return SurfaceView(surface=self, face_indices=face_indices, location=location)

    def get_random_location_on_face(
        self, face_index: int, **kwargs
    ) -> float | tuple[float, float] | ArrayLike:
        """
        Generate a random coordinate within a given face of an ungridtyped mesh.

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

        return get_random_location_on_face(
            self.uxgrid, face_index, rng=self.rng, **kwargs
        )

    def _save_data(
        self,
        ds: xr.Dataset | xr.DataArray,
        interval_number: int = 0,
        combine_data_files: bool = False,
    ) -> None:
        """
        Save the data to the specified directory. If `combine_data_files` is True, then all data variables are saved to a single NetCDF
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
                filename = _COMBINED_DATA_FILE_NAME.replace(
                    ".nc", f"{interval_number:06d}.nc"
                )

            data_file = self.data_dir / filename
            if data_file.exists():
                with xr.open_mfdataset(data_file) as ds_file:
                    ds_file = ds.merge(ds_file, compat="override")
                    ds_file.load()
            else:
                ds_file = ds

            temp_file = Path(temp_dir) / filename

            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in ds_file.data_vars}
            ds_file.to_netcdf(temp_file, encoding=encoding)
            ds_file.close()
            shutil.move(temp_file, data_file)

        return

    @property
    def uxds(self) -> UxDataset:
        """
        The data associated with the surface. This is an instance of UxDataset.
        """
        return self._uxds

    @property
    def data_dir(self):
        """
        Directory for data files.
        """
        return self.simdir / _SURFACE_DIR

    @property
    def grid_file(self):
        """
        Path to the grid file.
        """
        return self.simdir / _SURFACE_DIR / _GRID_FILE_NAME

    @property
    def area(self):
        """
        Total surface area of the target body.
        """
        if self._area is None:
            self._area = self.face_areas.sum()
        return self._area

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
    def radius(self):
        """
        Radius of the target body.
        """
        return self.target.radius

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
    def face_tree(self):
        if self._face_tree is None:
            self._face_tree = self.uxgrid.get_ball_tree(
                "face centers",
                distance_metric="haversine",
                coordinate_system="spherical",
                reconstruct=True,
            )

        return self._face_tree

    @abstractmethod
    def generate_face_distribution(
        self, **kwargs: Any
    ) -> tuple[NDArray, NDArray, NDArray]: ...

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

    @staticmethod
    def _distribute_points(
        distance: FloatLike,
        radius: FloatLike = 1.0,
        lon_range: PairOfFloats = (-180, 180),
        lat_range: PairOfFloats = (-90, 90),
    ) -> NDArray:
        """
        Distributes points on a sphere using Deserno's algorithm (Deserno 2004).

        Parameters
        ----------
        distance : float
            Approximate distance between points, used to determine the number of points, where n = 1/distance**2 when distributed over the whole sphere.
        radius : float, optional
            Radius of the sphere. Default is 1.0
        lon_range : tuple, optional
            Range of longitudes in degrees. Default is (-180,180).
        lat_range : tuple, optional
            Range of latitudes in degrees. Default is (-90,90).

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of cartesian points on the sphere.

        References
        ----------
        - Deserno, Markus., 2004. How to generate equidistributed points on the surface of a sphere. https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf

        """

        def _sph2cart(theta, phi, r):
            """
            Converts spherical coordinates to Cartesian coordinates.

            Parameters
            ----------
            theta : float
                Inclination angle in radians.
            phi : float
                Azimuthal angle in radians.
            r : float
                Radius.
            """
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            return x, y, z

        phi_range = np.deg2rad(lon_range) + np.pi
        theta_range = np.deg2rad(lat_range) + np.pi / 2
        points = []

        n = int(1 / distance**2)
        if n < 1:
            return

        a = 4 * np.pi / n
        d = np.sqrt(a)
        Mtheta = int(np.round(np.pi / d))
        dtheta = np.pi / Mtheta
        dphi = a / dtheta

        thetavals = np.pi * (np.arange(Mtheta) + 0.5) / Mtheta
        thetavals = thetavals[
            (thetavals >= theta_range[0]) & (thetavals < theta_range[1])
        ]

        for theta in thetavals:
            Mphi = int(np.round(2 * np.pi * np.sin(theta) / dphi))
            phivals = 2 * np.pi * np.arange(Mphi) / Mphi
            phivals = phivals[(phivals >= phi_range[0]) & (phivals < phi_range[1])]
            for phi in phivals:
                points.append(_sph2cart(theta, phi, radius))
        if len(points) == 0:
            return

        points = np.array(points, dtype=np.float64)
        points = points.T

        return points

    @property
    def uxgrid(self):
        """
        The grid object.
        """
        if self.uxds is not None:
            return self.uxds.uxgrid

    def _compute_pix_size(self, uxgrid: UxDataset | None = None) -> tuple[float, float]:
        """
        Compute the effective pixel size of the mesh based on the face areas.

        Parameters
        ----------
        uxgrid : UxDataset
            The grid object containing the mesh information. If not set, it will be retrieved self.uxds.uxgrid

        Returns
        -------
        tuple[float, float, float, float]
            The mean, standard deviation, minimum, and maximum of the pixel size in meters.
        """
        if uxgrid is None:
            if self.uxgrid is None:
                return None, None, None, None
            else:
                uxgrid = self.uxgrid
        face_areas = uxgrid.face_areas
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        pix_min = face_sizes.min().item() * self.radius
        pix_max = face_sizes.max().item() * self.radius
        return float(pix_mean), float(pix_std), float(pix_min), float(pix_max)

    @property
    def pix_mean(self):
        """
        The mean pixel size of the mesh.
        """
        if self._pix_mean is None and self.uxgrid is not None:
            self._pix_mean, self._pix_std, self._pix_min, self._pix_max = (
                self._compute_pix_size()
            )
        return self._pix_mean

    @property
    def pix_std(self):
        """
        The standard deviation of the pixel size of the mesh.
        """
        if self._pix_std is None and self.uxgrid is not None:
            self._pix_mean, self._pix_std, self._pix_min = self._compute_pix_size()
        return self._pix_std

    @property
    def pix_min(self):
        """
        The minimum pixel size of the mesh.
        """
        if self._pix_min is None and self.uxgrid is not None:
            self._pix_mean, self._pix_std, self._pix_min = self._compute_pix_size()
        return self._pix_min

    @property
    def pix_max(self):
        """
        The maximum pixel size of the mesh.
        """
        if self._pix_max is None and self.uxgrid is not None:
            self._pix_mean, self._pix_std, self._pix_min, self._pix_max = (
                self._compute_pix_size()
            )
        return self._pix_max

    @property
    def gridtype(self):
        """
        The name of the grid type.
        """
        return self._component_name

    @property
    def smallest_length(self):
        """
        The smallest length of the mesh.
        """
        if self._smallest_length is None:
            self._smallest_length = np.sqrt(np.min(self.face_areas)) * _SMALLFAC
        return self._smallest_length

    @property
    def n_face(self):
        """
        Total number of faces
        """
        return self.uxgrid.n_face

    @property
    def n_max_face_faces(self):
        """
        The maximum number of faces that surround a face.
        """
        return self.uxgrid.n_max_face_faces

    @property
    def face_areas(self):
        """
        The areas of each face.

        Notes
        -----
        Unlike uxarray.Grid.face_areas, this is in meters squared.

        """
        if self._face_areas is None:
            self._face_areas = self.uxgrid.face_areas.values * self.radius**2
        return self._face_areas

    @property
    def face_lat(self):
        """
        Latitude of the center of each face in degrees.
        """
        return self.uxgrid.face_lat.values

    @property
    def face_lon(self):
        """
        Longitude of the center of each face in degrees.
        """
        return self.uxgrid.face_lon.values

    @property
    def face_node_connectivity(self):
        """
        Indices of the nodes that make up each face.

        Dimensions: `(n_face, n_max_face_nodes)`

        Nodes are in counter-clockwise order.
        """
        return self.uxgrid.face_node_connectivity.values

    @property
    def face_face_connectivity(self):
        """
        Indices of the faces that surround each face.

        Dimensions: `(n_face, n_max_face_faces)`
        """
        return self.uxgrid.face_face_connectivity.values

    @property
    def face_x(self):
        """
        Cartesian x location of the center of each face in meters.
        """
        if self._face_x is None:
            self._face_x = self.uxgrid.face_x.values * self.radius
        return self._face_x

    @property
    def face_y(self):
        """
        Cartesian y location of the center of each face in meters.
        """
        if self._face_y is None:
            self._face_y = self.uxgrid.face_y.values * self.radius
        return self._face_y

    @property
    def face_z(self):
        """
        Cartesian z location of the center of each face in meters.
        """
        if self._face_z is None:
            self._face_z = self.uxgrid.face_z.values * self.radius
        return self._face_z

    @property
    def n_node(self):
        """
        Total number of nodes
        """
        return self.uxgrid.n_node

    @property
    def n_nodes_per_face(self):
        """
        The number of nodes that make up each face.

        Dimensions: `(n_node, )`
        """
        return self.uxgrid.n_nodes_per_face.values

    @property
    def node_lat(self):
        """
        Latitude of each node in degrees.
        """
        return self.uxgrid.node_lat.values

    @property
    def node_lon(self):
        """
        Longitude of each node in degrees.
        """
        return self.uxgrid.node_lon.values

    @property
    def node_x(self):
        """
        Cartesian x location of each node in meters.
        """
        if self._node_x is None:
            self._node_x = self.uxgrid.node_x.values * self.radius
        return self._node_x

    @property
    def node_y(self):
        """
        Cartesian y location of each node in meters.
        """
        if self._node_y is None:
            self._node_y = self.uxgrid.node_y.values * self.radius
        return self._node_y

    @property
    def node_z(self):
        """
        Cartesian z location of each node in meters.
        """
        if self._node_z is None:
            self._node_z = self.uxgrid.node_z.values * self.radius
        return self._node_z

    @property
    def node_face_connectivity(self):
        """
        Indices of the faces that surround each node.

        Dimensions: `(n_node, n_max_node_faces)`
        """
        return self.uxgrid.node_face_connectivity.values

    @property
    def node_elevation(self):
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
        self.uxds["node_elevation"][:] = value

    @property
    def face_elevation(self):
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
        self.uxds["face_elevation"][:] = value

    @property
    def face_indices(self) -> NDArray:
        """
        The indices of the faces of the surface.
        """
        return self.uxds.n_face.values

    @property
    def node_indices(self) -> NDArray:
        """
        The indices of the nodes of the surface.
        """
        return self.uxds.n_node.values


class SurfaceView:
    """
    This is used to generate a regional view of a subset of the surface mesh without making copies of any of the data.

    Parameters
    ----------
    surface : Surface
        The surface object that contains the mesh data.
    face_indices : NDArray | slice
        The indices of the faces to include in the view.
    node_indices : NDArray | slice | None, optional
        The indices of the nodes to include in the view. If None, all nodes connected to the faces are included.
    location : tuple[float, float] | None, optional
        The location of the center of the view in degrees. If this is set, then the view will contain `face_distance`, `node_distance`, `face_bearing`, and `node_bearing` arrays. Otherwise, these will be None.
    """

    def __init__(
        self,
        surface: Surface,
        face_indices: NDArray | slice,
        node_indices: NDArray | slice | None = None,
        location: tuple[float, float] | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_surface", None)
        object.__setattr__(self, "_area", None)
        object.__setattr__(self, "_n_face", None)
        object.__setattr__(self, "_n_node", None)
        object.__setattr__(self, "_face_distance", None)
        object.__setattr__(self, "_node_distance", None)
        object.__setattr__(self, "_face_bearing", None)
        object.__setattr__(self, "_node_bearing", None)
        object.__setattr__(self, "_location", None)
        object.__setattr__(self, "_face_indices", None)
        object.__setattr__(self, "_node_indices", None)

        self.surface = surface
        self._face_indices = face_indices
        if isinstance(face_indices, slice):
            self._n_face = self.surface.face_elevation[face_indices].size
        else:
            self._n_face = face_indices.size

        if node_indices is None:
            node_indices = np.unique(
                surface.uxds.uxgrid.face_node_connectivity.values[face_indices].ravel()
            )
            node_indices = node_indices[node_indices != INT_FILL_VALUE]

        self._node_indices = node_indices
        if isinstance(node_indices, slice):
            self._n_node = surface.node_elevation[node_indices].size
        else:
            self._n_node = node_indices.size

        if location is not None:
            self._location = validate_and_normalize_location(location)
            self._face_distance, self._node_distance = self.get_distance()
            self._face_bearing, self._node_bearing = self.get_initial_bearing()
        return

    def get_distance(
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
        return self.calculate_haversine_distance(
            lon1, lat1, face_lon2, face_lat2
        ), self.calculate_haversine_distance(lon1, lat1, node_lon2, node_lat2)

    def get_initial_bearing(
        self, location: tuple[float, float] | None = None
    ) -> tuple[NDArray, NDArray]:
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
            surface_functions.calculate_initial_bearing(
                lon1, lat1, face_lon2, face_lat2
            ),
            surface_functions.calculate_initial_bearing(
                lon1, lat1, node_lon2, node_lat2
            ),
        )

    def compute_volume(self, elevation: NDArray) -> NDArray:
        """
        Compute the volume of an array of elevation points

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
            raise ValueError(
                "elevation must be an array with the same size as the number of faces in the grid"
            )
        return np.sum(elevation * self.face_areas)

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
            face_areas=self.surface.face_areas,
            face_elevation=self.surface.face_elevation,
            node_face_connectivity=self.node_face_connectivity,
        )
        self.update_elevation(node_elevation, overwrite=True)
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
            raise ValueError(
                "kdiff must be a scalar or an array with the same size as the number of faces in the grid"
            )
        if np.any(kdiff < 0.0):
            raise ValueError("kdiff must be greater than 0.0")
        kdiffmax = np.max(kdiff)

        if abs(kdiffmax) < _VSMALL:
            return
        face_kappa = np.zeros(self.surface.n_face)
        face_kappa[self.face_indices] = kdiff
        if isinstance(self.face_indices, slice) and self.face_indices == slice(None):
            face_indices = np.arange(self.surface.n_face)
        else:
            face_indices = self.face_indices
        delta_face_elevation = surface_functions.apply_diffusion(
            face_areas=self.surface.face_areas,
            face_kappa=face_kappa,
            face_elevation=self.surface.face_elevation,
            face_face_connectivity=self.face_face_connectivity,
            face_indices=face_indices,
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
            raise ValueError(
                "critical_slope_angle must be between 0 and 90 degrees"
            ) from e

        if isinstance(self.face_indices, slice) and self.face_indices == slice(None):
            face_indices = np.arange(self.surface.n_face)
        else:
            face_indices = self.face_indices
        face_lon = np.deg2rad(self.surface.face_lon)
        face_lat = np.deg2rad(self.surface.face_lat)
        delta_face_elevation = surface_functions.slope_collapse(
            face_areas=self.surface.face_areas,
            face_elevation=self.surface.face_elevation,
            face_face_connectivity=self.face_face_connectivity,
            face_indices=face_indices,
            face_lon=face_lon,
            face_lat=face_lat,
            radius=self.surface.radius,
            critical_slope=critical_slope,
        )
        self.update_elevation(delta_face_elevation)
        self.add_data("ejecta_thickness", delta_face_elevation)
        self.interpolate_node_elevation_from_faces()

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
        mean = np.sum(noise[: self.n_face] * self.face_areas) / self.area
        noise -= mean
        self.update_elevation(noise)
        return

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
        Adds new data

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
            if isfacedata:
                n = self.n_face
            else:
                n = self.n_node
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
            raise ValueError(
                "data must be a scalar or an array with the same size as the number of faces or nodes in the grid"
            )

        if name not in self.surface.uxds.data_vars:
            self.surface._add_new_data(
                name, data=0.0, long_name=long_name, units=units, isfacedata=isfacedata
            )

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
                    raise ValueError(
                        "new_elev must be None, a scalar, or an array with the same size as the number of nodes, faces, or nodes+faces"
                    )
        except Exception as e:
            raise ValueError("new_elev must be None, a scalar, or an array") from e

        if update_face:
            self.add_data(
                name="face_elevation", data=new_face_elev, overwrite=overwrite
            )
        if update_node:
            self.add_data(
                name="node_elevation", data=new_node_elev, overwrite=overwrite
            )

        return

    def get_reference_surface(
        self, reference_radius: float, **kwargs: Any
    ) -> NDArray[np.float64]:
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

            region_surf = self.elevation_to_cartesian(region_coords, region_elevation)

            # Initial guess for the sphere center and radius
            guess_radius = 1.0 + region_elevation.mean()
            initial_guess = [0, 0, 0, guess_radius]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                try:
                    bounds = ([-1.0, -1.0, -1.0, 0.5], [1.0, 1.0, 1.0, 2.0])
                    fit_result, _ = curve_fit(
                        _sphere_function,
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
            A = f_vec[:, 0] ** 2 + f_vec[:, 1] ** 2 + f_vec[:, 2] ** 2
            B = -2 * (
                f_vec[:, 0] * reference_sphere_center[0]
                + f_vec[:, 1] * reference_sphere_center[1]
                + f_vec[:, 2] * reference_sphere_center[2]
            )
            C = (
                np.dot(reference_sphere_center, reference_sphere_center)
                - reference_sphere_radius**2
            )
            sqrt_term = B**2 - 4 * A * C
            valid = ~np.isnan(A) & (sqrt_term >= 0.0)

            # Initialize t with default value
            t = np.full_like(A, 1.0)

            # Calculate square root only for valid terms
            sqrt_valid_term = np.sqrt(np.where(valid, sqrt_term, 0))

            # Apply the formula only where valid
            t = np.where(valid, (-B + sqrt_valid_term) / (2 * A), t)
            if np.any(t[valid] < 0):
                t = np.where(valid & (t < 0), (-B - sqrt_valid_term) / (2 * A), t)

            elevations = self.surface.radius * (t * np.linalg.norm(f_vec, axis=1) - 1)
            return elevations

        # Find cells within the crater radius
        faces_within_region = self.face_distance <= reference_radius
        nodes_within_region = self.node_distance <= reference_radius
        points_within_region = np.concatenate(
            [faces_within_region, nodes_within_region]
        )

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
        reference_elevation[points_within_region] = _find_reference_elevations(
            region_coords, region_elevation
        )

        return reference_elevation

    def calculate_haversine_distance(
        self,
        lon1: FloatLike,
        lat1: FloatLike,
        lon2: FloatLike | ArrayLike,
        lat2: FloatLike | ArrayLike,
    ) -> NDArray[np.float64]:
        """
        Calculate the great circle distance between two points on a sphere.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians. Must be a single point
        lat1 : FloatLike
            Latitude of the first point in radians. Must be a single point
        lon2 : FloatLike, ArrayLike
            Longitude of the second point(s) in radians. Can be a single point or an array of points.
        lat2 : FloatLike, ArrayLike
            Latitude of the second point(s) in radians. Can be a single point or an array of points.

        Returns
        -------
        NDArray[np.float64]
            Great circle distance between the two points in meters.
        """
        return self.surface.calculate_haversine_distance(lon1, lat1, lon2, lat2)

    @staticmethod
    def elevation_to_cartesian(position: NDArray, elevation: NDArray) -> NDArray:
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

    @staticmethod
    def calculate_initial_bearing(
        lon1: FloatLike, lat1: FloatLike, lon2: FloatLike, lat2: FloatLike
    ) -> float:
        """
        Calculate the initial bearing from one point to another on the surface of a sphere.

        Parameters
        ----------
        lon1 : FloatLike
            Longitude of the first point in radians.
        lat1 : FloatLike
            Latitude of the first point in radians.
        lon2 : FloatLike
            Longitude of the second point in radians.
        lat2 : FloatLike
            Latitude of the second point in radians.

        Returns
        -------
        float
            Initial bearing from the first point to the second point in radians.
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

    @property
    def surface(self) -> Surface:
        """
        The surface object that contains the mesh data.
        """
        return self._surface

    @surface.setter
    def surface(self, value: Surface) -> None:
        """
        Set the surface object.

        Parameters
        ----------
        value : Surface
            The surface object to set.
        """
        self._surface = Surface.maker(value)

    @property
    def n_face(self) -> int:
        """
        The number of faces in the view.
        """
        return self._n_face

    @property
    def n_node(self) -> int:
        """
        The number of nodes in the view.
        """
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
        self.surface.face_elevation[self.face_indices] = value

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
        self.surface.node_elevation[self.node_indices] = value

    @property
    def location(self) -> tuple[float, float]:
        """
        The location of the center of the view.
        """
        return self._location

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
    def face_areas(self) -> NDArray:
        """
        The areas of the faces.
        """
        return self.surface.face_areas[self.face_indices]

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
    def face_node_connectivity(self) -> NDArray:
        """
        Indices of the nodes that make up the faces.

        Dimensions: `(n_face, n_max_face_nodes)`

        Nodes are in counter-clockwise order.
        """
        return self.surface.face_node_connectivity[self.face_indices, :]

    @property
    def face_face_connectivity(self) -> NDArray:
        """
        Indices of the faces that surround the faces.

        Dimensions: `(n_face, n_max_face_faces)`
        """
        return self.surface.face_face_connectivity[self.face_indices, :]

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
    def node_face_connectivity(self) -> NDArray:
        """
        Indices of the faces that surround the nodes.

        Dimensions: `(n_node, n_max_node_faces)`
        """
        return self.surface.node_face_connectivity[self.node_indices, :]

    @property
    def area(self) -> NDArray:
        """
        The total area of the faces in the view.
        """
        if self._area is None:
            self._area = self.face_areas.sum()
        return self._area

    @property
    def face_indices(self) -> NDArray:
        """
        The indices of the faces in the view.
        """
        return self._face_indices

    @property
    def node_indices(self) -> NDArray:
        """
        The indices of the nodes in the view.
        """
        return self._node_indices


import_components(__name__, __path__, ignore_private=True)
