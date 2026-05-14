from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
import warnings
from abc import abstractmethod
from contextlib import AbstractContextManager
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pyvista as pv
import rasterio
import uxarray as uxr
import xarray as xr
from affine import Affine
from matplotlib.axes import Axes
from numpy.typing import ArrayLike, NDArray
from pyproj import CRS, Transformer
from rasterio import DatasetReader, MemoryFile, windows
from rasterio.enums import Resampling
from rasterio.merge import merge
from rasterio.transform import rowcol
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window
from scipy.optimize import OptimizeWarning
from tqdm import tqdm
from uxarray import INT_FILL_VALUE, UxDataArray, UxDataset
from vtk import vtkUnstructuredGrid

from cratermaker.bindings import surface_bindings
from cratermaker.components.target import Target
from cratermaker.constants import _SMALLFAC, _VSMALL, FloatLike, PairOfFloats
from cratermaker.core.base import ComponentBase, CratermakerBase, import_components
from cratermaker.utils.general_utils import format_large_units, validate_and_normalize_location
from cratermaker.utils.montecarlo_utils import get_random_location_on_face

_N_TAG_LAYERS = 8


class Surface(ComponentBase):
    _registry: dict[str, type[Surface]] = {}

    """
    Used for handling surface-related data and operations in the cratermaker project.

    It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations.
    The Surface class extends UxDataset for the cratermaker project.
    """

    _registry: dict[str, type[Surface]] = {}

    def __init__(
        self,
        target: Target | str | None = None,
        simdir: str | Path | None = None,
        **kwargs,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        target : Target, optional
            The target body or name of a known target body for the impact simulation.
        reset : bool, optional
            Flag to indicate whether to reset the surface. Default is the value of `regrid`
        regrid : bool, optional
            Flag to indicate whether to regrid the surface. Default is False.
        simdir : str | Path
            |simdir|
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.components.target import Target

        super().__init__(simdir=simdir, **kwargs)

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
        object.__setattr__(self, "_crs", None)
        object.__setattr__(self, "_output_file_prefix", "surface")
        object.__setattr__(self, "_output_dir_name", "surface")
        object.__setattr__(self, "_grid_file_prefix", "grid")
        object.__setattr__(self, "_output_file_extension", "nc")
        object.__setattr__(self, "_is_new", None)

        self._output_file_pattern = [
            f"{self._output_file_prefix}*.{self._output_file_extension}",
            f"{self._grid_file_prefix}.{self._output_file_extension}",
        ]

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
        str_repr = super().__str__()
        str_repr += f"Target: {self.target.name}\n"
        str_repr += f"Grid File: {self.grid_file}\n"
        str_repr += f"Number of faces: {self.n_face}\n"
        str_repr += f"Number of nodes: {self.n_node}\n"
        return str_repr

    @classmethod
    def maker(
        cls: Surface,
        surface: str | Surface | None = None,
        target: Target | str | None = None,
        reset: bool = True,
        regrid: bool = False,
        simdir: str | Path | None = None,
        **kwargs,
    ) -> Surface:
        """
        Initialize a Surface model with the given name or instance.

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
            |simdir|
        **kwargs : Any
            |kwargs|

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

    def __getattr__(self, name: str):
        """
        Get attribute values from the uxds dataset.

        Parameters
        ----------
        name : str
            The name of the attribute to retrieve.
        """
        if not hasattr(self, "uxds"):
            raise AttributeError(f"{type(self).__name__!s} has no attribute {name!r}")
        uxds = object.__getattribute__(self, "uxds")
        if name in uxds:
            return uxds[name].data
        raise AttributeError(f"{type(self).__name__!s} has no attribute {name!r}")

    def saved_output_files(self, **kwargs: Any) -> list[Path]:
        """
        Check if the component has any output files in its output directory.

        Returns
        -------
        list[Path]
            A list of Path objects representing the files that would be removed during a reset operation. Returns an empty list if no files found
        """
        return self._full().saved_output_files(**kwargs)

    def reset(self, **kwargs: Any) -> None:
        """
        Reset the surface to its initial state.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|
        """
        super().reset(**kwargs)

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
        self.is_new = True

        return

    def extract_region(
        self,
        location: tuple[FloatLike, FloatLike],
        region_radius: FloatLike,
        at_least_one_face: bool = False,
        **kwargs: Any,
    ) -> LocalSurface | None:
        """
        Extract a regional grid based on a given location and radius.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.
        region_radius : float
            The radius of the region to extract in meters.
        at_least_one_face : bool, optional
            If True, ensure that at least one face is returned, even if the region radius is very small. Default is False.
        **kwargs : Any
            |kwargs|

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
                if at_least_one_face:
                    nearest_face = self.find_nearest_face(location)
                    face_indices = np.array([nearest_face])
                else:
                    return None

            # First select edges and nodes that are attached to these faces
            edge_indices = np.unique(self.face_edge_connectivity[face_indices].ravel())
            edge_indices = edge_indices[edge_indices != INT_FILL_VALUE]

            node_indices = np.unique(self.face_node_connectivity[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]

            # Now add in all faces that are connected to anything inside the region, so that the outermost border of the local region has a buffer of faces
            # These are needed for diffusion calculations
            neighbor_faces = self.face_face_connectivity[face_indices].ravel()
            node_faces = self.node_face_connectivity[node_indices].ravel()
            edge_faces = self.edge_face_connectivity[edge_indices].ravel()
            face_indices = np.unique(np.concatenate((face_indices, neighbor_faces, node_faces, edge_faces)))
            face_indices = face_indices[face_indices != INT_FILL_VALUE]

            # Add in all nodes that are attached to these buffer faces
            node_indices = np.unique(self.face_node_connectivity[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]

        else:  # This is the entire surface
            face_indices = self.face_indices
            edge_indices = self.edge_indices
            node_indices = self.node_indices

        return LocalSurface(
            surface=self,
            face_indices=face_indices,
            node_indices=node_indices,
            edge_indices=edge_indices,
            location=location,
            region_radius=region_radius,
            **vars(self.common_args),
        )

    def add_data(
        self,
        name: str,
        data: FloatLike | NDArray,
        long_name: str | None = None,
        units: str | None = None,
        isfacedata: bool = True,
        overwrite: bool = False,
        fill_value: float = 0.0,
        dtype=np.float64,
        positive_only: bool = False,
        **kwargs: Any,
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
        fill_value : float, optional
            The fill value to use for new data variables. Default is 0.0.
        dtype : data-type, optional
            The data type of the data variable. Default is np.float64.
        positive_only: bool, optional
            If True, only allow positive values on the data (data will be clipped at 0.0). Default is False.
        **kwargs : Any
            |kwargs|

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
            fill_value=fill_value,
            dtype=dtype,
            positive_only=positive_only,
            **kwargs,
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
            |kwargs|

        Notes
        -----
        When passing combined data, the first part of the array will be used for face elevation and the second part for node elevation.
        """
        return self._full().update_elevation(new_elevation=new_elevation, overwrite=overwrite, **kwargs)

    def add_tag(
        self,
        name: str,
        tag: int | None = None,
        long_name: str | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Adds an integer tag to the surface.

        Used primarily for tracking crater ids.

        Parameters
        ----------
        name : str
            The name of the tag variable to perform the add operation on.
        tag : int | None
            The integer to value of the tag. If None is provided, the tag layers will be reset to zero
        long_name : str | None
            The long name of the tag variable. This is only used if the tag variable doesn't exist in the surface dataset. If None, no long name will be added.
        **kwargs : Any
            |kwargs|
        """
        return self._full().add_tag(name=name, tag=tag, long_name=long_name, **kwargs)

    def remove_tag(self, name: str, tag: int, **kwargs: Any) -> None:
        """
        Removes an integer tag from the surface.

        Parameters
        ----------
        name : str
            The name of the tag variable to perform the remove operation on.
        tag : int
            The integer tag to be removed.
        **kwargs : Any
            |kwargs|
        """
        return self._full().remove_tag(name=name, tag=tag, **kwargs)

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

    def calculate_face_and_node_distances(self, location: tuple[float, float], validate: bool = True) -> tuple[NDArray, NDArray]:
        """
        Computes the distances from a given location to all faces and nodes.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.
        validate : bool, optional
            If the location should be validated and normalized. Only use if the validation is causing a performance bottleneck.

        Returns
        -------
        NDArray
            Array of distances for each face in meters.
        NDArray
            Array of distances for each node in meters.
        """
        return self._full().calculate_face_and_node_distances(location, validate=validate)

    def calculate_face_and_node_bearings(self, location: tuple[float, float]) -> tuple[NDArray, NDArray]:
        """
        Computes the initial bearing (relative to North) from a given location to all faces and nodes.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        NDArray
            Array of initial bearings for each face in degrees.
        NDArray
            Array of initial bearings for each node in degrees.

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
        return self._full().interpolate_node_elevation_from_faces()

    def get_random_location_on_face(self, face_index: int, **kwargs: Any) -> float | tuple[float, float] | ArrayLike:
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
        **kwargs : Any
            |kwargs|

        Returns
        -------
        (lon,lat) or ndarray[(lon,lat)] of given size
            A pair or array of pairs of longitude and latitude values in degrees.

        Notes
        -----
        This method is a wrapper for :py:meth:`cratermaker.utils.montecarlo_utils.get_random_location_on_face`.
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
        interval: int = 0,
        time_variables: dict | None = None,
        filename: str | None = None,
        **kwargs,
    ) -> None:
        """
        Save the surface data to the specified directory.

        Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the time dimension.

        Parameters
        ----------
        interval : int, optional
            Interval number to append to the data file name. Default is 0.
        time_variables : dict, optional
            Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the time dimension. Default is None.
        **kwargs : Any
            |kwargs|
        """
        return self._full().save(
            interval=interval,
            time_variables=time_variables,
            filename=filename,
            **kwargs,
        )

    def export(self, driver: str = "GPKG", interval: int | None = None, ask_overwrite: bool | None = None, **kwargs: Any) -> None:
        """
        Export the surface data to the specified format.

        Parameters
        ----------
        driver : str, optional
            The driver to use export the data to. Supported formats are 'VTK' or a driver supported by GeoPandas ('GPKG', 'ESRI Shapefile', etc.), and 'GeoTIFF'.
        interval : int | None, optional
            |interval_export|
        ask_overwrite : bool, optional
            |ask_overwrite_methods|
        **kwargs : Any
            |kwargs|
        """
        return self._full().export(
            driver=driver,
            interval=interval,
            ask_overwrite=ask_overwrite,
            **kwargs,
        )

    def to_vector_file(
        self,
        driver: str = "GPKG",
        interval: int | None = None,
        **kwargs,
    ) -> None:
        """
        Export the face-associated data from the surface view data to a vector file using GeoPandas.

        See `geopandas.GeoDataFrame.to_file <https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html>`_ for more detailed information on the available parameters.

        Parameters
        ----------
        driver : str, optional
            The file format driver to use for exporting. Default is 'GPKG'.
        interval : int | None, optional
            |interval_export|
        **kwargs : Any
            |kwargs|
        """
        return self._full().to_vector_file(driver=driver, interval=interval, **kwargs)

    def to_raster(
        self, variable_name: str = "face_elevation", **kwargs: Any
    ) -> tuple[NDArray[np.float32], tuple[float, float, float, float], Any, CRS]:
        """
        Rasterize a face-based variable into a 2D raster using rasterio.

        Parameters
        ----------
        variable_name : str, optional
            The name of the variable to rasterize. Default is "face_elevation".
        **kwargs : Any
            |kwargs|

        Returns
        -------
        raster : NDArray[np.float32]
            The rasterized variable as a 2D numpy array.
        extent : tuple[float, float, float, float]
            The extent of the raster in the format (xmin, xmax, ymin, ymax).
        transform : Affine
            The affine transform for the raster.
        crs : CRS
            The coordinate reference system of the raster.
        """
        return self._full().to_raster(variable_name, **kwargs)

    def get_raster_dims(self):
        """
        Get the dimensions of the raster based on the target radius and pixel size.

        Returns
        -------
        tuple[int, int]
            The width and height of the raster in pixels.
        """
        return self._full().get_raster_dims()

    def to_geotiff_file(
        self,
        interval: int | None = None,
        variable_name: str = "face_elevation",
        **kwargs,
    ) -> None:
        """
        Rasterize a face-based elevation variable into a GeoTIFF using rasterio.

        Parameters
        ----------
        interval : int | None, optional
            |interval_export|
        variable_name : str, optional
            The name of the variable to rasterize. Default is "face_elevation".
        **kwargs : Any
            |kwargs|
        """
        return self._full().to_geotiff(
            interval=interval,
            variable_name=variable_name,
            **kwargs,
        )

    def to_vtk_mesh(self, uxds: UxDataset | None = None, interval: int | None = None, **kwargs: Any) -> vtkUnstructuredGrid:
        """
        Exports the surface to a VTK PolyData object.

        Parameters
        ----------
        uxds : UxDataset, optional
            The dataset to export. If None, the method will use currently loaded data in the surface. Default is None.
        interval : int, optional
            The interval number to export. If provided, the method will either extract the interval number from uxds (if it has intervals saved), or, if no uxds is passed, load a saved interval from file.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        vtkUnstructuredGrid
            The VTK PolyData object representing the regional mesh.
        """
        return self._full().to_vtk_mesh(uxds=uxds, interval=interval, **kwargs)

    def to_vtk_file(
        self,
        interval: int | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Export the surface mesh to a VTK file and store it in the default export directory.

        Parameters
        ----------
        interval : int | None, optional
            |interval_export|
        **kwargs : Any
            |kwargs|
        """
        return self._full().to_vtk_file(
            interval=interval,
            **kwargs,
        )

    def plot(
        self,
        plot_style: Literal["map", "hillshade"] = "map",
        variable_name: str | None = None,
        interval: int | None = None,
        cmap: str | None = None,
        label: str | None = None,
        scalebar: bool | None = None,
        colorbar: bool = True,
        show: bool = True,
        save: bool = False,
        ax: Axes | None = None,
        close_when_done: bool | None = None,
        minimum_plot_width: float | None = 800.0,
        **kwargs: Any,
    ) -> Axes:
        """
        Plot an image of the surface.

        Parameters
        ----------
        plot_style : str, optional
            The style of the plot. Options are "map" and "hillshade". In "map" mode, the variable is displayed as a colored map. In "hillshade" mode, a hillshade image is generated using "face_elevation" data. If a different variable is passed to `variable`, then the hillshade will be overlayed with that variable's data. Default is "map".
        variable_name : str | None, optional
            The variable to plot. If None is provided then "face_elevation" is used in "map" mode.
        interval: int | None, optional
            The interval number of the surface to plot. If None, the currently loaded surface data will be used.
        cmap : str, optional
            The colormap to use for the plot. If None, a default colormap will be used ("cividis" by default and "grey" when plot_style=="hillshade" and variable=="face_elevation").
        label : str | None, optional
            A label for the plot. If None, no label will be added.
        scalebar : bool, optional
            If True, a scalebar will be added to the plot. Default is False for global surfaces.
        colorbar : bool, optional
            If True, a colorbar will be added to the plot when using "map" plot_style or "hillshade" with a variable overlay. Default is True.
        show : bool, optional
            If True, the plot will be displayed. Default for local surfaces is True.
        save : bool, optional
            If True, the plot will be saved to the default plot directory. Default is False.
        ax : matplotlib.axes.Axes, optional
            An existing Axes object to plot on. If None, a new figure and axes will be created.
        close_when_done : bool, optional
            If True, the figure will be closed after plotting. Default is True when save is True and show is False, and False otherwise.
        minimum_plot_width : float, optional
            Because the width of the plot is determined by the number of faces, small regions will generate small plots with labels that are hard to read. This parameter sets a lower limit to the width of the image that is generated by the plot. By default it is 800. Set to None to turn it off.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        matplotlib.axes.Axes
            A Matplotlib Axes object containing the plot.
        """
        if scalebar is None:
            scalebar = False
        return self._full().plot(
            variable_name=variable_name,
            plot_style=plot_style,
            cmap=cmap,
            interval=interval,
            label=label,
            scalebar=scalebar,
            colorbar=colorbar,
            show=show,
            save=save,
            ax=ax,
            minimum_plot_width=minimum_plot_width,
            **kwargs,
        )

    def pyvista_plotter(
        self,
        variable_name: str | None = None,
        variable: ArrayLike | None = None,
        interval: int | None = None,
        theme: str | None = None,
        transparent_background: bool | None = None,
        plotter: pv.Plotter | None = None,
        enable_interactive: bool = True,
        **kwargs: Any,
    ) -> pv.Plotter:
        """
        Show the surface region using an interactive 3D plot with PyVista.

        Parameters
        ----------
        variable_name : str, optional
            The name of the variable to plot. If the name of the variable is not already stored on the surface mesh, then the `variable` argument must also be passed. Default is None, which will plot a greyscale image of the surface.
        variable: (n_face) array, optional
            An array face values that will be used to color the surface mesh. This is required if `variable_name` is not stored on the mesh.
        interval : int | None, optional
            The interval number of the surface to show. If None, the currently loaded surface data will be displayed. Default is None.
        theme : str, optional
            The PyVista theme to use for the plot. If None, the default PyVista theme will be used.
        transparent_background : bool, optional
            If True, the background of the plot will be transparent. Default is False.
        plotter : pv.Plotter, optional
            An existing PyVista Plotter object to use for the plot. If None, a new Plotter object will be created. Default is None.
        enable_interactive : bool, optional
            If True, the key events for the plotter will be updated to include custom events for navigating between intervals. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        pyvista.Plotter
            The PyVista Plotter object for further customization.
        """
        return self._full().pyvista_plotter(
            variable_name=variable_name,
            variable=variable,
            interval=interval,
            theme=theme,
            transparent_background=transparent_background,
            plotter=plotter,
            enable_interactive=enable_interactive,
            **kwargs,
        )

    def compute_distances(
        self,
        locations: ArrayLike,
        reference_location: PairOfFloats,
    ) -> NDArray[np.float64]:
        """
        Calculate the great circle distance between one point and one or more other points in meters.

        Parameters
        ----------
        locations : FloatLike or ArrayLike
            Array of (lon, lat) locations of the second point or array of points in degrees.
        reference_location : PairOfFloats
            (lon, lat) location of the reference point in degrees.

        Returns
        -------
        NDArray
            Great circle distance between the two points in meters.

        Notes
        -----
        This is a wrapper for a compiled Rust function and is intended to be used as a helper to calculate_face_and_node_distances.
        """
        return self._full().compute_distances(locations=locations, reference_location=reference_location)

    def compute_bearings(
        self,
        locations: ArrayLike,
        reference_location: PairOfFloats,
    ) -> NDArray[np.float64]:
        """
        Calculate the initial bearing relative to North from one point to one or more other points in degrees.

        Parameters
        ----------
        locations : ArrayLike
            Longitude and latitude of the second point or array of points in degrees.
        reference_location : PairOfFloats
            Longitude and latitude of the first point in degrees.

        Returns
        -------
        NDArray
            Initial bearing from the first point to the second point or points in degrees.
        """
        return self._full().compute_bearings(locations=locations, reference_location=reference_location)

    def compute_location_from_distance_bearing(
        self,
        distance: FloatLike | ArrayLike,
        bearing: FloatLike | ArrayLike,
        reference_location: PairOfFloats,
    ) -> NDArray[np.float64]:
        """
        Calculate the longitude and latitude of one or more points given a reference point, initial bearings, and distances.

        Parameters
        ----------
        bearings : FloatLike or ArrayLike
            Initial bearing from the reference point to the target point or points in degrees.
        distances : FloatLike or ArrayLike
            Great circle distance from the reference point to the target point or points in meters.
        reference_location : PairOfFloats
            Longitude and latitude of the reference point in degrees.

        Returns
        -------
        NDArray
            Longitude and latitude of the target point or points in degrees.
        """
        return self._full().compute_location_from_distance_bearing(
            distance=distance,
            bearing=bearing,
            reference_location=reference_location,
        )

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

    def read_saved_output(self, interval: int | None = None, reset: bool = False, **kwargs: Any) -> UxDataset:
        """
        Load the grid and data files into a UxDataset object.

        Parameters
        ----------
        interval : int, optional
            Interval number to read from the data files. Default is None (all saved intervals)
        reset : bool, optional
            Flag to indicate whether to reset the surface. If True it reads in the grid but creates an empty dataset. Default is False.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        UxDataset
            An initialized UxDataset object containing the grid and data.
        """
        return self._full().read_saved_output(interval=interval, reset=reset, **kwargs)

    def _load_from_files(self, interval: int = -1, reset: bool = False, regrid: bool = False, **kwargs: Any) -> None:
        """
        Load the grid and data files into the surface object.

        This function loads the grid file and data files from the specified directory. If the grid file does not exist, it will attempt to create a new grid.  If the data files do not exist, it will create an empty dataset. If reset is True, it will delete all data files except the grid file.

        Parameters
        ----------
        interval : int, optional
            Interval number to read from the data files. Default is -1 (the last saved interval).
        reset : bool, optional
            Flag to indicate whether to reset the surface. Default is False.
        regrid : bool, optional
            Flag to indicate whether to regrid the surface. Default is False.
        """
        # Get the names of all data files in the data directory that are not the grid file
        regrid = self._regrid_if_needed(regrid=regrid, **kwargs)
        reset = reset or regrid
        self.is_new = reset
        ask_overwrite = self.ask_overwrite  # Store this in case of regridding. If regridding, we need it to be False for the reset
        if regrid:
            self.ask_overwrite = False

        # Read in only the last saved data file
        uxds = self.read_saved_output(interval=interval, reset=reset, **kwargs)
        if "interval" in uxds.dims:
            uxds = uxds.isel(interval=-1)
        object.__setattr__(self, "_uxds", uxds)

        if reset:
            self.reset(**kwargs)

        self.ask_overwrite = ask_overwrite  # Restore in case we had to set it to False for the reset during regridding

        return

    def _write_grid_file(self, uxgrid: uxr.Grid | None = None, grid_file: Path | str | None = None, **kwargs: Any) -> None:
        """
        Write the grid to a NetCDF file.

        Parameters
        ----------
        uxgrid : uxr.Grid, optional
            The grid to be written to the file. If None, the grid will be obtained from the surface object. Default is None.
        grid_file : Path, str, optional
            The path to the grid file. If None, the path will be obtained from the surface object. Default is None.
        **kwargs : Any
            |kwargs|
        """
        import uxarray.conventions.ugrid as ugrid

        if uxgrid is None:
            uxgrid = self.uxgrid
        if grid_file is None:
            grid_file = self.grid_file
        else:
            grid_file = Path(grid_file)

        grid_file.unlink(missing_ok=True)

        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            ds = uxgrid.to_xarray()

            # There is some issue that is causing a mismatch between the data variables and the grid_topology attributes that is causing problems when reading the file back in under some circumstances, such as deleting the grid file and reinitializing the surface. To get around this, we will check that all of the attributes are correct before writing the file.
            ugrid_names = ugrid.CONNECTIVITY_NAMES + ugrid.DIM_NAMES
            for name in ugrid_names:
                if name in ds.grid_topology.attrs and name not in ds:
                    del ds.grid_topology.attrs[name]
            if "edge_dimension" in ds.grid_topology.attrs and "n_edge" not in ds:
                del ds.grid_topology.attrs["edge_dimension"]
            if "face_coordinates" in ds.grid_topology.attrs:
                face_coords = ds.grid_topology.attrs["face_coordinates"].split()
                for coord in face_coords:
                    if coord not in ds:
                        del ds.grid_topology.attrs["face_coordinates"]
                        break

            ds.to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())

        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name, grid_file)

        return

    def _generate_grid(self, **kwargs: Any) -> None:
        """
        Generate a tessellated mesh of a sphere of based on the particular Surface component that is being used.
        """
        # Ignore divide by zero warnings that occur due to matmul on some systems (e.g. MacOS with M4 chips). These warnings do not affect the functionality of the code, so we can safely ignore them.
        orig_settings = np.seterr(divide="ignore", over="ignore", invalid="ignore")
        points = self._generate_face_distribution(**kwargs)

        threshold = min(10 ** np.floor(np.log10(self.pix / self.radius)), 1e-7)
        uxgrid = uxr.Grid.from_points(points, method="spherical_voronoi", threshold=threshold)
        uxgrid.attrs["_id"] = self._id
        self._write_grid_file(uxgrid=uxgrid)

        regrid = not self._is_same_grid
        if regrid:
            raise ValueError("Grid file does not match the expected parameters.")
        self._compute_face_size(uxgrid)
        np.seterr(**orig_settings)

        return

    @property
    def _is_same_grid(self):
        """Check if the existing grid matches the one defined by the current parameters and returns True if they match after regridding."""
        try:
            with xr.open_dataset(self.grid_file) as ds:
                ds.load()
                uxgrid = uxr.Grid.from_dataset(ds)
                old_id = uxgrid.attrs.get("_id")
                return old_id == self._id
        except Exception:
            # Failed to open an old file for whatever reason, so we'll need to regrid
            print("Generating a new grid.")
            return False

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
            A boolean indicating whether grid was regenerated.
        """
        # Find out if the file exists, if it does't we'll need to make a new grid
        self.output_dir.mkdir(parents=True, exist_ok=True)
        regrid = regrid or not Path(self.grid_file).exists()

        if not regrid:
            regrid = not self._is_same_grid

        if regrid:
            print("Creating a new grid")
            self._generate_grid(**kwargs)

        return regrid

    def _full(self) -> LocalSurface:
        return LocalSurface(
            self, face_indices=slice(None), node_indices=slice(None), edge_indices=slice(None), **vars(self.common_args)
        )

    def _save_data(
        self, ds: xr.Dataset | xr.DataArray, interval: int = 0, filename: Path | None = None, reset: bool = False, **kwargs: Any
    ) -> None:
        """
        Save the data to the specified directory.

        Parameters
        ----------
        ds : xr.Dataset or xr.DataArray
            The data to be saved.
        interval : int, Default is 0.
            Interval number to append to the data file name. Default is 0.
        filename : Path | None = None
            The name of the data file. If None, the name will be generated from the variable name and interval number.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        This function first saves to a temporary file and then moves that file to the final destination. This is done to avoid file
        locking issues with NetCDF files.
        """
        if isinstance(ds, xr.DataArray):
            ds = ds.to_dataset()

        if "interval" not in ds.dims:
            ds = ds.expand_dims(["interval"])
        if "interval" not in ds.coords:
            ds = ds.assign_coords({"interval": [interval]})

        ds.load()
        if filename is None:
            filename = self.output_filename(interval)
        else:
            filename = Path(filename)
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True, **kwargs) as temp_dir:
            data_file = self.output_dir / filename
            if data_file.exists():
                try:
                    with xr.open_mfdataset(data_file) as ds_old:
                        ds_old = ds_old.load()
                    ds_file = ds.merge(ds_old, compat="override")
                except Exception:
                    ds_file = ds
            else:
                ds_file = ds

            temp_file = Path(temp_dir) / filename

            comp = {"zlib": True, "complevel": 9}
            encoding = dict.fromkeys(ds_file.data_vars, comp)
            ds_file.to_netcdf(temp_file, encoding=encoding, **kwargs)
            ds_file.close()
            ds.close()
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
        interval: int = 0,
        dtype=np.float64,
        **kwargs: Any,
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
        interval : int, optional, default 0
            The interval number to use when saving the data to the data file.
        dtype : data-type, optional
            The data type of the data variable. Default is np.float64.
        **kwargs : Any
            |kwargs|
        """
        if self.uxgrid is None:
            with xr.open_dataset(self.grid_file, **kwargs) as ds:
                ds.load()
                uxgrid = uxr.Grid.from_dataset(ds)
        else:
            uxgrid = self.uxgrid
        if long_name is None and units is None:
            attrs = None
        else:
            attrs = {}
            if long_name:
                attrs["long_name"] = long_name
            if units:
                attrs["units"] = units

        if isfacedata:
            dims = "n_face"
            size = uxgrid.n_face
        else:
            dims = "n_node"
            size = uxgrid.n_node

        if data is None:
            data = np.zeros(size, dtype=dtype)
        elif np.isscalar(data):
            data = np.full(size, data, dtype=dtype)
        else:
            if data.size != size:
                raise ValueError("data must have the same size as the number of faces or nodes in the grid")
        uxda = UxDataArray(data=data, dims=dims, attrs=attrs, name=name, uxgrid=uxgrid, **kwargs)

        self._uxds[name] = uxda

        if save_to_file:
            self._save_data(uxda, interval=interval, **kwargs)
        return

    @abstractmethod
    def _generate_face_distribution(self, **kwargs: Any) -> NDArray: ...

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

    def get_location_extents(self, location: PairOfFloats, radius: FloatLike) -> tuple[float, float, float, float]:
        """
        Computes the longitude and latitude extents of a given region.

        Parameters
        ----------
        location : PairOfFloats
            The center of the region.
        radius : FloatLike
            The radius of the region.

        Returns
        -------
        lon_min : float
            Minimum longitude in degrees.
        lon_max : float
            Maximum longitude in degrees.
        lat_min : float
            Minimum latitude in degrees.
        lat_max : float
            Maximum latitude in degrees.

        """
        R = self.target.radius
        lon0, lat0 = location
        alpha_deg = np.degrees(radius / R)

        lat_min = max(-90.0, lat0 - alpha_deg)
        lat_max = min(90.0, lat0 + alpha_deg)

        # If we hit a pole, longitudes span everything
        if abs(lat0) + alpha_deg >= 90.0:
            lon_min, lon_max = -180.0, 180.0
        else:
            half_lon_span = alpha_deg / np.cos(np.radians(lat0))
            half_lon_span = min(half_lon_span, 180.0)

            lon_min = lon0 - half_lon_span
            lon_max = lon0 + half_lon_span
        return float(lon_min), float(lon_max), float(lat_min), float(lat_max)

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self.component_name, self.target.name, self.radius]

    @property
    def _id(self):
        """
        The hash id of the grid. This is used for determining if the grid needs to be regridded.
        """
        combined = ":".join(str(v) for v in self._hashvars)
        combined += ":v2026.2.0"  # Add version number because of API changes
        hash_object = hashlib.sha256(combined.encode())
        return hash_object.hexdigest()

    @property
    def uxds(self) -> UxDataset:
        """The data associated with the surface as an instance of UxDataset."""
        return self._uxds

    @property
    def uxgrid(self):
        """The grid object as an instance of a UxArray Grid."""
        if self.uxds is not None:
            return self.uxds.uxgrid

    @property
    def gridtype(self):
        """The name of the grid type."""
        return self._component_name

    @property
    def grid_file(self):
        """Path to the grid file."""
        return self.output_dir / f"{self._grid_file_prefix}.{self._output_file_extension}"

    @property
    def target(self):
        """The Target object associated with this Surface object."""
        return self._target

    @target.setter
    def target(self, value):
        from cratermaker.components.target import Target

        self._target = Target.maker(value)
        return

    @property
    def pix(self) -> float:
        """The effective pixel size of the mesh in meters."""
        if self._pix is None:
            self._pix = self.pix_mean
        return self._pix

    @property
    def pix_mean(self) -> float:
        """The mean pixel size of the mesh in meters."""
        if self._pix_mean is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_mean

    @property
    def pix_std(self) -> float:
        """The standard deviation of the pixel size of the mesh in meters."""
        if self._pix_std is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_std

    @property
    def pix_min(self) -> float:
        """The minimum pixel size of the mesh in meters."""
        if self._pix_min is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_min

    @property
    def pix_max(self) -> float:
        """The maximum pixel size of the mesh in meters."""
        if self._pix_max is None and self.uxgrid is not None:
            self._compute_face_size()
        return self._pix_max

    @property
    def radius(self) -> float:
        """Radius of the target body in meters."""
        return self.target.radius

    @property
    def area(self) -> float:
        """Total surface area of the target body in m²."""
        if self._area is None:
            self._area = float(self.face_area.sum())
        return self._area

    @property
    def face_area(self) -> NDArray[np.float64]:
        """An array of areas of each individual face in m²."""
        if self._face_area is None:
            self._compute_face_size()
        return self._face_area

    @property
    def face_size(self) -> NDArray[np.float64]:
        """The effective size of each face in meters, which s simply the square root of the face area, but is useful for certain comparisons and is equivalent to the `pix` variable from CTEM."""
        if self._face_size is None:
            self._compute_face_size()
        return self._face_size

    @property
    def smallest_length(self) -> float:
        """The smallest length of the mesh in meters."""
        if self._smallest_length is None:
            self._smallest_length = float(np.min(self.face_size) * _SMALLFAC)
        return self._smallest_length

    @property
    def face_lat(self) -> NDArray[np.float64]:
        """Latitude of the center of each face in degrees."""
        return self.uxgrid.face_lat.values

    @property
    def face_lon(self) -> NDArray[np.float64]:
        """Longitude of the center of each face in degrees."""
        return self.uxgrid.face_lon.values

    @property
    def face_x(self) -> NDArray[np.float64]:
        """Cartesian x location of the center of each face in meters."""
        if self._face_x is None:
            self._face_x = self.uxgrid.face_x.values * self.radius
        return self._face_x

    @property
    def face_y(self) -> NDArray[np.float64]:
        """Cartesian y location of the center of each face in meters."""
        if self._face_y is None:
            self._face_y = self.uxgrid.face_y.values * self.radius
        return self._face_y

    @property
    def face_z(self) -> NDArray[np.float64]:
        """Cartesian z location of the center of each face in meters."""
        if self._face_z is None:
            self._face_z = self.uxgrid.face_z.values * self.radius
        return self._face_z

    @property
    def node_lat(self) -> NDArray[np.float64]:
        """Latitude of each node in degrees."""
        return self.uxgrid.node_lat.values

    @property
    def node_lon(self) -> NDArray[np.float64]:
        """Longitude of each node in degrees."""
        return self.uxgrid.node_lon.values

    @property
    def node_x(self) -> NDArray[np.float64]:
        """Cartesian x location of each node in meters."""
        if self._node_x is None:
            self._node_x = self.uxgrid.node_x.values * self.radius
        return self._node_x

    @property
    def node_y(self) -> NDArray[np.float64]:
        """Cartesian y location of each node in meters."""
        if self._node_y is None:
            self._node_y = self.uxgrid.node_y.values * self.radius
        return self._node_y

    @property
    def node_z(self) -> NDArray[np.float64]:
        """Cartesian z location of each node in meters."""
        if self._node_z is None:
            self._node_z = self.uxgrid.node_z.values * self.radius
        return self._node_z

    @property
    def face_indices(self) -> NDArray[np.int64]:
        """The indices of the faces of the surface."""
        return self.uxds.n_face.values

    @property
    def node_indices(self) -> NDArray[np.int64]:
        """The indices of the nodes of the surface."""
        return self.uxds.n_node.values

    @property
    def edge_indices(self) -> NDArray[np.int64]:
        """The indices of the edges of the surface."""
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
        """The total area of all faces in each bin in m²."""
        if self._face_bin_area is None:
            self._compute_face_bins()

        return self._face_bin_area

    @property
    def face_bin_argmin(self) -> list[int]:
        """The index of the smallest face in each bin."""
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return self._face_bin_argmin

    @property
    def face_bin_argmax(self) -> list[int]:
        """The index of the largest face in each bin."""
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return self._face_bin_argmax

    @property
    def face_bin_min_areas(self) -> list[float]:
        """The area of the smallest face in each bin in m²."""
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return [float(self.face_area[face_index]) for face_index in self.face_bin_argmin]

    @property
    def face_bin_max_areas(self) -> list[float]:
        """The area of the largest face in each bin in m²."""
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return [float(self.face_area[face_index]) for face_index in self.face_bin_argmax]

    @property
    def face_bin_min_sizes(self) -> list[float]:
        """The effective size of the smallest face in each bin in meters."""
        if self._face_bin_argmin is None:
            self._compute_face_bins()

        return [float(self.face_size[face_index]) for face_index in self.face_bin_argmin]

    @property
    def face_bin_max_sizes(self) -> list[float]:
        """The effective size of the largest face in each bin in meters."""
        if self._face_bin_argmax is None:
            self._compute_face_bins()

        return [float(self.face_size[face_index]) for face_index in self.face_bin_argmax]

    @property
    def n_face(self) -> int:
        """Total number of faces."""
        return int(self.uxgrid.n_face)

    @property
    def n_node(self) -> int:
        """Total number of nodes."""
        return int(self.uxgrid.n_node)

    @property
    def n_edge(self) -> int:
        """Total number of edges."""
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
        """The maximum number of faces that surround a face."""
        return int(self.uxgrid.n_max_face_faces)

    @property
    def edge_face_distance(self) -> NDArray[np.float64]:
        """Distances between the centers of the faces that saddle each edge in meters."""
        if self._edge_face_distance is None:
            self._edge_face_distance = surface_bindings.compute_edge_distances(
                edge_connectivity=self.edge_face_connectivity,
                lon=np.radians(self.face_lon),
                lat=np.radians(self.face_lat),
                radius=self.radius,
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
            self._edge_lengths = surface_bindings.compute_edge_distances(
                edge_connectivity=self.edge_node_connectivity,
                lon=np.radians(self.node_lon),
                lat=np.radians(self.node_lat),
                radius=self.radius,
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
        """The BallTree for the face centers, which can be used for efficient spatial queries."""
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
        """The BallTree for the nodes, which can be used for efficient spatial queries."""
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
        """The BallTree for the edge centers, which can be used for efficient spatial queries."""
        if self._edge_tree is None:
            self._edge_tree = self.uxgrid.get_ball_tree(
                "edge centers",
                distance_metric="haversine",
                coordinate_system="spherical",
                reconstruct=True,
            )

        return self._edge_tree

    @staticmethod
    def get_crs(radius: float, name: str, location: tuple[float, float] | None = None) -> CRS:
        """
        Returns either a CRS for a global sphere or a local azimuthal equidistant projection centered at the given location if one is provided.

        Parameters
        ----------
        radius : float
            The radius of the sphere in meters.
        name : str
            The name of the target body.
        location : tuple[float, float] | None, optional
            The location of the center of the LAEA projection in degrees. If None, a global CRS is returned. Default is None.
        """
        if location is None:
            wkt = (
                f'GEOGCS["GCS_{name}_global",'
                f'    DATUM["D_{name}",'
                f'        ELLIPSOID["{name}",{radius:.6f},0,'
                f'            LENGTHUNIT["metre",1]]],'
                f'    PRIMEM["Reference_Meridian",0,'
                f'        ANGLEUNIT["degree",0.0174532925199433]]],'
                f"    CS[ellipsoidal,2],"
                f'        AXIS["geodetic latitude (Lat)",north,'
                f"            ORDER[1],"
                f'            ANGLEUNIT["degree",0.0174532925199433]],'
                f'        AXIS["geodetic longitude (Lon)",east,'
                f"            ORDER[2],"
                f'            ANGLEUNIT["degree",0.0174532925199433]],'
                f'REMARK["Created by Cratermaker"]]'
            )
        else:
            lon0, lat0 = location
            wkt = (
                f'PROJCRS["Proj_{name}_region_location_{lon0}_{lat0}_Lambert_Azimuthal_Equal_Area",'
                f'    BASEGEOGCRS["GCS_{name}_global",'
                f'        DATUM["D_{name}",'
                f'            ELLIPSOID["{name}",{radius:.6f},0,'
                f'                LENGTHUNIT["metre",1]]],'
                f'        PRIMEM["Reference_Meridian",0,'
                f'            ANGLEUNIT["degree",0.0174532925199433]]],'
                f'        CONVERSION["Lambert Azimuthal Equal Area",'
                f'            METHOD["Lambert Azimuthal Equal Area (Spherical)",'
                f'                ID["EPSG",1027]],'
                f'            PARAMETER["Latitude of natural origin",{lat0},'
                f'                ANGLEUNIT["degree",0.0174532925199433],'
                f'                ID["EPSG",8801]],'
                f'            PARAMETER["Longitude of natural origin",{lon0},'
                f'                ANGLEUNIT["degree",0.0174532925199433],'
                f'                ID["EPSG",8802]],'
                f'            PARAMETER["False easting",0,'
                f'                LENGTHUNIT["metre",1],'
                f'                ID["EPSG",8806]],'
                f'            PARAMETER["False northing",0,'
                f'                LENGTHUNIT["metre",1],'
                f'                ID["EPSG",8807]]],'
                f"        CS[Cartesian,2],"
                f'            AXIS["(E)",east,'
                f"                ORDER[1],"
                f'                LENGTHUNIT["metre",1]],'
                f'            AXIS["(N)",north,'
                f"                ORDER[2],"
                f'                LENGTHUNIT["metre",1]]],'
                f'REMARK["Created by Cratermaker"]'
            )

        return CRS.from_wkt(wkt)

    @property
    def crs(self) -> CRS:
        """Return a geographic CRS (lon/lat in degrees) on a sphere using the target radius from the surface. Axis order is Lon/East, Lat/North."""
        if self._crs is None:
            self._crs = self.get_crs(radius=self.radius, name=self.target.name)
        return self._crs

    @property
    def is_new(self) -> bool:
        """Whether this surface is newly created or has been loaded from disk and not reset."""
        return self._is_new

    @is_new.setter
    def is_new(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("is_new must be a boolean value.")
        self._is_new = value
        return

    def data_composer(self) -> DataComposer:
        """
        Creates a ``DataComposer`` tied to this ``Surface``. This should usually be called in a ``with`` statement.

        Returns
        -------
        DataComposer
        """
        return self._full().data_composer()

    @property
    def face_variables(self) -> list[str]:
        """
        Returns a list of the data variables currently being stored on the faces.
        """
        return self._full().face_variables

    @property
    def node_variables(self) -> list[str]:
        """
        Returns a list of the data variables currently being stored on the nodes.
        """
        return self._full().node_variables


class LocalSurface(CratermakerBase):
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
    reset : bool, optional
        Flag to indicate whether to remove any existing data files in the output directory. Default is True.
    **kwargs : Any
        |kwargs|
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
        super().__init__(**{**kwargs, "simdir": self.surface.simdir})

        object.__setattr__(self, "_uxds", None)
        object.__setattr__(self, "_uxgrid", None)
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
        object.__setattr__(self, "_edge_node_connectivity", None)
        object.__setattr__(self, "_face_edge_connectivity", None)
        object.__setattr__(self, "_face_node_connectivity", None)
        object.__setattr__(self, "_face_face_connectivity", None)
        object.__setattr__(self, "_node_face_connectivity", None)
        object.__setattr__(self, "_crs", None)
        object.__setattr__(self, "_output_file_extension", "nc")
        object.__setattr__(self, "_from_surface", None)
        object.__setattr__(self, "_to_surface", None)
        object.__setattr__(self, "_face_proj_x", None)
        object.__setattr__(self, "_face_proj_y", None)
        object.__setattr__(self, "_node_proj_x", None)
        object.__setattr__(self, "_node_proj_y", None)
        object.__setattr__(self, "_desloped_face_elevation", None)

        self._output_dir_name = self.surface._output_dir_name

        if location is not None:  # This is a true LocalSurface object
            self._location = validate_and_normalize_location(location)
            self._output_file_prefix = "local_surface"
            object.__setattr__(self, "_output_file_prefix", "local_surface")
            object.__setattr__(self, "_grid_file_prefix", "local_grid")
            self._output_file_pattern = [
                f"{self._output_file_prefix}*.{self._output_file_extension}",
                f"{self._grid_file_prefix}.{self._output_file_extension}",
            ]
        else:  # This is really a Surface object wearing a LocalSurface costume.
            object.__setattr__(self, "_output_file_prefix", self.surface._output_file_prefix)
            object.__setattr__(self, "_grid_file_prefix", self.surface._grid_file_prefix)
            object.__setattr__(self, "_region_radius", self.surface.radius)
            self._output_file_pattern = self.surface._output_file_pattern

        return

    def __str__(self) -> str:
        """
        String representation of the LocalSurface object.
        """
        str_repr = "<LocalSurface>\n"
        if self.is_local:
            str_repr += f"Location: {self.location[0]:.2f}°, {self.location[1]:.2f}°\n"

        if self.region_radius:
            str_repr += f"Region Radius: {format_large_units(self.region_radius, quantity='length')}\n"

        str_repr += f"Number of faces: {self.n_face}\n"
        str_repr += f"Number of nodes: {self.n_node}\n"
        return str_repr

    def __getattr__(self, name: str):
        """
        Get attribute from the surface's uxds.

        If the attribute is face-, node-, or edge-based, slice it to only include the faces, nodes, or edges in the local surface.

        Parameters
        ----------
        name : str
            Name of the attribute to get.
        """
        uxds = self.surface.uxds
        if name in uxds:
            arr = uxds[name].values
            if arr.ndim >= 1:
                if arr.shape[0] == self.surface.n_face:
                    return arr[self.face_indices, ...]
                elif arr.shape[0] == self.surface.n_node:
                    return arr[self.node_indices, ...]
                elif arr.shape[0] == self.surface.n_edge:
                    return arr[self.edge_indices, ...]
                else:
                    raise AttributeError(
                        f"Cannot match attribute {name!r} with shape {arr.shape} to face, node, or edge data of the surface"
                    )
            return arr
        raise AttributeError(f"{type(self).__name__!s} has no attribute {name!r}")

    def add_data(
        self,
        name: str,
        data: FloatLike | NDArray,
        long_name: str | None = None,
        units: str | None = None,
        isfacedata: bool = True,
        overwrite: bool = False,
        fill_value: float = 0.0,
        dtype=np.float64,
        positive_only: bool = False,
        **kwargs,
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
        fill_value : float, optional
            The fill value to use for new data variables. Default is 0.0.
        dtype : data-type, optional
            The data type of the data variable. Default is np.float64.
        positive_only: bool, optional
            If True, only allow positive values on the data (data will be clipped at 0.0). Default is False.
        **kwargs
            |kwargs|

        """
        # Check if the data is a scalar or an array
        if np.isscalar(data):
            n = self.n_face if isfacedata else self.n_node
            data = np.full(n, data, dtype=dtype)
        elif isinstance(data, list):
            data = np.array(data, dtype=dtype)
        else:
            data = np.asarray(data, dtype=dtype)
        if data.size == self.n_face:
            isfacedata = True
            indices = self.face_indices
        elif data.size == self.n_node:
            isfacedata = False
            indices = self.node_indices
        else:
            raise ValueError("data must be a scalar or an array with the same size as the number of faces or nodes in the grid")

        if name not in self.surface.uxds.data_vars:
            overwrite = True
            self.surface._add_new_data(
                name, data=fill_value, long_name=long_name, units=units, isfacedata=isfacedata, dtype=dtype, **kwargs
            )

        if overwrite:
            self.surface.uxds[name].data[indices] = data
        else:
            self.surface.uxds[name].data[indices] += data

        if positive_only:
            self.surface.uxds[name].data[indices] = np.maximum(0.0, self.surface.uxds[name].data[indices])

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
            |kwargs|

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
            raise ValueError(
                "new_elev must be None, a scalar, or an array with the same number of elements as either the faces or nodes of the surface mesh."
            ) from e

        if update_face:
            self.add_data(name="face_elevation", data=new_face_elev, overwrite=overwrite)
            if not update_node:
                self.interpolate_node_elevation_from_faces()
        if update_node:
            self.add_data(name="node_elevation", data=new_node_elev, overwrite=overwrite)

        return

    def add_tag(
        self,
        name: str,
        tag: int | None = None,
        long_name: str | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Adds an integer tag to the surface.

        Used primarily for tracking crater ids.

        Parameters
        ----------
        name : str
            The name of the tag to be added.
        tag : int | None
            The integer to apply the tag. If None is provided, the tag layers will be reset to zero
        long_name : str | None
            The long name of the tag to be added. If None, no long name will be added.
        **kwargs : Any
            |kwargs|
        """
        # Reset the tag layers if the tag is None or does not yet exist on the surface
        if name not in self.uxds or tag is None:
            dims = ("n_face", "layer")
            data = np.zeros((self.surface.n_face, _N_TAG_LAYERS), dtype=np.uint32)
            if long_name is not None:
                attrs = {"long_name": long_name}
            else:
                attrs = None

            uxda = uxr.UxDataArray(
                data=data,
                dims=dims,
                name=name,
                attrs=attrs,
                uxgrid=self.surface.uxgrid,
            )

            self.surface._uxds[name] = uxda
            if tag is None:
                return

        insert_layer = -1

        # Find the first layer that contains a contiguous set of empty spots
        for attempt in range(2):
            for i in range(_N_TAG_LAYERS):
                if np.all(self.uxds[name].data[:, i] == 0):
                    insert_layer = i
                    break
            if insert_layer == -1:  # If at first we don't succed, try sorting the tags and try again.
                if attempt == 0:
                    self.sort_tags(name)
                    break
                else:
                    # TODO: Implement an option to expand the number of layers if all layers are full
                    warnings.warn(f"All {name} layers are full. Overwriting the first layer.", stacklevel=2)
                    insert_layer = 0
        data = self.uxds[name].data
        data[:, insert_layer] = tag
        self.surface.uxds[name].data[self.face_indices, :] = data
        self._desloped_face_elevation = None  # Invalidate the desloped elevation cache since tags are often used for tracking craters, and we want to be able to get the original elevation of the crater faces if we need to.

        return

    def remove_tag(self, name: str, tag: int) -> None:
        """
        Removes an integer tag from the local surface.

        Parameters
        ----------
        tag : int
            The integer tag to be removed.
        """
        if name not in self.uxds:
            raise ValueError(f"Tag {name} does not exist on the surface")

        data = self.uxds[name].data
        if tag in data:
            data[data == tag] = 0
            self.surface.uxds[name].data[self.face_indices, :] = data

        return

    def sort_tags(self, name: str) -> None:
        """
        Sort the tag layers in descending order based on the tag values. This is used to ensure that the most recently added tags are in the first available layer, which is important for the add_tag method to function properly.

        Parameters
        ----------
        name : str
            The name of the tag to be sorted.
        """
        ds = self.surface.uxds[name].sel(n_face=self.face_indices)
        sorted_indices = np.argsort(ds.data, axis=1)[:, ::-1]
        ds = ds.isel(layer=xr.DataArray(sorted_indices, dims=["n_face", "layer"]))
        self.surface.uxds[name].data[self.face_indices, :] = ds.data
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

        delta_face_elevation = surface_bindings.apply_diffusion(face_kappa=kdiff, face_variable=self.face_elevation, region=self)
        self.update_elevation(delta_face_elevation)
        self.add_data("ejecta_thickness", long_name="ejecta thickness", units="m", data=delta_face_elevation, positive_only=True)
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

        delta_face_elevation = surface_bindings.slope_collapse(critical_slope=critical_slope, region=self)
        self.update_elevation(delta_face_elevation)
        self.add_data("ejecta_thickness", long_name="ejecta thickness", units="m", data=delta_face_elevation, positive_only=True)

    def compute_slope(self) -> NDArray[np.float64]:
        """
        Compute the slope of the surface.

        Returns
        -------
        NDArray[np.float64]
            The slope of all faces in degrees.
        """
        slope = surface_bindings.compute_slope(
            face_elevation=self.face_elevation,
            region=self,
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
            noise = surface_bindings.turbulence_noise(
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
        self, location: tuple[float, float] | None = None, validate: bool = True
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Computes the distances between nodes and faces and a given location.

        Parameters
        ----------
        location : tuple[float, float], option
            tuple containing the longitude and latitude of the location in degrees. If None, the location of the view center is used if it is set.
        validate : bool, optional
            If the location should be validated and normalized. Only use if the validation is causing a performance bottleneck.

        Returns
        -------
        NDArray
            Array of face distances in meters.
        NDArray
            Array of node distances in meters.
        """
        if location is None:
            if self.is_local:
                location = self.location
                validate = False
            else:
                raise ValueError("A value for location must be provided for global surfaces.")

        if len(location) == 1:
            location = location.item()
        if len(location) != 2:
            raise ValueError("location must be a single pair of (longitude, latitude).")
        if validate:
            location = validate_and_normalize_location(location)
        node_locations = np.vstack((self.node_lon, self.node_lat)).T
        face_locations = np.vstack((self.face_lon, self.face_lat)).T
        return self.compute_distances(
            locations=face_locations, reference_location=location, validate=False
        ), self.compute_distances(locations=node_locations, reference_location=location, validate=False)

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
            Array of initial bearings for each face in degrees.
        NDArray
            Array of initial bearings for each node in degrees.
        """
        if location is None:
            if self.is_local:
                location = self.location
            else:
                raise ValueError("A value for location must be provided for global surfaces.")

        if len(location) == 1:
            location = location.item()
        if len(location) != 2:
            raise ValueError("location must be a single pair of (longitude, latitude).")
        location = validate_and_normalize_location(location)
        face_bearing, node_bearing = (
            surface_bindings.compute_bearings(
                lon1=np.radians(location[0]),
                lat1=np.radians(location[1]),
                lon2=np.radians(self.face_lon),
                lat2=np.radians(self.face_lat),
            ),
            surface_bindings.compute_bearings(
                lon1=np.radians(location[0]),
                lat1=np.radians(location[1]),
                lon2=np.radians(self.node_lon),
                lat2=np.radians(self.node_lat),
            ),
        )
        return np.degrees(face_bearing), np.degrees(node_bearing)

    def interpolate_node_elevation_from_faces(self) -> None:
        """
        Update node elevations by area-weighted averaging of adjacent face elevations.

        For each node, the elevation is computed as the area-weighted average of the elevations
        of the surrounding faces.

        Returns
        -------
        None
        """
        node_elevation = surface_bindings.interpolate_node_elevation_from_faces(region=self)
        self.add_data(name="node_elevation", data=node_elevation, overwrite=True)
        return

    def get_reference_surface(
        self, reference_radius: float | None = None, only_faces: bool = False, **kwargs: Any
    ) -> NDArray[np.float64]:
        """
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.

        Parameters
        ----------
        reference_radius : float
            The radius of the reference region to compute the average over in meters. If None, the region_radius of the LocalSurface is used.
        only_faces : bool
            If True, only return the face elevation data of the reference surface. Default is False.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        NDArray[np.float64]
            The face and node elevation points of the reference sphere, or the original elevation points is the reference region is too small
        """

        def _find_reference_elevations(x: NDArray, y: NDArray, z: NDArray) -> NDArray:
            """
            Compute the mean plane that fits the points given by the projected x, projected y, and elevation array.

            Parameters
            ----------
            region_points : NDArray
                An array of shape (n, 3) where each row contains [x, y, elevation].

            Returns
            -------
            NDArray
                An array of elevation values corresponding to the mean plane for each projected x and y point.
            """
            # Create the design matrix for the plane fitting
            A = np.vstack([x, y, np.ones(len(x))]).T

            # Solve for the coefficients of the plane (Ax + By + C = z)
            coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)

            return coeffs[0] * x + coeffs[1] * y + coeffs[2]

        if reference_radius is None:
            reference_radius = self.region_radius
        # Find cells within the crater radius
        faces_within_region = self.face_distance <= reference_radius
        if only_faces:
            points_within_region = faces_within_region
            elevation = self.face_elevation
            n_reference = np.sum(faces_within_region)
        else:
            elevation = np.concatenate([self.face_elevation, self.node_elevation])
            nodes_within_region = self.node_distance <= reference_radius
            if not nodes_within_region.any():
                return elevation
            points_within_region = np.concatenate([faces_within_region, nodes_within_region])
            n_reference = np.sum(faces_within_region) + np.sum(nodes_within_region)

        # if there are not points within the radius, return the original elevations and use that as the reference
        if n_reference < 5:
            return elevation

        self.set_face_proj()
        if only_faces:
            x, y, z = self.face_proj_x, self.face_proj_y, self.face_elevation
        else:
            self.set_node_proj()
            x = np.concatenate([self.face_proj_x, self.node_proj_x])
            y = np.concatenate([self.face_proj_y, self.node_proj_y])
            z = np.concatenate([self.face_elevation, self.node_elevation])
        reference_elevation = elevation
        reference_elevation[points_within_region] = _find_reference_elevations(
            x[points_within_region], y[points_within_region], z[points_within_region]
        )

        return reference_elevation

    def extract_subregion(self, subregion_radius: FloatLike, **kwargs) -> LocalSurface | None:
        """
        Extract a subset of the LocalSurface region with a smaller radius than the original region.

        Parameters
        ----------
        subregion_radius : float
            The radius of the subregion to extract in meters.
        **kwargs
            |kwargs|

        Returns
        -------
        LocalSurface
            A LocalSurface object containing a view of the regional grid.

        """
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
            **vars(self.common_args),
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

    def compute_radial_gradient(self, variable: str | NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute the radial gradient of the local surface variable with respect to the local_location center.

        Parameters
        ----------
        variable : str or NDArray[np.float64]
            The variable to compute the radial gradient of. This can be a string (the name of a variable in the surface data) or an array of values the same size as the number of faces in the grid.

        Returns
        -------
        NDArray[np.float64]
            The radial gradient of all faces in meters per meter.
        """
        if isinstance(variable, str):
            if variable not in self.uxds:
                raise ValueError(f"Variable {variable} not found in the surface data.")
            variable = self.uxds[variable].data
        elif variable.size != self.n_face:
            raise ValueError("variable must be a string or an array with the same size as the number of faces in the grid")
        radial_gradient = surface_bindings.compute_radial_gradient(
            variable=variable,
            region=self,
        )

        return radial_gradient

    def export(
        self,
        driver: str = "GPKG",
        interval: int | None = None,
        ask_overwrite: bool | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Export the surface view data to a specified file format and stores it in the default export directory.

        If the format is "VTK," the data will be exported using the :py:meth:`Surface.to_vtk` method. Otherwise, it will use the :py:meth:`Surface.to_vector_file` to export the data.

        Parameters
        ----------
        driver : str, optional
            The driver to use export the data to. Supported formats are 'VTK' or a driver supported by GeoPandas ('GPKG', 'ESRI Shapefile', etc.), and 'GeoTIFF'.
        interval : int | None, optional
            |interval_export|
        ask_overwrite : bool | None, optional
            |ask_overwrite_methods|
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP

        # Temporarily set the ask_overwrite attribute for the duration of the export, but reset it to its original value afterwards.
        ask_overwrite_orig = self.ask_overwrite
        if ask_overwrite is not None:
            self.ask_overwrite = ask_overwrite

        if interval is not None:
            self.surface.save(interval=interval, **kwargs, skip_actions=True)
        if driver.upper() in ["VTK", "VTP"]:
            self.to_vtk_file(
                interval=interval,
                **kwargs,
            )
        elif driver.upper() in ["GEOTIFF", "GTIFF", "TIFF", "TIF"]:
            self.to_geotiff_file(
                interval=interval,
                **kwargs,
            )
        elif driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP:
            self.to_vector_file(
                driver=driver,
                interval=interval,
                **kwargs,
            )
        self.ask_overwrite = ask_overwrite_orig
        return

    def to_vector_file(
        self,
        driver: str = "GPKG",
        interval: int | None = None,
        to_file_kwargs: dict | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Export the face-associated data from the surface view data to a vector file using GeoPandas.

        See `geopandas.GeoDataFrame.to_file <https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html>`_ for more detailed information on the available parameters.

        Parameters
        ----------
        driver : str, optional
            The file format driver to use for exporting. Default is 'GPKG'.
        interval : int | None, optional
            |interval_export|
        to_file_kwargs : dict | None, optional
            Additional keyword arguments to pass to the GeoDataFrame.to_file method.
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP

        def _write_dataset(uxds, filename, layer_name, driver, **kwargs):
            if self.is_global:
                # Exclude periodic elements for global surfaces. For some reason, UxArray gets the windings wrong when using the split argument, so we have to exclude instead.
                gdfargs = {"engine": "geopandas", "periodic_elements": "exclude"}
            else:
                gdfargs = {"engine": "geopandas", "periodic_elements": "ignore"}

            # Convert to GeoDataFrame and set the CRS correctly for the type of surface
            gdf = uxds.uxgrid.to_geodataframe(**gdfargs).set_crs(self.surface.crs).to_crs(self.crs)
            variables = [v for v in uxds.data_vars if any(dim == "n_face" for dim in uxds[v].dims)]
            if not variables:
                raise ValueError("No face-based variables found to export to file.")

            # Shapefile can't store the crater id layers, so we remove variables with the pattern "idXX"
            if driver == "ESRI Shapefile":
                variables = [var for var in variables if not var.startswith("id")]

            for var in variables:
                if len(uxds[var].dims) == 1:
                    gds = uxds[var].to_geodataframe(**gdfargs)
                    gdf[var] = gds[var]
                elif "layer" in uxds[var].dims:
                    for layer in range(uxds.layer.size):
                        gds = uxds[var].isel(layer=layer).to_geodataframe(**gdfargs)
                        gdf[f"{var}{layer:02d}"] = gds[var]

            # Rename columns if exporting to Shapefile to comply with field name length limits
            if driver == "ESRI Shapefile":
                gdf = gdf.rename(columns={col: shp_key_fix(col) for col in gdf.columns})

            print(f"Exporting to {filename} using driver {driver}")
            if not self._overwrite_check(filename):
                return
            gdf.to_file(filename, layer=layer_name, driver=driver, **kwargs)
            return

        def shp_key_fix(key: str) -> str:
            """
            ESRI Shapefile format limits field names to 10 characters, so this function substitues longer names with shorter alternatives, truncates the results, and sets them to upper case.
            """
            alt_names = {
                "face_elevation": "faceelev",
                "face_distance": "facedis",
                "face_bearing": "facebearng",
                "face_indices": "faceindx",
            }
            for long, short in alt_names.items():
                if long in key:
                    key = key.replace(long, short)
            return key[:10].upper()

        # Common alias for Shapefile
        if driver.upper() == "SHP":
            driver = "ESRI Shapefile"

        if driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP:
            file_extension = VECTOR_DRIVER_TO_EXTENSION_MAP[driver.upper()]
        else:
            raise ValueError("Cannot infer file extension from driver {driver}.")

        uxds = self.read_saved_output(interval=interval, reset=False)
        interval_numbers = uxds.interval.values

        for interval in interval_numbers:
            uxdsi = uxds.sel(interval=interval).load()
            filename = self.export_dir / self.output_filename(interval).replace(self.output_file_extension, file_extension)
            _write_dataset(
                uxdsi,
                filename=filename,
                layer_name="face_data",
                driver=driver,
                **to_file_kwargs if to_file_kwargs is not None else {},
            )
        return

    def to_vtk_mesh(self, uxds: UxDataset | None = None, interval: int | None = None, **kwargs: Any) -> vtkUnstructuredGrid:
        """
        Exports the regional mesh to a VTK PolyData object.

        Parameters
        ----------
        uxds : UxDataset, optional
            The dataset to export. If None, the method will use currently loaded data in the surface. Default is None.
        interval : int, optional
            The interval number to export. If provided, the method will either extract the interval number from uxds (if it has intervals saved), or, if no uxds is passed, load a saved interval from file.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        vtkUnstructuredGrid
            The VTK PolyData object representing the regional mesh.
        """
        from vtk import (
            VTK_POLYGON,
            vtkPoints,
        )
        from vtkmodules.util.numpy_support import numpy_to_vtk
        from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
        from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter

        if uxds is None:
            if interval is None:
                uxds = self.uxds
            else:
                uxds = self.read_saved_output(interval=interval, reset=False).isel(interval=0)
                interval = None

        if interval is not None:
            if "interval" not in uxds:
                raise ValueError("uxds does not have an 'interval' variable")
            uxds = uxds.isel(interval=interval)

        node_xyz = np.c_[self.node_x, self.node_y, self.node_z]
        node_normals = node_xyz / self.surface.radius
        vtk_point_normals = numpy_to_vtk(node_normals.astype(np.float32), deep=True)
        vtk_point_normals.SetNumberOfComponents(3)
        vtk_point_normals.SetName("Normals")

        # Warp the mesh according to node elevation if it exists
        if "node_elevation" in uxds:
            warped_xyz = node_xyz + uxds.node_elevation.data[:, None] * node_normals
        else:
            warped_xyz = node_xyz

        node_x = warped_xyz[:, 0]
        node_y = warped_xyz[:, 1]
        node_z = warped_xyz[:, 2]

        vtk_data = vtkUnstructuredGrid()
        nodes = vtkPoints()
        for i in range(self.n_node):
            nodes.InsertNextPoint(node_x[i], node_y[i], node_z[i])
        vtk_data.SetPoints(nodes)
        vtk_data.Allocate(self.n_face)
        for i, n in enumerate(self.n_nodes_per_face):
            point_ids = self.face_node_connectivity[i][0:n]
            vtk_data.InsertNextCell(VTK_POLYGON, n, point_ids)

        grid = vtkUnstructuredGrid()
        grid.DeepCopy(vtk_data)

        for v in uxds.variables:
            array = numpy_to_vtk(uxds[v].values, deep=True)
            array.SetName(v)
            if "n_face" in uxds[v].dims:
                grid.GetCellData().AddArray(array)
            elif "n_node" in uxds[v].dims:
                grid.GetPointData().AddArray(array)
                if v == "node_elevation":
                    grid.GetPointData().SetActiveScalars(v)
            elif uxds[v].dims == ("time",) or uxds[v].size == 1:
                grid.GetFieldData().AddArray(array)

        geom_filter = vtkGeometryFilter()
        geom_filter.SetInputData(grid)
        geom_filter.Update()
        poly_data = geom_filter.GetOutput()

        poly_data.GetPointData().SetNormals(vtk_point_normals)
        normals_filter = vtkPolyDataNormals()
        normals_filter.SetInputData(poly_data)
        normals_filter.ComputeCellNormalsOn()
        normals_filter.Update()
        mesh = normals_filter.GetOutput()

        return mesh

    def to_vtk_file(
        self,
        interval: int | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Export the regional mesh to a VTK file and store it in the default export directory.

        Parameters
        ----------
        interval : int, optional
            |interval_export|
        **kwargs : Any
            |kwargs|
        """

        def _write_current_mesh(mesh, output_filename):
            from vtk import vtkXMLPolyDataWriter

            writer = vtkXMLPolyDataWriter()
            writer.SetDataModeToBinary()
            writer.SetCompressorTypeToZLib()
            if not self._overwrite_check(output_filename):
                return
            if output_filename.exists():
                output_filename.unlink()

            writer.SetFileName(output_filename)
            writer.SetInputData(mesh)
            print(f"Exporting to {output_filename}")
            writer.Write()
            return

        uxds = self.read_saved_output(interval=interval, reset=False)
        interval_numbers = uxds.interval.values
        if interval is not None:
            if interval < 0:
                interval = interval_numbers[interval]
            interval_numbers = [interval]

        # Check if we need to save the geometry file
        grid_filename = self.export_dir / f"{self._grid_file_prefix}.vtp"

        if not grid_filename.exists():
            grid = self.to_vtk_mesh(uxds=self.uxgrid.to_xarray())
            _write_current_mesh(grid, grid_filename)

        for interval in interval_numbers:
            uxdsi = uxds.sel(interval=interval).load()
            mesh = self.to_vtk_mesh(uxds=uxdsi)

            filename = self.export_dir / self.output_filename(interval=interval).replace(self.output_file_extension, "vtp")
            _write_current_mesh(mesh, filename)

        return

    def to_raster(
        self, uxda: UxDataArray | None = None, **kwargs: Any
    ) -> tuple[NDArray[np.float32], tuple[float, float, float, float], Any, CRS]:
        """
        Rasterize a face-based variable into a 2D raster using rasterio.

        Parameters
        ----------
        uxda : UxDataArray | None
            The UxDataArray containing the face-based variable to rasterize. If None, the method will use currently loaded data in the surface with "face_elevation" as the default variable. Default is None.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        raster : NDArray[np.float32]
            The rasterized variable as a 2D numpy array.
        extent : tuple[float, float, float, float]
            The extent of the raster in the format (xmin, xmax, ymin, ymax).
        transform : Affine
            The affine transform for the raster.
        crs : CRS
            The coordinate reference system of the raster.
        """
        # Check if rasterio is installed, and if not, just return without plotting
        try:
            from rasterio.features import rasterize
            from rasterio.transform import Affine, from_bounds
        except ImportError:
            warnings.warn(
                "rasterio is not installed. Cannot generate plot. On some platforms, you may need to install GDAL first before installing rasterio.",
                stacklevel=2,
            )
            return

        if uxda is None:
            uxda = self.uxds["face_elevation"].load()

        W, H = self.get_raster_dims()
        if self.is_global:
            # Splitting doesn't work well and makes a hash of the raster. So we'll just drop the periodic elements instead
            gdf = uxda.to_geodataframe(engine="geopandas", periodic_elements="exclude").set_crs(self.crs)
            xmin, xmax = -180.0, 180.0
            ymin, ymax = -90.0, 90.0
            # Get the approximate pixel size based on the number of faces
            npix = self.surface.n_face
            pix_area = 4.0 * np.pi * (self.surface.radius**2) / npix
            pix = np.sqrt(pix_area)
            deg_per_pix = 180.0 * pix / (np.pi * self.surface.radius)
            xres = yres = deg_per_pix

            transform = from_bounds(xmin, ymin, xmax, ymax, W, H)
        else:
            gdf = uxda.to_geodataframe(engine="geopandas", periodic_elements="ignore").set_crs(self.surface.crs).to_crs(self.crs)
            pix = self.pix
            R = (
                self.region_radius + 4 * pix
            )  # Add a 4 * pix border to enswure that the pixels along the border of the local region are plotted in full
            xmin, xmax = -R, R
            ymin, ymax = -R, R
            # pixel size (meters/pixel)
            xres = yres = (2.0 * R) / W
            # upper-left at (-R, +R); y increases downward in rasters
            transform = Affine.translation(-R, R) * Affine.scale(xres, -yres)

        vals = gdf[uxda.name].to_numpy()
        geoms = gdf.geometry.values
        shapes = [
            (geom, float(val))
            for geom, val in zip(geoms, vals, strict=False)
            if (geom is not None) and (not geom.is_empty) and np.isfinite(val)
        ]
        extent = (xmin, xmax, ymin, ymax)

        raster = rasterize(
            shapes=shapes,
            out_shape=(H, W),
            transform=transform,
            fill=np.nan,
            dtype=np.float32,
            all_touched=True,
        )
        return raster, extent, transform, gdf.crs

    def get_raster_dims(self) -> tuple[int, int]:
        """
        Get the dimensions of the raster based on the region radius (or target, if this is a global surface) and pixel size.

        Returns
        -------
        tuple[int, int]
            The width and height of the raster in pixels.
        """
        if self.is_global:
            # Splitting doesn't work well and makes a hash of the raster. So we'll just drop the periodic elements instead
            xmin, xmax = -180.0, 180.0
            ymin, ymax = -90.0, 90.0
            # Get the approximate pixel size based on the number of faces
            npix = self.surface.n_face
            pix_area = 4.0 * np.pi * (self.surface.radius**2) / npix
            pix = np.sqrt(pix_area)
            deg_per_pix = 180.0 * pix / (np.pi * self.surface.radius)
            xres = yres = deg_per_pix
            W = int(np.ceil((xmax - xmin) / xres))
            H = int(np.ceil((ymax - ymin) / yres))
        else:
            R = self.region_radius
            xmin, xmax = -R, R
            ymin, ymax = -R, R
            pix = self.pix
            W = H = int(max(1, np.ceil((2.0 * R) / pix)))
        return W, H

    def to_geotiff_file(
        self,
        interval: int | None = None,
        variable_name: str = "face_elevation",
        **kwargs,
    ) -> None:
        """
        Rasterize a face-based elevation variable into a GeoTIFF using rasterio.

        Parameters
        ----------
        interval : int, optional
            |interval_export|
        variable_name : str, optional
            The name of the variable to rasterize. Default is "face_elevation".
        """
        import matplotlib.pyplot as plt
        import rasterio as rio
        from cartopy import crs as ccrs

        def _write_dataset(uxda, output_filename, **kwargs):
            if self.is_global:
                projection = ccrs.PlateCarree()
            else:
                projection = ccrs.AzimuthalEquidistant(central_longitude=self.location[0], central_latitude=self.location[1])

            raster, extent, transform, crs = self.to_raster(uxda, **kwargs)

            H, W = raster.shape

            fig = plt.figure(figsize=(W, H), dpi=1)

            ax = fig.add_axes([0, 0, 1, 1], projection=projection)
            ax.set_axis_off()
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            ax.set_extent(extent, crs=projection)

            profile = {
                "driver": "GTiff",
                "height": H,
                "width": W,
                "count": 1,
                "dtype": raster.dtype,
                "crs": crs,
                "transform": transform,
                "nodata": np.nan,
            }
            print(f"Exporting to {output_filename}")
            with rio.open(output_filename, "w", **profile) as dst:
                dst.write(raster, 1)
            return

        # load data and select the face-based variables
        uxds = self.read_saved_output(interval=interval, reset=False)
        interval_numbers = uxds.interval.values
        file_extension = "tif"

        for interval in interval_numbers:
            uxda = uxds.sel(interval=interval)[variable_name].load()
            filename = self.export_dir / f"{self._output_file_prefix}{interval:06d}.{file_extension}"
            _write_dataset(
                uxda,
                filename,
                **kwargs,
            )

        return

    def save(
        self,
        interval: int = 0,
        time_variables: dict | None = None,
        filename: Path | str | None = None,
        **kwargs,
    ) -> None:
        """
        Save the region surface data to the specified directory.

        Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the time dimension.

        Parameters
        ----------
        interval : int, optional
            Interval number to append to the data file name. Default is 0.
        time_variables : dict, optional
            Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the time dimension. Default is None.
        filename : Path | str, optional
            The path to the output file. If None, the file will be saved to the default
        **kwargs : Any
            |kwargs|
        """
        self.surface.output_dir.mkdir(parents=True, exist_ok=True)
        if interval is None:
            interval = 0
        if time_variables is not None and not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")

        # Sort any layers that we might have before saving
        for var in self.uxds.data_vars:
            if "layer" in self.uxds[var].dims:
                self.sort_tags(name=var)

        self.uxds.close()

        ds = self.uxds.expand_dims(dim="interval").assign_coords({"interval": [interval]})
        if time_variables is not None:
            for k, v in time_variables.items():
                ds[k] = xr.DataArray(data=[v], name=k, dims=["interval"], coords={"interval": [interval]})

        self.surface._save_data(
            ds,
            interval=interval,
            filename=filename,
        )

        if self.is_local:  # Save the local grid if this is a local surface
            self._write_grid_file()

        save_args = {
            "interval": interval,
            "time_variables": time_variables,
            "filename": filename,
            **kwargs,
        }
        super().save(**save_args)

        return

    def _write_grid_file(self, uxgrid: uxr.Grid | None = None, grid_file: Path | str | None = None, **kwargs: Any) -> None:
        """
        Write the grid to a NetCDF file.

        Parameters
        ----------
        uxgrid : uxr.Grid, optional
            The grid to be written to the file. If None, the grid will be obtained from the surface object. Default is None.
        grid_file : Path, str, optional
            The path to the grid file. If None, the path will be obtained from the surface object. Default is None.
        **kwargs : Any
            |kwargs|
        """
        if uxgrid is None:
            uxgrid = self.uxgrid
        if grid_file is None:
            grid_file = self.surface.output_dir / f"{self._grid_file_prefix}.{self._output_file_extension}"
        else:
            grid_file = Path(grid_file)
        self.surface._write_grid_file(uxgrid=uxgrid, grid_file=grid_file, **kwargs)

        return

    def plot(
        self,
        filename: Path | str | None = None,
        plot_style: Literal["map", "hillshade"] = "map",
        variable_name: str | None = None,
        interval: int | None = None,
        cmap: str | None = None,
        label: str | None = None,
        scalebar: bool | None = None,
        colorbar: bool = True,
        show: bool = True,
        save: bool = False,
        ax: Axes | None = None,
        close_when_done: bool | None = None,
        minimum_plot_width: float | None = 800,
        **kwargs: Any,
    ) -> Axes:
        """
        Plot an image of the local region.

        Parameters
        ----------
        filename : Path | str | None, optional
            The path to save the plot to. If None, and save is True, the plot will be saved to the default plot directory with a filename based on the interval number. Default is None.
        plot_style : str, optional
            The style of the plot. Options are "map" and "hillshade". In "map" mode, the variable is displayed as a colored map. In "hillshade" mode, a hillshade image is generated using "face_elevation" data. If a different variable is passed to `variable`, then the hillshade will be overlayed with that variable's data. Default is "map".
        variable_name : str | None, optional
            The variable to plot. If None is provided then "face_elevation" is used in "map" mode.
        interval: int | None, optional
            The interval number of the surface to plot. If None, the currently loaded surface data will be used.
        cmap : str, optional
            The colormap to use for the plot. If None, a default colormap will be used ("cividis" by default and "grey" when plot_style=="hillshade" and variable=="face_elevation").
        label : str | None, optional
            A label for the plot. If None, no label will be added.
        scalebar : bool, optional
            If True, a scalebar will be added to the plot. Default is True.
        colorbar : bool, optional
            If True, a colorbar will be added to the plot when using "map" plot_style or "hillshade" with a variable overlay. Default is True.
        show : bool, optional
            If True, the plot will be displayed. Default for local surfaces is True.
        save : bool, optional
            If True, the plot will be saved to the default plot directory. Default is False unless filename is provided, in which case the default is True.
        ax : matplotlib.axes.Axes, optional
            An existing Axes object to plot on. If None, a new figure and axes will be created.
        close_when_done : bool, optional
            If True, the figure will be closed after plotting. Default is True when save is True and show is False, and False otherwise.
        minimum_plot_width : float, optional
            Because the width of the plot is determined by the number of faces, small regions will generate small plots with labels that are hard to read. This parameter sets a lower limit to the width of the image that is generated by the plot. By default it is 800. Set to None to turn it off.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        matplotlib.axes.Axes
            A Matplotlib Axes object containing the plot.
        """
        import matplotlib.pyplot as plt
        from matplotlib import rc_context, rcParams
        from matplotlib.colors import LightSource, Normalize
        from scipy.ndimage import gaussian_filter

        if scalebar is None:
            scalebar = True
        interval_numbers = self.surface.get_saved_interval_numbers()[0]
        if interval is None or interval == interval_numbers[-1]:
            uxds = self.uxds
        else:
            uxds = self.read_saved_output(interval=interval, reset=False)
            interval = uxds.interval.values.item()
            uxds = uxds.sel(interval=interval)

        if variable_name is not None and variable_name not in uxds:
            raise ValueError(f"Variable '{variable_name}' not found in the surface data.")
        if variable_name is None and plot_style == "map":
            variable_name = "face_elevation"
        if plot_style == "hillshade":
            if variable_name is None:
                do_overlay = False
            else:
                do_overlay = True
        if close_when_done is None:
            close_when_done = save and not show
        colorbar = colorbar and (plot_style == "map" or do_overlay)

        if variable_name is not None:
            ret = self.to_raster(uxds[variable_name].load())
            variable_long_name = uxds[variable_name].attrs.get("long_name", variable_name)
            variable_units = uxds[variable_name].attrs.get("units", "")
            variable_raster = ret[0]
            extent = ret[1]
            H, W = variable_raster.shape
            vmin = kwargs.pop("vmin", np.nanmin(variable_raster))
            vmax = kwargs.pop("vmax", np.nanmax(variable_raster))
            norm = Normalize(vmin=vmin, vmax=vmax)
        else:
            variable_long_name = ""
            variable_units = ""
            vmin = kwargs.pop("vmin", 0.0)
            vmax = kwargs.pop("vmax", 1.0)

        if plot_style == "hillshade":
            hill_args = {"dx": self.pix, "dy": self.pix, "fraction": 1.0}
            ret = self.to_raster(uxds["face_elevation"].load())
            elevation = ret[0]
            extent = ret[1]
            elevation = gaussian_filter(elevation, sigma=2, mode="constant", cval=np.nan)
            H, W = elevation.shape
            azimuth = 300.0
            solar_angle = 20.0
            ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
            if do_overlay:
                if cmap is None:
                    cmap = "cividis"
                cmap = plt.get_cmap(cmap)
                variable_raster = np.clip((variable_raster - vmin) / (vmax - vmin), 0.0, 1.0)
                rgb = cmap(variable_raster)
                blended = ls.shade_rgb(rgb, elevation, blend_mode="overlay", **hill_args)
                if np.any(np.isnan(variable_raster[~np.isnan(elevation)])):
                    hillshade = ls.hillshade(elevation, **hill_args)
                    graymap = plt.get_cmap("gray")
                    hillshade = graymap(hillshade)
                    cvals = np.where(blended == 0, hillshade, blended)
                else:
                    cvals = blended
            else:
                if cmap is None:
                    cmap = "gray"
                cvals = ls.hillshade(elevation, **hill_args)
            interpolation = kwargs.pop("interpolation", "lanczos")
        elif plot_style == "map":
            if cmap is None:
                cmap = "cividis"
            cvals = variable_raster
            interpolation = kwargs.pop("interpolation", "bicubic")
        else:
            raise ValueError("plot_style must be either 'map' or 'hillshade'")

        if minimum_plot_width is not None:
            dpi = max(minimum_plot_width, W)

        # Plot with (1, 1) inch figure and dpi=resolution for exact pixel size
        if ax is None:
            _, ax = plt.subplots(figsize=(1, 1), dpi=dpi, frameon=False, layout="constrained")

        # This helps scale plot elements that are sized in points, because we're using a somewhat non-standard dpi and figure size to plot the rasterio-genrated surface data
        scale_fraction = 72 / dpi
        scaled_params = {
            k: v * scale_fraction
            for k, v in rcParams.items()
            if v is not None and type(v) is float and ("pad" in k or "width" in k or "size" in k)
        }
        scaled_params["font.size"] = 24 * scale_fraction

        with rc_context(scaled_params):
            ax.imshow(cvals, interpolation=interpolation, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal", extent=extent)
            ax.axis("off")
            xmin, xmax, ymin, ymax = extent
            # Add scale bar before saving/showing image
            if scalebar and self.is_local:
                # Determine max physical size for the scale bar
                max_physical_size = xmax / 2 / np.sqrt(2)

                # Choose "nice" scale bar length
                nice_values = np.array([1, 10, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000])  # in meters
                scale_length = nice_values[nice_values <= max_physical_size].max()
                bar_height = 0.2 * scale_length
                scale_text = f"{int(scale_length)} m" if scale_length < 1000 else f"{int(scale_length / 1000)} km"

                # Position in lower right corner
                x_start = xmax - scale_length + xmin * 0.05
                y_start = -(ymax - bar_height + ymin * 0.05)

                rect = plt.Rectangle((x_start, y_start), scale_length, bar_height, color="black")
                ax.add_patch(rect)
                # Label above the scale bar
                ax.text(
                    x_start + scale_length / 2,
                    y_start + 2 * bar_height,
                    scale_text,
                    color="black",
                    ha="center",
                    va="bottom",
                    fontweight="bold",
                )
            if label:
                # Label in the upper left corner
                x_start = xmin
                y_start = ymax
                if self.is_local:
                    va = "top"  # Local surfaces are a circle with empty corners, so we can put the label in the upper left corner without it overlapping the surface image.
                else:
                    va = "bottom"  # Global surfaces fill most of the entire image, so we need to place the label above the plotting area to avoid overlap.
                ax.text(
                    x_start,
                    y_start,
                    label,
                    color="black",
                    ha="left",
                    va=va,
                    fontweight="bold",
                )
            if colorbar:
                norm = Normalize(vmin=vmin, vmax=vmax)
                sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
                cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.5)
                cbar.set_label(label=f"{variable_long_name} [{variable_units}]")

            if save:
                file_prefix = f"{self.output_file_prefix}_{plot_style}"
                if variable_name is not None:
                    file_prefix += f"_{variable_name}"
                if interval is None:
                    uxds = self.uxds
                    filename = self.plot_dir / f"{file_prefix}.{self.output_image_file_extension}"
                else:
                    uxds = self.read_saved_output(interval=interval, reset=False)
                    interval = uxds.interval.values.item()
                    uxds = uxds.sel(interval=interval)
                    filename = self.plot_dir / f"{file_prefix}{interval:06d}.{self.output_image_file_extension}"

                plt.savefig(filename, pad_inches=0, dpi=ax.figure.get_dpi())
            if show:
                plt.show()
            if close_when_done and (save or show):
                plt.close()
        return ax

    def pyvista_plotter(
        self,
        variable_name: str | None = None,
        variable: ArrayLike | None = None,
        interval: int | None = None,
        theme: str | None = None,
        transparent_background: bool | None = None,
        plotter: pv.Plotter | None = None,
        enable_interactive: bool = True,
        **kwargs: Any,
    ) -> pv.Plotter:
        """
        Show the local surface region using an interactive 3D plot with PyVista.

        Parameters
        ----------
        variable_name : str | None, optional
            The name of the variable to plot. If the name of the variable is not already stored on the surface mesh, then the `variable` argument must also be passed. Default is None, which will plot a greyscale image of the surface.
        variable : (n_face) array, optional
            An array face values that will be used to color the surface mesh. This is required if `variable_name` is not a face variable that is already saved in the the uxds dataset. Default is None.
        interval : int | None, optional
            The interval number of the surface to plot. If None, the currently loaded surface data will be used. Default is None.
        theme : str, optional
            The PyVista plot theme to use. If None, the default PyVista theme will be used. Default is None.
        transparent_background : bool, optional
            If True, the background of the plot will be transparent. Default is None, which will use the default background setting for the chosen plot theme.
        plotter : pyvista.Plotter, optional
            A pre-existing Plotter object to use. If None, then a new one will be created and returned. Default is None.
        enable_interactive : bool, optional
            If True, the default PyVista key events will be updated to include custom events for toggling scalar visibility, changing the camera view, and showing a help message. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        pyvista.Plotter
            The PyVista Plotter object for further customization.
        """
        from cratermaker.constants import PYVISTA_ADD_MESH_KWARGS
        from cratermaker.utils.general_utils import toggle_pyvista_actor, update_pyvista_help_message

        def _set_title(uxds, variable_name):
            if variable_name in uxds and "long_name" in uxds[variable_name].attrs:
                if "units" in uxds[variable_name].attrs:
                    title = f"{uxds[variable_name].attrs['long_name']} ({uxds[variable_name].attrs['units']})"
                else:
                    title = uxds[variable_name].attrs["long_name"]
            else:
                title = variable_name
            return title

        def update_scalars(plotter, cmap):
            camera_orig = plotter.camera
            mesh_actor = plotter.actors[self.target.name]
            if not plotter.scalar_bars:
                scalar_bar_actor = None
            else:
                old_title = next(iter(plotter.scalar_bars.keys()))
                scalar_bar_actor = plotter.scalar_bars[old_title]
            scalar_names = [""]  # Add entry for no scalars
            for v in mesh_actor.mapper.dataset.array_names:
                if v in mesh_actor.mapper.dataset.cell_data:
                    scalar_names.append(v)
            scalar_names = list(set(scalar_names))  # Remove duplicates
            scalar_names.sort()
            if mesh_actor.mapper.GetScalarVisibility():
                idx = scalar_names.index(mesh_actor.mapper.array_name)
                idx += 1
                if idx == len(scalar_names):
                    idx = 0
            else:  # We are currently not showing scalars, so start over from the first non-empty slot
                idx = 1

            next_scalar_name = scalar_names[idx]
            if next_scalar_name == "":
                mesh_actor.mapper.SetScalarVisibility(False)
                if scalar_bar_actor is not None:
                    scalar_bar_actor.SetVisibility(False)
            else:
                mesh_actor.mapper.SetScalarVisibility(True)
                mesh_actor.mapper.array_name = next_scalar_name
                mesh_actor.mapper.dataset.active_scalars_name = next_scalar_name
                mesh_actor.mapper.scalar_range = mesh_actor.mapper.dataset.get_data_range(arr_var=next_scalar_name)
                mesh_actor.mapper.SetColorModeToMapScalars()
                mesh_actor.mapper.lookup_table.cmap = cmap
                mesh_actor.mapper.SetScalarVisibility(True)
                title = _set_title(self.uxds, next_scalar_name)
                if scalar_bar_actor is None:
                    scalar_bar_actor = plotter.add_scalar_bar(title=title, mapper=mesh_actor.mapper)
                else:
                    scalar_bar_actor.SetVisibility(True)
                    scalar_bar_actor.SetTitle(title)
            plotter.camera = camera_orig
            plotter.update()
            return

        def _help_message(help_message: str | None = None):
            help_message += "\nv: Isometric view"
            help_message += "\nUp/Down: Zoom in/out"
            help_message += "\n+/-: Increase/decrease point size"
            help_message += "\nw: Wireframe view"
            help_message += "\ns: Shaded view"
            help_message += "\nC: Enable cell picking"
            help_message += "\nh: Toggle this help message"
            help_message += "\nq: Quit"
            help_actor = pv.CornerAnnotation(0, help_message, name="help")
            help_actor.SetVisibility(False)
            return help_actor

        def reset_view(plotter):
            if self.is_global:
                plotter.reset_camera()
                plotter.view_yz()
            else:
                # Compute camera position so that it sits over the local region and the region fills the frame
                center_face_ind = self.surface.find_nearest_face(self.location)
                local_center = np.array(
                    [
                        self.surface.face_x[center_face_ind],
                        self.surface.face_y[center_face_ind],
                        self.surface.face_z[center_face_ind],
                    ]
                )
                # Add a 20% buffer to keep the region fully in view
                d = 1.2 * self.radius / np.tan(np.radians(plotter.camera.view_angle / 2))
                distance_multiplier = 1.0 + d / self.target.radius
                plotter.camera_position = local_center * distance_multiplier
                plotter.camera.focal_point = local_center
                plotter.camera.clipping_range = (0.35 * plotter.camera.distance, 2.0 * plotter.camera.distance)
            return

        if theme is not None:
            pv.set_plot_theme(theme)
        if transparent_background is not None:
            pv.global_theme.transparent_background = transparent_background

        new_plotter = plotter is None
        if new_plotter:
            plotter = pv.Plotter()
            plotter.enable_hidden_line_removal()

        if interval is None:
            uxds = self.uxds
        else:
            uxds = self.read_saved_output(interval=interval).isel(interval=0)
            interval = None
        mesh = self.to_vtk_mesh(uxds)

        if new_plotter:
            reset_view(plotter)

        face_variables = []
        component_variables = []
        for v in uxds.data_vars:
            if uxds[v].ndim > 0 and uxds[v].shape[0] == self.n_face:
                mesh.cell_data[v] = uxds[v].data
                face_variables.append(v)
                if uxds[v].ndim == 2:
                    component_variables.append(v)

        component = kwargs.pop("component", None)
        if variable_name is not None and variable_name in face_variables:
            title = _set_title(uxds, variable_name)
            if variable_name in component_variables:
                component = 0 if component is None else component
                title += f" (layer {component})"
        else:
            if variable is not None:
                variable = np.asarray(variable, dtype=np.float64)
                if variable.shape[0] != self.n_face:
                    raise ValueError("variable must be a string or an array with the same size as the number of faces in the grid")
                mesh.cell_data[variable_name] = variable
                title = variable_name
                if variable.ndim == 2:
                    component = 0 if component is None else component
                    title += f" (layer {component})"
            elif variable_name is not None:
                raise ValueError(
                    f"Variable '{variable_name}' not found in the surface data. The 'variable' argument must be provided with scalar values for the faces."
                )

        if variable_name is None:
            scalars = "face_elevation"  # Add a default scalar to enable the cmap to be applied properly if we decide to turn them back on later
        else:
            scalars = variable_name

        cmap = kwargs.pop("cmap", "cividis")
        add_mesh_kwargs = {k: v for k, v in kwargs.items() if k in PYVISTA_ADD_MESH_KWARGS}
        add_mesh_kwargs = {
            "name": self.target.name,
            "show_edges": False,
            "show_scalar_bar": False,
            "color": "grey",
            "pbr": True,
            **add_mesh_kwargs,
        }
        mesh_actor = plotter.add_mesh(mesh, scalars=scalars, component=component, cmap=cmap, **add_mesh_kwargs)
        if new_plotter and self.is_global:
            plotter.view_yz()

        if variable_name is None:
            mesh_actor.mapper.SetScalarVisibility(False)
        else:
            plotter.add_scalar_bar(title=title, mapper=mesh_actor.mapper)
        if enable_interactive and new_plotter:
            plotter = update_pyvista_help_message(plotter, new_message="j: Cycle through scalar face variables")
            plotter.add_key_event("j", lambda plotter=plotter, cmap=cmap: update_scalars(plotter, cmap=cmap))
            plotter.add_key_event("r", lambda plotter=plotter: reset_view(plotter))
        return plotter

    def export_region_polygon(self, driver: str = "GPKG", **kwargs: Any) -> None:
        """
        Export the local surface region as a polygon to a vector file.

        This will create a polygon that can be used for the `OpenCraterTool <https://github.com/thomasheyer/OpenCraterTool>`__ plugin in QGIS.

        Parameters
        ----------
        driver : str, optional
            The OGR driver to use for exporting the polygon. Default is "GPKG"
        **kwargs : Any
            |kwargs|
        """
        import geopandas as gpd
        import pandas as pd

        from cratermaker.components.crater import Crater
        from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP

        if driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP:
            file_extension = VECTOR_DRIVER_TO_EXTENSION_MAP[driver.upper()]
        else:
            raise ValueError("Cannot infer file extension from driver {driver}.")

        # If this is a local surface, we will use the Crater functionality to generate a polygon representation of it
        if self.is_local:
            region = Crater.maker(radius=self.region_radius, location=self.location)
            geoms = region.to_geoseries(surface=self.surface, split_antimeridian=False, use_measured_properties=False).to_crs(
                self.crs
            )
            df = pd.DataFrame(
                {"area": [self.area], "area_name": ["local_region"]},
                index=[0, 1],
            )
            gdf = gpd.GeoDataFrame(data=df, geometry=geoms, crs=self.crs)
            filename = self.export_dir / f"local_surface_AREA.{file_extension}"

            if file_extension == "shp":
                gdf.to_file(filename, driver=driver, **kwargs)
            else:
                gdf.to_file(filename, layer="local_surface_AREA", driver=driver, **kwargs)

        return

    def compute_distances(
        self,
        locations: ArrayLike,
        reference_location: PairOfFloats | None = None,
        validate: bool = True,
    ) -> NDArray[np.float64]:
        """
        Calculate the great circle distance between one point and one or more other points in meters.

        Parameters
        ----------
        locations : ArrayLike
            Longitude and latitude of the second point or array of points in degrees.
        reference_location : PairOfFloats | None, optional
            Longitude and latitude of the first point in degrees. If None, then the center of the local region will be used. Default is None.
        validate : bool, optional
            If the locations should be validated and normalized. Only use if the validation is causing a performance bottleneck.

        Returns
        -------
        NDArray
            Great circle distance between the two points in meters.

        Notes
        -----
        This is a wrapper for a compiled Rust function and is intended to be used as a helper to calculate_face_and_node_distances.
        """
        if reference_location is None:
            if self.is_local:
                reference_location = self.location
            else:
                raise ValueError("reference_location must be provided for global surfaces")

        if validate:
            locations = validate_and_normalize_location(locations)
            reference_location = validate_and_normalize_location(reference_location)
        locations = np.radians(locations)
        lon1, lat1 = np.radians(reference_location)
        if locations.ndim == 1 and locations.size == 2:
            locations = np.expand_dims(locations, axis=0)

        return surface_bindings.compute_distances(
            lon1=lon1, lat1=lat1, lon2=locations[:, 0], lat2=locations[:, 1], radius=self.surface.radius
        )

    def compute_bearings(self, locations: ArrayLike, reference_location: PairOfFloats | None = None) -> NDArray[np.float64]:
        """
        Calculate the initial bearing from one point to one or more other points in degrees.

        Parameters
        ----------
        locations : ArrayLike
            Longitude and latitude of the second point or array of points in degrees.
        reference_location : PairOfFloats | None, optional
            Longitude and latitude of the first point in degrees. If None, then the center of the local region will be used. Default is None.

        Returns
        -------
        NDArray
            Initial bearing from the first point to the second point or points in degrees.

        Notes
        -----
        This is intended to be used as a helper to calculate_face_and_node_bearings.
        """
        if reference_location is None:
            if self.is_local:
                reference_location = self.location
            else:
                raise ValueError("reference_location must be provided for global surfaces")
        locations = np.asarray(locations, dtype=np.float64)
        if locations.ndim == 1 and locations.size == 2:
            locations = locations.reshape((1, 2))

        lon1, lat1 = np.radians(validate_and_normalize_location(reference_location))
        lon2, lat2 = np.radians(locations[:, 0]), np.radians(locations[:, 1])
        bearing = surface_bindings.compute_bearings(lon1=lon1, lat1=lat1, lon2=lon2, lat2=lat2)
        return np.degrees(bearing)

    def compute_location_from_distance_bearing(
        self,
        distance: FloatLike | ArrayLike,
        bearing: FloatLike | ArrayLike,
        reference_location: PairOfFloats | None = None,
    ) -> NDArray[np.float64]:
        """
        Calculate the longitude and latitude of one or more points given a reference point, initial bearings, and distances.

        Parameters
        ----------
        bearing : FloatLike | ArrayLike
            One or more initial bearing values from the reference point to the target point or points in degrees.
        distance : FloatLike | ArrayLike
            One or more great circle distance from the reference point to the target point or points in meters.
        reference_location : PairOfFloats | None, optional
            One pair of longitude and latitude of the reference point in degrees. If None, then the center of the local region will be used. Default is None.

        Returns
        -------
        tuple or list of tuples the same length as bearing and distance.
            longitude and latitude as a tuple of floats in degrees.
        """
        if reference_location is None:
            if self.is_local:
                reference_location = self.location
            else:
                raise ValueError("reference_location must be provided for global surfaces")
        lon1, lat1 = np.radians(validate_and_normalize_location(reference_location))
        bearing = np.atleast_1d(np.radians(bearing))
        distance = np.atleast_1d(distance).astype(np.float64)
        if bearing.shape != distance.shape:
            raise ValueError("bearing and distance must have the same shape")
        if bearing.ndim > 1:
            raise ValueError("bearing and distance must have the same number of elements")
        lonlat2 = surface_bindings.compute_location_from_distance_bearing(
            lon1=lon1, lat1=lat1, distances=distance, bearings=bearing, radius=self.surface.radius
        )
        return validate_and_normalize_location(np.degrees(lonlat2))

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

    def saved_output_files(self, **kwargs: Any) -> list[Path]:
        """
        Check if the component has any output files in its output directory.

        Returns
        -------
        list[Path]
            A list of Path objects representing the files that would be removed during a reset operation. Returns an empty list if no files found
        """
        output_files = super().saved_output_files(**kwargs)
        # remove grid file from the list, as it is not removed during reset
        if self.grid_file in output_files:
            output_files.remove(self.grid_file)

        return output_files

    def get_location_extents(self) -> tuple[float, float, float, float]:
        """
        Computes the longitude and latitude extents of the local region.

        Returns
        -------
        lon_min : float
            Minimum longitude in degrees.
        lon_max : float
            Maximum longitude in degrees.
        lat_min : float
            Minimum latitude in degrees.
        lat_max : float
            Maximum latitude in degrees.

        """
        if self.is_global:
            raise ValueError("Cannot get the location extents of a global region.")

        return self.surface.get_location_extents(self.location, self.radius)

    @property
    def surface(self) -> Surface:
        """The Surface object that contains the mesh data for the global surface."""
        return self._surface

    @property
    def target(self):
        """The Target object associated with the LocalSurface."""
        return self.surface.target

    @property
    def n_edge(self) -> int:
        """The number of edges in the view."""
        if self._n_edge is None:
            if isinstance(self._edge_indices, slice):
                if np.isscalar(self._surface._uxds.uxgrid.n_edge):
                    self._n_edge = int(self._surface._uxds.uxgrid.n_edge)
                else:
                    self._n_edge = int(self._surface._uxds.uxgrid.n_edge[self._edge_indices].size)
            else:
                self._n_edge = self._edge_indices.size
        return self._n_edge

    @property
    def n_face(self) -> int:
        """The number of faces in the view."""
        if self._n_face is None:
            if isinstance(self.face_indices, slice):
                if np.isscalar(self.surface.face_elevation.size):
                    self._n_face = int(self.surface.face_elevation.size)
                else:
                    self._n_face = int(self.surface.face_elevation[self.face_indices].size)
            elif isinstance(self.face_indices, np.ndarray):
                self._n_face = int(self.face_indices.size)

        return self._n_face

    @property
    def n_node(self) -> int:
        """The number of nodes in the view."""
        if self._n_node is None:
            if isinstance(self.node_indices, slice):
                if np.isscalar(self.surface.node_elevation.size):
                    self._n_node = int(self.surface.node_elevation.size)
                else:
                    self._n_node = int(self.surface.node_elevation[self.node_indices].size)
            elif isinstance(self.node_indices, np.ndarray):
                self._n_node = int(self.node_indices.size)

        return self._n_node

    @property
    def n_nodes_per_face(self) -> NDArray:
        """The number of nodes per face in the view."""
        return self.surface.n_nodes_per_face[self.face_indices]

    @property
    def face_size(self) -> NDArray:
        """The effective pixel size of faces in the view in meters."""
        return self.surface.face_size[self.face_indices]

    @property
    def location(self) -> tuple[float, float]:
        """The (longitude, latitude) location of the center of the view in degrees."""
        return self._location

    @property
    def region_radius(self) -> FloatLike:
        """The radius of the region to include in the view in meters."""
        return self._region_radius

    @property
    def face_bearing(self) -> NDArray:
        """The initial bearing from the location to the faces measured clockwise relative to North in degrees."""
        if self._location is not None and self._face_bearing is None:
            self._face_bearing, self._node_bearing = self.calculate_face_and_node_bearings()
        return self._face_bearing

    @property
    def node_bearing(self) -> NDArray:
        """The initial bearing from the location to the nodes measured clockwise relative to North in degrees."""
        if self._location is not None and self._node_bearing is None:
            self._face_bearing, self._node_bearing = self.calculate_face_and_node_bearings()
        return self._node_bearing

    @property
    def face_distance(self) -> NDArray:
        """The distance from the location to the faces in meters."""
        if self._location is not None and self._face_distance is None:
            self._face_distance, self._node_distance = self.calculate_face_and_node_distances()
        return self._face_distance

    @property
    def node_distance(self) -> NDArray:
        """The distance from the location to the nodes in meters."""
        if self._location is not None and self._node_distance is None:
            self._face_distance, self._node_distance = self.calculate_face_and_node_distances()
        return self._node_distance

    @property
    def face_area(self) -> NDArray:
        """The areas of the faces in the view in m²."""
        return self.surface.face_area[self.face_indices]

    @property
    def face_lat(self) -> NDArray:
        """Latitude of the center of the faces in the view in degrees."""
        return self.surface.face_lat[self.face_indices]

    @property
    def face_lon(self) -> NDArray:
        """Longitude of the center of the faces in the viewin degrees."""
        return self.surface.face_lon[self.face_indices]

    @property
    def face_x(self) -> NDArray:
        """Cartesian x location of the center of the faces in the view in meters."""
        return self.surface.face_x[self.face_indices]

    @property
    def face_y(self) -> NDArray:
        """Cartesian y location of the center of the faces in the view in meters."""
        return self.surface.face_y[self.face_indices]

    @property
    def face_z(self) -> NDArray:
        """Cartesian z location of the center of the faces in the view in meters."""
        return self.surface.face_z[self.face_indices]

    @property
    def edge_face_distance(self) -> NDArray:
        """
        Distances between the edges and the faces in the view in meters.

        Dimensions: `(n_edge)`
        """
        return self.surface.edge_face_distance[self.edge_indices]

    @property
    def edge_length(self) -> NDArray:
        """
        Lengths of the edges in the view in meters.

        Dimensions: `(n_edge)`
        """
        return self.surface.edge_length[self.edge_indices]

    @property
    def node_lat(self) -> NDArray:
        """Latitude of the nodes in the view in degrees."""
        return self.surface.node_lat[self.node_indices]

    @property
    def node_lon(self) -> NDArray:
        """Longitude of the nodes in the view in degrees."""
        return self.surface.node_lon[self.node_indices]

    @property
    def node_x(self) -> NDArray:
        """Cartesian x location of the nodes in the view in meters."""
        return self.surface.node_x[self.node_indices]

    @property
    def node_y(self) -> NDArray:
        """Cartesian y location of the nodes in the view in meters."""
        return self.surface.node_y[self.node_indices]

    @property
    def node_z(self) -> NDArray:
        """Cartesian z location of the nodes in the view in meters."""
        return self.surface.node_z[self.node_indices]

    @property
    def area(self) -> float:
        """The total area of the faces in the view in m²."""
        if self._area is None:
            self._area = float(self.face_area.sum())
        return self._area

    @property
    def face_indices(self) -> NDArray:
        """The indices of the faces in the view."""
        if self._face_indices is None:
            raise ValueError("face_indices must be set to use this object.")
        return self._face_indices

    @property
    def node_indices(self) -> NDArray:
        """The indices of the nodes in the view."""
        if self._node_indices is None:
            self._node_indices = np.unique(self.surface.face_node_connectivity[self.face_indices].ravel())
            self._node_indices = self._node_indices[self._node_indices != INT_FILL_VALUE]

        return self._node_indices

    @property
    def edge_indices(self) -> NDArray:
        """The indices of the edges in the view."""
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

    @property
    def crs(self) -> CRS:
        """
        Return a local Azimuthal Equidistant (LAEA) CRS centered on `self.location` (in meters).

        If `self.location` is not set, fall back to the parent `Surface` geographic CRS.
        The LAEA CRS uses the target body's spherical radius.
        """
        if self._crs is None:
            self._crs = self.surface.get_crs(
                radius=self.surface.radius,
                name=self.surface.target.name,
                location=self.location,
            )
        return self._crs

    def read_saved_output(self, interval: int | None = None, reset: bool = False, **kwargs: Any) -> uxr.UxDataset:
        """
        Read the saved local surface data from disk for a given interval number and return it as a UxDataset.

        Parameters
        ----------
        interval : int | None, optional
            Interval number of data file to read. Default is None (all intervals are read).
        reset : bool, optional
            Flag to indicate whether to reset the surface. If True it reads in the grid but creates an empty dataset and saves it to interval 0. Default is False.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        UxDataset
            An initialized UxDataset object containing the grid and data.
        Dataset
            The xarray Dataset containing only the local surface data.
        """
        if self.is_local:
            uxds_global = self.surface.read_saved_output(interval=interval, **kwargs)
            if not reset:
                return uxr.UxDataset(uxds_global.sel(n_face=self.face_indices, n_node=self.node_indices), uxgrid=self.uxgrid)
            else:
                return uxr.UxDataset(uxgrid=self.uxgrid)

        if reset:
            interval = 0
            filename = self.output_dir / self.output_filename(interval)
            filename.unlink(missing_ok=True)
        files = super().read_saved_output(interval=interval, **kwargs)
        if isinstance(files, tuple):
            ds, grid = files
        else:
            grid = files
            ds = None
        if grid is None:
            raise ValueError("No grid file found.")
        uxgrid = uxr.Grid.from_dataset(grid)

        if ds is None:
            reset = True
        if reset:
            uxds = uxr.UxDataset(uxgrid=uxgrid)
        else:
            uxds = uxr.UxDataset.from_xarray(ds, uxgrid=uxgrid)
            for v in uxds.data_vars:
                uxds[v].attrs = ds[v].attrs.copy()
                # Ensure that the tags are the correct data type
                if "layer" in ds[v].dims:
                    uxds[v] = uxds[v].astype(np.uint32)
        return uxds

    @property
    def uxgrid(self) -> uxr.Grid:
        """The Uxarrray Grid representation of the local surface."""
        if self._uxgrid is None:
            if self.is_global:
                self._uxgrid = self.surface.uxgrid
            else:
                grid_ds = xr.Dataset(
                    data_vars={
                        "face_x": (("n_face",), self.face_x),
                        "face_y": (("n_face",), self.face_y),
                        "face_z": (("n_face",), self.face_z),
                        "face_lat": (("n_face",), self.face_lat),
                        "face_lon": (("n_face",), self.face_lon),
                        "node_x": (("n_node",), self.node_x),
                        "node_y": (("n_node",), self.node_y),
                        "node_z": (("n_node",), self.node_z),
                        "node_lat": (("n_node",), self.node_lat),
                        "node_lon": (("n_node",), self.node_lon),
                        "edge_face_distance": (("n_edge",), self.edge_face_distance),
                        "edge_length": (("n_edge",), self.edge_length),
                        "edge_face_connectivity": (
                            ("n_edge", "two"),
                            self.edge_face_connectivity,
                        ),
                        "face_edge_connectivity": (
                            ("n_face", "n_max_face_edges"),
                            self.face_edge_connectivity,
                        ),
                        "face_node_connectivity": (
                            ("n_face", "n_max_face_nodes"),
                            self.face_node_connectivity,
                        ),
                        "face_face_connectivity": (
                            ("n_face", "n_max_face_faces"),
                            self.face_face_connectivity,
                        ),
                        "node_face_connectivity": (
                            ("n_node", "n_max_node_faces"),
                            self.node_face_connectivity,
                        ),
                        "edge_node_connectivity": (
                            ("n_edge", "two"),
                            self.edge_node_connectivity,
                        ),
                    },
                )
                grid_ds["grid_topology"] = self.surface.uxgrid._ds.grid_topology
                self._uxgrid = uxr.Grid.from_dataset(grid_ds)
        return self._uxgrid

    @property
    def uxds(self) -> UxDataset:
        """The UxDataset representation of the local surface."""
        if self.is_global:
            return self.surface.uxds
        return uxr.UxDataset(self.surface.uxds.sel(n_face=self.face_indices, n_node=self.node_indices), uxgrid=self.uxgrid)

    @property
    def grid_file(self):
        """Path to the grid file."""
        return self.output_dir / f"{self._grid_file_prefix}.{self._output_file_extension}"

    @property
    def pix(self) -> float:
        """The effective pixel size of the base surface in meters."""
        return self.surface.pix

    @property
    def radius(self) -> float:
        """The radius of the local region in meters."""
        return self.region_radius

    @property
    def plot_dir(self) -> Path:
        """The directory to save plots to."""
        return self.surface.plot_dir

    @property
    def from_surface(self) -> Transformer:
        """A pyproj Transformer object to convert from the surface CRS to the local CRS."""
        if self._from_surface is None:
            self._from_surface = Transformer.from_crs(
                self.surface.crs,
                self.crs,
                always_xy=True,
            )
        return self._from_surface

    @property
    def to_surface(self) -> Transformer:
        """A pyproj Transformer object to convert from the local CRS to the surface CRS."""
        if self._to_surface is None:
            self._to_surface = Transformer.from_crs(
                self.crs,
                self.surface.crs,
                always_xy=True,
            )
        return self._to_surface

    def set_face_proj(self):
        """
        Set the projected x and y coordinates of the faces relative to the LocalSurface center.
        """
        if self.face_distance is not None and self.face_bearing is not None and self._face_proj_x is None:
            self._face_proj_x = self.face_distance * np.sin(np.radians(self.face_bearing))
            self._face_proj_y = self.face_distance * np.cos(np.radians(self.face_bearing))
        return

    def set_node_proj(self):
        """
        Set the projected x and y coordinates of the nodes relative to the LocalSurface center.
        """
        if self.node_distance is not None and self.node_bearing is not None and self._node_proj_x is None:
            self._node_proj_x = self.node_distance * np.sin(np.radians(self.node_bearing))
            self._node_proj_y = self.node_distance * np.cos(np.radians(self.node_bearing))
        return

    @property
    def face_proj_x(self) -> NDArray:
        """The projected x coordinates of the faces relative to the LocalSurface center in meters."""
        return self._face_proj_x

    @property
    def face_proj_y(self) -> NDArray:
        """The projected y coordinates of the faces relative to the LocalSurface center in meters."""
        return self._face_proj_y

    @property
    def node_proj_x(self) -> NDArray:
        """The projected x coordinates of the nodes relative to the LocalSurface center in meters."""
        return self._node_proj_x

    @property
    def node_proj_y(self) -> NDArray:
        """The projected y coordinates of the nodes relative to the LocalSurface center in meters."""
        return self._node_proj_y

    @property
    def desloped_face_elevation(self) -> NDArray:
        """The face elevations with the mean slope of the region removed in meters."""
        return self._desloped_face_elevation

    def compute_desloped_face_elevation(self):
        """Compute the face elevations with the mean slope of the region removed."""
        if self.is_local and self._desloped_face_elevation is None:
            reference_elevation = self.get_reference_surface(only_faces=True)
            self._desloped_face_elevation = self.uxds.face_elevation.data - reference_elevation
        return

    @property
    def is_local(self) -> bool:
        """Whether this surface is a local region or the full global surface."""
        return self.location is not None

    @property
    def is_global(self) -> bool:
        """Whether this surface is the full global surface."""
        return self.location is None

    def data_composer(self) -> DataComposer:
        """
        Creates a ``DataComposer`` tied to this ``LocalSurface``. This should usually be called in a ``with`` statement.

        Returns
        -------
        DataComposer
        """
        return DataComposer(self)

    @property
    def face_variables(self) -> list[str]:
        """
        Returns a list of the data variables currently being stored on the faces.
        """
        return [v for v in self.uxds.variables if "n_face" in self.uxds[v].dims]

    @property
    def node_variables(self) -> list[str]:
        """
        Returns a list of the data variables currently being stored on the nodes.
        """
        return [v for v in self.uxds.variables if "n_node" in self.uxds[v].dims]


class DataComposer(AbstractContextManager):
    """
    Created using ``Surface.data_composer`` or ``LocalSurface.data_composer``. Contains some static methods for loading LOLA data.
    """

    def __init__(self, surface: LocalSurface):
        """
        **Warning:** This object should not be instantiated directly. Instead, use ``Surface.data_composer()`` or ``LocalSurface.data_composer()``.
        """
        self._localsurface = surface
        self._data_list: list[DatasetReader] = []
        self._finished = False

    def add_data(
        self,
        data: DatasetReader | str | list[DatasetReader | str],
    ):
        """
        Adds data to eventually be applied to the surface.

        The data isn't applied until ``finish`` is called or, if appplicable, ``self`` exits the ``with`` context.

        Parameters
        ----------
        dataset : DatasetReader | str | list[DatasetReader | str]
            A single dataset or a list of multiple datasets to be merged.
            Calls ``rasterio.open`` when necessary.
        """
        if self._finished:
            raise ValueError(f"{type(self).__name__} is already finished or cancelled.")

        if not isinstance(data, list):
            data = [data]

        for i, src in enumerate(data):
            if not isinstance(src, DatasetReader):
                data[i] = rasterio.open(src)

            self._data_list.append(data[i])

    def cancel(self):
        """
        Cancels the execution and closes all open datasets.
        """
        for dataset in self._data_list:
            dataset.close()
        self._data_list = []
        self._finished = True

    def finish(self):
        """
        Apply all the datasets to the mesh and close them. This is implicitly called when exiting a with context.
        """
        if self._finished:
            raise ValueError(f"{type(self).__name__} is already finished or cancelled.")

        def _orderable_distance(lon1, lat1, lon2, lat2):
            _nnon1 = np.deg2rad(lon1)
            lat1 = np.deg2rad(lat1)
            lon2 = np.deg2rad(lon2)
            lat2 = np.deg2rad(lat2)

            dlon = (lon2 - lon1 + np.pi) % (2 * np.pi)
            dlat = lat2 - lat1
            return np.pow(np.sin(dlat / 2), 2) + np.cos(lat1) * np.cos(lat2) * np.pow(np.sin(dlon / 2), 2)

        def _read(dataset: DatasetReader, window: Window | None = None) -> tuple[NDArray, Window]:
            block_windows = []
            for _, block in dataset.block_windows(1):
                if window is None or windows.intersect(window, block):
                    block_windows.append(block)

            res_window: Window = windows.union(block_windows)

            print(f"    Reading {dataset.name} ({res_window.width}x{res_window.height}px of {dataset.width}x{dataset.height}px)")

            res = np.empty((res_window.height, res_window.width))

            for block in tqdm(block_windows):
                out = res[
                    block.row_off - res_window.row_off : block.row_off - res_window.row_off + block.height,
                    block.col_off - res_window.col_off : block.col_off - res_window.col_off + block.width,
                ]
                dataset.read(1, window=block, out=out)

            return res, res_window

        print(f"Applying {len(self._data_list)} dataset{'' if len(self._data_list) == 1 else 's'} to the mesh")

        lons = np.concatenate((self._localsurface.face_lon, self._localsurface.node_lon))
        lats = np.concatenate((self._localsurface.face_lat, self._localsurface.node_lat))

        read_data: list[NDArray] = []
        data_windows: list[Window] = []

        pix_dist = np.full((len(self._data_list), len(lons)), np.inf)
        for n, dataset in enumerate(self._data_list):
            to_data = Transformer.from_crs(self._localsurface.surface.crs, dataset.crs)
            from_data = Transformer.from_crs(dataset.crs, self._localsurface.surface.crs)

            x_coords, y_coords = to_data.transform(lons, lats)

            mask1 = np.isfinite(x_coords) & np.isfinite(y_coords)
            rows, cols = rowcol(dataset.transform, x_coords[mask1], y_coords[mask1])
            mask2 = (cols >= 0) & (cols < dataset.width) & (rows >= 0) & (rows < dataset.height)
            mask1[mask1] = mask2
            rows = rows[mask2]
            cols = cols[mask2]

            if self._localsurface.is_local:
                cmin, rmin, cmax, rmax = np.min(cols), np.min(rows), np.max(cols), np.max(rows)
                data, window = _read(dataset, window=Window(cmin, rmin, cmax - cmin + 1, rmax - rmin + 1))
            else:
                data, window = _read(dataset)

            read_data.append(data)
            data_windows.append(window)

            values = data[rows - window.row_off, cols - window.col_off]

            print("        Calculating pixel distances")

            mask3 = np.isfinite(values) & (values != dataset.nodata)
            mask1[mask1] = mask3
            x_pix, y_pix = from_data.transform(*dataset.xy(rows[mask3], cols[mask3]))
            pix_dist[n, mask1] = _orderable_distance(lons[mask1], lats[mask1], x_pix, y_pix)

        idx = np.argmin(pix_dist, axis=0)
        global_mask = np.isfinite(pix_dist[idx, np.arange(len(idx))])
        elevation = np.zeros_like(idx, dtype=np.float32)

        print("    Setting elevation data")
        for n, (dataset, data, window) in enumerate(zip(self._data_list, read_data, data_windows, strict=True)):
            to_data = Transformer.from_crs(self._localsurface.surface.crs, dataset.crs)
            mask = global_mask & (idx == n)

            if "KILOMETER" in dataset.units:
                scale_factor = 1000.0
            else:
                scale_factor = 1.0

            x_coords, y_coords = to_data.transform(lons[mask], lats[mask])
            rows, cols = rowcol(dataset.transform, x_coords, y_coords)
            elevation[mask] = data[rows - window.row_off, cols - window.col_off] * scale_factor
            dataset.close()

        self._localsurface.update_elevation(elevation)
        self._data_list = []
        self._finished = True
        print("    Done")

    def __enter__(self) -> DataComposer:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.finished:
            return
        if exc_type is None:
            self.finish()
        else:
            self.cancel()

    @staticmethod
    def _get_lola_cylindrical_url_from_pds(
        pds_file_resolution: int, location: PairOfFloats, boundary_offset: tuple[int, int] = (0, 0)
    ) -> str:
        """
        Retrieve the appropriate LOLA DEM file url for a given location and resolution from the PDS.

        Parameters
        ----------
        pds_file_resolution : int
            DEM resolution in meters per pixel.
        location : PairOfFloats,
            The longitude and latitude of the location in degrees.
        boundary_offset : tuple[int, int], optional
            Offset to apply to tile index to access neighboring tiles. Default is (0, 0). This is used when the local region crosses one or more tile boundaries.

        Returns
        -------
        url : str
            Link to the DEM file

        Notes
        -----
        This is meant to be used for mid-latitude regions on the Moon between -60 and 60 degrees latitude.
        For resolutions of 128, 256, and 512 pix/deg, this will return the SLDEM2015 datasets.

        """
        lola_cylindrical_url = (
            "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/float_img/"
        )
        sldem_global_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/sldem2015/global/float_img/"
        sldem_tile_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/sldem2015/tiles/float_img/"

        def _get_pds_file_suffix(pds_file_resolution, location, boundary_offset):
            if location[0] < 0:
                location = (location[0] + 360, location[1])
            if pds_file_resolution == 512:
                dlon = 45
                dlat = 30
            elif pds_file_resolution == 256:
                dlon = 120
                dlat = 60

            latdir = "n" if location[1] + boundary_offset[1] * dlat > 0 else "s"
            latval = np.abs(location[1] / dlat) + boundary_offset[1]
            if latdir == "n":
                latlo = int(np.ceil(latval - 1) * dlat)
                lathi = int(latlo + dlat)
            else:
                lathi = int(np.floor(latval + 1) * dlat)
                latlo = int(lathi - dlat)

            lonval = np.abs(location[0] / dlon) + boundary_offset[0]

            # Pole crossing, so flip the longitudes by 180 degrees
            if lathi > 90 or latlo > 90:
                lathi = 90
                latlo = lathi - dlat
                lonval += 180 / dlon

            lonlo = int(np.floor(lonval) * dlon) % 360
            lonhi = int(lonlo + dlon)

            # Check to make sure the lo and hi values are in the right direction
            if (latdir.upper() == "N" and latlo > lathi) or (latdir.upper() == "S" and latlo < lathi):
                tmp = latlo
                latlo = lathi
                lathi = tmp

            if pds_file_resolution == 512:
                return f"{latlo:02d}{latdir.lower()}_{lathi:02d}{latdir.lower()}_{lonlo:03d}_{lonhi:03d}"
            elif pds_file_resolution == 256:
                return f"{latlo:d}{latdir.lower()}_{lathi:d}{latdir.lower()}_{lonlo:03d}_{lonhi:03d}"

        if pds_file_resolution >= 256:
            url = f"{sldem_tile_url}sldem2015_{pds_file_resolution:d}_{_get_pds_file_suffix(pds_file_resolution, location, boundary_offset)}_float.xml"
        elif pds_file_resolution == 128:
            url = f"{sldem_global_url}sldem2015_{pds_file_resolution:d}_60s_60n_000_360_float.xml"
        else:
            url = f"{lola_cylindrical_url}ldem_{pds_file_resolution:d}_float.xml"

        return url

    @staticmethod
    def get_lola_cylindrical_files_from_pds(
        resolution: FloatLike, lat_range: PairOfFloats, lon_range: PairOfFloats
    ) -> tuple[list[str], int]:
        """
        Retrieve the appropriate cylindrically projected LOLA DEM file url or list of urls for a given location and resolution from the PDS.

        This will determine which, if any, boundaries are crossed and will return a list of files to cover a region.

        Parameters
        ----------
        resolution : FloatLike
            Requested resolution in degrees per pixel. The closest available resolution will be used.
        lat_range : tuple of float, optional
            The (min_lat, max_lat) in degrees of the local region.
        lon_range : tuple of float, optional
            The (min_lon, max_lon) in degrees of the local region.
        """
        AVAILABLE_RESOLUTIONS = [4, 16, 64, 128, 256, 512]  #  pix / deg
        diffs = [abs(resolution - res) for res in AVAILABLE_RESOLUTIONS]
        pds_file_resolution = AVAILABLE_RESOLUTIONS[np.argmin(diffs)]

        lat_min, lat_max = lat_range
        lon_min, lon_max = lon_range
        center = ((lon_min + lon_max) / 2.0, (lat_min + lat_max) / 2.0)

        # First, retrive the file for the centerpoint:
        filelist = [DataComposer._get_lola_cylindrical_url_from_pds(pds_file_resolution, center)]
        if pds_file_resolution < 256:
            return (
                filelist,
                pds_file_resolution,
            )  # These files cover the entire globe, no need to determine if boundaries are crossed

        combo = [(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_min), (lon_max, lat_max)]

        for loc in combo:
            f = DataComposer._get_lola_cylindrical_url_from_pds(pds_file_resolution, loc)
            if f not in filelist:
                filelist.append(f)

        return filelist, pds_file_resolution

    @staticmethod
    def get_lola_polar_files_from_pds(resolution: FloatLike, lat_range: PairOfFloats) -> tuple[list[str], int]:
        """
        Retrieve the appropriate polar projected LOLA DEM file url for a given location and resolution from the PDS.

        Parameters
        ----------
        resolution : FloatLike
            Requested resolution in meters per pixel. The closest available resolution will be used.
        lat_range : tuple of float
            The (min_lat, max_lat) in degrees of the local region.

        """
        import rasterio

        src_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/polar/float_img/"

        AVAILABLE_RESOLUTIONS = [400, 240, 200, 120, 100, 80, 60, 40, 30, 20, 10, 5]
        MIN_LAT = [45, 60, 45, 60, 45, 80, 60, 80, 75, 80, 85, 87.5]

        if lat_range[0] > 0:
            pole = "n"
            min_lat_index = 0  # in the northern hemisphere, the minimum latitude defines coverage
        else:
            pole = "s"
            min_lat_index = 1  # in the southern hemisphere, the maximum latitude defines coverage
        combo = [
            (res, minlat)
            for res, minlat in zip(AVAILABLE_RESOLUTIONS, MIN_LAT, strict=True)
            if abs(lat_range[min_lat_index]) >= minlat
        ]
        if len(combo) == 0:
            raise ValueError("No available LOLA polar DEM files cover the requested latitude range.")

        valid_resolutions, valid_min_lat = zip(*combo, strict=True)

        diffs = [abs(resolution - res) for res in valid_resolutions]
        pds_file_resolution = valid_resolutions[np.argmin(diffs)]
        pds_lat_min = valid_min_lat[np.argmin(diffs)]
        filename = f"ldem_{pds_lat_min:d}{pole}_{pds_file_resolution:d}m"
        url = f"{src_url}{filename}_float.xml"

        return [url], pds_file_resolution

    @staticmethod
    def get_lola_dem_file_list(
        pix: FloatLike,
        lat_range: PairOfFloats,
        lon_range: PairOfFloats,
    ) -> tuple[list[str], int]:
        """
        Determine the set of LOLA DEM files needed to cover the local region.

        Parameters
        ----------
        pix : FloatLike
            The approximate resolution in meters per pixel to use to select the DEM files.
        lat_range : PairOfFloats, optional
            The (min_lat, max_lat) in degrees of the local region.
        lon_range : PairOfFloats, optional
            The (min_lon, max_lon) in degrees of the local region.
        """
        target_pds_resolution = np.pi / 180.0 * 1737.53e3 / pix  # The moon's radius
        if target_pds_resolution > 10 and (
            lat_range[0] > 60 or lat_range[1] < -60
        ):  # Use polar files high latitude, high resolution regions.
            return DataComposer.get_lola_polar_files_from_pds(pix, lat_range=lat_range)
        else:  # Use cylindrical for all other cases
            return DataComposer.get_lola_cylindrical_files_from_pds(
                resolution=target_pds_resolution, lat_range=lat_range, lon_range=lon_range
            )

    def populate_with_lola_data(self, pix: FloatLike | None = None):
        """
        Loads the DEM files from LOLA needed to cover the surface.

        Parameters
        ----------
        pix : FloatLike
            The approximate resolution in meters per pixel to use to select the DEM files.
            Defaults to the resolution of the surface used.
        """
        if pix is None:
            pix = self._localsurface.pix

        if self._localsurface.is_local:
            lon_min, lon_max, lat_min, lat_max = self._localsurface.get_location_extents()
            self.add_data(DataComposer.get_lola_dem_file_list(pix, lat_range=(lat_min, lat_max), lon_range=(lon_min, lon_max))[0])
        else:
            self.add_data(DataComposer.get_lola_dem_file_list(pix, lat_range=(-60, 60), lon_range=(-180, 180))[0])
            self.add_data(DataComposer.get_lola_polar_files_from_pds(pix, lat_range=(60, -60))[0])
            self.add_data(DataComposer.get_lola_polar_files_from_pds(pix, lat_range=(-60, 60))[0])

    @property
    def localsurface(self) -> LocalSurface:
        """
        The ``LocalSurface`` used to create the ``DataComposer``.
        """
        return self._localsurface

    @property
    def surface(self) -> Surface:
        """
        The ``Surface`` used to create the ``DataComposer``.
        """
        return self._localsurface.surface

    @property
    def finished(self):
        """
        Whether ``finish()`` or ``cancel()`` has been called or not.
        """
        return self._finished


import_components(__name__, __path__)
