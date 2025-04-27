from __future__ import annotations
import xarray as xr
from xarray import DataArray, Dataset
import uxarray as uxr
from uxarray import INT_FILL_VALUE, Grid, UxDataArray, UxDataset
import os
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import shutil
from pathlib import Path
import tempfile
from typing import List, Union
from numpy.typing import NDArray, ArrayLike
from numpy.random import Generator
from .target import Target
import warnings
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils.montecarlo import get_random_location_on_face
from cratermaker.utils.general_utils import parameter
from cratermaker.components.grid import Grid
from cratermaker.core.base import _to_config
from cratermaker.constants import _DATA_DIR, _GRID_FILE_NAME, _SMALLFAC, _COMBINED_DATA_FILE_NAME
from cratermaker._cratermaker import surface_functions
from .base import CratermakerBase, _rng_init, _simdir_init, CommonArgs
from typing import Any

class Surface:
    """
    This class is used for handling surface-related data and operations in the cratermaker project. It provides methods for 
    setting elevation data, calculating distances and bearings, and other surface-related computations.
    
    The Surface class extends UxDataset for the cratermaker project.
    
    Parameters
    ----------
    uxds : UxDataset
        The UxDataset object that contains the surface data.
    grid : str, optional
        The name of the grid used for the surface. Default is "icosphere".
    target : Target, optional
        The target body or name of a known target body for the impact simulation. 
    compute_face_areas : bool, optional
        Flag to indicate whether to compute face areas. Default is False.    
    simdir : str | Path
        The main project simulation directory. Defaults to the current working directory if None.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """ 

    def __init__(self, 
                 uxds: UxDataset, 
                 grid: Grid | str | None = None,
                 target: Target | str | None = None,
                 compute_face_areas: bool = False,
                 simdir: os.PathLike = Path.cwd(),
                 rng: Generator | None = None, 
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs):

        argproc = CratermakerBase(simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state)
        self._rng = argproc.rng
        self._rng_seed = argproc.rng_seed
        self._rng_state = argproc.rng_state
        self._simdir = argproc.simdir

        # Additional initialization for Surface
        self.grid = grid
        self._target = Target.make(target, **kwargs)
        self._compute_face_areas = compute_face_areas
        self._name = "Surface"
        self._description = "Surface class for cratermaker"
        self._area = None
        self._smallest_length = None
        self._node_tree = None
        self._face_tree = None
        self._uxds = uxds

        return
    
    def load_from_data(self, compute_face_areas):
        if compute_face_areas: 
            # Compute face area needed in the non-normalized units for future calculations
            self.face_areas = self.uxds.uxgrid.face_areas.values * self.target.radius**2
            self.smallest_length = np.sqrt(self.face_areas.min()) * _SMALLFAC   
        
        self.face_lat = self.uxds.uxgrid.face_lat.values
        self.face_lon = self.uxds.uxgrid.face_lon.values
        self.node_lat = self.uxds.uxgrid.node_lat.values
        self.node_lon = self.uxds.uxgrid.node_lon.values
        self.node_elevation = self.uxds["node_elevation"].values
        self.face_elevation = self.uxds["face_elevation"].values
        self.ejecta_thickness = self.uxds["ejecta_thickness"].values
        self.ray_intensity = self.uxds["ray_intensity"].values
    
    def save_to_data(self):
        self.uxds["node_elevation"].values = self.node_elevation
        self.uxds["face_elevation"].values = self.face_elevation
        self.uxds["ejecta_thickness"].values = self.ejecta_thickness
        self.uxds["ray_intensity"].values = self.ray_intensity
    
    def full_view(self):
        return SurfaceView(self, slice(None), slice(None))

    def to_config(self, remove_common_args: bool = False, **kwargs: Any) -> dict[str, Any]:
        """
        Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file. 

        Parameters
        ----------
        obj : Any
            The object whose attributes will be stored.  It must have a _user_defined attribute.
        remove_common_args : bool, optional
            If True, remove the set of common arguments that are shared among all components of the project from the configuration. Defaults to False.
        **kwargs : Any
            Additional keyword arguments for subclasses.

        Returns
        -------
        dict[str, Any]
            A dictionary of the object's attributes that can be serialized to YAML.
        Notes
        -----
        - The function will ignore any attributes that are not serializable to human-readable YAML. Therefore, it will ignore anything that cannot be converted into a str, int, float, or bool.
        - The function will convert Numpy types to their native Python types.
        """
        return _to_config(self, remove_common_args=remove_common_args, **kwargs)

    @classmethod
    def make(cls: Surface, 
             grid: str | None = None,
             target: Target | None = None, 
             reset_surface: bool = True, 
             regrid: bool = False,
             simdir: str | None = None,
             rng: Generator | None = None,
             rng_seed: int | None = None,
             rng_state: dict | None = None,
             **kwargs) -> Surface:
        """
        Factory method to create a Surface instance from a grid file.

        Parameters
        ----------
        grid : str, optional
            The name of the grid used for the surface. Default is "icosphere".
        target : Target, optional
            The target body or name of a known target body for the impact simulation. 
        reset_surface : bool, optional
            Flag to indicate whether to reset the surface. Default is True.
        regrid : bool, optional
            Flag to indicate whether to regrid the surface. Default is False.
        simdir : str | Path
            The main project simulation directory. Defaults to the current working directory if None.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
        Additional keyword arguments.

        Returns
        -------
        Surface
            An initialized Surface object.
        """
        argproc = CratermakerBase(simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state)
        rng = argproc.rng
        rng_seed = argproc.rng_seed
        rng_state = argproc.rng_state
        simdir = argproc.simdir

        target = Target.make(target, **kwargs)

        kwargs = {**kwargs, **vars(argproc.common_args)}
        grid = Grid.make(grid=grid, target=target, regrid=regrid, **kwargs) 

        # Get the names of all data files in the data directory that are not the grid file
        data_dir = grid.file.parent
        data_file_list = list(data_dir.glob("*.nc"))
        if grid.file in data_file_list:
            data_file_list.remove(grid.file)
            
        # Generate a new surface if either it is explicitly requested via parameter or a data file doesn't yet exist 
        reset_surface = reset_surface or not data_file_list or grid.regrid
        
        # If reset_surface is True, delete all data files except the grid file 
        if reset_surface or grid.regrid:
            for f in data_file_list:
                f.unlink()  
            data_file_list = []
        
        # Initialize UxDataset with the loaded data
        try:
            if data_file_list:
                surf = uxr.open_mfdataset(grid.file, data_file_list, use_dual=False).isel(time=-1)
                surf.uxgrid = uxr.open_grid(grid.file, use_dual=False)
            else:
                surf = uxr.UxDataset()
                surf.uxgrid = uxr.open_grid(grid.file, use_dual=False)
        except:
            raise ValueError("Error loading grid and data files")
        
        surf = cls(surf,
                   target = target,
                   simdir = simdir,
                   rng = rng,
                   rng_seed = rng_seed,
                   rng_state = rng_state,
                   compute_face_areas = True,
                   grid = grid) 
        
        if reset_surface:
            surf.generate_data(data=0.0,
                               name="ejecta_thickness",
                               long_name="ejecta thickness",
                               units= "m",
                               save_to_file=True)     
            surf.generate_data(data=0.0,
                               name="ray_intensity",
                               long_name="ray intensity value",
                               units= "",
                               save_to_file=True)                         
            surf.set_elevation(0.0,save_to_file=True)

        surf.grid_config = grid.to_config(remove_common_args=True)
        surf.grid_config.pop("radius", None) # Radius is determined by the target when the grid is associated with a Surface, so this is redundant 
        surf.load_from_data(compute_face_areas=True)
        
        return surf
    
    def _calculate_binary_op(self, *args, **kwargs):
        """Override to make the result a complete instance of ``cratermaker.Surface``."""
        ds = super()._calculate_binary_op(*args, **kwargs)

        if isinstance(ds, Surface):
            ds._grid = self._grid
            ds._description = self._description
            ds._simdir = self._simdir
            ds._target = self._target
            ds._rng = self._rng
            ds._rng_seed = self._rng_seed
            ds._rng_state = self._rng_state
            ds._smallest_length = self._smallest_length
            ds._compute_face_areas = False
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         simdir=self.simdir,
                         rng=self.rng,
                         rng_seed=self.rng_seed,
                         rng_state=self.rng_state,
                         compute_face_areas=False,
                         )
        return ds
    
    @classmethod
    def _construct_direct(cls, *args, **kwargs):
        """Override to make the result a ``cratermaker.Surface`` class."""

        return cls(uxr.UxDataset._construct_direct(*args, **kwargs))    

    def _copy(self, **kwargs):
        """Override to make the result a complete instance of ``cratermaker.Surface``."""
        copied = super()._copy(**kwargs)

        copied._grid = self._grid
        copied._description = self._description
        copied._simdir = self._simdir
        copied._smallest_length = self._smallest_length
        copied._area = self._area
        copied._target = self._target
        copied._rng = self._rng
        copied._rng_seed = self._rng_seed
        copied._rng_state = self._rng_state
        copied._compute_face_areas = False
        
        return copied    
  
    def _replace(self, *args, **kwargs):
        """Override to make the result a complete instance of ``cratermaker.Surface``."""
        ds = super()._replace(*args, **kwargs)

        if isinstance(ds, Surface):
            ds._grid = self._grid
            ds._description = self._description
            ds._simdir = self._simdir
            ds._target = self._target
            ds._rng = self._rng
            ds._rng_seed = self._rng_seed
            ds._rng_state = self._rng_state
            ds._smallest_length = self._smallest_length
            ds._area = self._area
            ds._compute_face_areas = False
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         simdir=self.simdir,
                         compute_face_areas=False,
                         rng=self.rng,
                         rng_seed=self.rng_seed,
                         rng_state=self.rng_state)
            ds._smallest_length = self._smallest_length
            ds._area = self._area
        return ds   

    @property
    def uxds(self) -> UxDataset:
        """
        The data associated with the surface. This is an instance of UxDataset.
        """
        return self._uxds

    @uxds.setter
    def uxds(self, value: UxDataset):
        if not isinstance(value, UxDataset):
            raise TypeError("data must be an instance of UxDataset")
        self._uxds = value

    @property
    def data_dir(self):
        """
        Directory for data files.
        """
        return self.simdir / _DATA_DIR

    @property
    def grid_file(self):
        """
        Path to the grid file.
        """
        return self.simdir / _DATA_DIR / _GRID_FILE_NAME

    @property
    def grid(self):
        """
        The type of grid used for the surface.
        """
        return self._grid
    
    @grid.setter
    def grid(self, value):
        if not isinstance(value, (Grid, str)):
            raise TypeError("grid must be a string or Grid object")
        self._grid = value
    
    @parameter
    def smallest_length(self):
        """
        Smallest length value that is directly modeled on the grid. This is used to determine the maximum distance of ejecta to 
        consider, for instance
        """
        return self._smallest_length

    @smallest_length.setter
    def smallest_length(self, value):
        if not isinstance(value, float):
            raise TypeError("smallest_length must be a float")
        self._smallest_length = value

    @parameter
    def grid_config(self):
        """
        The grid configuration used for the surface.
        """
        return self._grid_config
    
    @grid_config.setter
    def grid_config(self, value):
        if not isinstance(value, dict):
            raise TypeError("grid_config must be a dictionary")
        self._grid_config = value

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
        self._target = Target.make(value)
        return 
    
    def to_config(self, remove_common_args: bool = False, **kwargs: Any) -> dict[str, Any]:
        """
        Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file. 

        Parameters
        ----------
        obj : Any
            The object whose attributes will be stored.  It must have a _user_defined attribute.
        remove_common_args : bool, optional
            If True, remove the set of common arguments that are shared among all components of the project from the configuration. Defaults to False.
        **kwargs : Any
            Additional keyword arguments for subclasses.

        Returns
        -------
        dict[str, Any]
            A dictionary of the object's attributes that can be serialized to YAML.
        Notes
        -----
        - The function will ignore any attributes that are not serializable to human-readable YAML. Therefore, it will ignore anything that cannot be converted into a str, int, float, or bool.
        - The function will convert Numpy types to their native Python types.
        """
        return _to_config(self, remove_common_args=remove_common_args, **kwargs)

    @property
    def simdir(self):
        """
        The main project simulation directory.

        Returns
        -------
        Path
            The initialized simulation directory as a Path object. Will be a relative path if possible, otherwise will be absolute. If it doesn't exist, it will be created.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value):
        self._simdir = _simdir_init(value)
        
    @property
    def rng_seed(self):
        """
        The random rng_seed for the simulation RNG.

        Returns
        -------
        int or None
            The integer rng_seed used to initialize the RNG, or None if not set.
        """
        return self._rng_seed

    @rng_seed.setter
    def rng_seed(self, value):
        if value is not None:
            if not isinstance(value, int) or np.isnan(value) or np.isinf(value) or value < 0:
                raise TypeError("rng_seed must be a positive integer")
            self._rng_seed = int(value)
        else:
            self._rng_seed = None

    @property
    def rng(self):
        """
        The random number generator used for stochastic elements of the simulation.

        Returns
        -------
        numpy.random.Generator or None
            The RNG instance, or None if not initialized.
        """
        return self._rng

    @property
    def node_tree(self):
        if self._node_tree is None:
            self._node_tree = self.uxds.uxgrid.get_ball_tree("nodes", distance_metric="haversine", coordinate_system="spherical", reconstruct=True)
        
        return self._node_tree

    @property
    def face_tree(self):
        if self._face_tree is None:
            self._face_tree = self.uxds.uxgrid.get_ball_tree("face centers",  distance_metric="haversine", coordinate_system="spherical", reconstruct=True)
        
        return self._face_tree
    
    @rng.setter
    def rng(self, value):
        self._rng, _ = _rng_init(rng=value, rng_seed=self.rng_seed, rng_state=self.rng_state)

    @property 
    def rng_state(self):
        """
        The state of the random number generator.

        Returns
        -------
        dict or None
            A dictionary representing the RNG state, or None if the RNG is not initialized.
        """
        return self.rng.bit_generator.state if self.rng is not None else None
    
    @rng_state.setter
    def rng_state(self, value):
        _, self._rng_state = _rng_init(rng=self.rng, rng_seed=self.rng_seed, rng_state=value)

    @property
    def common_args(self) -> CommonArgs:
        return CommonArgs(simdir=self.simdir, rng=self.rng, rng_seed=self.rng_seed, rng_state=self.rng_state)    
        
    def generate_data(self,
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
        uxgrid = uxr.open_grid(self.grid_file,use_dual=False)
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
            data = np.zeros(size,dtype=np.float64) 
        elif np.isscalar(data):
            data = np.full(size,data)
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
         
        self.uxds[name] = uxda
        
        if save_to_file:
            self._save_data(uxda, self.data_dir, interval_number, combine_data_files)
        return 
        
    def set_elevation(self, 
                    new_elev: NDArray[np.float64] | List[FloatLike] | None = None,
                    save_to_file: bool = False, 
                    combine_data_files: bool = False,
                    interval_number: int = 0,
                    ) -> None:
        """
        Set elevation data for the target's surface mesh.

        Parameters
        ----------
        new_elev : array_like, optional
            New elevation data to be set. If None, the elevation is set to zero. 
        save_to_file : bool, default False
            If True, save the elevation data to the elevation file.
        combine_data_files : bool, default False
            If True, combine the current data with the existing data for previous intervals in the data file.
        interval_number : int, default 0
            The interval number to use when saving the elevation data to the elevation file.
        """
        if new_elev is None or np.isscalar(new_elev):
            gen_node = True
            gen_face = True 
        elif new_elev.size == self.uxds.uxgrid.n_node:
            gen_node = True
            gen_face = False
        elif new_elev.size == self.uxds.uxgrid.n_face:
            gen_node = False
            gen_face = True
        else:
            gen_node = False
            gen_face = False
            raise ValueError("new_elev must be None, a scalar, or an array with the same size as the number of nodes in the grid")
    
        try:
            self.save_to_data()
        except AttributeError:
            pass
            
        if gen_node:
            self.generate_data(data=new_elev, 
                               name="node_elevation",
                               long_name="elevation of nodes",
                               units= "m",
                               isfacedata=False,
                               save_to_file=save_to_file,
                               combine_data_files=combine_data_files,
                               interval_number=interval_number
                              )   
        if gen_face:
            self.generate_data(data=new_elev, 
                               name="face_elevation",
                               long_name="elevation of faces",
                               units="m",
                               isfacedata=True,
                               save_to_file=save_to_file,
                               combine_data_files=combine_data_files,
                               interval_number=interval_number
                              )
        self.load_from_data(compute_face_areas=False)
        return 

    @staticmethod
    def calculate_haversine_distance(lon1: FloatLike, 
                    lat1: FloatLike, 
                    lon2: FloatLike, 
                    lat2: FloatLike,
                    radius: FloatLike = 1.0) -> float:
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
        radius : FloatLike
            Radius of the sphere in meters.

        Returns
        -------
        float
            Great circle distance between the two points in meters.
        """
        # Calculate differences in coordinates
        dlon = lon2 - lon1
        dlat = lat2 - lat1

        # Haversine formula
        a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
        c = 2 * np.arcsin(np.sqrt(a))
        return radius * c
    
    def get_distance(self, 
                     view: 'SurfaceView',
                     location: tuple[float, float]) -> UxDataArray:
        """
        Computes the distances between nodes and faces and a given location.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        UxDataArray
            DataArray of distances for each node in meters.
        UxDataArray
            DataArray of distances for each face in meters.
        """
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.node_lon[view.node_indices])
        node_lat2 = np.deg2rad(self.node_lat[view.node_indices])
        face_lon2 = np.deg2rad(self.face_lon[view.face_indices])
        face_lat2 = np.deg2rad(self.face_lat[view.face_indices])
        return self.calculate_haversine_distance(lon1,lat1,node_lon2,node_lat2,self.target.radius), self.calculate_haversine_distance(lon1,lat1,face_lon2,face_lat2,self.target.radius) 

    @staticmethod
    def calculate_initial_bearing(lon1: FloatLike, 
                                lat1: FloatLike, 
                                lon2: FloatLike, 
                                lat2: FloatLike) -> float:
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
    
    def get_initial_bearing(self, view: 'SurfaceView', location: tuple[float, float]) -> tuple[NDArray, NDArray]:
        """
        Computes the initial bearing between nodes and faces and a given location.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        DataArray
            DataArray of initial bearings for each node in radians.
        DataArray
            DataArray of initial bearings for each face in radians.
        """
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.node_lon[view.node_indices])
        node_lat2 = np.deg2rad(self.node_lat[view.node_indices])
        face_lon2 = np.deg2rad(self.face_lon[view.face_indices])
        face_lat2 = np.deg2rad(self.face_lat[view.face_indices])
        return (surface_functions.calculate_initial_bearing(lon1,lat1,node_lon2,node_lat2), 
                surface_functions.calculate_initial_bearing(lon1,lat1,face_lon2,face_lat2))
    
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
            The index of the nearest node in the grid to the given point.
        int
            The index of the nearest face in the grid to the given point.

        Notes
        -----
        The method uses the ball tree query method that is included in the UxArray.Grid class.
        """          
        
        coords = np.asarray(location)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", Warning)
            
            node_ind = self.node_tree.query(coords=coords, k=1, return_distance=False)
        
            face_ind = self.face_tree.query(coords=coords, k=1, return_distance=False)
        return node_ind.item(), face_ind.item()

    def get_reference_surface(self, 
                              view: 'SurfaceView',
                              face_crater_distance: NDArray, 
                              node_crater_distance: NDArray,
                              location: tuple[FloatLike, FloatLike], 
                              region_radius: float) -> NDArray[np.float64]:
        """
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.

        Parameters
        ----------
        location : tuple[float, float]
            tuple containing the longitude and latitude of the reference location in degrees.
        region_radius : float
            The radius of the region to compute the average over in meters.

        Returns
        -------
        center_vector : ndarray
            The vector pointing to the center of the cap from the sphere's center.
        """
        def sphere_function(coords, x_c, y_c, z_c, r):
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
            return (x - x_c)**2 + (y - y_c)**2 + (z - z_c)**2 - r**2
   
        # Find cells within the crater radius
        faces_within_radius = face_crater_distance <= region_radius
        nodes_within_radius = node_crater_distance <= region_radius

        face_elevation = self.face_elevation[view.face_indices]
        node_elevation = self.node_elevation[view.node_indices]
        
        if np.sum(faces_within_radius) < 5 or not nodes_within_radius.any():
            return face_elevation, node_elevation
       
        face_grid = np.column_stack((self.uxds.uxgrid.face_x.values[view.face_indices], 
                               self.uxds.uxgrid.face_y.values[view.face_indices], 
                               self.uxds.uxgrid.face_z.values[view.face_indices]))
        region_faces = face_grid[faces_within_radius]
        region_elevation = face_elevation[faces_within_radius] / self.target.radius
        region_surf = self.elevation_to_cartesian(region_faces, region_elevation) 

        x, y, z = region_surf.T

        # Initial guess for the sphere center and radius
        guess_radius = 1.0 + region_elevation.mean()
        initial_guess = [0, 0, 0, guess_radius]  

        # Perform the curve fitting
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            try:
                bounds =([-1.0, -1.0, -1.0, 0.5],[1.0, 1.0, 1.0, 2.0]) 
                popt, _  = curve_fit(sphere_function, region_surf[faces_within_radius], np.zeros_like(x[faces_within_radius]), p0=initial_guess, bounds=bounds)
            except:
                popt = initial_guess

        # Extract the fitted sphere center and radius
        reference_sphere_center = popt[:3]
        reference_sphere_radius = popt[3]
        
        def find_reference_elevations(coords):
            # Find the point along the original vector that intersects the sphere 
            f_vec = coords / self.target.radius      
            A = f_vec[:,0]**2 + f_vec[:,1]**2 + f_vec[:,2]**2
            B = -2 * (f_vec[:,0] * reference_sphere_center[0] + f_vec[:,1] * reference_sphere_center[1] + f_vec[:,2] * reference_sphere_center[2])
            C = np.dot(reference_sphere_center, reference_sphere_center) - reference_sphere_radius**2
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
                
            elevations = self.target.radius * (t * np.linalg.norm(f_vec, axis=1)  - 1)
            return elevations
        
        # Calculate the distances between the original face points and the intersection points
        reference_face_elevation = face_elevation
        reference_face_elevation[faces_within_radius] = find_reference_elevations(region_faces)
        
        # Now do the same thing to compute the nodal values 
        node_grid = np.column_stack((self.uxds.uxgrid.node_x.values[view.node_indices], 
                               self.uxds.uxgrid.node_y.values[view.node_indices], 
                               self.uxds.uxgrid.node_z.values[view.node_indices]))
        region_nodes = node_grid[nodes_within_radius]
        
        reference_node_elevation = node_elevation
        reference_node_elevation[nodes_within_radius] = find_reference_elevations(region_nodes)
        
        return reference_face_elevation, reference_node_elevation

    @staticmethod
    def elevation_to_cartesian(position: NDArray, 
                            elevation: NDArray
                            ) -> NDArray:
        
        
        runit = position / np.linalg.norm(position, axis=1, keepdims=True)
        
        return position + elevation[:, np.newaxis] * runit
   
    def extract_region(self,
                       location: tuple[FloatLike, FloatLike],
                       region_radius: FloatLike): 
        
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
        Surface
            A new Surface object containing the regional grid. 
            
        Notes
        -----
        Though not well documented, the mapping from the region grid back to the original grid can be done using the variables
        "subgrid_face_indices" and "subgrid_node_indices" in the region grid. For example, the following will extract a region
        surface, modify the face elevation, then insert the modified values back onto the original surface:
        
        .. code-block:: python    
              
            import xarray as xr
                
            region_surf = surf.extract_region(location, region_radius)
            region_surf['face_elevation'] = xr.full_like(region_surf['face_elevation'], 1.0)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] = region_surf['face_elevation']
            
        """ 
        
        region_angle = np.rad2deg(region_radius / self.target.radius)
        coords = np.asarray(location)

        ind = self.face_tree.query_radius(coords, region_angle)
        if len(ind) == 0:
            return None
        
        return SurfaceView(self, ind)
    
    def get_random_location_on_face(self, 
                                    face_index: int, 
                                    size: int = 1,
                                    **kwargs
                                    ) -> Union[float, tuple[float, float], ArrayLike]:
        """
        Generate a random coordinate within a given face of an unstructured mesh.

        Parameters
        ----------
        grid : uxarray.Grid
            The grid object containing the mesh information.
        face_index : int
            The index of the face within the grid to obtain the random sample.
        size : int or tuple of ints, optional
            The number of samples to generate. If size is None (the default), a single tuple is returned. If size is greater than 1, 
            then a structured array with fields 'lon' and 'lat' is returned.
                
        Returns
        -------
        (lon,lat) or ndarray[(lon,lat)] of given size
            A pair or array of pairs of longitude and latitude values in degrees.

        Notes
        -----
        This method is a wrapper for :func:`cratermaker.utils.montecarlo.get_random_location_on_face`. 
        """
        
        return get_random_location_on_face(self.uxds.uxgrid, face_index, size,**kwargs)

    @staticmethod
    def _save_data(ds: xr.Dataset | xr.DataArray,
                out_dir: os.PathLike,
                interval_number: int = 0,
                combine_data_files: bool = False
                ) -> None:
        """
        Save the data to the specified directory. If `combine_data_files` is True, then all data variables are saved to a single NetCDF
        file. If False, then only the data variables for the current interval are saved to a NetCDF file with the interval number
        appended. 
        
        Parameters
        ----------
        ds : xr.Dataset or xr.DataArray
            The data to be saved.
        out_dir : PathLike
            Directory to save the data.
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
            ds = ds.assign_coords({"time":[interval_number]})      
            
            
        with tempfile.TemporaryDirectory() as temp_dir:
            if combine_data_files:
                filename = _COMBINED_DATA_FILE_NAME
            else:
                filename = _COMBINED_DATA_FILE_NAME.replace(".nc", f"{interval_number:06d}.nc")
                
            data_file = os.path.join(out_dir, filename)
            if os.path.exists(data_file):
                ds_file = xr.open_mfdataset(data_file)
                ds_file = ds.merge(ds_file, compat="override")
            else:
                ds_file = ds    
                
            temp_file = os.path.join(temp_dir, filename)
            
            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in ds_file.data_vars}
            ds_file.to_netcdf(temp_file, encoding=encoding)
            ds_file.close()     
            shutil.move(temp_file, data_file)

        return



class SurfaceView:
    def __init__(self, 
                 surf: Surface, 
                 face_indices: NDArray | slice, 
                 node_indices: NDArray | slice | None = None):
        self.surf = surf
        self.face_indices = face_indices
        if node_indices is None:
            node_indices = np.unique(surf.uxds.uxgrid.face_node_connectivity.values[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]
        
        self.node_indices = node_indices


def _save_surface(surf: Surface, 
         out_dir: os.PathLike | None = None,
         combine_data_files: bool = False,
         interval_number: int = 0,
         time_variables: dict | None = None,
         *args, **kwargs, 
         ) -> None:
    """
    Save the surface data to the specified directory. Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the time dimension. If 'interval_number' is included as a key in `time_variables`, then this will be appended to the data file name.

    Parameters
    ----------
    surface : Surface
        The surface object to be saved. 
    out_dir : str, optional
        Directory to save the surface data. If None, the data is saved to the current working directory.
    combine_data_files : bool, optional
        If True, combine all data variables into a single NetCDF file, otherwise each variable will be saved to its own NetCDF file. Default is False.
    interval_number : int, optional
        Interval number to append to the data file name. Default is 0.
    time_variables : dict, optional
        Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the time dimension. Default is None.
    """
    do_not_save = ["face_areas"]
    if out_dir is None:
        out_dir = surf.data_dir
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)         
      
    if time_variables is None:
        time_variables = {"elapsed_time":float(interval_number)}  
    else:
        if not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")
        
    # Variables that we do not want to save as they are computed at runtime     
    
    surf.save_to_data()
    surf.uxds.close()
    
    ds = surf.uxds.expand_dims(dim="time").assign_coords({"time":[interval_number]})
    for k, v in time_variables.items():
        ds[k] = xr.DataArray(data=[v], name=k, dims=["time"], coords={"time":[interval_number]})
                
    drop_vars = [k for k in ds.data_vars if k in do_not_save]
    if len(drop_vars) > 0:
        ds = ds.drop_vars(drop_vars)
        
    surf._save_data(ds, out_dir, interval_number, combine_data_files)

    return

