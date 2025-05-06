from __future__ import annotations
from pathlib import Path
from abc import abstractmethod
import os
import shutil
import tempfile
import numpy as np
import uxarray as uxr
from uxarray import INT_FILL_VALUE, UxDataArray, UxDataset
import xarray as xr
from scipy.optimize import curve_fit, OptimizeWarning
from numpy.typing import NDArray, ArrayLike
from typing import Any
import hashlib
from cratermaker.constants import _GRID_FILE_NAME, _DATA_DIR, _SMALLFAC, _COMBINED_DATA_FILE_NAME 
from cratermaker.components.target import Target
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from cratermaker.utils.montecarlo_utils import get_random_location_on_face
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker._cratermaker import surface_functions
import warnings

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

    def __init__(self, 
                 target: Target | str | None = None,
                 simdir: str | Path | None = None,
                 **kwargs):

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
        object.__setattr__(self, "_smallest_length", None)
        super().__init__(simdir=simdir, **kwargs)
                       
        self._data_variable_init = {
            "node_elevation": {
                "units" : "m", 
                "long_name": "elevation of nodes", 
                "initial_value" : 0.0,
                "isfacedata" : False 
            },
            "face_elevation": {
                "units" : "m",
                "long_name": "elevation of faces",
                "initial_value" : 0.0,
                "isfacedata" : True
            },
            "ejecta_thickness": {
                "units" : "m",
                "long_name": "ejecta thickness",
                "initial_value" : 0.0,
                "isfacedata" : True
            },
            "ray_intensity": {
                "units" : "",
                "long_name": "ray intensity value",
                "initial_value" : 0.0,
                "isfacedata" : True
            },
        }

        self.target = Target.maker(target, **kwargs)

        return

    def __repr__(self) -> str:
        base = super().__repr__()
        return (
            f"{base}\n"
            f"Target: {self.target.name}\n"
            f"Grid File: {self.grid_file}"
        )

    @classmethod
    def maker(cls: Surface, 
             surface: str | Surface | None = None,
             target: Target | str | None = None,
             reset: bool = False,
             regrid: bool = False, 
             simdir: str | Path | None = None,
             **kwargs) -> Surface:
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

        surface = super().maker(component=surface, 
                             target=target, 
                             reset = reset,
                             regrid = regrid,
                             simdir=simdir, 
                             **kwargs)
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
                if reset: # Create an empty dataset 
                    self._uxds = uxr.UxDataset()
                else: # Read data from from existing datafiles
                    with uxr.open_mfdataset(uxgrid, data_file_list, use_dual=False) as ds:
                        self._uxds = ds.isel(time=-1).load()
                self._uxds.uxgrid = uxr.Grid.from_dataset(uxgrid)
                self._uxgrid = uxgrid
        except:
            raise ValueError("Error loading grid and data files")
        
        if reset:
            self.reset(**kwargs)
        else:
            self.load_from_data()

        return

    def save_to_files(self, 
                      combine_data_files: bool = False,
                      interval_number: int = 0,
                      time_variables: dict | None = None,
                      *args, **kwargs) -> None:
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
            time_variables = {"elapsed_time":float(interval_number)}  
        else:
            if not isinstance(time_variables, dict):
                raise TypeError("time_variables must be a dictionary")
            
        # Variables that we do not want to save as they are computed at runtime     
        
        self.save_to_data()
        self.uxds.close()
        
        ds = self.uxds.expand_dims(dim="time").assign_coords({"time":[interval_number]})
        for k, v in time_variables.items():
            ds[k] = xr.DataArray(data=[v], name=k, dims=["time"], coords={"time":[interval_number]})
                    
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
        shutil.move(temp_file.name,self.grid_file)      
        self._uxgrid = uxgrid

        regrid = self.regrid_if_needed(**kwargs)
        assert(not regrid)

        self._pix_mean, self._pix_std, self._pix_min, self._pix_max = self._compute_pix_size(uxgrid)

        return

    def regrid_if_needed(self, force : bool = False, **kwargs: Any) -> bool:
        """
        Check if the existing grid matches the desired parameters determine if regridding is necessary.

        This function checks if a grid file exists and matches the specified parameters based on a unique hash generated from these 
        parameters. If the grid does not exist or does not match the parameters it generates a new grid and returns True. 

        Parameters
        ----------
        force : bool, optional
            Flag to force regridding even if the grid file exists. Default is False.

        Returns
        -------
        bool
            A boolean indicating whether the grid should be regenerated. 
        """
    
        # Find out if the file exists, if it does't we'll need to make a new grid
        self.data_dir.mkdir(parents=True, exist_ok=True)
        regrid = force or not Path(self.grid_file).exists()

        if not regrid:
            try: 
                with xr.open_dataset(self.grid_file) as ds:
                    ds.load()
                    uxgrid = uxr.Grid.from_dataset(ds)
                    old_id = uxgrid.attrs.get("_id")
                    regrid = old_id != self._id
            except:
                regrid = True

        if regrid:
            print("Creating a new grid")
            try:
                self.generate_grid(**kwargs)
            except:
                raise RuntimeError("Failed to create a new grid")

        return regrid

    def reset(self, **kwargs: Any) -> None:
        """
        Reset the surface to its initial state.
        """
        for name, entry in self._data_variable_init.items():
            self.generate_data(name=name,
                               data=entry["initial_value"],
                               long_name=entry["long_name"],
                               units=entry["units"],
                               isfacedata=entry["isfacedata"],
                               save_to_file=True)
                               

        self.load_from_data()

    def load_from_data(self):
        self.face_lat = self.uxgrid.face_lat.values
        self.face_lon = self.uxgrid.face_lon.values
        self.node_lat = self.uxgrid.node_lat.values
        self.node_lon = self.uxgrid.node_lon.values
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
        return self.simdir / _DATA_DIR

    @property
    def grid_file(self):
        """
        Path to the grid file.
        """
        return self.simdir / _DATA_DIR / _GRID_FILE_NAME

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
            self._node_tree = self.uxgrid.get_ball_tree("nodes", distance_metric="haversine", coordinate_system="spherical", reconstruct=True)
        
        return self._node_tree

    @property
    def face_tree(self):
        if self._face_tree is None:
            self._face_tree = self.uxgrid.get_ball_tree("face centers",  distance_metric="haversine", coordinate_system="spherical", reconstruct=True)
        
        return self._face_tree
    
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
         
        self._uxds[name] = uxda
        
        if save_to_file:
            self._save_data(uxda, interval_number, combine_data_files)
        return 
        
    def set_elevation(self, 
                    new_elev: NDArray[np.float64] | list[FloatLike] | None = None,
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
        self.load_from_data()
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
                     location: tuple[float, float]) -> tuple[float, float]:
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
        if len(location) == 1:
            location = location.item()
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.node_lon[view.node_indices])
        node_lat2 = np.deg2rad(self.node_lat[view.node_indices])
        face_lon2 = np.deg2rad(self.face_lon[view.face_indices])
        face_lat2 = np.deg2rad(self.face_lat[view.face_indices])
        return self.calculate_haversine_distance(lon1,lat1,node_lon2,node_lat2,self.radius), self.calculate_haversine_distance(lon1,lat1,face_lon2,face_lat2,self.radius) 

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
        if len(location) == 1:
            location = location.item()
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
        
        if len(location) == 1:
            location = location.item()
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
        region_elevation = face_elevation[faces_within_radius] / self.radius
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
            f_vec = coords / self.radius      
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
                
            elevations = self.radius * (t * np.linalg.norm(f_vec, axis=1)  - 1)
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
                
            region_surf = surface.extract_region(location, region_radius)
            region_surf['face_elevation'] = xr.full_like(region_surf['face_elevation'], 1.0)
            surface['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] = region_surf['face_elevation']
            
        """ 
        
        region_angle = np.rad2deg(region_radius / self.radius)
        if len(location) == 1:
            location = location.item()
        coords = np.asarray(location)  

        ind = self.face_tree.query_radius(coords, region_angle)
        if len(ind) == 0:
            return None
        
        return SurfaceView(self, ind)
    
    def get_random_location_on_face(self, 
                                    face_index: int, 
                                    **kwargs
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
        
        return get_random_location_on_face(self.uxgrid, face_index, rng=self.rng, **kwargs)

    def _save_data(self,
                   ds: xr.Dataset | xr.DataArray,
                   interval_number: int = 0,
                   combine_data_files: bool = False) -> None:
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
            ds = ds.assign_coords({"time":[interval_number]})      
            
            
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
            
            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in ds_file.data_vars}
            ds_file.to_netcdf(temp_file, encoding=encoding)
            ds_file.close()     
            shutil.move(temp_file, data_file)

        return
    
    def __del__(self):
        try:
            if hasattr(self, "_uxds") and hasattr(self._uxds, "close"):
                if hasattr(self._uxds, "uxgrid") and hasattr(self._uxds.uxgrid, "_ds"):
                    self._uxds.uxgrid._ds.close()
                self._uxds.close()
        except Exception:
            pass
        try:
            if hasattr(self, "_grid") and hasattr(self._grid, "uxgrid") and hasattr(self._grid.uxgrid, "_ds"):
                self._grid.uxgrid._ds.close()
        except Exception:
            pass

    @abstractmethod
    def generate_face_distribution(self, **kwargs: Any) -> tuple[NDArray,NDArray,NDArray]: ...

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
    def _distribute_points(distance: FloatLike,
                           radius: FloatLike=1.0, 
                           lon_range: PairOfFloats=(-180,180), 
                           lat_range: PairOfFloats=(-90,90)) -> NDArray:
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
            
        .. rubric:: References
        
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
        theta_range = np.deg2rad(lat_range) + np.pi/2 
        points = []
        
        n = int(1/distance**2)
        if n < 1: 
            return 
            
        a = 4 * np.pi / n
        d = np.sqrt(a)
        Mtheta = int(np.round(np.pi / d))
        dtheta = np.pi / Mtheta
        dphi = a / dtheta
        
        thetavals = np.pi * (np.arange(Mtheta) + 0.5) / Mtheta
        thetavals = thetavals[(thetavals >= theta_range[0]) & (thetavals < theta_range[1])]
        
        for theta in thetavals:
            Mphi = int(np.round(2 * np.pi * np.sin(theta) / dphi))
            phivals = 2 * np.pi * np.arange(Mphi) / Mphi
            phivals = phivals[(phivals >= phi_range[0]) & (phivals < phi_range[1])]
            for phi in phivals:
                points.append(_sph2cart(theta, phi, radius))
        if len(points) == 0:
            return
        
        points = np.array(points,dtype=np.float64)
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
            self._pix_mean, self._pix_std, self._pix_min, self._pix_max = self._compute_pix_size()
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
            self._pix_mean, self._pix_std, self._pix_min, self._pix_max = self._compute_pix_size()
        return self._pix_max
    
    @property
    def gridtype(self):
        """
        The name of the grid type.
        """
        return self._component_name
    
    @property
    def face_areas(self):
        """
        The face areas of the mesh.
        """
        if self._face_areas is None:
            self._face_areas = self.uxgrid.face_areas.values * self.radius**2
        return self._face_areas
    
    @property
    def smallest_length(self):
        """
        The smallest length of the mesh.
        """
        if self._smallest_length is None:
            self._smallest_length = np.sqrt(np.min(self.face_areas)) * _SMALLFAC  
        return self._smallest_length


class SurfaceView:
    def __init__(self, 
                 surface: Surface, 
                 face_indices: NDArray | slice, 
                 node_indices: NDArray | slice | None = None):
        self.surface = surface
        self.face_indices = face_indices
        if node_indices is None:
            node_indices = np.unique(surface.uxds.uxgrid.face_node_connectivity.values[face_indices].ravel())
            node_indices = node_indices[node_indices != INT_FILL_VALUE]
        
        self.node_indices = node_indices






import_components(__name__, __path__, ignore_private=True)
