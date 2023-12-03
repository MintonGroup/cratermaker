import xarray as xr
import uxarray as uxr
from uxarray import UxDataset
from glob import glob
import os
import numpy as np
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
from pathlib import Path
from ..utils.general_utils import float_like
from .target import Target
from typing import Tuple, List
from numpy.typing import NDArray


class Surface(UxDataset):
    __slots__ = UxDataset.__slots__ + ('_name', '_description', 'target_name', 'pix', 'grid_type')    
    """Surface class for cratermaker"""
    def __init__(self, *args, **kwargs):

        # Call the super class constructor with the dataset
        super().__init__(*args, **kwargs)
        
        # Additional initialization for Surface
        self._name = "Surface"
        self._description = "Surface class for cratermaker"
        

    def set_elevation(self, new_elev: NDArray[np.float64] | List[float_like] | None = None) -> uxr.UxDataArray:
        """
        Set elevation data for the target's surface mesh.

        Parameters
        ----------
        nCells: int, optional
            The number of cells in the mesh. Must be passed if new_elev is not
        new_elev : array_like, optional
            New elevation data to be set. If None, the elevation is set to zero. Must be passed if nCells is not
        """        
       
        if new_elev is None:
            new_elev = np.zeros(self.uxgrid.n_face,dtype=np.float64)
        else:
            if new_elev.size != self.uxgrid.n_face:
                raise ValueError("new_elev must have the same size as the number of faces in the grid")
        uxda = uxr.UxDataArray(
                data=new_elev,
                dims=["n_face"],
                attrs={"long_name":"elevation of cells"},
                uxgrid=self.uxgrid
                )
        
        return uxda


    @staticmethod
    def calculate_haversine_distance(lon1: float_like, 
                    lat1: float_like, 
                    lon2: float_like, 
                    lat2: float_like,
                    radius: float_like) -> np.float64:
        """
        Calculate the great circle distance between two points on a sphere.

        Parameters
        ----------
        lon1 : float_like
            Longitude of the first point in radians.
        lat1 : float_like
            Latitude of the first point in radians.
        lon2 : float_like
            Longitude of the second point in radians.
        lat2 : float_like
            Latitude of the second point in radians.
        radius : float_like
            Radius of the sphere in meters.

        Returns
        -------
        np.float64
            Great circle distance between the two points in meters.
        """
        # Calculate differences in coordinates
        dlon = lon2 - lon1
        dlat = lat2 - lat1

        # Haversine formula
        a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
        c = 2 * np.arcsin(np.sqrt(a))
        return radius * c
    

    def get_cell_distance(self, 
                        location: Tuple[np.float64, np.float64],
                        target_radius: np.float64) -> uxr.UxDataArray:
        """
        Computes the distances between cell centers and a given location.

        Parameters
        ----------
        location : Tuple[np.float64, np.float64]
            Tuple containing the longitude and latitude of the location in radians.

        Returns
        -------
        UxArray.UxDataArray
            DataArray of distances for each cell in meters.
        """
        return self.calculate_haversine_distance(location[0],location[1],self.uxgrid.face_lon,self.uxgrid.face_lat,target_radius)
    

    @staticmethod
    def calculate_initial_bearing(lon1: float_like, 
                                lat1: float_like, 
                                lon2: float_like, 
                                lat2: float_like) -> np.float64:
        """
        Calculate the initial bearing from one point to another on the surface of a sphere.

        Parameters
        ----------
        lon1 : float_like
            Longitude of the first point in radians.
        lat1 : float_like
            Latitude of the first point in radians.
        lon2 : float_like
            Longitude of the second point in radians.
        lat2 : float_like
            Latitude of the second point in radians.

        Returns
        -------
        np.float64
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


    def get_cell_initial_bearing(self, location: Tuple[np.float64, np.float64]) -> uxr.UxDataArray:
        """
        Computes the initial bearing between cell centers and a given location.

        Parameters
        ----------
        location : Tuple[np.float64, np.float64]
            Tuple containing the longitude and latitude of the location in radians.

        Returns
        -------
        xarray.DataArray
            DataArray of initial bearings for each cell in radians.
        """
        return self.calculate_initial_bearing(location[0], location[1], self.uxgrid.face_lon, self.uxgrid.face_lat)


    def get_average_surface(self,
                            location: Tuple[float_like, float_like], 
                            radius: np.float64) -> Tuple[np.float64, np.float64]:
        """
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.

        Parameters
        ----------
        location : Tuple[float, float]
            Tuple containing the longitude and latitude of the reference location in radians.
        radius : float
            The reference radius of to compute the average over in meters.

        Returns
        -------
        center_vector : ndarray
            The vector pointing to the center of the cap from the sphere's center.
        radius : float
            The radius of the cap.
        """

        # Find cells within the crater radius
        # cells_within_radius = data['crater_distance'] <= radius

        # bearings = data['crater_bearing'].where(cells_within_radius, drop=True)
        # distances = data['crater_distance'].where(cells_within_radius, drop=True)

        # # Convert bearings to vector components
        # # Bearing is angle from north, positive clockwise, but we need standard mathematical angle, positive counter-clockwise
        # angles = np.deg2rad(90) - bearings  # Convert bearing to angle in radians
        # x_components = np.cos(angles) * distances
        # y_components = np.sin(angles) * distances

        # # Calculate the weighted average vector components
        # # Weight by the area of each cell to give more importance to larger cells
        # cell_areas = mesh['areaCell'].where(cells_within_radius, drop=True)
        # weighted_x = (x_components * cell_areas).sum() / cell_areas.sum()
        # weighted_y = (y_components * cell_areas).sum() / cell_areas.sum()

        # # Calculate the weighted mean elevation to get the z-component
        # elevation_values = data['elevation'].where(cells_within_radius, drop=True)
        # weighted_z = (elevation_values * cell_areas).sum() / cell_areas.sum()

        # # Combine components to form the cap center vector
        # center_vector = -np.array([weighted_x.item(), weighted_y.item(), weighted_z.item()])

        # # The radius of the cap is the length of the cap center vector
        # radius = np.linalg.norm(center_vector)
        center_vector = None
        return center_vector 

    

def initialize_surface(grid_file: str = "surface_grid.nc",
         dem_file: str = "surface_dem.nc",
         data_files: str = "surface_*.nc",
         make_new_grid: bool = False,
         reset_surface: bool = True,
         grid_temp_dir: str | Path | os.PathLike = ".mesh",
         pix: float_like | None = None,
         target: Target | str | None = None,
         *args, **kwargs):

    if not target:
        target = Target("Moon")
    elif isinstance(target, str):
        try:
            target = Target(target)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a ")
    
    # Load the grid and data files
    data_file_list = glob(data_files)
    if grid_file in data_file_list:
        data_file_list.remove(grid_file)
    if dem_file not in data_file_list:
        data_file_list.append(dem_file)        
    
    if not os.path.exists(grid_temp_dir):
        os.mkdir(grid_temp_dir)
    
    grid_temp_dir = Path(grid_temp_dir)
    grid_file = Path(grid_file)
    if not grid_file.is_absolute():
        grid_file = Path.cwd() / grid_file
        
    dem_file = Path(dem_file)
    if not dem_file.is_absolute():
        dem_file = Path.cwd() / dem_file   
    grid_file = str(grid_file)
    dem_file = str(dem_file)
        
    # Generate a new surface if either it is explicitly requested via parameter or a data file doesn't yet exist 
    make_new_grid = make_new_grid or not os.path.exists(grid_file)
    reset_surface = reset_surface or not os.path.exists(dem_file) or make_new_grid

    if make_new_grid:
        generate_grid(target.radius,pix,grid_file,dem_file,grid_temp_dir, *args, **kwargs)
       
    if reset_surface:
        generate_surface_dem(grid_file,dem_file)
         
    # Initialize UxDataset with the loaded data
    surf = uxr.open_mfdataset(grid_file, data_file_list, latlon=True, use_dual=False)
    
    return Surface(surf,uxgrid=surf.uxgrid)
        
        
def generate_grid(target_radius: float_like, 
                cell_size: float_like, 
                grid_file: str | Path | os.PathLike = "surface_mesh.nc", 
                dem_file: str | Path | os.PathLike = "surface_dem.nc", 
                grid_temp_dir:  str | Path | os.PathLike = ".mesh",
                *args, **kwargs) -> Surface:
    """
    Generate a tessellated mesh of a sphere using the jigsaw-based mesh builder in MPAS-tools.

    This function generates temporary files in the `grid_temp_dir` directory and saves the final mesh to `grid_file`

    Returns
    -------
    xarray Dataset
    An MPAS-style dataset with the mesh data and variables 
    """
    def make_uniform_cell_size(cell_size: float_like) -> Tuple[NDArray,NDArray,NDArray]:
        """
        Create cell width array for this mesh on a regular latitude-longitude grid.
        Returns
        -------
        cellWidth : ndarray
            m x n array of cell width in km
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            longitude in degrees (length m and between -90 and 90)
        """
        dlat = 10
        dlon = 10
        constantCellWidth = cell_size * 1e-3 # build_spherical_mesh assumes units of km, so must be converted

        nlat = int(180/dlat) + 1
        nlon = int(360/dlon) + 1

        lat = np.linspace(-90., 90., nlat)
        lon = np.linspace(-180., 180., nlon)

        cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
        return cellWidth, lon, lat

    cellWidth, lon, lat = make_uniform_cell_size(cell_size)
    orig_dir = os.getcwd()
    os.chdir(grid_temp_dir)
    build_spherical_mesh(cellWidth, lon, lat, out_filename=str(grid_file), earth_radius=target_radius, plot_cellWidth=False)
    os.chdir(orig_dir)
  
    return 

def generate_surface_dem(grid_file: os.PathLike, dem_file: os.PathLike):
    uxgrid = uxr.open_grid(grid_file,latlon=True,use_dual=False)
    new_elev = np.zeros(uxgrid.n_face,dtype=np.float64)
    ds = xr.DataArray(
            data=new_elev,
            dims=["n_face"],
            attrs={"long_name":"elevation of cells"},
            name="elevation"
            ) 
    ds = ds.to_dataset()
    ds.to_netcdf(dem_file) 
    ds.close()
    return
    