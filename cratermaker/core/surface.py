import xarray as xr
from xarray import DataArray, Dataset
import uxarray as uxr
from uxarray import UxDataArray, UxDataset
from glob import glob
import os
import sys
import numpy as np
import shutil
import tempfile
from typing import Tuple, List
from numpy.typing import NDArray
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
import logging
from .target import Target
from ..utils.custom_types import FloatLike, PairOfFloats
import warnings

# Default file names and directories
_DATA_DIR = "surface_data"
_COMBINED_DATA_FILE_NAME = "surface_data.nc"
_GRID_FILE_NAME = "grid.nc"
_GRID_TEMP_DIR = ".grid"

# Mapping from MPAS to UGRID dimension names
_DIM_MAP = {"n_node": "nVertices", 
            "n_face": "nCells",
            "n_edge": "nEdges",
            }

class Surface(UxDataset):
    """

    This class is used for handling surface-related data and operations in the 
    cratermaker project. It provides methods for setting elevation data, 
    calculating distances and bearings, and other surface-related computations.
    
    The Surface class extends UxDataset for the cratermaker project.

    """   
    __slots__ = UxDataset.__slots__ + ('_name', '_description', '_grid_temp_dir', '_data_dir', '_grid_file', '_target_radius', '_pix', '_grid_type', '_average_region_center', '_average_region_vector')    
    
    """Surface class for cratermaker"""
    def __init__(self, *args, **kwargs):
        self._grid_temp_dir = kwargs.pop('grid_temp_dir', None)
        self._data_dir = kwargs.pop('data_dir', None)
        self._grid_file = kwargs.pop('grid_file', None)
        self._target_radius = kwargs.pop('target_radius', None)
        self._pix = kwargs.pop('pix', None)
        self._grid_type = kwargs.pop('grid_type', None)
        
        # Call the super class constructor with the UxDataset
        super().__init__(*args, **kwargs)
        
        # Additional initialization for Surface
        self._name = "Surface"
        self._description = "Surface class for cratermaker"
        self._average_region_center = np.zeros(3)
        self._average_region_vector = np.array([0.0, 0.0, self.target_radius])
        

    def generate_data(self,
                      name: str,
                      long_name: str | None = None,
                      data: FloatLike | NDArray | None = None,
                      isfacedata: bool = True,
                      save_to_file: bool = False,
                    ) -> None:
        """
        Generate a NetCDF data file using the provided grid file and save it to the specified path.


        Parameters
        ----------
        name : str
            Name of the data variable. This will also be used as the data file name.
        long_name : str, optional
            Long name of the data variable that will be saved as an attribute.
        data : scalar or array-like
            Data file to be saved. If data is a scalar, then the data file will be filled with that value. If data is an array, then the data file will be filled with the array values. The data array must have the same size as the number of faces or nodes in the grid.
        isfacedata : bool, optional
            Flag to indicate whether the data is face data or node data. Default is True.
        save_to_file: bool, optional
            Specify whether the data should be saved to a file. Default is False.
        Returns
        -------
        None
        """    
        uxgrid = uxr.open_grid(self.grid_file,latlon=True,use_dual=False)
        # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
        # The API change does not affect the functionality of the code, so we can safely ignore the warning        
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", FutureWarning)        
            if isfacedata: 
                dims = ["n_face"]
                size = uxgrid.n_face
                dim_map = {'n_face': 'nCells'}
            else:
                dims = ["n_node"]
                size = uxgrid.n_node
                dim_map = {'n_node': 'nVertices'}
        
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
                    attrs=None if long_name is None else {"long_name": long_name},
                    name=name,
                    uxgrid=uxgrid
                    ) 
            if save_to_file:
                data_file = os.path.join(self.data_dir, f"{name}.nc")
                uxda.rename(dim_map).to_netcdf(data_file) 
                uxda.close()
                
            self[name] = uxda
        return 

        
    def set_elevation(self, 
                    new_elev: NDArray[np.float64] | List[FloatLike] | None = None,
                    save_to_file: bool = False, 
                    ) -> None:
        """
        Set elevation data for the target's surface mesh.

        Parameters
        ----------
        new_elev : array_like, optional
            New elevation data to be set. If None, the elevation is set to zero. 
        save_to_file : bool, default False
            If True, save the elevation data to the elevation file.
        """
        if new_elev is None or np.isscalar(new_elev):
            gen_node = True
            gen_face = True 
        elif new_elev.size == self.uxgrid.n_node:
            gen_node = True
            gen_face = False
        elif new_elev.size == self.uxgrid.n_face:
            gen_node = False
            gen_face = True
        else:
            gen_node = False
            gen_face = False
            raise ValueError("new_elev must be None, a scalar, or an array with the same size as the number of nodes in the grid")
    
        if gen_node:
            self.generate_data(data=new_elev, 
                               name="node_elevation",
                               long_name="elevation of nodes",
                               isfacedata=False,
                               save_to_file=save_to_file
                              )   
        if gen_face:
            self.generate_data(data=new_elev, 
                               name="face_elevation",
                               long_name="elevation of faces",
                               isfacedata=True,
                               save_to_file=save_to_file
                              )
        return 


    @staticmethod
    def calculate_haversine_distance(lon1: FloatLike, 
                    lat1: FloatLike, 
                    lon2: FloatLike, 
                    lat2: FloatLike,
                    radius: FloatLike = 1.0) -> np.float64:
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
    

    def get_distance(self, 
                     location: Tuple[np.float64, np.float64]) -> UxDataArray:
        """
        Computes the distances between nodes and faces and a given location.

        Parameters
        ----------
        location : Tuple[np.float64, np.float64]
            Tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        UxDataArray
            DataArray of distances for each node in meters.
        UxDataArray
            DataArray of distances for each face in meters.
        """
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.uxgrid.node_lon)
        node_lat2 = np.deg2rad(self.uxgrid.node_lat)        
        face_lon2 = np.deg2rad(self.uxgrid.face_lon)
        face_lat2 = np.deg2rad(self.uxgrid.face_lat)
        return self.calculate_haversine_distance(lon1,lat1,node_lon2,node_lat2,self.target_radius), self.calculate_haversine_distance(lon1,lat1,face_lon2,face_lat2,self.target_radius) 
    

    @staticmethod
    def calculate_initial_bearing(lon1: FloatLike, 
                                lat1: FloatLike, 
                                lon2: FloatLike, 
                                lat2: FloatLike) -> np.float64:
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

    
    def get_initial_bearing(self, location: Tuple[np.float64, np.float64]) -> UxDataArray:
        """
        Computes the initial bearing between nodes and faces and a given location.

        Parameters
        ----------
        location : Tuple[np.float64, np.float64]
            Tuple containing the longitude and latitude of the location in degrees.

        Returns
        -------
        DataArray
            DataArray of initial bearings for each node in radians.
        DataArray
            DataArray of initial bearings for each face in radians.
        """
        lon1 = np.deg2rad(location[0])
        lat1 = np.deg2rad(location[1])
        node_lon2 = np.deg2rad(self.uxgrid.node_lon)
        node_lat2 = np.deg2rad(self.uxgrid.node_lat)
        face_lon2 = np.deg2rad(self.uxgrid.face_lon)
        face_lat2 = np.deg2rad(self.uxgrid.face_lat)       
        return self.calculate_initial_bearing(lon1,lat1,node_lon2,node_lat2), self.calculate_initial_bearing(lon1,lat1,face_lon2,face_lat2)
    
    
    def find_nearest_index(self,point):
        """
        Find the index of the nearest node and face to a given point.

        This method calculates the Haversine distance from the given point to each face in the grid,
        and returns the index of the face with the minimum distance.

        Parameters
        ----------
        point : tuple
            A tuple containing two elements: (longitude, latitude) in degrees.

        Returns
        -------
        int
            The index of the nearest node in the grid to the given point.
        int
            The index of the nearest face in the grid to the given point.

        Notes
        -----
        The method converts the longitude and latitude values from degrees to radians before
        calculating distances. The Haversine formula is used to compute the distances on the
        surface of a sphere with a radius of 1.0 unit. 
        """        
        lon1 = np.deg2rad(point[0])
        lat1 = np.deg2rad(point[1])
        node_lon2 = np.deg2rad(self.uxgrid.node_lon)
        node_lat2 = np.deg2rad(self.uxgrid.node_lat)          
        face_lon2 = np.deg2rad(self.uxgrid.face_lon)
        face_lat2 = np.deg2rad(self.uxgrid.face_lat)        
        node_distances, face_distances = self.calculate_haversine_distance(lon1,lat1,node_lon2,node_lat2,radius=1.0), self.calculate_haversine_distance(lon1,lat1,face_lon2,face_lat2,radius=1.0)
        return np.argmin(node_distances.data), np.argmin(face_distances.data)


    def get_average_surface(self,
                            location: Tuple[FloatLike, FloatLike], 
                            region_radius: np.float64) -> NDArray[np.float64]:
        """
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.

        Parameters
        ----------
        location : Tuple[float, float]
            Tuple containing the longitude and latitude of the reference location in degrees.
        region_radius : float
            The radius of the region to compute the average over in meters.

        Returns
        -------
        center_vector : ndarray
            The vector pointing to the center of the cap from the sphere's center.
        """

        # Find cells within the crater radius
        if 'face_crater_distance' in self:
            cells_within_radius = self['face_crater_distance'] <= region_radius
        else:
            _, cells_within_radius = self.get_distance(location)
            cells_within_radius = cells_within_radius <= region_radius
       
        vert_vars = ['face_x', 'face_y', 'face_z'] 
        region_mesh = self.uxgrid._ds[vert_vars].where(cells_within_radius, drop=True)
        region_elevation = self['face_elevation'].where(cells_within_radius, drop=True)
        region_surf = elevation_to_cartesian(region_mesh, region_elevation)

        # Fetch x, y, z values of the mesh within the region
        region_delta = region_surf - region_mesh

        # Fetch the areas of the cells within the region
        cell_areas = self['face_areas'].where(cells_within_radius, drop=True)

        # Calculate the weighted average of the x, y, and z coordinates to get the average surface vector
        weighted_x = (region_delta['face_x'] * cell_areas).sum() / cell_areas.sum()
        weighted_y = (region_delta['face_y'] * cell_areas).sum() / cell_areas.sum()
        weighted_z = (region_delta['face_z'] * cell_areas).sum() / cell_areas.sum()
        average_region_elevation = ((region_elevation * cell_areas).sum() / cell_areas.sum()).item()
        average_region_center = np.array([weighted_x.item(), weighted_y.item(), weighted_z.item()])
       
        # This will compute the vector that describes the orientation of the center of the center of the original mesh surface. 
        _, center_face_index, = self.find_nearest_index(location)
        surface_vector = np.array([self.uxgrid.face_x[center_face_index], self.uxgrid.face_y[center_face_index], self.uxgrid.face_z[center_face_index]])
        
        # Normalize to get the unit vector
        surface_vector = surface_vector / np.linalg.norm(surface_vector)
        
        # Sum the target radius and the average elevation to get the end point of the cap vector
        surface_vector *= (self.target_radius + average_region_elevation)
        
        # Now we can get the vector pointing to the center of the cap from the generating sphere center
        average_region_vector = surface_vector - average_region_center 

        self.average_region_center = average_region_center
        self.average_region_vector = average_region_vector

        return


    @property
    def grid_temp_dir(self):
        """
        Directory for temporary grid files.
        """
        return self._grid_temp_dir

    @grid_temp_dir.setter
    def grid_temp_dir(self, value):
        self._grid_temp_dir = value
        if not os.path.exists(self._grid_temp_dir):
            os.makedirs(self._grid_temp_dir)

    @property
    def data_dir(self):
        """
        Directory for data files.
        """
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value):
        self._data_dir = value
        if not os.path.exists(self._data_dir):
            os.makedirs(self._data_dir)

    @property
    def grid_file(self):
        """
        Path to the grid file.
        """
        return self._grid_file

    @grid_file.setter
    def grid_file(self, value):
        self._grid_file = value

    @property
    def target_radius(self):
        """
        Radius of the target body.
        """
        return self._target_radius

    @target_radius.setter
    def target_radius(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("target_radius must be a float or an integer")
        self._target_radius = value

    @property
    def pix(self):
        """
        Approximate pixel size or resolution used to generate the mesh.
        """
        return self._pix

    @pix.setter
    def pix(self, value):
        if not isinstance(value, FloatLike):
            raise TypeError("pix must be of type FloatLike")
        self._pix = value

    @property
    def grid_type(self):
        """
        Type of the grid used.
        """
        return self._grid_type

    @grid_type.setter
    def grid_type(self, value):
        self._grid_type = value
        

    @property
    def average_region_center(self):
        """
        Center of the average region (used as a reference surface for crater morphology calculations).
        """
        return self._average_region_center

    @average_region_center.setter
    def average_region_center(self, value):
        self._average_region_center = value
        
    @property
    def average_region_vector(self):
        """
        Vector pointing from the center of the generating sphere to the center of the average region.
        """
        return self._average_region_vector
    
    @average_region_vector.setter
    def average_region_vector(self, value):
        self._average_region_vector = value
    
        
def initialize_surface(make_new_grid: bool = False,
         reset_surface: bool = True,
         pix: FloatLike | None = None,
         target: Target | str | None = None,
         simdir: os.PathLike | None = None,
         *args, **kwargs) -> Surface:
    """
    Initialize a Surface object with specified parameters and directory structure.

    This function creates necessary directories, generates grid and surface DEM if required,
    and initializes a Surface object with the loaded data.

    Parameters
    ----------
    make_new_grid : bool, default False
        If True, generate a new grid.
    reset_surface : bool, default True
        If True, reset the surface data.
    pix : FloatLike | None, optional
        Pixel size or resolution of the grid.
    target : Target | str | None, optional
        The target body for the surface, either as a Target object or a string name.
    *args
        Variable length argument list for additional parameters.
    **kwargs
        Arbitrary keyword arguments.

    Returns
    -------
    Surface
        An initialized Surface object.

    Raises
    ------
    ValueError
        If the provided target name is invalid.
    TypeError
        If the target is neither a Target instance nor a valid string name.
    """
    
    if simdir is None:
        simdir = os.getcwd()
    if not target:
        target = Target("Moon")
    elif isinstance(target, str):
        try:
            target = Target(target)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a valid name of a target body")
    
    # Verify directory structure exists and create it if not
    grid_temp_dir_path = os.path.join(simdir, _GRID_TEMP_DIR) 
    if not os.path.exists(grid_temp_dir_path):
        os.mkdir(grid_temp_dir_path)
    
    data_dir_path = os.path.join(simdir, _DATA_DIR)     
    if not os.path.exists(data_dir_path):
        os.mkdir(data_dir_path)
        
    grid_file_path = os.path.join(data_dir_path,_GRID_FILE_NAME)
    
    # Check to see if the grid is correct for this particular set of parameters. If not, then delete it and regrid
    make_new_grid = make_new_grid or not os.path.exists(grid_file_path)
    if not make_new_grid:
        uxgrid = uxr.open_grid(grid_file_path) 
        if "pix" not in uxgrid.parsed_attrs or uxgrid.parsed_attrs["pix"] != pix:
            make_new_grid = True
        elif "grid_type" not in uxgrid.parsed_attrs or uxgrid.parsed_attrs["grid_type"] != "uniform": # this will need to be updated when other grid types are added
            make_new_grid = True
    
    if make_new_grid:
        reset_surface = True
        generate_grid(target=target,
                      pix=pix,
                      grid_file=grid_file_path,
                      grid_temp_dir=grid_temp_dir_path)
    
    # Get the names of all data files in the data directory that are not the grid file
    data_file_list = glob(os.path.join(data_dir_path, "*.nc"))
    if grid_file_path in data_file_list:
        data_file_list.remove(grid_file_path)
        
    # Generate a new surface if either it is explicitly requested via parameter or a data file doesn't yet exist 
    reset_surface = reset_surface or make_new_grid or not data_file_list
    
    # If reset_surface is True, delete all data files except the grid file 
    if reset_surface:
        for f in data_file_list:
            os.remove(f)
        data_file_list = []        
    
    # Initialize UxDataset with the loaded data
    try:
        if data_file_list:
            surf = uxr.open_mfdataset(grid_file_path, data_file_list, latlon=True, use_dual=False)
        else:
            surf = uxr.UxDataset()
            surf.uxgrid = uxr.open_grid(grid_file_path, latlon=True, use_dual=False)
    except:
        raise ValueError("Error loading grid and data files")
    surf = Surface(surf,
                   uxgrid=surf.uxgrid,
                   source_datasets=surf.source_datasets,
                   grid_temp_dir = grid_temp_dir_path,
                   data_dir = data_dir_path,
                   grid_file = grid_file_path,
                   target_radius = target.radius ,
                   ) 
    
    if reset_surface:
        surf.set_elevation(0.0,save_to_file=True)
        
    # Compute face area needed future calculations
    # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
    # The API change does not affect the functionality of the code, so we can safely ignore the warning
    # if self._face_areas is not None: it allows for using the cached result
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", FutureWarning)         
        surf['face_areas'] = uxr.UxDataArray(surf.uxgrid.face_areas, dims=('n_face',), name='face_areas', attrs={'long_name': 'area of faces'})
    
    return surf


def _make_uniform_face_size(cell_size: FloatLike) -> Tuple[NDArray,NDArray,NDArray]:
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


def generate_grid(target: Target | str, 
                pix: FloatLike, 
                grid_file: os.PathLike,
                grid_temp_dir: os.PathLike)  -> Surface:
    """
    Generate a tessellated mesh of a sphere using the jigsaw-based mesh builder in MPAS-tools.

    This function generates temporary files in the `grid_temp_dir` directory and saves the final mesh to `grid_file`.

    Parameters
    ----------
    target : str or Target
        Name of target body or a Target object
    pix : FloatLike
        Desired cell size for the mesh.
    grid_file : os.PathLike
        Path where the grid file will be saved.
    grid_temp_dir : os.PathLike
        Path to the directory for storing temporary grid files.

    Returns
    -------
    A cratermaker Surface object with the generated grid as the uxgrid attribute and with an elevation variable set to zero.
    """
    from matplotlib._api.deprecation import MatplotlibDeprecationWarning
    if isinstance(target, str):
        try:
            target = Target(target)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a valid name of a target body")
    
    cellWidth, lon, lat = _make_uniform_face_size(pix)
    orig_dir = os.getcwd()
    os.chdir(grid_temp_dir)
    # Configure logger to suppress output
    logger = logging.getLogger("mpas_logger")
    file_handler = logging.FileHandler('mesh.log')
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)   
    
    # We can't rely on the jigsaw executable being in the PATH, but the executable will be bundled with the cratermaker project, so 
    # we'll look there first.
    
    # Store the original PATH
    original_path = os.environ.get('PATH', '')  

    try:
        # Add the directory containing the jigsaw binary to PATH
        jigsaw_bin_dir = os.path.join(sys.prefix, 'site-packages', 'cratermaker', 'bin')
        os.environ['PATH'] = jigsaw_bin_dir + os.pathsep + original_path
        
        print("Building grid with jigsaw...")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)
            build_spherical_mesh(cellWidth, lon, lat, out_filename=str(grid_file), earth_radius=target.radius,logger=logger)
    except:
        print("Error building grid with jigsaw. See mesh.log for details.")
        raise
    os.chdir(orig_dir)
    print("Done")
    
    # Create the attribute dictionary that will enable the grid to be identified in case it needs to be regridded 
    with xr.open_dataset(grid_file) as ds:
        ds = ds.assign_attrs(pix=pix, grid_type="uniform") 
    
    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    # Write to the temporary file
    ds.to_netcdf(temp_file.name)

    # Replace the original file only if writing succeeded
    shutil.move(temp_file.name,grid_file)    
  
    return 


def save_surface(surf: Surface, 
                 out_dir: os.PathLike | None = None,
                 combine_data_files: bool = False,
                 *args, **kwargs, 
                 ) -> None:
    """
    Save the surface data to the specified directory. Each data variable is saved to a separate NetCDF file.

    Parameters
    ----------
    surface : Surface
        The surface object to be saved. 
    out_dir : str, optional
        Directory to save the surface data. If None, the data is saved to the current working directory.
    """
    if out_dir is None:
        out_dir = surf.data_dir
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)         
        
    surf.close()
    
    # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
    # The API change does not affect the functionality of the code, so we can safely ignore the warning
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", FutureWarning)
    
        if combine_data_files:
            with tempfile.TemporaryDirectory() as temp_dir:
                outpath = os.path.join(temp_dir, _COMBINED_DATA_FILE_NAME)
                dim_map = {k: _DIM_MAP[k] for k in surf.dims if k in _DIM_MAP}
                surf.rename(dim_map).to_netcdf(outpath)
                shutil.move(outpath,os.path.join(out_dir, _COMBINED_DATA_FILE_NAME)) 
        else: 
            with tempfile.TemporaryDirectory() as temp_dir:
                for var in surf.data_vars:
                    dim_map = {k: _DIM_MAP[k] for k in surf[var].dims if k in _DIM_MAP}  # only map dimensions that are in the variable
                    outname = var + ".nc" 
                    outpath =os.path.join(temp_dir, outname)
                    surf[var].rename(dim_map).to_netcdf(outpath)
                    shutil.move(outpath,os.path.join(out_dir, outname))
        
    return


def elevation_to_cartesian(position: Dataset, 
                           elevation: DataArray
                           ) -> Dataset:
    
    vars = list(position.data_vars)
    if len(vars) != 3:
        raise ValueError("Dataset must contain exactly three coordinate variables")
    
    # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
    # The API change does not affect the functionality of the code, so we can safely ignore the warning
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", FutureWarning)    
        dim_var = list(position.dims)[0]

    rvec = np.column_stack((position[vars[0]], position[vars[1]], position[vars[2]]))
    runit = rvec / np.linalg.norm(rvec, axis=1, keepdims=True)
    
    ds_new = Dataset(
                     {
                      vars[0]: ((dim_var,), rvec[:,0] + elevation.values * runit[:,0]),
                      vars[1]: ((dim_var,), rvec[:,1] + elevation.values * runit[:,1]),
                      vars[2]: ((dim_var,), rvec[:,2] + elevation.values * runit[:,2]),
                     }
                    )
    return ds_new

