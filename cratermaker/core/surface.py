import xarray as xr
import uxarray as uxr
from uxarray import UxDataset
from glob import glob
import os
import numpy as np
import shutil
import tempfile
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
from mpas_tools.viz.paraview_extractor import extract_vtk
import logging
from ..utils.general_utils import float_like
from .target import Target
from typing import Tuple, List
from numpy.typing import NDArray

# Default file names and directories
_DATA_DIR = "surface_data"
_GRID_FILE_NAME = "grid.nc"
_NODE_ELEVATION_FILE_NAME = "elevation_node.nc"
_FACE_ELEVATION_FILE_NAME = "elevation_face.nc"
_GRID_TEMP_DIR = ".grid"

# Mapping from MPAS to UGRID dimension names
_DIM_MAP = {"n_node": "nVertices", 
            "n_face": "nCells",
            "n_edge": "nEdges",
            }
class Surface(UxDataset):
    """
    Surface class that extends UxDataset for cratermaker project.

    This class is used for handling surface-related data and operations in the 
    cratermaker project. It provides functionalities for setting elevation data, 
    calculating distances and bearings, and other surface-related computations.

    Attributes
    ----------
    grid_temp_dir : str
        Directory for temporary grid files.
    data_dir : str
        Directory for data files.
    grid_file : str
        Path to the grid file.
    node_elevation_file : str
        Path to the node elevation file.
    face_elevation_file : str
        Path to the face elevation file.
    target_name : str
        Name of the target body.
    pix : float_like
        Pixel size or resolution of the grid.
    grid_type : str
        Type of the grid used.

    Methods
    -------
    set_elevation(new_elev=None)
        Set elevation data for the target's surface mesh.
    calculate_haversine_distance(lon1, lat1, lon2, lat2, radius)
        Calculate the great circle distance between two points on a sphere.
    get_cell_distance(location, target_radius)
        Computes the distances between cell centers and a given location.
    calculate_initial_bearing(lon1, lat1, lon2, lat2)
        Calculate the initial bearing from one point to another on the surface of a sphere.
    get_cell_initial_bearing(location)
        Computes the initial bearing between cell centers and a given location.
    get_average_surface(location, radius)
        Calculate the orientation of a hemispherical cap that represents the average surface within a given region.
    """   
    __slots__ = UxDataset.__slots__ + ('_name', '_description','grid_temp_dir','data_dir','grid_file','node_elevation_file','face_elevation_file','target_name', 'pix', 'grid_type')
    
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
        new_elev : array_like, optional
            New elevation data to be set. If None, the elevation is set to zero. Must be passed if nCells is not
        """
        if new_elev is None:     
            self['elevation_face'] = generate_data(grid_file=self.grid_file,data_file=self.face_elevation_file,data=new_elev, name="elevation",long_name="elevation of faces", isfacedata=True, save_to_file=False)
            self['elevation_node'] = generate_data(grid_file=self.grid_file,data_file=self.node_elevation_file,data=new_elev, name="elevation",long_name="elevation of nodes", isfacedata=False, save_to_file=False)
        elif new_elev.size == self.uxgrid.n_face:
            self['elevation_face'] = generate_data(grid_file=self.grid_file,data_file=self.face_elevation_file,data=new_elev, name="elevation",long_name="elevation of faces", isfacedata=True, save_to_file=False)
        elif new_elev.size == self.uxgrid.n_node:
            self['elevation_node'] = generate_data(grid_file=self.grid_file,data_file=self.node_elevation_file,data=new_elev, name="elevation",long_name="elevation of nodes", isfacedata=False, save_to_file=False)
        else:
            raise ValueError("new_elev must be None or an array with the same size as the number of faces or nodes in the grid")
          
        return 


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
        lon = np.deg2rad(location[0])
        lat = np.deg2rad(location[1])
        return self.calculate_haversine_distance(lon,lat,self.uxgrid.face_lon,self.uxgrid.face_lat,target_radius)
    

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


def initialize_surface(make_new_grid: bool = False,
         reset_surface: bool = True,
         pix: float_like | None = None,
         target: Target | str | None = None,
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
    pix : float_like | None, optional
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
    grid_temp_dir_path = os.path.join(os.getcwd(), _GRID_TEMP_DIR) 
    if not os.path.exists(grid_temp_dir_path):
        os.mkdir(grid_temp_dir_path)
    
    data_dir_path = os.path.join(os.getcwd(), _DATA_DIR)     
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
    
    # Now redo the elevation data files if necessary 
    node_elevation_file_path = os.path.join(data_dir_path,_NODE_ELEVATION_FILE_NAME)
    face_elevation_file_path = os.path.join(data_dir_path,_FACE_ELEVATION_FILE_NAME)
    
    # Load the grid and data files
    data_file_list = glob(os.path.join(data_dir_path, "*.nc"))
    if grid_file_path in data_file_list:
        data_file_list.remove(grid_file_path)
    
    # Generate a new surface if either it is explicitly requested via parameter or a data file doesn't yet exist 
    reset_surface = reset_surface or not os.path.exists(node_elevation_file_path) or not os.path.exists(face_elevation_file_path) or make_new_grid  
    
    # If reset_surface is True, delete all data files except the grid file 
    if reset_surface:
        for f in data_file_list:
            os.remove(f)
        data_file_list = []
        generate_data(grid_file=grid_file_path,
                      data_file=node_elevation_file_path,
                      name="elevation_node",
                      long_name="elevation of nodes",
                      save_to_file = True,
                      isfacedata=False,
                      )
        generate_data(grid_file=grid_file_path,
                      data_file=face_elevation_file_path,
                      name="elevation_face",
                      long_name="elevation of faces",
                      save_to_file = True,
                      isfacedata=True,
                      )        
        

    if node_elevation_file_path not in data_file_list:
        data_file_list.append(node_elevation_file_path)
    if face_elevation_file_path not in data_file_list:
        data_file_list.append(face_elevation_file_path)
        
    # Initialize UxDataset with the loaded data
    try:
        surf = uxr.open_mfdataset(grid_file_path, data_file_list, latlon=True, use_dual=False)
    except:
        raise ValueError("Error loading grid and data files")
    surf = Surface(surf,uxgrid=surf.uxgrid,source_datasets=surf.source_datasets) 
    
    surf.grid_temp_dir = grid_temp_dir_path
    surf.data_dir = data_dir_path
    surf.grid_file = grid_file_path
    surf.node_elevation_file = node_elevation_file_path
    surf.face_elevation_file = face_elevation_file_path
    
    return surf


def _make_uniform_cell_size(cell_size: float_like) -> Tuple[NDArray,NDArray,NDArray]:
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
                pix: float_like, 
                grid_file: os.PathLike,
                grid_temp_dir: os.PathLike)  -> Surface:
    """
    Generate a tessellated mesh of a sphere using the jigsaw-based mesh builder in MPAS-tools.

    This function generates temporary files in the `grid_temp_dir` directory and saves the final mesh to `grid_file`.

    Parameters
    ----------
    target : str or Target
        Name of target body or a Target object
    pix : float_like
        Desired cell size for the mesh.
    grid_file : os.PathLike
        Path where the grid file will be saved.
    grid_temp_dir : os.PathLike
        Path to the directory for storing temporary grid files.

    Returns
    -------
    A cratermaker Surface object with the generated grid as the uxgrid attribute and with an elevation variable set to zero.
    """
    if isinstance(target, str):
        try:
            target = Target(target)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a valid name of a target body")
    
    cellWidth, lon, lat = _make_uniform_cell_size(pix)
    orig_dir = os.getcwd()
    os.chdir(grid_temp_dir)
    # Configure logger to suppress output
    logger = logging.getLogger("mpas_logger")
    file_handler = logging.FileHandler('mesh.log')
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)     

    print("Building grid with jigsaw...")
    try:
        build_spherical_mesh(cellWidth, lon, lat, out_filename=str(grid_file), earth_radius=target.radius, plot_cellWidth=False, logger=logger)
    except:
        print("Error building grid with jigsaw. See mesh.log for details.")
        raise
    os.chdir(orig_dir)
    
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


def generate_data(grid_file: os.PathLike,
                  data_file: os.PathLike,
                  name: str,
                  long_name: str | None = None,
                  data: float_like | NDArray | None = None,
                  isfacedata: bool = True,
                  save_to_file: bool = False) -> None:
    """
    Generate a NetCDF data file using the provided grid file and save it to the specified path.


    Parameters
    ----------
    data_file : os.PathLike
        Path where the grid file will be saved.
    grid_file : os.PathLike
        Path where the grid file can be found.
    name : str
        Name of the data variable.
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
    uxgrid = uxr.open_grid(grid_file,latlon=True,use_dual=False)
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
    uxda = uxr.UxDataArray(
            data=data,
            dims=dims,
            attrs=None if long_name is None else {"long_name": long_name},
            name=name,
            uxgrid=uxgrid
            ) 
    if save_to_file:
        uxds = uxda.to_dataset()
        uxds.rename(dim_map).to_netcdf(data_file) 
        uxds.close()
    return uxda 


def save_surface(surf: Surface, 
                 out_dir: os.PathLike | None = None,
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
    surf.close()
    for var in surf.data_vars:
        dim_map = {k: _DIM_MAP[k] for k in surf[var].dims if k in _DIM_MAP}  # only map dimensions that are in the variable
        surf[var].rename(dim_map).to_netcdf(os.path.join(out_dir, var + ".nc")) 
        
    return


def export_vtk(surf: Surface, 
               out_dir: os.PathLike = "vtk_files",
               *args, **kwargs) -> None:
    """
    Export the surface mesh to a VTK file.

    Parameters
    ----------
    out_dir : str, Default "vtk_files"
        Directory to store the VTK files.
    """
    ignore_time = "time" not in surf.dims
    
    # Save the current state of the surface before exporting 
    save_surface(surf)
    
    # Combine the grid and data into one file
    extract_vtk(
        filename_pattern=os.path.join(surf.data_dir, '*.nc'),
        mesh_filename=surf.grid_file,
        variable_list=['allOnVertices', 'allOnCells'], 
        dimension_list=['maxEdges=','vertexDegree='], 
        combine=True,
        include_mesh_vars=True,
        out_dir=out_dir, 
        ignore_time=ignore_time, 
        *args, **kwargs)
    
    return
