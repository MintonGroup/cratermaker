import xarray as xr
from xarray import DataArray, Dataset
import uxarray as uxr
from uxarray import UxDataArray, UxDataset
from glob import glob
import os
import sys
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import shutil
import tempfile
from typing import Tuple, List, Literal, get_args, Any
from typing_extensions import Type
import hashlib
from numpy.typing import NDArray
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
import logging
from .target import Target
from ..utils.custom_types import FloatLike, PairOfFloats
import warnings
from ..cython.perlin import apply_noise
from abc import ABC, abstractmethod

# Define valid grid types
GridType = Literal["uniform", "hires_local"]

# Derive valid_grid_types list from GridType
valid_grid_types: Type[List[str]] = list(get_args(GridType))

# Default file names and directories
_DATA_DIR = "surface_data"
_COMBINED_DATA_FILE_NAME = "surf.nc"
_GRID_FILE_NAME = "grid.nc"
_GRID_TEMP_DIR = ".grid"

# This is a factor used to determine the smallest length scale in the grid
_SMALLFAC = 1.0e-5

# Mapping from MPAS to UGRID dimension names
_DIM_MAP = {"n_node": "nVertices", 
            "n_face": "nCells",
            "n_edge": "nEdges",
            }
class GridStrategy(ABC):
    @abstractmethod
    def generate_face_distribution(self) -> Tuple[NDArray,NDArray,NDArray]:
        pass
    
    def generate_grid(self, 
                      grid_file: os.PathLike,
                      grid_temp_dir: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> Tuple[os.PathLike, os.PathLike]: 
        """
        Generate a tessellated mesh of a sphere using the jigsaw-based mesh builder in MPAS-tools.

        This function generates temporary files in the `grid_temp_dir` directory and saves the final mesh to `grid_file`.

        Parameters
        ----------
        grid_file : os.PathLike
            The file path to the grid file.
        grid_temp_dir : os.PathLike
            The directory for temporary grid files. 
        grid_hash : str, optional
            Hash of the grid parameters. Default is None, which will generate a new hash.
            
        Notes
        -----
        The grid configuration is determined by the `grid_strategy` attribute of the Surface object. The `grid_strategy` attribute
        determines the type of grid to be generated and its associated parameters. For detailed information on the parameters specific
        to each grid type, refer to the documentation of the respective grid parameter classes (`UniformGrid`, `HiResLocalGrid`, etc.).
        
        """
        
        from matplotlib._api.deprecation import MatplotlibDeprecationWarning
        if not hasattr(self,"radius" ):
            raise ValueError("radius must be set in the grid strategy")
            
        # Verify directory structure exists and create it if not
        orig_dir = os.getcwd()
        cellWidth, lon, lat = self.generate_face_distribution()
            
        try:
            # verify directory structure, and create anything that's missing.
            if not os.path.exists(grid_temp_dir):
                os.mkdir(grid_temp_dir)
            
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
                    build_spherical_mesh(cellWidth, lon, lat, out_filename=str(grid_file), earth_radius=self.radius,logger=logger)
            except:
                print("Error building grid with jigsaw. See mesh.log for details.")
                raise
            os.chdir(orig_dir)
            print("Done")
            
            if not grid_hash:
                grid_hash = self.generate_hash() 
            
            with xr.open_dataset(grid_file) as ds:
                # Create a temporary file
                with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                    # Write to the temporary file
                    ds.assign_attrs({"grid_hash":grid_hash}).to_netcdf(temp_file.name) 
                    temp_file.flush()
                    os.fsync(temp_file.fileno())
    
            # Replace the original file only if writing succeeded
            shutil.move(temp_file.name,grid_file)    
        finally:
            os.chdir(orig_dir)
    
        return 

    def generate_hash(self) -> str:
        """
        Generate a hash of the grid parameters.

        Returns
        -------
        str
            Hash of the grid parameters.
        """
        attribute_pairs = []
        for attr, value in vars(self).items():
            attribute_pairs.append(f"{attr}:{value}")
        combined = ":".join(attribute_pairs) 
        hash_object = hashlib.sha256(combined.encode())
        return hash_object.hexdigest()
    
    def check_and_regrid(self,
                         grid_file: os.PathLike,
                         grid_temp_dir: os.PathLike,
                         **kwargs: Any,
                        ) -> bool:
        """
        Check if the existing grid matches the desired parameters and regrid if necessary.

        This function checks if a grid file exists and matches the specified parameters based on a unique hash generated from these 
        parameters. If the grid does not exist or does not match the parameters, a new grid is generated. The hash of the grid 
        parameters is also returned so that it can be stored in the new grid

        Parameters
        ----------
        grid_file : PathLike
            The file path where the grid is saved or will be saved. 
        grid_file : os.PathLike
            The file path to the grid file. 
        grid_temp_dir : os.PathLike
            The directory for temporary grid files. 

        Returns
        -------
        bool
            A boolean indicating whether a new grid was generated. `True` if regridding was necessary, `False` otherwise.

        """
    
        # Generate the hash for the current parameters
        grid_hash = self.generate_hash()

        # Find out if the file exists, if it does't we'll need to make a new grid
        make_new_grid = not os.path.exists(grid_file)
        
        if not make_new_grid:
            uxgrid = uxr.open_grid(grid_file)
            try: 
                old_hash = uxgrid.parsed_attrs.get("grid_hash")
                make_new_grid = old_hash != grid_hash
            except:
                make_new_grid = True
            
        if make_new_grid:
            if os.path.exists(grid_file):
                os.remove(grid_file)
            self.generate_grid(grid_file=grid_file, grid_temp_dir=grid_temp_dir, grid_hash=grid_hash, **kwargs) 
            
            # Check to make sure we can open the grid file, then store the hash in the metadata
            uxgrid = uxr.open_grid(grid_file)
            new_hash = uxgrid.parsed_attrs.get("grid_hash")
            assert(new_hash == grid_hash)

        return make_new_grid
    
class UniformGrid(GridStrategy):
    """
    Create a uniform grid configuration with the given pixel size.
    
    Parameters
    ----------
    pix : float
        The approximate face size for the mesh in meters.
    radius: FloatLike
        The radius of the target body in meters.
        
    Returns
    -------
    UniformGrid
        An instance of the UniformGrid class initialized with the given pixel size. 
    """    
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike,
                 **kwargs: Any):

        self.pix = pix
        self.radius = radius

    def generate_face_distribution(self) -> Tuple[NDArray,NDArray,NDArray]:
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
        constantCellWidth = self.pix * 1e-3 # build_spherical_mesh assumes units of km, so must be converted

        nlat = int(180/dlat) + 1
        nlon = int(360/dlon) + 1

        lat = np.linspace(-90., 90., nlat)
        lon = np.linspace(-180., 180., nlon)

        cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
        return cellWidth, lon, lat

class HiResLocalGrid(GridStrategy):
    """
    Create a uniform grid configuration with the given pixel size.
    
    Parameters
    ----------
    pix : FloatLike
        The approximate face size inside the local region in meters.
    radius: FloatLike
        The radius of the target body in meters.
    local_radius : FloatLike
        The radius of the local region in meters.
    local_location : PairOfFloats
        The longitude and latitude of the location in degrees.
    superdomain_scale_factor : FloatLike
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters 
        that are modeled outside the local region are those whose ejecta could just reach the boundary.
        
    Returns
    -------
    UniformGrid
        An instance of the UniformGrid class initialized with the given pixel size. 
    """        
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike, 
                 local_radius: FloatLike, 
                 local_location: PairOfFloats,
                 superdomain_scale_factor: FloatLike,
                 **kwargs: Any):
        self.pix = pix
        self.radius = radius
        self.local_radius = local_radius
        self.local_location = local_location
        self.superdomain_scale_factor = superdomain_scale_factor

    def generate_face_distribution(self) -> Tuple[NDArray, NDArray, NDArray]:
        """
        Create a cell width array for a high-resolution local mesh around a given location, with adaptive cell sizes outside the local region.

            
        Returns
        -------
        cellWidth : ndarray
            m x n array of cell width in km
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            latitude in degrees (length m and between -90 and 90)
        """
        dlat = 10
        dlon = 10

        nlat = int(180 / dlat) + 1
        nlon = int(360 / dlon) + 1

        lat = np.linspace(-90., 90., nlat)
        lon = np.linspace(-180., 180., nlon)
        cellWidth = np.zeros((nlat, nlon))

        # Convert location to radians for distance calculation
        loc_lon_rad = np.radians(self.local_location[0])
        loc_lat_rad = np.radians(self.local_location[1])

        for i in range(nlat):
            for j in range(nlon):
                # Convert grid point to radians
                lon_rad = np.radians(lon[j])
                lat_rad = np.radians(lat[i])

                # Calculate distance from the location to the grid point
                distance = Surface.calculate_haversine_distance(loc_lon_rad, loc_lat_rad, lon_rad, lat_rad)

                if distance <= self.local_radius:
                    # Inside the local region
                    cellWidth[i, j] = self.pix * 1e-3  # Convert to km
                else:
                    # Outside the local region, adapt cell size based on distance from boundary
                    boundary_distance = distance - self.local_radius
                    cellWidth[i, j] = (boundary_distance / self.superdomain_scale_factor + self.pix) * 1e-3  # Convert to km

        return cellWidth, lon, lat

class Surface(UxDataset):
    """
    This class is used for handling surface-related data and operations in the cratermaker project. It provides methods for 
    setting elevation data, calculating distances and bearings, and other surface-related computations.
    
    The Surface class extends UxDataset for the cratermaker project.
    
    Parameters
    ----------
    *args
        Variable length argument list for additional parameters to pass to the ``uxarray.UxDataset`` class.
    target : Target, optional
        The target body or name of a known target body for the impact simulation. 
    data_dir : os.PathLike, optional
        The directory for data files.
    grid_file : os.PathLike, optional
        The file path to the grid file.
    compute_face_areas : bool, optional
        Flag to indicate whether to compute face areas. Default is False.    
    **kwargs
        This is used to pass additional keyword arguments to pass to the ``uxarray.UxDataset`` class.
    """
    __slots__ = UxDataset.__slots__ + ('_name', '_description', '_data_dir', '_grid_file', '_smallest_length', '_area', '_target', '_grid_hash')    

    def __init__(self, 
                 *args, 
                 target: Target | None = None, 
                 data_dir: os.PathLike | None = None,
                 grid_file: os.PathLike | None = None,
                 compute_face_areas: bool = False,
                 **kwargs):

        # Call the super class constructor with the UxDataset
        super().__init__(*args, **kwargs)
       
        # Additional initialization for Surface
        self._name = "Surface"
        self._description = "Surface class for cratermaker"
        self._data_dir = data_dir
        self._grid_file = grid_file
        self._target = target
        self._area = None
        self._grid_hash = None
        self._smallest_length = None
       
        if compute_face_areas: 
            # Compute face area needed future calculations
            # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
            # The API change does not affect the functionality of the code, so we can safely ignore the warning
            # if self._face_areas is not None: it allows for using the cached result
            with warnings.catch_warnings(): 
                warnings.simplefilter("ignore", FutureWarning)
                if 'face_areas' not in self:
                    self['face_areas'] = uxr.UxDataArray(self.uxgrid.face_areas * self.target.radius**2, dims=('n_face',), name='face_areas', attrs={'long_name': 'area of faces', 'units': 'm^2'})
                    self.smallest_length = np.sqrt(self['face_areas'].min().item()) * _SMALLFAC        

    @classmethod
    def initialize(cls, 
                   target: Target | None,
                   data_dir: os.PathLike | None = None,
                   grid_file: os.PathLike | None = None,
                   grid_temp_dir: os.PathLike | None = None,
                   reset_surface: bool = True, 
                   grid_type: GridType = valid_grid_types[0], 
                   **kwargs):
        """
        Factory method to create a Surface instance from a grid file.

        Parameters
        ----------
        target : Target, optional
            Target object or name of known body for the simulation. Default is Target("Moon")
        data_dir : os.PathLike
            The directory for data files. Default is set to `${PWD}/surface_data`.
        grid_file : os.PathLike
            The file path to the grid file. Default is set to be from the current workding directory, to `{data_dir}/grid.nc'.
        grid_temp_dir : os.PathLike
            The directory for temporary grid files. Default is set to `${PWD}/.grid`.
        reset_surface : bool, optional
            Flag to indicate whether to reset the surface. Default is True.
        grid_type : ["uniform", "hires_local"], optional
            The type of grid to be generated. Default is "uniform".            
        **kwargs : dict
            Additional keyword arguments for initializing the Surface instance based on the specific grid_type.

        Returns
        -------
        Surface
            An initialized Surface object.
                   
        Notes
        -----
        The grid configuration is determined by `grid_type` and its associated parameters. For detailed information on the 
        parameters specific to each grid type, refer to the documentation of the respective grid parameter classes (`UniformGrid`, 
        `HiResLocalGrid`, etc.).

        See Also
        --------
        cratermaker.core.surface.UniformGrid : Parameters for a uniform grid configuration.
        cratermaker.core.surface.HiResLocalGrid : Parameters for a high-resolution local grid configuration.  
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
        
        if grid_type not in valid_grid_types:
            raise ValueError(f"Invalid grid_type {grid_type}. Valid options are {valid_grid_types}")
         
        pix = kwargs.pop('pix', None) 
        if pix is not None:
            pix = np.float64(pix)
        else:    
            pix = np.sqrt(4 * np.pi * target.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid

        # Verify directory structure exists and create it if not
        if not data_dir:
            data_dir = os.path.join(os.getcwd(), _DATA_DIR)
            
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
            reset_surface = True
      
        if not grid_temp_dir:   
            grid_temp_dir = os.path.join(os.getcwd(), _GRID_TEMP_DIR) 
            
        if not os.path.exists(grid_temp_dir):
            os.mkdir(grid_temp_dir)
       
        if not grid_file: 
            grid_file = os.path.join(data_dir,_GRID_FILE_NAME)
            
        # Process the grid parameters from the arguments and build the strategy object 
        if grid_type == "uniform":
            grid_strategy = UniformGrid(pix=pix, radius=target.radius)
        elif grid_type == "hires_local":
            grid_strategy = HiResLocalGrid(pix=pix, radius=target.radius, **kwargs)       
   
        # Check if a grid file exists and matches the specified parameters based on a unique hash generated from these parameters. 
        make_new_grid = grid_strategy.check_and_regrid(grid_file=grid_file, grid_temp_dir=grid_temp_dir, **kwargs)
        
        # Get the names of all data files in the data directory that are not the grid file
        data_file_list = glob(os.path.join(data_dir, "*.nc"))
        if grid_file in data_file_list:
            data_file_list.remove(grid_file)
            
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
                surf = uxr.open_mfdataset(grid_file, data_file_list, latlon=True, use_dual=False).isel(Time=-1)
                surf.uxgrid = uxr.open_grid(grid_file, latlon=True, use_dual=False)
            else:
                surf = uxr.UxDataset()
                surf.uxgrid = uxr.open_grid(grid_file, latlon=True, use_dual=False)
        except:
            raise ValueError("Error loading grid and data files")
        
        surf = cls(surf,
                   uxgrid=surf.uxgrid,
                   source_datasets=surf.source_datasets,
                   target = target,
                   data_dir = data_dir,
                   grid_file = grid_file,
                   compute_face_areas = True,
                    ) 
        
        if reset_surface:
            surf.set_elevation(0.0,save_to_file=True)
        
        return surf        
        
    def __getitem__(self, key):
        """Override to make sure the result is an instance of ``cratermaker.Surface``"""

        value = super().__getitem__(key)

        if isinstance(value, uxr.UxDataset):
            value = Surface(value,
                            uxgrid=self.uxgrid,
                            source_datasets=self.source_datasets,
                            target=self.target,
                            grid_file=self.grid_file,
                            data_dir=self.data_dir,
                            compute_face_areas=False,
                            )
        return value
    
    def _calculate_binary_op(self, *args, **kwargs):
        """Override to make the result a complete instance of ``cratermaker.Surface``."""
        ds = super()._calculate_binary_op(*args, **kwargs)

        if isinstance(ds, Surface):
            ds._name = self._name
            ds._description = self._description
            ds._data_dir = self._data_dir
            ds._grid_file = self._grid_file
            ds._target = self._target
            ds._smallest_length = self._smallest_length
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         grid_file=self.grid_file,
                         data_dir=self.data_dir,
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

        deep = kwargs.get('deep', None)

        if deep == True:
            # Deep copy all of the properties
            copied._name = self._name.copy()
            copied._description = self._description()
            copied._data_dir = self._data_dir.copy()
            copied._grid_file = self._grid_file.copy()
            copied._target = self._target.copy()
            copied._smallest_length = self._smallest_length.copy()
            copied._area = self._area.copy()
            copied._grid_hash = self._grid_hash.copy()
        else:
            # Point to the existing properties
            copied._name = self._name
            copied._description = self._description
            copied._data_dir = self._data_dir
            copied._grid_file = self._grid_file
            copied._target = self._target
            copied._smallest_length = self._smallest_length
            copied._area = self._area
            copied._grid_hash = self._grid_hash
        return copied    
  
    def _replace(self, *args, **kwargs):
        """Override to make the result a complete instance of ``cratermaker.Surface``."""
        ds = super()._replace(*args, **kwargs)

        if isinstance(ds, Surface):
            ds._name = self._name
            ds._description = self._description
            ds._data_dir = self._data_dir
            ds._grid_file = self._grid_file
            ds._target = self._target
            ds._smallest_length = self._smallest_length
            ds._area = self._area
            ds._grid_hash = self._grid_hash
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         grid_file=self.grid_file,
                         data_dir=self.data_dir,
                         compute_face_areas=False,
                         )
            ds._smallest_length = self._smallest_length
            ds._area = self._area
            ds._grid_hash = self._grid_hash
        return ds   

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
        
    @property
    def area(self):
        """
        Total surface area of the target body.
        """
        if self._area is None:
            self._area = self['face_areas'].sum().assign_attrs({'long_name': 'total surface area', 'units': 'm^2'})
        return self._area

    @property
    def target(self):
        """
        The target body for the impact simulation. Set during initialization.
        """
        return self._target

    @target.setter
    def target(self, value):
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value    
        
    @property
    def grid_hash(self):
        """
        The hash of the grid parameters.
        """
        return self._grid_hash
    
    @grid_hash.setter
    def grid_hash(self, value):
        if not isinstance(value, str):
            raise TypeError("grid_hash must be a string")
        self._grid_hash = value
    
    @property
    def grid_strategy(self):
        """
        The hash of the grid parameters.
        """
        return self._grid_strategy
    
    @grid_strategy.setter
    def grid_strategy(self, value):
        if not isinstance(value, GridStrategy):
            raise TypeError("grid_strategy must be of type GridStrategy")
        self._grid_strategy = value
        
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
        uxgrid = uxr.open_grid(self.grid_file,latlon=True,use_dual=False)
        if long_name is None and units is None:
            attrs = None
        else:
            attrs = {}
            if long_name:
                attrs["long_name"] = long_name
            if units:
                attrs["units"] = units
        # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
        # The API change does not affect the functionality of the code, so we can safely ignore the warning        
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", FutureWarning)        
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
                    uxgrid=uxgrid
                    ) 
            self[name] = uxda
            
            if save_to_file:
                _save_data(uxda, self.data_dir, interval_number, combine_data_files)
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
        return self.calculate_haversine_distance(lon1,lat1,node_lon2,node_lat2,self.target.radius), self.calculate_haversine_distance(lon1,lat1,face_lon2,face_lat2,self.target.radius) 

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

    def get_reference_surface(self,
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
        if 'face_crater_distance' and 'node_crater_distance' in self:
            faces_within_radius = self['face_crater_distance'] <= region_radius
            nodes_within_radius = self['node_crater_distance'] <= region_radius
        else:
            nodes_within_radius, faces_within_radius = self.get_distance(location)
            faces_within_radius = faces_within_radius <= region_radius
            nodes_within_radius = nodes_within_radius <= region_radius
        
        if np.sum(faces_within_radius) < 5 or not nodes_within_radius.any():
            self['reference_face_elevation'] = self['face_elevation']
            self['reference_node_elevation'] = self['node_elevation']
            return
       
        inc_face = self.n_face
        face_vars = ['face_x', 'face_y', 'face_z'] 
        face_grid = self.n_face.uxgrid._ds[face_vars].sel(n_face=inc_face)
        region_faces = face_grid.where(faces_within_radius, drop=False)
        region_elevation = self['face_elevation'].where(faces_within_radius, drop=False)
        region_surf = self.elevation_to_cartesian(region_faces, region_elevation) / self.target.radius

        x, y, z = region_surf['face_x'], region_surf['face_y'], region_surf['face_z']
        region_vectors = np.vstack((x, y, z)).T

        # Initial guess for the sphere center and radius
        guess_radius = 1.0 + region_elevation.mean().values.item() / self.target.radius
        initial_guess = [0, 0, 0, guess_radius]  

        # Perform the curve fitting
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            try:
                bounds =([-1.0, -1.0, -1.0, 0.5],[1.0, 1.0, 1.0, 2.0]) 
                popt, _  = curve_fit(sphere_function, region_vectors[faces_within_radius], np.zeros_like(x[faces_within_radius]), p0=initial_guess, bounds=bounds)
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
        face_vectors  = np.vstack((region_faces.face_x, region_faces.face_y, region_faces.face_z)).T
        self['reference_face_elevation'] = xr.where(faces_within_radius, find_reference_elevations(face_vectors), self['face_elevation'])
        
        # Now do the same thing to compute the nodal values 
        inc_node = self.n_node
        node_vars = ['node_x', 'node_y', 'node_z'] 
        node_grid = self.n_node.uxgrid._ds[node_vars].sel(n_node=inc_node)
        region_nodes = node_grid.where(nodes_within_radius, drop=False) 
        node_vectors  = np.vstack((region_nodes.node_x, region_nodes.node_y, region_nodes.node_z)).T
        
        self['reference_node_elevation'] = xr.where(nodes_within_radius, find_reference_elevations(node_vectors), self['node_elevation'])
        
        return

    @staticmethod
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
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", FutureWarning)
        dim_map = {k: _DIM_MAP[k] for k in ds.dims if k in _DIM_MAP}
        ds = ds.rename(dim_map)
        if isinstance(ds, xr.DataArray):
            ds = ds.to_dataset()
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
                
            if "Time" not in ds_file.dims:
                ds_file = ds_file.expand_dims(["Time"])
            if "Time" not in ds_file.coords:
                ds_file = ds_file.assign_coords({"Time":[interval_number]})                   
                
            temp_file = os.path.join(temp_dir, filename)
            ds_file.to_netcdf(temp_file) 
            ds_file.close()     
            shutil.move(temp_file, data_file)

    return
    

def save(surf: Surface, 
         out_dir: os.PathLike | None = None,
         combine_data_files: bool = False,
         interval_number: int = 0,
         time_variables: dict | None = None,
         *args, **kwargs, 
         ) -> None:
    """
    Save the surface data to the specified directory. Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the Time dimension. If 'interval_number' is included as a key in `time_variables`, then this will be appended to the data file name.

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
        Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the Time dimension. Default is None.
    """
    do_not_save = ["face_areas"]
    if out_dir is None:
        out_dir = surf.data_dir
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)         
      
    if time_variables is None:
        time_variables = {"elapsed_time":np.float64(interval_number)}  
    else:
        if not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")
        
    # Variables that we do not want to save as they are computed at runtime     
        
    surf.close()
    
    # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
    # The API change does not affect the functionality of the code, so we can safely ignore the warning
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore", FutureWarning)
        ds = surf.expand_dims(dim="Time").assign_coords({"Time":[interval_number]})
        for k, v in time_variables.items():
            ds[k] = xr.DataArray(data=[v], name=k, dims=["Time"], coords={"Time":[interval_number]})
                    
        drop_vars = [k for k in ds.data_vars if k in do_not_save]
        if len(drop_vars) > 0:
            ds = ds.drop_vars(drop_vars)
            
        _save_data(ds, out_dir, interval_number, combine_data_files)

        return




