import xarray as xr
from xarray import DataArray, Dataset
import uxarray as uxr
from uxarray import UxDataArray, UxDataset
from glob import glob
import os
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import shutil
import tempfile
from abc import ABC, abstractmethod
from typing import Tuple, List, Literal, get_args, Any, Union
from typing_extensions import Type
import hashlib
from numpy.typing import NDArray, ArrayLike
from numpy.random import Generator
from .target import Target
from ..utils.general_utils import validate_and_convert_location
from ..utils.custom_types import FloatLike, PairOfFloats
from ..utils.montecarlo import get_random_location_on_face
import warnings

# Define valid grid types
GridType = Literal["icosphere","arbitrary", "hires local"]

# Derive valid_grid_types list from GridType
valid_grid_types: Type[List[str]] = list(get_args(GridType))

# Default file names and directories
_DATA_DIR = "surface_data"
_COMBINED_DATA_FILE_NAME = "surf.nc"
_GRID_FILE_NAME = "grid.nc"

# This is a factor used to determine the smallest length scale in the grid
_SMALLFAC = 1.0e-5

class GridStrategy(ABC):
    def __init__(self, **kwargs: Any):
        self._grid = None

    @abstractmethod
    def generate_face_distribution(self) -> Tuple[NDArray,NDArray,NDArray]:
        pass
    

    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> Tuple[os.PathLike, os.PathLike]:                       
        """
        Generate a tessellated mesh of a sphere of evenly distributed points


        Parameters
        ----------
        grid_file : os.PathLike
            The file path to the grid file.
        grid_hash : str, optional
            Hash of the grid parameters. Default is None, which will generate a new hash.
            
        Notes
        -----
        The grid configuration is determined by the `grid_strategy` attribute of the Surface object. The `grid_strategy` attribute
        determines the type of grid to be generated and its associated parameters. For detailed information on the parameters specific
        to each grid type, refer to the documentation of the respective grid parameter classes (`UnifromIcosphereGrid`, 
        `ArbitraryResolutionGrid`, `HiResLocalGrid`, etc.).
        
        """       

        points = self.generate_face_distribution() 
        grid = uxr.Grid.from_points(points, method="spherical_voronoi")
        if not grid_hash:
            grid_hash = self.generate_hash() 
        grid.attrs["grid_hash"] = grid_hash
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            grid.to_xarray().to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())
            
        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name,grid_file)            
        print("Mesh generation complete")
        self.grid = grid 
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


    def check_if_regrid(self,
                         grid_file: os.PathLike,
                         **kwargs: Any,
                        ) -> bool:
        """
        Check if the existing grid matches the desired parameters determine if regridding is necessary.

        This function checks if a grid file exists and matches the specified parameters based on a unique hash generated from these 
        parameters. If the grid does not exist or does not match the parameters it returns True. 

        Parameters
        ----------
        grid_file : PathLike
            The file path where the grid is saved or will be saved. 
        grid_file : os.PathLike
            The file path to the grid file. 

        Returns
        -------
        bool
            A boolean indicating whether the grid should be regenerated. 
        """
    
        # Generate the hash for the current parameters
        grid_hash = self.generate_hash()

        # Find out if the file exists, if it does't we'll need to make a new grid
        make_new_grid = not os.path.exists(grid_file)
        
        if not make_new_grid:
            uxgrid = uxr.open_grid(grid_file)
            try: 
                old_hash = uxgrid.attrs.get("grid_hash")
                make_new_grid = old_hash != grid_hash
            except:
                make_new_grid = True
                
        return make_new_grid

    
    def create_grid(self,
                     grid_file: os.PathLike,
                     **kwargs: Any,
                     ) -> bool:
        """

        Creates a new grid file based on the grid parameters and stores the new grid as the grid property of the object. 

        Parameters
        ----------
        grid_file : PathLike
            The file path where the grid will be saved. 
        grid_file : os.PathLike
            The file path to the grid file. 
        """
    
        # Generate the hash for the current parameters
        grid_hash = self.generate_hash()

        if os.path.exists(grid_file):
            os.remove(grid_file)
        self.generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs) 
        
        # Check to make sure we can open the grid file, then store the hash in the metadata
        uxgrid = uxr.open_grid(grid_file)
        new_hash = uxgrid.attrs.get("grid_hash")
        assert(new_hash == grid_hash)

        return 
    
    @property
    def pix(self):
        """
        The approximate face size for the mesh inside the local region in meters.
        """
        return self._pix
    
    @pix.setter
    def pix(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("pix must be a positive float")
        self._pix = value
        
    @property
    def radius(self):
        """
        The radius of the target body in meters.
        """
        return self._radius
    
    @radius.setter
    def radius(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("radius must be a positive float")
        self._radius = value
        
    @property
    def grid(self):
        """
        The grid object.
        """
        return self._grid
    
    @grid.setter
    def grid(self, value: uxr.Grid):
        if not isinstance(value, uxr.Grid):
            raise TypeError("grid must be an instance of uxarray.Grid")
        self._grid = value
    
class IcosphereGrid(GridStrategy):    
    """
    Create a uniform grid configuration using an icosphere. This is the most accurate and efficient way to create a uniform grid, but is limited to a few resolutions.
    
    Parameters
    ----------
    gridlevel : float
        The subdivision level of the icosphere. The number of faces is 20 * 4**level. The default level is 8.
    radius: FloatLike
        The radius of the target body in meters.
        
    Returns
    -------
    IcosphereGrid
        An instance of the IcosphereGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 gridlevel: int = 8, 
                 radius: FloatLike = 1.0, 
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.gridlevel = gridlevel
        self.radius = radius
        
        
    def generate_face_distribution(self) -> NDArray:
        """
        Creates the points that define the mesh centers.
           
        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """ 
        from trimesh.creation import icosphere
       
        print(f"Generating a mesh with icosphere level {self.gridlevel}.")  
        mesh = icosphere(self.gridlevel)
        points = mesh.vertices.T
        return points
   
    
    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> Tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        print(f"Effective pixel size: {pix_mean:.2f} +/- {pix_std:.2f} m")
        return    
    
class ArbitraryResolutionGrid(GridStrategy):
    """
    Create a uniform grid configuration with an arbitrary user-defined pixel size. This will not be as nice as the regular IcosphereGrid, but can be any resolution desired.
    
    Parameters
    ----------
    pix : float
        The approximate face size for the mesh in meters.
    radius: FloatLike
        The radius of the target body in meters.
        
    Returns
    -------
    ArbitraryResolutionGrid
        An instance of the ArbitraryResolutionGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike, 
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.pix = pix
        self.radius = radius
    

    def generate_face_distribution(self) -> NDArray:
        """
        Creates the points that define the mesh centers.
           
        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        
        """                

        print(f"Generating a mesh with uniformly distributed faces of size ~{self.pix} m.")
        points = distribute_points(distance=self.pix/self.radius) 
        points[:,0] = np.array([0,0,1])
        points[:,-1] = np.array([0,0,-1])
        return points
   
    
    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> Tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        print(f"Effective pixel size: {pix_mean:.2f} +/- {pix_std:.2f} m")
        return


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
    HiResLocalGrid
        An instance of the HiResLocalGrid clas initialized with the given point distribution
    """        
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike, 
                 local_radius: FloatLike, 
                 local_location: PairOfFloats,
                 superdomain_scale_factor: FloatLike,
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.pix = pix
        self.radius = radius        
        self.local_radius = local_radius
        self.local_location = local_location
        self.superdomain_scale_factor = superdomain_scale_factor
        
    
    def _generate_variable_size_array(self) -> Tuple[NDArray, NDArray, NDArray]:
        """
        Create an array of target pixel sizes for pairs of longitude and latitude values for a high-resolution local mesh around a 
        given location, with adaptive cell sizes outside the local region.
        
            
        Returns
        -------
        pix_array: ndarray
            m x n array of pixel sizes
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            latitude in degrees (length m and between -90 and 90)
        """
        
        def _pix_func(lon,lat):
            lon_rad = np.radians(lon)
            lat_rad = np.radians(lat)
            loc_lon_rad = 0.0 #np.radians(self.local_location[0])
            loc_lat_rad = 0.0 #np.radians(self.local_location[1])

            # Calculate distance from the location to the grid point
            distance = Surface.calculate_haversine_distance(loc_lon_rad, loc_lat_rad, lon_rad, lat_rad, self.radius)
            ans = np.where(distance <= self.local_radius, self.pix, (distance - self.local_radius) / self.superdomain_scale_factor + self.pix)
            return ans
        
        from scipy.interpolate import interp1d

        # Suppose we know pix(lat, lon)
        # Step 1: Construct a fine preliminary grid to estimate integrals
       
        Lat = np.linspace(-90., 90., 733)
        Lon = np.linspace(-180., 180., 733)
        LAT, LON = np.meshgrid(Lat, Lon, indexing='ij')

        pix_values = _pix_func(LON, LAT)   # Evaluate pix on this fine grid
        w = 1.0 / pix_values

        # Step 2: Integrate w over longitude to get W_lat(lat)
        W_lat_vals = np.trapezoid(w, x=Lon, axis=1)  # integrate along lon dimension
        W_lat_cumulative = np.cumsum(W_lat_vals)
        W_lat_cumulative -= W_lat_cumulative[0] # normalize from 0 to 1
        W_lat_cumulative /= W_lat_cumulative[-1]  

        # Create a function to invert W_lat(lat)
        f_lat = interp1d(W_lat_cumulative, Lat, bounds_error=False, fill_value='extrapolate')
        
        M = int(2*np.pi*self.radius/pix_values.max()) - 1
        while M > 0:
            badval=False

            # Step 3: For each lat interval, choose lon lines similarly
            N = M + 1  # number of lon lines
            lat_lines = f_lat(np.linspace(0, 1, N))
            lon_lines = np.zeros((M, N))

            for i in range(M):
                lat_low, lat_high = lat_lines[i], lat_lines[i+1]
                # Extract w in this lat band
                mask = (LAT >= lat_low) & (LAT <= lat_high)
                w_band = w[mask].reshape(-1, len(Lon))  
                # Integrate this band over lat
                w_band_vals = np.trapezoid(w_band, x=Lat[(Lat>=lat_low)&(Lat<=lat_high)], axis=0)
                W_lon_band_cumulative = np.cumsum(w_band_vals)
                if W_lon_band_cumulative[-1] > 0:
                    W_lon_band_cumulative -= W_lon_band_cumulative[0]
                    W_lon_band_cumulative /= W_lon_band_cumulative[-1]
                    f_lon = interp1d(W_lon_band_cumulative, Lon, bounds_error=False, fill_value='extrapolate')
                    lon_lines[i,:] = f_lon(np.linspace(0, 1, N))
                else:
                    badval=True
                    M = M // 2
                    break
                
            if not badval:
                break
            
        if badval:
            raise ValueError("Could not generate a grid with the given parameters. Please try again with different parameters.")
            
        LAT = np.zeros((M, N))
        LON = np.zeros((M, N))
        
        for i in range(M):
            LAT[i, :] = lat_lines[i]  # every point on this horizontal line has the same lat

        LON = lon_lines
            
        pix_array = _pix_func(LON, LAT) 

        return pix_array, LON, LAT

    def _rotate_point_cloud(self,points):
        """
        Rotate a point cloud so that the point at [-r,0,0] moves to (lon,lat)
        using the convention:
        - Longitude [-180, 180] degrees, increasing eastward
        - Latitude [-90, 90] degrees, increasing northward
        - (0,0) lon/lat corresponds to [r,0,0] in Cartesian space
        
        Parameters:
            points (np.ndarray): Nx3 array of (x,y,z) points.
            lon (float): Target longitude in degrees.
            lat (float): Target latitude in degrees.
            r (float): Radius of the sphere (default=1).
        
        Returns:
            np.ndarray: Rotated Nx3 point cloud.
        """

        from scipy.spatial.transform import Rotation as R
        # Convert target lon, lat to radians
        lon_rad, lat_rad = np.radians(self.local_location)
        
        # Compute target unit vector (correcting for lon,lat convention)
        target = np.array([
            self.radius * np.cos(lat_rad) * np.cos(lon_rad),  
            self.radius * np.cos(lat_rad) * np.sin(lon_rad),  
            self.radius * np.sin(lat_rad)                     
        ])

        # Original vector (the point we want to move)
        original = np.array([-self.radius, 0, 0])  # Starts at [-r, 0, 0]
        if np.isclose(lon_rad, 0.0) and np.isclose(lat_rad, 0.0):
            rotation = R.from_euler('z', 180, degrees=True)  # 180-degree rotation around z-axis
            return rotation.apply(points)

        # Compute the axis of rotation (cross product)
        axis = np.cross(original, target)
        axis_norm = np.linalg.norm(axis)
        
        # If the axis is zero (no rotation needed), return the original points
        if axis_norm < 1e-10:
            return points

        axis /= axis_norm  # Normalize the axis

        # Compute the rotation angle (dot product)
        angle = np.arccos(np.clip(np.dot(original / self.radius, target / self.radius), -1.0, 1.0))  # Normalize for dot product

        # Create rotation object
        rotvec = axis * angle  # Convert to rotation vector
        rotation = R.from_rotvec(rotvec)  # Create rotation from axis-angle

        return rotation.apply(points)

    def generate_face_distribution(self) -> NDArray:
        """
        Creates the points that define the mesh centers.
           
        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        
        """                
        
        print(f"Generating a mesh with variable resolution faces")
        print(f"Center of local region: {self.local_location}")
        print(f"Size of local region: {self.local_radius:.2f} m")
        print(f"Hires region pixel size: {self.pix:.2f} m")
        print(f"Lores region pixel size: {self.pix * self.superdomain_scale_factor:.2f} m")
        pix_array, lon, lat = self._generate_variable_size_array()
        
        points = []
        n = lon.shape[0]
        m = lat.shape[1]
       
        for i in range(n-1):
            for j in range(m-1):
                lon_range = (lon[i,j], lon[i,j+1])
                lat_range = (lat[i,j], lat[i+1,j])
                p = distribute_points(distance=pix_array[i,j]/self.radius, lon_range=lon_range, lat_range=lat_range) 
                if p is not None:
                    points.append(p)
                
        points = np.concatenate(points, axis=1)
        points = self._rotate_point_cloud(points.T).T
        
        return points


    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> Tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_min = face_sizes.min().item() * self.radius
        pix_max = face_sizes.max().item() * self.radius
        print(f"Generated {self.grid.n_face} faces")
        print(f"Effective pixel size range: {pix_min:.2f},{pix_max:.2f} m")
        return

    @property
    def local_radius(self):
        """
        The radius of the local region in meters.
        """
        return self._local_radius

    @local_radius.setter
    def local_radius(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("local_radius must be a positive float")
        if value > 2 * np.pi * self.radius:
            raise ValueError("local_radius must be less than 2 * pi * radius of the target body")
        self._local_radius = value
        
    @property
    def local_location(self):
        """
        The longitude and latitude of the location in degrees.
        """
        return self._local_location
    
    @local_location.setter
    def local_location(self, value: PairOfFloats):
        if not isinstance(value, tuple) or len(value) != 2:
            raise TypeError("local_location must be a tuple of two floats")
        self._local_location = validate_and_convert_location(value)
        
    @property
    def superdomain_scale_factor(self):
        """
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters that are 
        modeled outside the local region are those whose ejecta could just reach the boundary.
        """
        return self._superdomain_scale_factor
    
    @superdomain_scale_factor.setter
    def superdomain_scale_factor(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value < 1.0:
            raise TypeError("superdomain_scale_factor must be a positive float greater than or equal to 1")
        self._superdomain_scale_factor = value
   
    
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
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.        
    **kwargs
        This is used to pass additional keyword arguments to pass to the ``uxarray.UxDataset`` class.
    """
    __slots__ = UxDataset.__slots__ + ('_name', '_description', '_data_dir', '_grid_file', '_smallest_length', '_area', '_target', '_rng')    

    def __init__(self, 
                 *args, 
                 target: Target | None = None, 
                 data_dir: os.PathLike | None = None,
                 grid_file: os.PathLike | None = None,
                 compute_face_areas: bool = False,
                 rng: Generator | None = None, 
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
        self._smallest_length = None
        self.rng = rng
       
        if compute_face_areas: 
            # Compute face area needed in the non-normalized units for future calculations
            self['face_areas'] = self.uxgrid.face_areas.assign_attrs(units='m^2') * self.target.radius**2
            self.smallest_length = np.sqrt(self['face_areas'].min().item()) * _SMALLFAC     
        
        return   

    @classmethod
    def initialize(cls, 
                   target: Target | None,
                   data_dir: os.PathLike | None = None,
                   grid_file: os.PathLike | None = None,
                   reset_surface: bool = True, 
                   grid_type: GridType = valid_grid_types[0], 
                   rng: Generator | None = None,
                   regrid: bool = False,
                   **kwargs):
        """
        Factory method to create a Surface instance from a grid file.

        Parameters
        ----------
        target : Target, optional
            Target object or name of known body for the simulation. Default is Target("Moon")
        data_dir : os.PathLike
            The directory for data files. Default is set to ``${PWD}/surface_data``.
        grid_file : os.PathLike
            The file path to the grid file. Default is set to be from the current working directory, to ``{data_dir}/grid.nc``.
        reset_surface : bool, optional
            Flag to indicate whether to reset the surface. Default is True.
        grid_type : ["icosphere", "arbitrary", "hires local"], optional
            The type of grid to be generated. Default is "icosphere".            
        rng : Generator, optional
            A random number generator instance. If not provided, the default numpy RNG will be used. 
        **kwargs : dict
            Additional keyword arguments for initializing the Surface instance based on the specific grid_type.

        Returns
        -------
        Surface
            An initialized Surface object.
                   
        Notes
        -----
        The grid configuration is determined by `grid_type` and its associated parameters. For detailed information on the 
        parameters specific to each grid type, refer to the documentation of the respective grid parameter classes (`IcosphereGrid`,`ArbitraryResolutionGrid`, 
        `HiResLocalGrid`, etc.).

        See Also
        --------
        cratermaker.core.surface.IcosphereGrid : Parameters for an icosphere grid configuration.
        cratermaker.core.surface.ArbitraryResolutionGrid : Parameters for a uniform arbitrary resolution grid configuration.
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
        
        gridlevel = kwargs.pop('gridlevel', 8)
         
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
      
        if not grid_file: 
            grid_file = os.path.join(data_dir,_GRID_FILE_NAME)
            
        # Process the grid parameters from the arguments and build the strategy object 
        if grid_type == "icosphere":
            grid_strategy = IcosphereGrid(gridlevel=gridlevel, radius=target.radius)
        if grid_type == "arbitrary":
            grid_strategy = ArbitraryResolutionGrid(pix=pix, radius=target.radius)
        elif grid_type == "hires local":
            grid_strategy = HiResLocalGrid(pix=pix, radius=target.radius, **kwargs)       
   
        # Check if a grid file exists and matches the specified parameters based on a unique hash generated from these parameters. 
        if not regrid: 
            make_new_grid = grid_strategy.check_if_regrid(grid_file=grid_file, **kwargs)
        else:
            make_new_grid = True
        
        if make_new_grid:
            print("Creating a new grid")
            grid_strategy.create_grid(grid_file=grid_file, **kwargs)
        else:
            print("Using existing grid")
        
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
                surf = uxr.open_mfdataset(grid_file, data_file_list, use_dual=False).isel(time=-1)
                surf.uxgrid = uxr.open_grid(grid_file, use_dual=False)
            else:
                surf = uxr.UxDataset()
                surf.uxgrid = uxr.open_grid(grid_file, use_dual=False)
        except:
            raise ValueError("Error loading grid and data files")
        
        surf = cls(surf,
                   uxgrid=surf.uxgrid,
                   source_datasets=surf.source_datasets,
                   target = target,
                   data_dir = data_dir,
                   grid_file = grid_file,
                   rng = rng,
                   compute_face_areas = True,
                  ) 
        
        if reset_surface:
            surf.generate_data(data=0.0,
                               name="ejecta_thickness",
                               long_name="ejecta thickness",
                               units= "m",
                               save_to_file=True
                              )     
            surf.generate_data(data=0.0,
                               name="ray_intensity",
                               long_name="ray intensity value",
                               units= "",
                               save_to_file=True
                              )                         
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
                            rng=self.rng,
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
            ds._rng = self._rng
            ds._smallest_length = self._smallest_length
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         grid_file=self.grid_file,
                         data_dir=self.data_dir,
                         rng=self.rng,
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

        copied._name = self._name
        copied._description = self._description
        copied._data_dir = self._data_dir
        copied._grid_file = self._grid_file
        copied._smallest_length = self._smallest_length
        copied._area = self._area
        copied._target = self._target
        copied._rng = self._rng
        
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
            ds._rng = self._rng
            ds._smallest_length = self._smallest_length
            ds._area = self._area
        else:
            ds = Surface(ds,
                         uxgrid=self.uxgrid,
                         source_datasets=self.source_datasets,
                         target=self.target,
                         grid_file=self.grid_file,
                         data_dir=self.data_dir,
                         compute_face_areas=False,
                         rng=self.rng,
                         )
            ds._smallest_length = self._smallest_length
            ds._area = self._area
        return ds   

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
    def rng(self):
        """
        A random number generator instance.
        
        Returns
        -------
        Generator
        """ 
        return self._rng
    
    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()           
        
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
            node_tree = self.uxgrid.get_ball_tree("nodes", distance_metric="haversine", coordinate_system="spherical")
            node_ind = node_tree.query(coords=coords, k=1, return_distance=False)
        
            face_tree = self.uxgrid.get_ball_tree("face centers",  distance_metric="haversine", coordinate_system="spherical")
            face_ind = face_tree.query(coords=coords, k=1, return_distance=False)
        return node_ind.item(), face_ind.item()

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
        region_elevation = self['face_elevation'].where(faces_within_radius, drop=False) / self.target.radius
        region_surf = self.elevation_to_cartesian(region_faces, region_elevation) 

        x, y, z = region_surf['face_x'], region_surf['face_y'], region_surf['face_z']
        region_vectors = np.vstack((x, y, z)).T

        # Initial guess for the sphere center and radius
        guess_radius = 1.0 + region_elevation.mean().values.item() 
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
   
    def extract_region(self,
                       location: Tuple[FloatLike, FloatLike],
                       region_radius: FloatLike): 
        
        """
        Extract a regional grid based on a given location and radius.
        
        Parameters
        ----------
        location : Tuple[float, float]
            Tuple containing the longitude and latitude of the location in degrees.
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
        try:
            region_grid = self.uxgrid.subset.bounding_circle(center_coord=location, r=region_angle,element="face centers")
        except ValueError:
            return None
        
        region_surf = self.isel(n_face=region_grid._ds["subgrid_face_indices"], n_node=region_grid._ds["subgrid_node_indices"])
        region_surf.uxgrid = region_grid
      
        return region_surf
    
    def get_random_location_on_face(self, 
                                    face_index: int, 
                                    size: int = 1,
                                    **kwargs
                                    ) -> Union[np.float64, Tuple[np.float64, np.float64], ArrayLike]:
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
        
        return get_random_location_on_face(self.uxgrid, face_index, size,**kwargs)


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


def distribute_points(distance: FloatLike,
                      radius: FloatLike=1.0, 
                      lon_range: PairOfFloats=(-180,180), 
                      lat_range: PairOfFloats=(-90,90)):
    """
    Distributes points on a sphere using Deserno's algorithm [1]_.
        
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
    .. [1] Deserno, Markus., 2004. How to generate equidistributed points on the surface of a sphere. https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    
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


def save(surf: Surface, 
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
        time_variables = {"elapsed_time":np.float64(interval_number)}  
    else:
        if not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")
        
    # Variables that we do not want to save as they are computed at runtime     
        
    surf.close()
    
    ds = surf.expand_dims(dim="time").assign_coords({"time":[interval_number]})
    for k, v in time_variables.items():
        ds[k] = xr.DataArray(data=[v], name=k, dims=["time"], coords={"time":[interval_number]})
                
    drop_vars = [k for k in ds.data_vars if k in do_not_save]
    if len(drop_vars) > 0:
        ds = ds.drop_vars(drop_vars)
        
    _save_data(ds, out_dir, interval_number, combine_data_files)

    return




