import importlib
import pkgutil
from abc import ABC, abstractmethod
import os
import shutil
import tempfile
import numpy as np
import uxarray as uxr
from numpy.typing import NDArray
from typing import Any
import hashlib
from cratermaker.utils.custom_types import FloatLike, PairOfFloats

class GridMaker(ABC):
    def __init__(self, **kwargs: Any):
        object.__setattr__(self, "_user_defined", set())
        self._user_defined.add("gridtype")
        self._grid = None

    def to_config(self) -> dict:
        """
        Only include those parameters the user actually set.
        """
        return {name: getattr(self, name) for name in self._user_defined}
    
    @abstractmethod
    def generate_face_distribution(self) -> tuple[NDArray,NDArray,NDArray]: ...

    def generate_grid(self,
                      grid_file: os.PathLike,
                      grid_hash: str | None = None,
                      **kwargs: Any) -> tuple[os.PathLike, os.PathLike]:                       
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
        The grid configuration is determined by the `gridtype` attribute of the Surface object. The `gridtype` attribute
        determines the type of grid to be generated and its associated parameters. For detailed information on the parameters specific to each grid type, refer to the documentation of the respective grid parameter classes (`UnifromIcosphereGrid`, 
        `ArbitraryResolutionGrid`, `HiResLocalGrid`, etc.).
        
        """       

        points = self.generate_face_distribution(**kwargs) 
        grid = uxr.Grid.from_points(points, method="spherical_voronoi")
        if not grid_hash:
            grid_hash = self.generate_hash(**kwargs) 
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
    
    def generate_hash(self, **kwargs: Any) -> str:
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

    @staticmethod 
    def _distribute_points(distance: FloatLike,
                           radius: FloatLike=1.0, 
                           lon_range: PairOfFloats=(-180,180), 
                           lat_range: PairOfFloats=(-90,90)) -> NDArray:
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

    @property
    def pix(self):
        """
        The approximate face size for a cell of the mesh.
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

    @property
    def gridtype(self):
        """
        The registered name of this scaling gridtype set by the @register_scaling_gridtype decorator.
        """ 
        return self._gridtype

_registry: dict[str, GridMaker] = {}

def register_grid_type(name: str):
    """
    Class decorator to register a grid component under the given key.
    """
    def decorator(cls):
        cls._gridtype = name 
        _registry[name] = cls
        return cls
    return decorator

def available_grid_types() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_grid_type(name: str):
    """Return the component instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_grid_type) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")