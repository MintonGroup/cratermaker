import importlib
import pkgutil
from pathlib import Path
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
from cratermaker.utils.general_utils import _to_config, parameter
from cratermaker.constants import _GRID_FILE_NAME, _DATA_DIR
from cratermaker.core.base import CratermakerBase

class GridMaker(CratermakerBase, ABC):
    def __init__(self, 
                 simdir: str | Path = Path.cwd(),
                 radius: FloatLike = 1.0, 
                 **kwargs: Any):
        super().__init__(simdir=simdir, **kwargs)
        object.__setattr__(self, "_grid", None)
        object.__setattr__(self, "_radius", None)
        object.__setattr__(self, "_grid_file", None)
        self.radius = radius

    
    @abstractmethod
    def generate_face_distribution(self, **kwargs: Any) -> tuple[NDArray,NDArray,NDArray]: ...

    def generate_grid(self, **kwargs: Any) -> tuple[os.PathLike, os.PathLike]:                       
        """
        Generate a tessellated mesh of a sphere of evenly distributed points


        Notes
        -----
        The grid configuration is determined by the `gridtype` attribute of the Surface object. The `gridtype` attribute
        determines the type of grid to be generated and its associated parameters. For detailed information on the parameters specific to each grid type, refer to the documentation of the respective grid parameter classes (`UnifromIcosphereGrid`, 
        `ArbitraryResolutionGrid`, `HiResLocalGrid`, etc.).
        
        """       

        points = self.generate_face_distribution(**kwargs) 
        grid = uxr.Grid.from_points(points, method="spherical_voronoi")
        grid.attrs["_id"] = self._id
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            grid.to_xarray().to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())
            
        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name,self.file)            
        print("Mesh generation complete")
        self.grid = grid 
        return         

    def check_if_regrid(self, **kwargs: Any) -> bool:
        """
        Check if the existing grid matches the desired parameters determine if regridding is necessary.

        This function checks if a grid file exists and matches the specified parameters based on a unique hash generated from these 
        parameters. If the grid does not exist or does not match the parameters it returns True. 

        Returns
        -------
        bool
            A boolean indicating whether the grid should be regenerated. 
        """
    
        # Find out if the file exists, if it does't we'll need to make a new grid
        make_new_grid = not os.path.exists(self.file)
        
        if not make_new_grid:
            uxgrid = uxr.open_grid(self.file)
            try: 
                old_id = uxgrid.attrs.get("_id")
                make_new_grid = old_id != self._id
            except:
                make_new_grid = True
                
        return make_new_grid

    def create_grid(self, **kwargs: Any):
        """

        Creates a new grid file based on the grid parameters and stores the new grid as the grid property of the object. 

        """
    
        # Generate the hash for the current parameters
        self.file.unlink(missing_ok=True)
        self.generate_grid(**kwargs) 
        
        # Check to make sure we can open the grid file and that the hash matches
        regrid = self.check_if_regrid(**kwargs)
        assert(not regrid)

        return 

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._gridtype, self._radius]
    
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

    @parameter
    def radius(self):
        """
        The radius of the target body in meters.
        """
        return self._radius
    
    @radius.setter
    def radius(self, value: FloatLike):
        if not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("radius must be a positive float")
        self._radius = float(value)
        
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

    @parameter
    def gridtype(self):
        """
        The registered name of this scaling gridtype set by the @register_scaling_gridtype decorator.
        """ 
        return self._gridtype
    
    @property
    def file(self):
        """
        The grid file path.
        """
        return self._simdir / _DATA_DIR / _GRID_FILE_NAME


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

def _init_grid(grid: str | GridMaker | None = None, 
               **kwargs: Any) -> GridMaker:
    """
    Initialize a grid object based on the given name or class.
    
    Parameters
    ----------
    grid : str or GridMaker, optional
        The name of the grid type or an instance of a grid class. If None, a default grid type will be used.
    **kwargs: Any
        Additional keyword arguments to pass to the grid constructor.
    
    Returns
    -------
    GridMaker
        An instance of the specified grid type.
    """

    if grid is None:
        grid = "icosphere"
    if isinstance(grid, str):
        if grid not in available_grid_types():
            raise KeyError(f"Unknown grid model: {grid}. Available models: {available_grid_types()}")
        return get_grid_type(grid)(**kwargs)
    elif isinstance(grid, type) and issubclass(grid, GridMaker):
        return grid(**kwargs)
    elif isinstance(grid, GridMaker):
        return grid
    else:
        raise TypeError(f"grid must be a string or a subclass of Gridmaker, not {type(grid)}")
