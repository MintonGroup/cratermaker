from __future__ import annotations
from pathlib import Path
from abc import abstractmethod
import os
import shutil
import tempfile
import numpy as np
import uxarray as uxr
import xarray as xr
from numpy.typing import NDArray
from typing import Any
import hashlib
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import  parameter
from cratermaker.constants import _GRID_FILE_NAME, _DATA_DIR
from cratermaker.core.target import Target
from cratermaker.utils.component_utils import ComponentBase, import_components

class Grid(ComponentBase):
    _registry: dict[str, type[Grid]] = {}

    def __init__(self, radius: FloatLike = 1.0, regrid: bool = False, **kwargs: Any):
        super().__init__(**kwargs)
        object.__setattr__(self, "_grid", None)
        object.__setattr__(self, "_radius", None)
        object.__setattr__(self, "_file", None)
        object.__setattr__(self, "_regrid", regrid)
        self.radius = radius


    @classmethod
    def make(cls, 
             grid : str | type[Grid] | Grid| None = None, 
             target: Target | None = None,
             radius: FloatLike = 1.0, 
             simdir: str | Path = Path.cwd(),
             regrid: bool = False,
             **kwargs: Any) -> Grid:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        grid : str or Grid or None
            The name of the grid construction model to use, or an instance of Grid. If None, the default "icosphere" is used.
        target : Target, optional
            The target object. If this is passed, the radius of the target will be used as the radius of the grid. 
        radius : float, optional
            The radius of the body to make the mesh grid with in meters. Default is 1.0. Ignored if target is passed.
        regrid : bool, optional
            If True, the grid will be regridded even if a grid file already exists. Default is False.
        kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        Grid
            An instance of the specified Grid model.

        Raises
        ------
        KeyError
            If the specified morphology model name is not found in the registry.
        TypeError
            If the specified morphology model is not a string or a subclass of Grid.
        """

        # Call the base class version of make and pass the morphology argument as the component argument
        if target is not None and isinstance(target, Target):
            radius = target.radius or radius

        # Verify directory structure exists and create it if not
        data_dir = simdir / _DATA_DIR

        data_dir.mkdir(parents=True, exist_ok=True)

        grid_file = data_dir / _GRID_FILE_NAME 
        regrid = regrid or not grid_file.exists()

        # Create the grid making object
        if grid is None:
            grid = "icosphere"
        grid = super().make(component=grid, simdir=simdir, radius=radius, regrid=regrid, **kwargs)

        # Check if a grid file exists and matches the specified parameters based on a unique hash generated from these parameters. 
        if not regrid: 
            make_new_grid = grid.check_if_regrid(**kwargs)
        else:
            make_new_grid = True
        
        if make_new_grid:
            print("Creating a new grid")
            grid.create_grid(**kwargs)
            grid._regrid = True
        else:
            print("Using existing grid")

        return grid


    @abstractmethod
    def generate_face_distribution(self, **kwargs: Any) -> tuple[NDArray,NDArray,NDArray]: ...

    def generate_grid(self, **kwargs: Any) -> None: 
        """
        Generate a tessellated mesh of a sphere of evenly distributed points

        Notes
        -----
        The grid configuration is determined by the `name` attribute of the Surface object. The `name` attribute
        determines the type of grid to be generated and its associated parameters. For detailed information on the parameters specific to each grid type, refer to the documentation of the respective grid parameter classes (`UnifromIcosphereGrid`, 
        `ArbitraryResolutionGrid`, `HiResLocalGrid`, etc.).
        
        """       

        points = self.generate_face_distribution(**kwargs) 
        uxgrid = uxr.Grid.from_points(points, method="spherical_voronoi")
        uxgrid.attrs["_id"] = self._id
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            uxgrid.to_xarray().to_netcdf(temp_file.name)
            temp_file.flush()
            os.fsync(temp_file.fileno())
            
        # Replace the original file only if writing succeeded
        shutil.move(temp_file.name,self.file)      
        self.uxgrid = uxgrid
        print("Mesh generation complete")

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
        make_new_grid = not Path(self.file).exists()
        
        if not make_new_grid:
            with xr.open_dataset(self.file) as ds:
                uxgrid = uxr.Grid.from_dataset(ds)
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
    def uxgrid(self):
        """
        The grid object.
        """
        return self._uxgrid
    
    @uxgrid.setter
    def uxgrid(self, value: uxr.Grid):
        if not isinstance(value, uxr.Grid):
            raise TypeError("grid must be an instance of uxarray.Grid")
        self._uxgrid = value

    @property
    def file(self):
        """
        The grid file path.
        """
        return self._simdir / _DATA_DIR / _GRID_FILE_NAME

    @property
    def regrid(self):
        """
        Whether the grid was rebuilt when the class was constructed.
        """
        return self._regrid


import_components(__name__, __path__, ignore_private=True)

