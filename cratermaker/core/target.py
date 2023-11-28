import numpy as np
import os
import xarray as xr
from pathlib import Path
from dataclasses import dataclass
from typing import Union, Optional
import os
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
from .material import Material
from ..utils.general_utils import set_properties, create_catalogue, check_properties

@dataclass
class Target:
    """
    Represents the target body in the crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.

    Attributes
    ----------
    material : Material
        The material composition of the target.
    size : float
        The size of the target, in relevant units. 
    """
       
    # Set up instance variables
    name: Optional[str]
    radius: Optional[Union[np.float64,float,int]]
    gravity: Optional[Union[np.float64,float,int]]
    material_name: Optional[str]
    material: Optional[Material]
    mean_impact_velocity: Optional[Union[np.float64,float,int]]
    pix: Optional[Union[np.float64,float,int]]
    transition_scale_type: Optional[str] # Options are silicate and ice
    cachedir: Optional[Union[str,Path]] 
    ds_file: Optional[Union[str,Path]] 
    ds: Optional[xr.Dataset]
    
    def __getitem__(self,key):
        return self.ds[key]
    
    def __setitem__(self, key, value):
        self.ds[key] = value
        return
    
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __post_init__(self):

        if self.cachedir is None:
            self.cachedir = Path.cwd() / ".cache" 
        elif not isinstance(self.cachedir, Path):
            self.cachedir = Path(self.cachedir)
            
        if self.ds_file is None:
            self.ds_file = Path.cwd() / "surface_dem.nc"
        elif not isinstance(self.ds_file, Path):
            self.ds_file = Path(self.ds_file)
            
        if not os.path.exists(self.cachedir):
            os.mkdir(self.cachedir)            
        
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = np.float64(9.80665) # 1 g in SI units
        
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", 41100.0),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", 29100.0),
            ("Earth",   6371.01e3, 1.000 * gEarth, "Wet Soil" , 24600.0),
            ("Moon",    1737.53e3, 0.1657* gEarth, "Soft Rock", 22100.0),
            ("Mars",    3389.92e3, 0.379 * gEarth, "Soft Rock", 10700.0),
            ("Ceres",   469.7e3,   0.029 * gEarth, "Ice"      , 5300.0),
            ("Vesta",   262.7e3,   0.025 * gEarth, "Soft Rock", 5300.0),
        ]      
        # Mean velocities for terrestrial planets based on analysis of simulations from Minton & Malhotra (2010) of main belt-derived asteroid
        # Mean velocities for the asteroids are from Bottke et al. (1994)
        
        self.catalogue = create_catalogue(body_properties, body_values)
        
        # Set properties for the Target object based on the arguments passed to the function
        if self.name:
            self.material = "TEMP" 
            self.set_properties(catalogue=self.catalogue, key=self.name)
            self.material = Material(name=self.material_name)
        else: 
            raise ValueError('No target defined!')
        
        
        if self.radius is not None:
            self.radius = np.float64(self.radius)
        if self.gravity is not None:
            self.gravity = np.float64(self.gravity)
        if self.mean_impact_velocity is not None:
            self.mean_impact_velocity = np.float64(self.mean_impact_velocity)
        if self.transition_scale_type is None:
            self.transition_scale_type = "silicate"
        if self.pix is not None:
            self.pix = np.float64(self.pix)
        else:    
            self.pix = np.sqrt(4 * np.pi * self.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid        

        if os.path.exists(self.ds_file):
            self.ds = xr.open_dataset(self.ds_file)
        else:
            self.make_new_surface()
            
        # Check to make sure all required properties are set 
        check_properties(self)
        
        return
    

    def make_new_surface(self):
        """
        Generate a tessellated mesh of a sphere using the jigsawpy library and convert it to xarray format.

        This function sets up Jigsaw mesh files and a mesh body, defines a basic sphere using an ellipsoid mesh model,
        sets mesh options, generates a tessellated mesh, and converts the mesh to an xarray Dataset. 

        Returns
        -------
        None
            The function does not return a value but updates the `self.mesh` attribute with the generated uxarray object.
        """
        def cellWidthVsLatLon(pix):
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
            constantCellWidth = pix * 1e-3 # build_spherical_mesh assumes units of km, so must be converted

            nlat = int(180/dlat) + 1
            nlon = int(360/dlon) + 1

            lat = np.linspace(-90., 90., nlat)
            lon = np.linspace(-180., 180., nlon)

            cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
            return cellWidth, lon, lat
        
        cellWidth, lon, lat = cellWidthVsLatLon(self.pix)
        orig_dir = os.getcwd()
        os.chdir(self.cachedir)
        build_spherical_mesh(cellWidth, lon, lat, out_filename=str(self.ds_file), earth_radius=self.radius, plot_cellWidth=False)
        os.chdir(orig_dir)
    
        self.ds = xr.load_dataset(self.ds_file)
        self.set_elevation()
        self.ds.to_netcdf(self.ds_file)
        
        return 
        
        
    def set_elevation(self,new_elev=None):
        if new_elev is None:
            new_elev = np.zeros(self.ds.nCells.size,dtype=np.float64)
            
        dem = xr.DataArray(
            data=new_elev,
            dims=["nCells"],
            attrs={"long_name":"elevation of cells"}
            )
        self.ds['elevation'] = dem 
        return
    
    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """         
        set_properties(self,**kwargs)
        return

    
    @property
    def escape_velocity(self):
        return np.sqrt(2 * self.radius * self.gravity)    