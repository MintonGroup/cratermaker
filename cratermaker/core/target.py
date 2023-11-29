import numpy as np
import os
import xarray as xr
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Tuple
from operator import sub
import os
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
from .material import Material
from ..utils.general_utils import set_properties, create_catalogue, check_properties, float_like

@dataclass
class Target:
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.

    Attributes
    ----------
    name : str or None
        Name of the target body.
    radius : float_like or None
        Radius of the target body in meters.
    diameter : float_like or None
        Diameter of the target body in meters.
    gravity : float_like or None
        Surface gravity of the target body in m/s^2.
    material_name : str or None
        Name of the material composition of the target body.
    material : Material or None
        Material composition of the target body.
    mean_impact_velocity : float_like or None
        Mean impact velocity in m/s.
    pix : float_like or None
        Pixel resolution for the mesh.
    transition_scale_type : str or None
        Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
    cachedir : str, Path, or os.PathLike or None
        Directory path for caching data.
    ds_file : str, Path, or os.PathLike or None
        File path for the target dataset file.
    catalogue : dict or None
        Optional input of catalogue solar system targets to replace the built-in catalogue.
    ds : xarray.Dataset
        xarray Dataset representing the surface mesh and associated data.
    """
       
    # Set up instance variables
    name: str | None = None
    radius: float_like | None = None
    diameter: float_like | None = None
    gravity: float_like | None = None
    material_name: str | None = None
    material: Material | None = None
    mean_impact_velocity: float_like | None = None
    pix: float_like | None = None
    transition_scale_type: str | None = None
    cachedir: str | Path | os.PathLike | None = None
    ds_file: str | Path | os.PathLike | None = None
    catalogue: dict | None = None
    ds: xr.Dataset = field(default_factory=xr.Dataset)
    
    def __getitem__(self,key):
        return self.ds[key]
    
    def __setitem__(self, key, value):
        self.ds[key] = value
        return
    
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        """
        Initialize the target object, setting properties from the provided arguments,
        and creating a catalogue of known solar system targets if not provided.
        """    
        
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = np.float64(9.80665) # 1 g in SI units
        
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity", "transition_scale_type"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", 41100.0, "silicate"),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", 29100.0, "silicate"),
            ("Earth",   6371.01e3, 1.000 * gEarth, "Wet Soil" , 24600.0, "silicate"),
            ("Moon",    1737.53e3, 0.1657* gEarth, "Soft Rock", 22100.0, "silicate"),
            ("Mars",    3389.92e3, 0.379 * gEarth, "Soft Rock", 10700.0, "silicate"),
            ("Ceres",   469.7e3,   0.029 * gEarth, "Ice"      , 5300.0,  "ice"),
            ("Vesta",   262.7e3,   0.025 * gEarth, "Soft Rock", 5300.0,  "silicate"),
        ]      
        # Mean velocities for terrestrial planets based on analysis of simulations from Minton & Malhotra (2010) of main belt-derived asteroid
        # Mean velocities for the asteroids are from Bottke et al. (1994)
       
        if self.catalogue is None: 
            self.catalogue = create_catalogue(body_properties, body_values)

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
        
        # ensure that only either diamter of radius is passed
        values_set = sum(x is not None for x in [self.diameter, self.radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius may be set")
        elif values_set == 1:
            # Be sure to perform the conversion here before the catalogue gets evaluated, in case of potential overrides (e.g. passing diameter as an argument to override a catalogue radius value)
            if self.diameter is not None:
                self.diameter = np.float64(self.diameter)
                self.radius = self.diameter / 2
            elif self.radius is not None:
                self.radius = np.float64(self.radius)
                self.diameter = self.radius * 2 

        # Set properties for the Target object based on the arguments passed to the function
        self.set_properties(**asdict(self))        
        self.material = Material(name=self.material_name)
        
        # Check to make sure diameter and radius conversion happens when catalogue values are used
        if self.diameter is not None:
            self.diameter = np.float64(self.diameter)
            self.radius = self.diameter / 2
        elif self.radius is not None:
            self.radius = np.float64(self.radius)
            self.diameter = self.radius * 2       
            
        if self.radius is not None:
            self.radius = np.float64(self.radius)
        if self.gravity is not None:
            self.gravity = np.float64(self.gravity)
        if self.mean_impact_velocity is not None:
            self.mean_impact_velocity = np.float64(self.mean_impact_velocity)
        if self.pix is not None:
            self.pix = np.float64(self.pix)
        else:    
            self.pix = np.sqrt(4 * np.pi * self.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid
            
        if os.path.exists(self.ds_file):
            self.ds = xr.open_dataset(self.ds_file)
            
        valid_transition_scale_types = ["silicate", "ice"]
        if self.transition_scale_type is not None:
            if not isinstance(self.transition_scale_type, str):
                raise ValueError(f"Transition scale type must be a string and one of {valid_transition_scale_types}")
            self.transition_scale_type = self.transition_scale_type.lower()
            if self.transition_scale_type not in valid_transition_scale_types:
                raise ValueError(f"{self.transition_scale_type} is not a valid transition_scale_type. Must be one of {valid_transition_scale_types}")
            
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
        """
        Set elevation data for the target's surface mesh.

        Parameters
        ----------
        new_elev : np.ndarray, optional
            New elevation data to be set. If None, the elevation is set to zero.
        """        
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
        """
        Calculate the escape velocity for the target body.

        Returns
        -------
        np.float64
            Escape velocity in m/s.
        """        
        return np.sqrt(2 * self.radius * self.gravity)

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
    
    
    def get_cell_distance(self, location: Tuple[np.float64, np.float64]) -> xr.DataArray:
        """
        Computes the distances between cell centers and a given location.

        Parameters
        ----------
        location : Tuple[np.float64, np.float64]
            Tuple containing the longitude and latitude of the location in radians.

        Returns
        -------
        xarray.DataArray
            DataArray of distances for each cell in meters.
        """
        return self.calculate_haversine_distance(location[0],location[1],self.ds.lonCell,self.ds.latCell,self.radius)
    
    
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

    
    def get_cell_initial_bearing(self, location: Tuple[np.float64, np.float64]) -> xr.DataArray:
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
        return self.calculate_initial_bearing(location[0], location[1], self.ds.lonCell, self.ds.latCell)
    

    def get_average_surface(self, location: Tuple[np.float64, np.float64], radius: np.float64) -> Tuple[np.float64, np.float64]:
        """
        Calculate the orientation and radius of the hemispherical cap.

        Parameters
        ----------
        location : Tuple[float, float]
            Tuple containing the longitude and latitude of the reference location in radians.
        radius : float
            The reference radius of to compute the average over in meters.

        Returns
        -------
        cap_center_vector : ndarray
            The vector pointing to the center of the cap from the sphere's center.
        cap_radius : float
            The radius of the cap.
        """
        # Calculate distances from the crater center to each cell
        distances = self.get_cell_distance(location)

        # Find cells within the crater radius
        cells_within_radius = distances <= radius

        # Get the bearings and distances for cells within the crater radius
        bearings = self.ds['bearing'].where(cells_within_radius, drop=True)
        distances = distances.where(cells_within_radius, drop=True)

        # Convert bearings to vector components
        # Bearing is angle from north, positive clockwise, but we need standard mathematical angle, positive counter-clockwise
        angles = np.deg2rad(90) - bearings  # Convert bearing to angle in radians
        x_components = np.cos(angles) * distances
        y_components = np.sin(angles) * distances

        # Calculate the weighted average vector components
        # Weight by the area of each cell to give more importance to larger cells
        cell_areas = self.ds['areaCell'].where(cells_within_radius, drop=True)
        weighted_x = (x_components * cell_areas).sum() / cell_areas.sum()
        weighted_y = (y_components * cell_areas).sum() / cell_areas.sum()

        # Calculate the weighted mean elevation to get the z-component
        elevation_values = self.ds['elevation'].where(cells_within_radius, drop=True)
        weighted_z = (elevation_values * cell_areas).sum() / cell_areas.sum()

        # Combine components to form the cap center vector
        cap_center_vector = np.array([weighted_x.item(), weighted_y.item(), weighted_z.item()])

        # The radius of the cap is the length of the cap center vector
        cap_radius = np.linalg.norm(cap_center_vector)

        return cap_center_vector, cap_radius

# Example usage:
#cap_center_vector, cap_radius = sim.target.calculate_cap_orientation(sim.crater.location, sim.crater.final_radius)
    

    