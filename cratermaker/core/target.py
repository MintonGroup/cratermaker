import numpy as np
import os
import xarray as xr
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Tuple
import os
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
    transition_scale_type : str or None
        Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
    catalogue : dict or None
        Optional input of catalogue solar system targets to replace the built-in catalogue.
    """
       
    # Set up instance variables
    name: str | None = None
    radius: float_like | None = None
    diameter: float_like | None = None
    gravity: float_like | None = None
    material_name: str | None = None
    material: Material | None = None
    mean_impact_velocity: float_like | None = None
    transition_scale_type: str | None = None
    catalogue: dict | None = None
    
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
