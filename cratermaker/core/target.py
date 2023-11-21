import numpy as np
from dataclasses import dataclass
from .material import Material
from ..utils.general_utils import set_properties, create_catalogue
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
    name: str = None
    radius: np.float64 = None
    gravity: np.float64 = None
    material_name: str = None
    material: Material = None
    mean_impact_velocity: np.float64 = None
    transition_scale_type: str = "silicate" # Options are silicate and ice
    
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = 9.80665 # 1 g in SI units
        
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", 41100.0),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", 29100.0),
            ("Earth",   6371.01e3, 1.0   * gEarth, "Wet Soil" , 24600.0),
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
        
        return
    
    @property
    def escape_velocity(self):
        return np.sqrt(2 * self.radius * self.gravity)
    
    
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
    