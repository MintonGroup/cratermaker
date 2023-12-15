import numpy as np
import os
import xarray as xr
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Tuple, Dict, Sequence
import os
from ..utils.general_utils import set_properties, create_catalogue, check_properties
from ..utils.custom_types import FloatLike
import warnings


@dataclass
class Material:
    """
    Represents the material properties relevant to the crater simulation.

    This class defines various physical properties of the material involved in the cratering process.
    
    """

    # Define all valid properties for the Target object
    name: str | None = None
    K1: FloatLike | None = None
    mu: FloatLike | None = None
    Ybar: FloatLike | None = None
    density: FloatLike | None = None
    catalogue: Dict[str, Dict[str, FloatLike]] | None = None

    config_ignore = ['catalogue']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        """
        Initialize the target object, setting properties from the provided arguments,
        and creating a catalogue of known solar system targets if not provided.
        
        Parameters
        ----------
        name : str
            The name of the material. If the material is matched to one that is present in the catalogue, the rest of the properties 
            will be retrieved for it unless specified. If the name is not known from the catalogue, then all other properties must 
            be supplied and in order to build a custom material.
        K1 : FloatLike
            Variable used in crater scaling (see _[1])
        mu : FloatLike
            Variable used in crater scaling (see _[1])
        Ybar : FloatLike
            The strength of the material, (Pa)
        density : FloatLike
            Volumentric density of material, (kg/m^2)
        catalogue : Dict[str, Dict[str, FloatLike]]
            An optional dictionary containing a catalogue of known materials to use for the simulation. The catalogue should be 
            constructed using a nested dictionary, where the first level of keys are the names of the materials, and the second level
            are the corresponding property names (K1, mu, Ybar, density). 
            If not provided, a default catalogue will be used.
            
        Notes
        -----
        The material properties defined here include crater scaling relationship values that were used in CTEM from Richardson (2009) [1]_.
        These values used are from Holsapple (1993) [2]_ and Kraus et al. (2011) [3]_.
        
        References
        ----------
        .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029
        .. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333–373. https://doi.org/10.1146/annurev.ea.21.050193.002001
        .. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724–738. https://doi.org/10.1016/j.icarus.2011.05.016

        """ 
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "density" 
        ]
        material_values = [
            ("Water",     2.30,     0.55,   0.0,        1000.0),
            ("Sand",      0.24,     0.41,   0.0,        1750.0),
            ("Dry Soil",  0.24,     0.41,   0.18e6,     1500.0),
            ("Wet Soil",  0.20,     0.55,   1.14e6,     2000.0),
            ("Soft Rock", 0.20,     0.55,   7.60e6,     2250.0),
            ("Hard Rock", 0.20,     0.55,   18.0e6,     2500.0),
            ("Ice",       15.625,   0.48,   0.0,        900.0), 
        ]        
       
        if self.catalogue is None: 
            self.catalogue = create_catalogue(material_properties, material_values)
        
        # Set properties for the Material object based on the catalogue value)
        self.set_properties(**asdict(self))
        
        # Check to make sure all required properties are set 
        check_properties(self)
        
        # Ensure types are cast correctly
        self.K1 = np.float64(self.K1)
        self.mu = np.float64(self.mu)
        self.Ybar = np.float64(self.Ybar)
        self.density = np.float64(self.density)
        
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


@dataclass
class Target:
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.
    """
       
    # Set up instance variables
    name: str | None = None
    radius: FloatLike | None = None
    diameter: FloatLike | None = None
    gravity: FloatLike | None = None
    material_name: str | None = None
    material: Material | None = None
    mean_impact_velocity: FloatLike | None = None
    transition_scale_type: str | None = None
    catalogue: dict | None = None
    
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        """
        Initialize the target object, setting properties from the provided arguments,
        and creating a catalogue of known solar system targets if not provided.

        Attributes
        ----------
        name : str or None
            Name of the target body.
        radius : FloatLike or None
            Radius of the target body in meters.
        diameter : FloatLike or None
            Diameter of the target body in meters.
        gravity : FloatLike or None
            Surface gravity of the target body in m/s^2.
        material_name : str or None
            Name of the material composition of the target body.
        material : Material or None
            Material composition of the target body.
        mean_impact_velocity : FloatLike or None
            Mean impact velocity in m/s.
        transition_scale_type : str or None
            Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
        catalogue : dict or None
            Optional input of catalogue solar system targets to replace the built-in catalogue.
        
        """    
        
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = np.float64(9.80665) # 1 g in SI units
        
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "transition_scale_type"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", "silicate"),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", "silicate"),
            ("Earth",   6371.01e3, 1.000 * gEarth, "Wet Soil" , "silicate"),
            ("Moon",    1737.53e3, 0.1657* gEarth, "Soft Rock", "silicate"),
            ("Mars",    3389.92e3, 0.379 * gEarth, "Soft Rock", "silicate"),
            ("Ceres",   469.7e3,   0.029 * gEarth, "Ice"      , "ice"),
            ("Vesta",   262.7e3,   0.025 * gEarth, "Soft Rock", "silicate"),
        ]      

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

        valid_transition_scale_types = ["silicate", "ice"]
        if self.transition_scale_type is not None:
            if not isinstance(self.transition_scale_type, str):
                raise ValueError(f"Transition scale type must be a string and one of {valid_transition_scale_types}")
            self.transition_scale_type = self.transition_scale_type.lower()
            if self.transition_scale_type not in valid_transition_scale_types:
                raise ValueError(f"{self.transition_scale_type} is not a valid transition_scale_type. Must be one of {valid_transition_scale_types}")
           
        if self.mean_impact_velocity is not None:
            self.mean_impact_velocity = np.float64(self.mean_impact_velocity)
            
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
    