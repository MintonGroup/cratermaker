import numpy as np
import os
import xarray as xr
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Tuple, Dict, Sequence, Optional, Any
import os
from ..utils.general_utils import set_properties, create_catalogue, check_properties
from ..utils.custom_types import FloatLike
import warnings

class Material:
    """
    Represents the material properties relevant to the crater simulation.

    This class defines various physical properties of the material involved in the cratering process.

    """
    config_ignore = ['catalogue']  # Instance variables to ignore when saving to file
    def __init__(self,
                 name: str | None = None,
                 K1: FloatLike | None = None,
                 mu: FloatLike | None = None,
                 Ybar: FloatLike | None = None,
                 density: FloatLike | None = None,
                 catalogue: Dict[str, Dict[str, FloatLike]] | None = None,
                 **kwargs: Any,
                 ):
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
            Volumentric density of material, (kg/m^3)
        catalogue : Dict[str, Dict[str, FloatLike]]
            An optional dictionary containing a catalogue of known materials to use for the simulation. The catalogue should be 
            constructed using a nested dictionary, where the first level of keys are the names of the materials, and the second level
            are the corresponding property names (K1, mu, Ybar, density). 
            If not provided, a default catalogue will be used.
        **kwargs : Any
            Additional keyword argumments that could be set by the user.
            
        Notes
        -----
        The material properties defined here include crater scaling relationship values that were used in CTEM from Richardson (2009) [1]_.  These values used are from Holsapple (1993) [2]_ and Kraus et al. (2011) [3]_.
        
        References
        ----------
        .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029
        .. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333–373. https://doi.org/10.1146/annurev.ea.21.050193.002001
        .. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724–738. https://doi.org/10.1016/j.icarus.2011.05.016        

        """ 
        # Set the attributes for this class
        self._name = name
        self._K1 = None
        self._mu = None
        self._Ybar = None
        self._density = None
        self._catalogue = None
       
        self.catalogue = catalogue
        
        # Set properties for the Material object based on the arguments passed to the function
        self.set_properties(name=name,
                            K1=K1,
                            mu=mu,
                            Ybar=Ybar,
                            density=density,
                            catalogue=self.catalogue,
                            **kwargs) 
        
        # Check to make sure all required properties are set 
        check_properties(self)
        
        return    

    @property
    def catalogue(self):
        """
        A nested dictionary containing a catalogue of known materials to use for the simulation. 
        
        Returns
        -------
        Dict[str, Dict[str, FloatLike]] or None
            A catalogue of known materials to use for the simulation.
        """
        return self._catalogue
    
    @catalogue.setter
    def catalogue(self, value):
                      
        if not isinstance(value, dict) and value is not None:
            raise TypeError("catalogue must be a dict or None") 
             
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
       
        if value is None: 
            self._catalogue = create_catalogue(material_properties, material_values)
        else:
            self._catalogue = value    

    @property
    def name(self):
        """
        The name of the material.
        
        Returns
        -------
        str 
            Name of the material.
        """
        return self._name
    
    @name.setter
    def name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("name must be a string or None")
        self._name = value
        
    @property
    def K1(self):
        """
        K1 crater scaling relationship term. 
        
        Returns
        -------
        np.float64 
        
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._K1
    
    @K1.setter
    def K1(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("K1 must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("K1 must be a positive number")
        self._K1 = np.float64(value)
        
    @property
    def mu(self):
        """
        mu crater scaling relationship term.
        
        Returns
        -------
        np.float64 
        
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._mu
    
    @mu.setter
    def mu(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("mu must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("mu must be a positive number")
        self._mu = np.float64(value)
        
    @property
    def Ybar(self):
        """
        The strength of the material in Pa.
        
        Returns
        -------
        np.float64 
            
                    
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._Ybar
    
    @Ybar.setter
    def Ybar(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("Ybar must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("Ybar must be a positive number")
        self._Ybar = np.float64(value)
        
    @property
    def density(self):
        """
            Volumentric density of material in kg/m^3.
        
        Returns
        -------
        np.float64 
        """
        return self._density
    
    @density.setter
    def density(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("density must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("density must be a positive number")
        self._density = np.float64(value)
        
    
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
    

class Target:
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.
    """
       
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __init__(self, 
                 name: str, 
                 radius: FloatLike | None = None, 
                 diameter: FloatLike | None = None, 
                 gravity: FloatLike | None = None, 
                 material_name: str | None = None, 
                 material: Material | None = None,
                 mean_impact_velocity: FloatLike | None = None, 
                 transition_scale_type: str | None = None,
                 catalogue: Dict[str, Dict[str, FloatLike]] | None = None,
                 **kwargs: Any,
                 ):
        """
        Initialize the target object, setting properties from the provided arguments,
        and creating a catalogue of known solar system targets if not provided.

        Parameters
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
        **kwargs : Any
            Additional keyword argumments that could be set by the user.
        """    

        # Set the attributes for this class
        self._name = name
        self._gravity = None
        self._material_name = None
        self._mean_impact_velocity = None
        self._transition_scale_type = None
        self._material = None
        
        # ensure that only either diamter of radius is passed
        values_set = sum(x is not None for x in [diameter, radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius may be set")        
        
        self._radius = None
        self._diameter = None
        self._catalogue = None
       
        self.catalogue = catalogue 
            
        # Set properties for the Target object based on the arguments passed to the function
        self.set_properties(name=name, 
                            radius=radius, 
                            diameter=diameter, 
                            gravity=gravity, 
                            material_name=material_name, 
                            mean_impact_velocity=mean_impact_velocity, 
                            transition_scale_type=transition_scale_type, 
                            catalogue=self.catalogue,
                            **kwargs)
        self.material = Material(name=self.material_name)
        
        return
    
    @property
    def catalogue(self):
        """
        A nested dictionary containing a catalogue of known solar system targets to use for the simulation. 
        
        Returns
        -------
        Dict[str, Dict[str, FloatLike]] or None
            A catalogue of known solar system targets to use for the simulation.
        """
        return self._catalogue

    @catalogue.setter
    def catalogue(self, value):
        
        if not isinstance(value, dict) and value is not None:
            raise TypeError("catalogue must be a dict or None")
        
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
        
        if value is None: 
            self._catalogue = create_catalogue(body_properties, body_values)
        else:
            self._catalogue = value    
    

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

    
    @property
    def name(self):
        """
        The name of the target body.
        
        Returns
        -------
        str 
        """
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("name must be a string or None")
        self._name = value

       
    @property
    def radius(self) -> Optional[float]:
        """
        Radius of the target body in m.
        
        Returns
        -------
        np.float64 
        """        
        return self._radius

    @radius.setter
    def radius(self, value: float):
        if not isinstance(value, float) or value <= 0:
            raise ValueError("Radius must be a positive number.")
        self._radius = np.float64(value)
        self._diameter = np.float64(value) * 2


    @property
    def diameter(self) -> Optional[float]:
        """
        Diameter of the target body in m.
        
        Returns
        -------
        np.float64 
        """                
        return self._diameter

    @diameter.setter
    def diameter(self, value: float):
        if not isinstance(value, float) or value <= 0:
            raise ValueError("Diameter must be a positive number.")
        self._diameter = np.float64(value)
        self._radius = np.float64(value) / 2

        
    @property
    def gravity(self):
        """
        Surface gravity of the target body in m/s^2.
        
        Returns
        -------
        np.float64
        """
        return self._gravity

    @gravity.setter
    def gravity(self, value):
        if value is not None and not isinstance(value, FloatLike):
            raise TypeError("gravity must be a numeric value or None")
        self._gravity = np.float64(value)


    @property
    def material_name(self):
        """
        The name of the material composition of the target body.
        
        Returns
        -------
        str 
        """
        return self._material_name

    @material_name.setter
    def material_name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("material_name must be a string or None")
        self._material_name = value
        
        
    @property
    def material(self):
        """
        The material properties associated with the target body.
        
        Returns
        -------
        Material
        """
        return self._material
    
    @material.setter
    def material(self, value):
        if not isinstance(value, Material) and value is not None:
            raise TypeError("material must be a Material object or None")
        self._material = value


    @property
    def mean_impact_velocity(self):
        """
        Mean impact velocity in m/s^2.
        
        Returns
        -------
        np.float64
        """
        return self._mean_impact_velocity

    @mean_impact_velocity.setter
    def mean_impact_velocity(self, value):
        if value is not None and not isinstance(value, FloatLike):
            raise TypeError("mean_impact_velocity must be a numeric value or None")
        self._mean_impact_velocity = np.float64(value)


    @property
    def transition_scale_type(self):
        """
        The type of simple-to-complex transition scaling used for the surface, either 'silicate' or 'ice'.
        
        Returns
        -------
        str 
        """
        return self._transition_scale_type

    @transition_scale_type.setter
    def transition_scale_type(self, value):
        valid_types = ["silicate", "ice"]
        if value not in valid_types and value is not None:
            raise ValueError(f"Invalid transition_scale_type: {value}. Must be one of {valid_types} or None")
        self._transition_scale_type = value
