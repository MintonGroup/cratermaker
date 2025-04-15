import numpy as np
from typing import Dict, Optional, Any
from ..utils.general_utils import set_properties, create_catalogue, check_properties, group, ParameterGroups, to_config
from ..utils.custom_types import FloatLike
import inspect 
from astropy.constants import G


class Target:
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.
    """
       
    config_ignore = ['_catalogue','_material']  # Instance variables to ignore when saving to file
    def __init__(self, 
                 name: str, 
                 radius: FloatLike | None = None, 
                 diameter: FloatLike | None = None, 
                 gravity: FloatLike | None = None, 
                 transition_scale_type: str | None = None,
                 material_name: str | None = None,                  
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
        self._transition_scale_type = None
        self._material_name = None
        
        # ensure that only either diamter of radius is passed
        values_set = sum(x is not None for x in [diameter, radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius may be set")        
        
        self._radius = None
        self._diameter = None
        self._bulk_density = None
        self._catalogue = None
        self.catalogue = catalogue 
            
        # Set properties for the Target object based on the arguments passed to the function
        self.set_properties(name=name, 
                            radius=radius, 
                            diameter=diameter, 
                            gravity=gravity, 
                            material_name=material_name,
                            transition_scale_type=transition_scale_type, 
                            catalogue=self.catalogue,
                            **kwargs)

        check_properties(self) 
        return
    
    def to_config(self) -> Dict[str, Any]:
        """
        Serialize the user-defined parameters of the object into a Dict.

        Returns
        -------
        dict
            A dictionary containing the serializable attributes of the object.

        """   
        config = to_config(self, required_counts={"size": 1, "bulk_properties": 2})
        return config

    
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

    @group("size", "bulk_properties")
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
        if self._gravity is None and self._bulk_density is not None:
            self.gravity = 4 * np.pi * G.value / 3.0 * self._bulk_density * self._radius
        elif self._bulk_density is None and self._gravity is not None:
            self.bulk_density = 3 * self._gravity / (4 * np.pi * G.value) * self._radius        

    @group("size", "bulk_properties")
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
        self.radius = np.float64(value) / 2
        if self._bulk_density is not None:
            self._gravity = 4 * np.pi * G.value / 3.0 * self._bulk_density * self._radius
        elif self._gravity is not None:
            self._bulk_density = 3 * self._gravity / (4 * np.pi * G.value) * self._radius

    @group("bulk_properties")
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
        if self._radius is not None:
            self._bulk_density = 3 * self._gravity / (4 * np.pi * G.value) * self._radius
        elif self._bulk_density is not None:
            self._radius = 3 * self._gravity / (4 * np.pi * G.value) * self._bulk_density

    @group("bulk_properties")
    def bulk_density(self) -> Optional[float]:
        """
        Bulk density of the target body in kg/m^3.
        
        Returns
        -------
        np.float64 
        """        
        return self._bulk_density

    @bulk_density.setter
    def bulk_density(self, value: float):
        if not isinstance(value, float) or value <= 0:
            raise ValueError("Bulk density must be a positive number.")
        self._bulk_density = np.float64(value)
        if self._radius is not None:
            self._gravity = 4 * np.pi * G.value / 3.0 * self._bulk_density * self._radius
        elif self._gravity is not None:
            self._radius = 3 * self._gravity / (4 * np.pi * G.value) * self._bulk_density

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
