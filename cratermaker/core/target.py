import numpy as np
from typing import Any
from ..utils.general_utils import _set_properties, _check_properties
from ..utils.custom_types import FloatLike
from astropy.constants import G
from ..plugins.target_catalogue import get_target_catalogue


class Target:
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.
    """
       
    def __init__(self, 
                 name: str, 
                 radius: FloatLike | None = None, 
                 diameter: FloatLike | None = None,
                 mass: FloatLike | None = None, 
                 gravity: FloatLike | None = None, 
                 bulk_density: FloatLike | None = None,
                 transition_scale_type: str | None = None,
                 material_name: str | None = None,      
                 catalogue_name: str = "default",            
                 **kwargs: Any,
                 ):
        """
        Initialize the target object, setting properties from the provided arguments.

        Parameters
        ----------
        name : str or None
            Name of the target body.
        radius : FloatLike or None
            Radius of the target body in km.
        diameter : FloatLike or None
            Diameter of the target body in km.
        mass : FloatLike or None
            Mass of the target body in kg.
        gravity : FloatLike or None
            Surface gravity of the target body in m/s^2.
        bulk_density : FloatLike or None
            Bulk density of the target body in kg/m^3.
        transition_scale_type : str or None
            Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
        material_name : str or None
            Name of the material composition of the target body.            
        catalogue_name : str
            Name of the target catalogue to use. Default is "default".
        **kwargs : Any
            Additional keyword argumments that could be set by the user.
        """    

        object.__setattr__(self, "_name",       None)
        object.__setattr__(self, "_radius",       None)
        object.__setattr__(self, "_diameter",     None)
        object.__setattr__(self, "_mass",      None)
        object.__setattr__(self, "_gravity",      None)
        object.__setattr__(self, "_bulk_density", None)
        object.__setattr__(self, "_transition_scale_type", None)
        object.__setattr__(self, "_material_name", None)
        object.__setattr__(self, "_catalogue_name", None)
        object.__setattr__(self, "_user_defined", set())   # which public props were set by user
        object.__setattr__(self, "_updating",     False)   # guard against recursive updates
        object.__setattr__(self, "_bulk_property_names", ("radius", "diameter", "gravity", "mass", "bulk_density"))

        # ensure that only either diamter of radius is passed
        size_values_set = sum(x is not None for x in [diameter, radius])
        if size_values_set > 1:
            raise ValueError("Only one of diameter or radius may be set")

        bulk_values_set = size_values_set + sum(x is not None for x in [gravity, bulk_density])
        if bulk_values_set > 2:
            raise ValueError("Only two of diameter/radius, gravity, mass, and bulk_density may be set")

        catalogue = get_target_catalogue(catalogue_name).get_targets()
        # Set properties for the Target object based on the arguments passed to the function
        self._set_properties(name=name, 
                            radius=radius, 
                            diameter=diameter,
                            mass=mass, 
                            gravity=gravity, 
                            bulk_density=bulk_density,
                            material_name=material_name,
                            catalogue_name=catalogue_name,
                            catalogue=catalogue,
                            transition_scale_type=transition_scale_type, 
                            **kwargs)

        _check_properties(self) 
        return
    
    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)

        # If this attribute is not being updated internally, mark it as a user-defined parameter
        if not self._updating:
            # mark that the *public* name was user‐defined
            public_name = name.lstrip("_")
            self._user_defined.add(public_name)

            # If this is a bulk property attribute, recompute all the others to keep them consistent
            if name in self._bulk_property_names:
                object.__setattr__(self, "_updating", True)
                self._update_bulk_properties()
                object.__setattr__(self, "_updating", False)

    def _update_bulk_properties(self):
        """Given any two of _radius/_diameter/_gravity/_bulk_density, compute the others."""
        r = self._radius
        d = self._diameter
        Gval = G.value

        # keep diameter and radius in sync
        if r is not None and d != 2*r:
            object.__setattr__(self, "_diameter", 2*r)
            object.__setattr__(self, "_user_defined", self._user_defined - {"diameter"})
        elif d is not None and r != d/2:
            object.__setattr__(self, "_radius", d/2)
            object.__setattr__(self, "_user_defined", self._user_defined - {"radius"})

        # Check to see if we have currently defined enough bulk property values to compute the others. Otherwise, we'll skip this
        # computation for now. 
        provided = {
            "radius":       self._radius is not None,
            "gravity":      self._mass is not None,
            "gravity":      self._gravity is not None,
            "bulk_density": self._bulk_density is not None,
        }
        n = sum(provided.values())
        if n < 2:
            return

        # compute gravity from radius and density
        if provided["radius"] and provided["bulk_density"] and not provided["gravity"]:
            computed = 4*np.pi*Gval*self._radius*self._bulk_density/3
            object.__setattr__(self, "_gravity", computed)
            object.__setattr__(self, "_user_defined",
                               self._user_defined - {"gravity"})

        # compute density from radius and gravity 
        elif provided["radius"] and provided["gravity"] and not provided["bulk_density"]:
            computed = 3*self._gravity/(4*np.pi*Gval*self._radius)
            object.__setattr__(self, "_bulk_density", computed)
            object.__setattr__(self, "_user_defined",
                               self._user_defined - {"bulk_density"})

        # compute radius & diameter from gravity and density 
        elif provided["gravity"] and provided["bulk_density"] and not provided["radius"]:
            computed_r = 3*self._gravity/(4*np.pi*Gval*self._bulk_density)
            object.__setattr__(self, "_radius", computed_r)
            object.__setattr__(self, "_user_defined",
                               self._user_defined - {"radius"})
            # diameter will sync via the first block on next call

    @property
    def radius(self) -> float | None:
        return self._radius

    @radius.setter
    def radius(self, value: float):
        if value is not None and value <= 0:
            raise ValueError("Radius must be positive")
        # assign via __setattr__
        setattr(self, "_radius", float(value))

    @property
    def diameter(self) -> float | None:
        return self._diameter

    @diameter.setter
    def diameter(self, value: float):
        if value is not None and value <= 0:
            raise ValueError("Diameter must be positive")
        setattr(self, "_diameter", float(value))

    @property
    def gravity(self) -> float | None:
        return self._gravity

    @gravity.setter
    def gravity(self, value: float):
        if value is not None and value < 0:
            raise ValueError("Gravity must be non-negative")
        setattr(self, "_gravity", float(value))

    @property
    def bulk_density(self) -> float | None:
        return self._bulk_density

    @bulk_density.setter
    def bulk_density(self, value: float):
        if value is not None and value <= 0:
            raise ValueError("Bulk density must be positive")
        setattr(self, "_bulk_density", float(value))

    def to_config(self) -> dict:
        """
        Only include those parameters the user actually set.
        """
        return {name: getattr(self, name) for name in self._user_defined}

    def _set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        parameter setting is handled by the `util._set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.
        """         
        
        return _set_properties(self,**kwargs)

    
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
    def catalogue_name(self):
        """
        The name of the target catalogue to use 
        
        Returns
        -------
        str 
        """
        return self._catalogue_name

    @catalogue_name.setter
    def catalogue_name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("catalogue_name must be a string or None")
        self._catalogue_name = value
        
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
