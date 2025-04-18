import numpy as np
from typing import Any
from ..utils.general_utils import _set_properties, _check_properties, _to_config
from ..utils.custom_types import FloatLike
from astropy.constants import G
from ..components.target_catalogue import get_target_catalogue

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
        transition_scale_type : str or None
            Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
        material_name : str or None
            Name of the material composition of the target body.            
        catalogue_name : str
            Name of the target catalogue to use. Default is "default".
        **kwargs : Any
            Additional keyword argumments that could be set by the user.

        Notes
        -----
        - The `radius` and `diameter` parameters are mutually exclusive. Only one of them should be provided.
        - Parameters set explicitly using keyword arguments will override those drawn from the catalogue.
        """    

        object.__setattr__(self, "_name",       None)
        object.__setattr__(self, "_radius",       None)
        object.__setattr__(self, "_diameter",     None)
        object.__setattr__(self, "_mass",      None)
        object.__setattr__(self, "_transition_scale_type", None)
        object.__setattr__(self, "_material_name", None)
        object.__setattr__(self, "_catalogue_name", None)
        object.__setattr__(self, "_user_defined", set())   # which public props were set by user
        object.__setattr__(self, "_updating",     False)   # guard against recursive updates

        # ensure that only either diamter of radius is passed
        size_values_set = sum(x is not None for x in [diameter, radius])
        if size_values_set > 1:
            raise ValueError("Only one of diameter or radius may be set")

        catalogue = get_target_catalogue(catalogue_name)().get_targets()
        # Set properties for the Target object based on the arguments passed to the function
        _set_properties(self, 
                        name=name, 
                        radius=radius, 
                        diameter=diameter,
                        mass=mass, 
                        material_name=material_name,
                        catalogue_name=catalogue_name,
                        catalogue=catalogue,
                        transition_scale_type=transition_scale_type, 
                        **kwargs
                    )

        _check_properties(self) 
        return
    
    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)

        # If this attribute is not being updated internally, mark it as a user-defined parameter
        if not self._updating:
            # mark that the *public* name was user‐defined
            public_name = name.lstrip("_")
            self._user_defined.add(public_name)
            if name in ("radius", "diameter"):
                object.__setattr__(self, "_updating", True)
                r = self._radius
                d = self._diameter
                if r is not None and d != 2*r:
                    object.__setattr__(self, "_diameter", 2*r)
                    object.__setattr__(self, "_user_defined", self._user_defined - {"diameter"})
                elif d is not None and r != d/2:
                    object.__setattr__(self, "_radius", d/2)
                    object.__setattr__(self, "_user_defined", self._user_defined - {"radius"})
                object.__setattr__(self, "_updating", False)

    @property
    def radius(self) -> np.float64 | None:
        return self._radius

    @radius.setter
    def radius(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Radius must be positive")
        setattr(self, "_radius", np.float64(value))

    @property
    def diameter(self) -> np.float64 | None:
        return self._diameter

    @diameter.setter
    def diameter(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Diameter must be positive")
        setattr(self, "_diameter", np.float64(value))

    @property
    def mass(self) -> np.float64 | None:
        return self._mass

    @mass.setter
    def mass(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Mass must be positive")
        setattr(self, "_mass", np.float64(value))

    def to_config(self, **kwargs: Any) -> dict:
        return _to_config(self)

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


    # The following are computed properties based on radius and mass
    @property
    def escape_velocity(self) -> np.float64:
        """
        Calculate the escape velocity for the target body in SI units.

        Returns
        -------
        np.float64
            Escape velocity in m/s.
        """        
        return np.sqrt(2 * self.radius*1e-3 * self.gravity)
    
    @property
    def gravity(self) -> np.float64 | None:
        """
        Calculate the gravitational acceleration at the surface of the target body in SI units.

        Returns
        -------
        np.float64
            Gravitational acceleration in m/s^2. 
        """
        return 4*np.pi*G.value*(self.radius*1e-3)*self.bulk_density/3

    @property
    def bulk_density(self) -> np.float64 | None:
        """
        Calculate the bulk density of the target body in SI units.

        Returns
        -------
        np.float64
            Bulk density in kg/m^3.
        """
        return self.mass / (4/3 * np.pi * (self.radius*1e-3)**3)  # in kg/m^3

