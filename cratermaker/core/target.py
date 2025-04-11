import numpy as np
from typing import Any
from ..utils.general_utils import _set_properties, _to_config, _create_catalogue
from ..utils.custom_types import FloatLike
from astropy.constants import G

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
        object.__setattr__(self, "_catalogue",  None)
        object.__setattr__(self, "_user_defined", set())   # which public props were set by user
        object.__setattr__(self, "_updating",     False)   # guard against recursive updates


        # ensure that only either diamter of radius is passed
        size_values_set = sum(x is not None for x in [diameter, radius])
        if size_values_set > 1:
            raise ValueError("Only one of diameter or radius may be set")

        catalogue = kwargs.pop("catalogue", self.catalogue)

        # Set properties for the Target object based on the arguments passed to the function
        _set_properties(self, 
                        name=name, 
                        radius=radius, 
                        diameter=diameter,
                        mass=mass, 
                        material_name=material_name,
                        catalogue=catalogue,
                        transition_scale_type=transition_scale_type, 
                        **kwargs
                    )
        arg_check = sum(x is None for x in [self.name, self.diameter, self.mass, self.transition_scale_type])
        if arg_check > 0:
            raise ValueError("Invalid Target")
        return
    
    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        # If this attribute is not being updated internally, mark it as a user-defined parameter
        if not self._updating:
            public_name = name.lstrip("_")
            cls = type(self)
            param = getattr(cls, public_name, None)
            if isinstance(param, property) and getattr(param, 'fset', None) is not None:
                # mark that the *public* name was user‐defined
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
    def radius(self) -> float | None:
        return self._radius

    @radius.setter
    def radius(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Radius must be positive")
        setattr(self, "_radius", float(value))

    @property
    def diameter(self) -> float | None:
        return self._diameter

    @diameter.setter
    def diameter(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Diameter must be positive")
        setattr(self, "_diameter", float(value))

    @property
    def mass(self) -> float | None:
        return self._mass

    @mass.setter
    def mass(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Mass must be positive")
        setattr(self, "_mass", float(value))

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
    def catalogue(self):
        """
        The target catalogue used for the target body.
        
        Returns
        -------
        dict 
        """
        def _make_target_catalogue():
            # Define some built-in catalogue values for known solar system targets of interest
            target_properties = [
                "name",    "radius",   "mass",      "material_name", "transition_scale_type"
            ]
            # The catalogue was created with Swiftest
            target_values = [
                ("Mercury", 2439.40e3, 3.301001e+23, "Soft Rock", "silicate"),
                ("Venus", 6051.84e3, 4.867306e+24, "Hard Rock", "silicate"),
                ("Earth", 6371.01e3, 5.972168e+24, "Wet Soil", "silicate"),
                ("Moon", 1737.53e3, 7.345789e+22, "Soft Rock", "silicate"),
                ("Mars", 3389.92e3, 6.416909e+23, "Soft Rock", "silicate"),
                ("Phobos", 11.17e3, 1.080000e+16, "Soft Rock", "silicate"),
                ("Deimos", 6.30e3, 1.800000e+15, "Soft Rock", "silicate"),
                ("Ceres", 469.70e3, 9.383516e+20, "Ice", "ice"),
                ("Vesta", 262.70e3, 2.590270e+20, "Soft Rock", "silicate"),
                ("Io", 1821.49e3, 8.929649e+22, "Hard Rock", "silicate"),
                ("Europa", 1560.80e3, 4.798574e+22, "Ice", "ice"),
                ("Ganymede", 2631.20e3, 1.481479e+23, "Ice", "ice"),
                ("Callisto", 2410.30e3, 1.075661e+23, "Ice", "ice"),
                ("Titan", 2575.50e3, 1.345181e+23, "Ice", "ice"),
                ("Rhea", 764.50e3, 2.306459e+21, "Ice", "ice"),
                ("Dione", 562.50e3, 1.095486e+21, "Ice", "ice"),
                ("Tethys", 536.30e3, 6.174430e+20, "Ice", "ice"),
                ("Enceladus", 252.30e3, 1.080318e+20, "Ice", "ice"),
                ("Mimas", 198.80e3, 3.750939e+19, "Ice", "ice"),
                ("Ariel", 578.90e3, 1.250019e+21, "Ice", "ice"),
                ("Umbriel", 584.70e3, 1.279535e+21, "Ice", "ice"),
                ("Titania", 788.90e3, 3.338178e+21, "Ice", "ice"),
                ("Oberon", 761.40e3, 3.076577e+21, "Ice", "ice"),
                ("Miranda", 235.70e3, 6.442623e+19, "Ice", "ice"),
                ("Triton", 1352.60e3, 2.140292e+22, "Ice", "ice"),
                ("Charon", 606.00e3, 1.589680e+21, "Ice", "ice"),
                ("Pluto", 1188.30e3, 1.302498e+22, "Ice", "ice"),
                ("Arrokoth", 9.13e3, 7.485000e+14, "Ice", "ice"),
            ]
            
            return _create_catalogue(target_properties, target_values)

        if self._catalogue is None:
            self._catalogue = _make_target_catalogue()
        return self._catalogue 

    @property
    def catalogue_key(self):
        """
        The key used to identify the property used as the key in a catalogue.
        """
        return "name"

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
    def escape_velocity(self) -> float:
        """
        Calculate the escape velocity for the target body in SI units.

        Returns
        -------
        float
            Escape velocity in m/s.
        """        
        return np.sqrt(2 * self.radius * self.gravity)
    
    @property
    def gravity(self) -> float | None:
        """
        Calculate the gravitational acceleration at the surface of the target body in SI units.

        Returns
        -------
        float
            Gravitational acceleration in m/s^2. 
        """
        return 4*np.pi*G.value*(self.radius)*self.bulk_density/3

    @property
    def bulk_density(self) -> float | None:
        """
        Calculate the bulk density of the target body in SI units.

        Returns
        -------
        float
            Bulk density in kg/m^3.
        """
        return self.mass / (4/3 * np.pi * (self.radius)**3)  # in kg/m^3

