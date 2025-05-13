from __future__ import annotations

from typing import Any

import numpy as np
from astropy.constants import G

from cratermaker.constants import FloatLike
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import (
    _create_catalogue,
    _set_properties,
    format_large_units,
    parameter,
)


class Target(ComponentBase):
    """
    Represents the target body in a crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.
    """

    _catalogue_header = [
        "name",
        "radius",
        "mass",
        "material",
        "transition_scale_type",
    ]
    # The catalogue was created with Swiftest
    _catalogue_values = [
        ("Mercury", 2439.40e3, 3.301001e23, "Soft Rock", "silicate"),
        ("Venus", 6051.84e3, 4.867306e24, "Hard Rock", "silicate"),
        ("Earth", 6371.01e3, 5.972168e24, "Wet Soil", "silicate"),
        ("Moon", 1737.53e3, 7.345789e22, "Soft Rock", "silicate"),
        ("Mars", 3389.92e3, 6.416909e23, "Soft Rock", "silicate"),
        ("Phobos", 11.17e3, 1.080000e16, "Soft Rock", "silicate"),
        ("Deimos", 6.30e3, 1.800000e15, "Soft Rock", "silicate"),
        ("Ceres", 469.70e3, 9.383516e20, "Soft Rock", "ice"),
        ("Vesta", 262.70e3, 2.590270e20, "Soft Rock", "silicate"),
        ("Io", 1821.49e3, 8.929649e22, "Hard Rock", "silicate"),
        ("Europa", 1560.80e3, 4.798574e22, "Ice", "ice"),
        ("Ganymede", 2631.20e3, 1.481479e23, "Ice", "ice"),
        ("Callisto", 2410.30e3, 1.075661e23, "Ice", "ice"),
        ("Titan", 2575.50e3, 1.345181e23, "Ice", "ice"),
        ("Rhea", 764.50e3, 2.306459e21, "Ice", "ice"),
        ("Dione", 562.50e3, 1.095486e21, "Ice", "ice"),
        ("Tethys", 536.30e3, 6.174430e20, "Ice", "ice"),
        ("Enceladus", 252.30e3, 1.080318e20, "Ice", "ice"),
        ("Mimas", 198.80e3, 3.750939e19, "Ice", "ice"),
        ("Ariel", 578.90e3, 1.250019e21, "Ice", "ice"),
        ("Umbriel", 584.70e3, 1.279535e21, "Ice", "ice"),
        ("Titania", 788.90e3, 3.338178e21, "Ice", "ice"),
        ("Oberon", 761.40e3, 3.076577e21, "Ice", "ice"),
        ("Miranda", 235.70e3, 6.442623e19, "Ice", "ice"),
        ("Triton", 1352.60e3, 2.140292e22, "Ice", "ice"),
        ("Charon", 606.00e3, 1.589680e21, "Ice", "ice"),
        ("Pluto", 1188.30e3, 1.302498e22, "Ice", "ice"),
        ("Arrokoth", 9.13e3, 7.485000e14, "Ice", "ice"),
    ]
    _catalogue = _create_catalogue(_catalogue_header, _catalogue_values)

    _density_catalogue = {
        "Water": 1000.0,
        "Sand": 1750.0,
        "Dry Soil": 1500.0,
        "Wet Soil": 2000.0,
        "Soft Rock": 2250.0,
        "Hard Rock": 2500.0,
        "Ice": 900.0,
    }

    def __init__(
        self,
        name: str,
        radius: FloatLike | None = None,
        diameter: FloatLike | None = None,
        mass: FloatLike | None = None,
        transition_scale_type: str | None = None,
        material: str | None = None,
        density: FloatLike | None = None,
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
        material : str or None
            Name of the material composition of the target body.
        density : FloatLike or None
            Volumetric density of the surface of the target body in kg/m^3.

        **kwargs : Any
            Additional keyword argumments that could be set by the user.

        Notes
        -----
        - The `radius` and `diameter` parameters are mutually exclusive. Only one of them should be provided.
        - Parameters set explicitly using keyword arguments will override those drawn from the catalogue.
        """
        object.__setattr__(self, "_updating", False)  # guard against recursive updates
        super().__init__(**kwargs)
        object.__setattr__(self, "_name", None)
        object.__setattr__(self, "_radius", None)
        object.__setattr__(self, "_mass", None)
        object.__setattr__(self, "_transition_scale_type", None)
        object.__setattr__(self, "_material", None)
        object.__setattr__(self, "_density", None)

        # ensure that only either diamter of radius is passed
        size_values_set = sum(x is not None for x in [diameter, radius])
        if size_values_set > 1:
            raise ValueError("Only one of diameter or radius may be set")
        if diameter is not None:
            radius = diameter / 2.0

        catalogue = kwargs.pop("catalogue", self.__class__._catalogue)

        # Set properties for the Target object based on the arguments passed to the function
        _set_properties(
            self,
            name=name,
            radius=radius,
            mass=mass,
            material=material,
            catalogue=catalogue,
            density=density,
            transition_scale_type=transition_scale_type,
            **kwargs,
        )
        arg_check = sum(
            x is None
            for x in [self.name, self.diameter, self.mass, self.transition_scale_type]
        )
        if arg_check > 0:
            raise ValueError("Invalid Target")

        if self._density is None:
            if self.material in self._density_catalogue:
                self._density = self._density_catalogue[self.material]

    def __str__(self) -> str:
        diameter = format_large_units(self.diameter, quantity="length")
        escape_velocity = format_large_units(self.escape_velocity, quantity="velocity")
        return (
            f"<Target: {self.name}>\n"
            f"Material: {self.material}\n"
            f"Diameter: {diameter}\n"
            f"Mass: {self.mass:.2e} kg\n"
            f"Surface density: {self.density:.1f} kg/m³\n"
            f"Transition Type: {self.transition_scale_type}\n"
            f"Escape Velocity: {escape_velocity}\n"
            f"Gravity: {self.gravity:.3f} m/s²"
        )

    @parameter
    def radius(self) -> float | None:
        return self._radius

    @radius.setter
    def radius(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Radius must be positive")
        setattr(self, "_radius", float(value))

    @property
    def diameter(self) -> float | None:
        if self._radius is not None:
            return 2 * self._radius

    @property
    def mass(self) -> float | None:
        return self._mass

    @mass.setter
    def mass(self, value: FloatLike):
        if value is not None and value <= 0:
            raise ValueError("Mass must be positive")
        setattr(self, "_mass", float(value))

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

    @parameter
    def material(self):
        """
        The name of the material composition of the target body.

        Returns
        -------
        str
        """
        return self._material

    @material.setter
    def material(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("material must be a string or None")
        self._material = value

    @parameter
    def density(self):
        """
        The volumetric density of the surface of the target body in kg/m^3.

        Returns
        -------
        float
        """
        return self._density

    @density.setter
    def density(self, value):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("density must be a numeric value or None")
            if value <= 0:
                raise ValueError("density must be a positive number")
            self._density = float(value)
        else:
            self._density = None
        return

    @property
    def catalogue(self) -> str:
        from cratermaker.utils.general_utils import format_large_units

        if self.__class__._catalogue is None:
            return "This Projectile component does not have a catalogue."
        lines = []
        # Use the self.__class__._catalogue_header list to make the header
        header = "|".join([f"{h:<11}" for h in self.__class__._catalogue_header])
        lines.append(header)
        lines.append(len(header) * "-")
        for name, entry in self.__class__._catalogue.items():
            line = ""
            for k, v in entry.items():
                if k == "radius":
                    val = format_large_units(v, quantity="length")
                elif k == "mass":
                    val = f"{v:.2e} kg"
                else:
                    val = v
                line += f"|{val:<11}"
            lines.append(f"{name:<11}{line}")
        return "\n".join(lines)

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
            raise ValueError(
                f"Invalid transition_scale_type: {value}. Must be one of {valid_types} or None"
            )
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
        return G.value * self.mass / (self.radius**2)

    @classmethod
    def maker(
        cls: type[Target],
        target: Target | str = "Moon",
        radius: FloatLike | None = None,
        diameter: FloatLike | None = None,
        mass: FloatLike | None = None,
        transition_scale_type: str | None = None,
        material: str | None = None,
        density: FloatLike | None = None,
        **kwargs: Any,
    ) -> Target:
        """
        Initialize the target object, setting properties from the provided arguments.

        Parameters
        ----------
        target : str, Target, or None
            Name of the target body or a Target object.
        radius : FloatLike or None
            Radius of the target body in km.
        diameter : FloatLike or None
            Diameter of the target body in km.
        mass : FloatLike or None
            Mass of the target body in kg.
        transition_scale_type : str or None
            Simple-to-complex transition scaling to use for the surface (either "silicate" or "ice").
        material : str or None
            Name of the material composition of the target body.
        density : FloatLike or None
            Volumetric density of the surface of the target body in kg/m^3.
        **kwargs : Any
            Additional keyword argumments that could be set by the user.

        Notes
        -----
        - The `radius` and `diameter` parameters are mutually exclusive. Only one of them should be provided.
        - Parameters set explicitly using keyword arguments will override those drawn from the catalogue.
        """

        if target is None:
            try:
                target = cls(
                    name="Moon",
                    radius=radius,
                    diameter=diameter,
                    mass=mass,
                    transition_scale_type=transition_scale_type,
                    material=material,
                    density=density,
                    **kwargs,
                )
            except Exception as e:
                raise RuntimeError("Error initializing target.") from e
        elif isinstance(target, str):
            try:
                target = cls(
                    name=target,
                    radius=radius,
                    diameter=diameter,
                    mass=mass,
                    transition_scale_type=transition_scale_type,
                    material=material,
                    density=density,
                    **kwargs,
                )
            except KeyError:
                raise ValueError(
                    f"Target '{target}' not found in the catalogue. Please provide a valid target name."
                )
        elif not isinstance(target, Target):
            raise TypeError("target must be a string or a Target object")
        target._component_name = target.name

        return target


import_components(__name__, __path__, ignore_private=True)
