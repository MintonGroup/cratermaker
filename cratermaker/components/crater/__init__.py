from __future__ import annotations

import hashlib
import math
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
from geopandas import GeoSeries
from numpy.random import Generator

from cratermaker.core.base import CratermakerBase
from cratermaker.utils.general_utils import format_large_units, validate_and_normalize_location

if TYPE_CHECKING:
    from cratermaker.components.projectile import Projectile
    from cratermaker.components.scaling import Scaling
    from cratermaker.components.surface import Surface
    from cratermaker.components.target import Target


@dataclass(frozen=True, slots=True)
class CraterFixed:
    id: np.uint32 = None
    semimajor_axis: float | None = None
    semiminor_axis: float | None = None
    orientation: float | None = None
    transient_diameter: float | None = None
    projectile_diameter: float | None = None
    projectile_velocity: float | None = None
    projectile_angle: float | None = None
    projectile_mass: float | None = None
    location: tuple[float, float] | None = None
    morphology_type: str | None = None
    time: float | None = None
    radius: float | None = field(default=None, init=False)
    diameter: float | None = field(default=None, init=False)
    transient_radius: float | None = field(default=None, init=False)
    projectile_radius: float | None = field(default=None, init=False)
    projectile_density: float | None = field(default=None, init=False)
    projectile_vertical_velocity: float | None = field(default=None, init=False)
    projectile_direction: float | None = field(default=None, init=False)

    def __post_init__(self):
        # Compute derived attributes that depend on other fields and cannot be set directly
        if self.semimajor_axis is not None and self.semiminor_axis is not None:
            object.__setattr__(self, "radius", math.sqrt(self.semimajor_axis * self.semiminor_axis))
            object.__setattr__(self, "diameter", self.radius * 2)
        if self.transient_diameter is not None:
            object.__setattr__(self, "transient_radius", self.transient_diameter / 2.0)
        if self.projectile_diameter is not None:
            object.__setattr__(self, "projectile_radius", self.projectile_diameter / 2.0)
            if self.projectile_mass is not None:
                object.__setattr__(
                    self, "projectile_density", (3.0 * self.projectile_mass) / (4.0 * math.pi * self.projectile_radius**3)
                )
        if self.projectile_velocity is not None and self.projectile_angle is not None:
            object.__setattr__(
                self, "projectile_vertical_velocity", self.projectile_velocity * math.sin(math.radians(self.projectile_angle))
            )
        if self.orientation is not None:
            object.__setattr__(self, "projectile_direction", self.orientation)


class CraterVariable:
    def __init__(
        self,
        measured_diameter: float | None = None,
        measured_radius: float | None = None,
        measured_semimajor_axis: float | None = None,
        measured_semiminor_axis: float | None = None,
        measured_orientation: float | None = None,
        measured_location: tuple[float, float] | None = None,
        measured_rim_height: float | None = None,
        measured_floor_depth: float | None = None,
        degradation_state: float | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_measured_semimajor_axis", None)
        object.__setattr__(self, "_measured_semiminor_axis", None)
        object.__setattr__(self, "_measured_radius", None)
        object.__setattr__(self, "_measured_orientation", None)
        object.__setattr__(self, "_measured_location", None)
        object.__setattr__(self, "_measured_rim_height", None)
        object.__setattr__(self, "_measured_floor_depth", None)
        object.__setattr__(self, "_degradation_state", None)

        if measured_diameter is not None:
            self.measured_diameter = measured_diameter
        if measured_radius is not None:
            self.measured_radius = measured_radius
        if measured_semimajor_axis is not None:
            self.measured_semimajor_axis = measured_semimajor_axis
        if measured_semiminor_axis is not None:
            self.measured_semiminor_axis = measured_semiminor_axis
        if measured_orientation is not None:
            self.measured_orientation = measured_orientation
        if measured_location is not None:
            self.measured_location = measured_location
        if measured_rim_height is not None:
            self.measured_rim_height = measured_rim_height
        if measured_floor_depth is not None:
            self.measured_floor_depth = measured_floor_depth
        if degradation_state is not None:
            self.degradation_state = degradation_state
        return

    def __repr__(self):
        return (
            f"CraterVariable(measured_diameter={self.measured_diameter}, "
            f"measured_radius={self.measured_radius}, "
            f"measured_semimajor_axis={self.measured_semimajor_axis}, "
            f"measured_semiminor_axis={self.measured_semiminor_axis}, "
            f"measured_orientation={self.measured_orientation}, "
            f"measured_location={self.measured_location}, "
            f"measured_rim_height={self.measured_rim_height}, "
            f"measured_floor_depth={self.measured_floor_depth}, "
            f"degradation_state={self.degradation_state})"
        )

    def as_dict(self):
        dict_repr = {
            "measured_orientation": self.measured_orientation,
            "measured_location": self.measured_location,
            "measured_rim_height": self.measured_rim_height,
            "measured_floor_depth": self.measured_floor_depth,
            "degradation_state": self.degradation_state,
        }
        if self.measured_semimajor_axis is not None:
            if self.measured_radius is not None and self.measured_radius == self.measured_semimajor_axis:
                dict_repr["measured_diameter"] = self.measured_diameter
            else:
                dict_repr["measured_semimajor_axis"] = self.measured_semimajor_axis
                dict_repr["measured_semiminor_axis"] = self.measured_semiminor_axis
        return dict_repr

    @property
    def measured_diameter(self) -> float | None:
        """Final diameter of the crater in meters."""
        return self.measured_radius * 2 if self.measured_radius is not None else None

    @measured_diameter.setter
    def measured_diameter(self, value: float | None):
        if value is None:
            self.measured_radius = None
        else:
            self.measured_radius = value / 2.0
        return

    @property
    def measured_radius(self) -> float | None:
        """Final radius of the crater in meters."""
        if self.measured_semimajor_axis is not None and self.measured_semiminor_axis is not None:
            return math.sqrt(self.measured_semimajor_axis * self.measured_semiminor_axis)
        return None

    @measured_radius.setter
    def measured_radius(self, value: float | None):
        if value is None:
            self.measured_semimajor_axis = None
            self.measured_semiminor_axis = None
        else:
            self.measured_semimajor_axis = value
            self.measured_semiminor_axis = value
        return

    @property
    def measured_semimajor_axis(self) -> float | None:
        return self._measured_semimajor_axis

    @measured_semimajor_axis.setter
    def measured_semimajor_axis(self, value: float | None):
        if value is None:
            self._measured_semimajor_axis = None
            return
        if value <= 0:
            raise ValueError("measured_semimajor_axis must be positive.")
        if self._measured_semiminor_axis is not None and value < self._measured_semiminor_axis:
            self._measured_semimajor_axis = self._measured_semiminor_axis
            self._measured_semiminor_axis = float(value)
        else:
            self._measured_semimajor_axis = float(value)
        return

    @property
    def measured_semiminor_axis(self) -> float | None:
        return self._measured_semiminor_axis

    @measured_semiminor_axis.setter
    def measured_semiminor_axis(self, value: float | None):
        if value is None:
            self._measured_semiminor_axis = None
            return
        if value <= 0:
            raise ValueError("measured_semiminor_axis must be positive.")
        if self._measured_semimajor_axis is not None and value > self._measured_semimajor_axis:
            self._measured_semiminor_axis = self._measured_semimajor_axis
            self._measured_semimajor_axis = float(value)
        else:
            self._measured_semiminor_axis = float(value)
        return

    @property
    def measured_orientation(self) -> float | None:
        return self._measured_orientation

    @measured_orientation.setter
    def measured_orientation(self, value: float | None):
        if value is None:
            self._measured_orientation = None
            return
        self._measured_orientation = float(value) % 360
        return

    @property
    def measured_location(self) -> tuple[float, float] | None:
        return self._measured_location

    @measured_location.setter
    def measured_location(self, value: tuple[float, float] | None):
        if value is None:
            self._measured_location = None
            return
        self._measured_location = validate_and_normalize_location(value)
        return

    @property
    def measured_rim_height(self) -> float | None:
        return self._measured_rim_height

    @measured_rim_height.setter
    def measured_rim_height(self, value: float | None):
        if value is None:
            self._measured_rim_height = None
            return
        self._measured_rim_height = float(value)
        return

    @property
    def measured_floor_depth(self) -> float | None:
        return self._measured_floor_depth

    @measured_floor_depth.setter
    def measured_floor_depth(self, value: float | None):
        if value is None:
            self._measured_floor_depth = None
            return
        self._measured_floor_depth = float(value)
        return

    @property
    def degradation_state(self) -> float | None:
        return self._degradation_state

    @degradation_state.setter
    def degradation_state(self, value: float | None):
        if value is None:
            self._degradation_state = None
            return
        self._degradation_state = float(value)
        return


class Crater:
    def __init__(self, crater: Crater | None = None, fixed_cls=CraterFixed, var_cls=CraterVariable, **kwargs):
        if crater is not None:
            if not isinstance(crater, Crater):
                raise TypeError("crater must be an instance of Crater or None")
            fixed_fields = asdict(crater._fixed)
            for f in fields(crater._fixed):
                if not f.init:
                    fixed_fields.pop(f.name)
            var_fields = crater._var.as_dict()
            # Be sure to scrub any fields that are redundant with any potential kwargs
            kwargs = self._scrub_measured_input_size(**kwargs)
            if "measured_semimajor_axis" in kwargs or "measured_semiminor_axis" in kwargs:
                var_fields.pop("measured_diameter", None)
                var_fields.pop("measured_radius", None)
        else:
            fixed_fields = {}
            var_fields = {}

        for f in fields(fixed_cls):
            if f.name in kwargs and f.init and (f.name not in fixed_fields or fixed_fields[f.name] is None):
                fixed_fields[f.name] = kwargs.pop(f.name)
            if f.name not in fixed_fields and f.name in kwargs and f.init:
                fixed_fields[f.name] = kwargs.pop(f.name)

        var_fields.update({k: v for k, v in kwargs.items() if k not in fixed_fields})
        self._fixed = fixed_cls(**fixed_fields)
        self._var = var_cls(**var_fields)
        return

    def __getattr__(self, name: str):
        """
        Delegates attribute access to the fixed and variable fields.

        If an attribute is not found in the Crater class, it will first check the fixed fields, then the variable fields.
        """
        if hasattr(self._fixed, name):
            return getattr(self._fixed, name)
        if hasattr(self._var, name):
            return getattr(self._var, name)
        raise AttributeError(f"{type(self).__name__} has no attribute {name!r}")

    def __setattr__(self, name: str, value: Any):
        """
        Delegates attribute setting to the variable fields.

        Only attributes in the variable fields can be set. Fixed fields are immutable.
        """
        if name in ("_fixed", "_var", "_fixed_fields", "_var_fields"):
            object.__setattr__(self, name, value)
            return

        if name in CraterFixed.__dataclass_fields__:
            raise AttributeError(f"{name!r} is a fixed attribute and cannot be modified.")

        # Call the setters in CraterVariable if they exist
        if hasattr(self._var, name):
            setattr(self._var, name, value)
            return

    def as_dict(self, ignore_keys: list[str] = [], **kwargs) -> dict:
        """
        Expose delegated attributes for serialization.
        """
        base = {}
        base.update(asdict(self._fixed))
        base.update(self._var.as_dict())
        if ignore_keys is not None:
            for key in ignore_keys:
                base.pop(key, None)
        return base

    def __dir__(self):
        """
        Expose delegated attributes for autocomplete/introspection.
        """
        base = set(super().__dir__())
        base.update(dir(self.fixed))
        base.update(dir(self.var))
        return sorted(base)

    def __str__(self):
        if self.time is None:
            timetext = "Not set"
        else:
            timetext = f"{format_large_units(self.time, quantity='time')}"
        if self.semimajor_axis != self.semiminor_axis:
            size_text = f"Elliptical size {format_large_units(self.semimajor_axis, quantity='length')} x {format_large_units(self.semiminor_axis, quantity='length')}\nMean Diameter: {format_large_units(self.diameter, quantity='length')}"
        else:
            size_text = f"Diameter: {format_large_units(self.diameter, quantity='length')}"
        if self.measured_semimajor_axis is not None and self.measured_semimajor_axis != self.semimajor_axis:
            if self.measured_semimajor_axis != self.measured_semiminor_axis:
                size_text += f"\nMeasured elliptical size {format_large_units(self.measured_semimajor_axis, quantity='length')} x {format_large_units(self.measured_semiminor_axis, quantity='length')}\nMeasured Mean Diameter: {format_large_units(self.measured_diameter, quantity='length')}"
            else:
                size_text += f"\nMeasured diameter: {format_large_units(self.measured_diameter, quantity='length')}"
        if self.measured_orientation is not None and self.measured_orientation != self.orientation:
            size_text += f"\nMeasured orientation: {self.measured_orientation:.1f}°"
        if self.measured_location is not None and self.measured_location != self.location:
            size_text += f"\nMeasured location (lon, lat): ({self.measured_location[0]:.4f}°, {self.measured_location[1]:.4f}°)\n"

        return (
            f"ID: {self.id}\n"
            f"{size_text}\n"
            f"transient_diameter: {format_large_units(self.transient_diameter, quantity='length')}\n"
            f"projectile_diameter: {format_large_units(self.projectile_diameter, quantity='length')}\n"
            f"projectile_mass: {self.projectile_mass:.4e} kg\n"
            f"projectile_density: {self.projectile_density:.0f} kg/m^3\n"
            f"projectile_velocity: {format_large_units(self.projectile_velocity, quantity='velocity')}\n"
            f"projectile_angle: {self.projectile_angle:.1f}°\n"
            f"projectile_direction: {self.projectile_direction:.1f}°\n"
            f"location (lon,lat): ({self.location[0]:.4f}°, {self.location[1]:.4f}°)\n"
            f"morphology_type: {self.morphology_type}\n"
            f"time: {timetext}"
        )

    def _scrub_measured_input_size(self, **kwargs: Any):
        if "measured_diameter" in kwargs:
            measured_diameter = kwargs.pop("measured_diameter")
            kwargs["measured_semimajor_axis"] = measured_diameter / 2.0
            kwargs["measured_semiminor_axis"] = measured_diameter / 2.0
        elif "measured_radius" in kwargs:
            measured_radius = kwargs.pop("measured_radius")
            kwargs["measured_semimajor_axis"] = measured_radius
            kwargs["measured_semiminor_axis"] = measured_radius
        return kwargs

    @classmethod
    def maker(
        cls: type[Crater],
        crater: Crater | None = None,
        scaling: str | Scaling | None = None,
        target: str | Target | None = None,
        projectile: str | Projectile | None = None,
        diameter: float | None = None,
        radius: float | None = None,
        final_diameter: float | None = None,
        final_radius: float | None = None,
        transient_diameter: float | None = None,
        transient_radius: float | None = None,
        semimajor_axis: float | None = None,
        semiminor_axis: float | None = None,
        orientation: float | None = None,
        projectile_diameter: float | None = None,
        projectile_radius: float | None = None,
        projectile_mass: float | None = None,
        projectile_density=None,
        projectile_velocity: float | None = None,
        projectile_mean_velocity: float | None = None,
        projectile_vertical_velocity=None,
        projectile_angle: float | None = None,
        projectile_direction: float | None = None,
        projectile_location: tuple[float, float] | None = None,
        location: tuple[float, float] | None = None,
        time: float | None = None,
        measured_semimajor_axis: float | None = None,
        measured_semiminor_axis: float | None = None,
        measured_orientation: float | None = None,
        measured_diameter: float | None = None,
        measured_radius: float | None = None,
        measured_location: tuple[float, float] | None = None,
        measured_rim_height: float | None = None,
        measured_floor_depth: float | None = None,
        degradation_state: float | None = None,
        simdir: str | Path | None = None,
        rng: Generator = None,
        rng_seed: str | int | None = None,
        rng_state: dict | None = None,
        check_redundant_inputs: bool = True,
        **kwargs: Any,
    ):
        """
        Create a Crater object with the given parameters.

        Parameters
        ----------
        crater : Crater, optional
            A Crater object to copy parameters from.
        scaling : str or Scaling, optional
            A string key or instance of a scaling model. If none provided, a default will be used
        target : str or Target, optional
            A string key or instance of a target model. If none provided, a default will be used.
        projectile : str or Projectile, optional
            A string key or instance of a projectile model. If none provided, a default will be used.
        diameter : float, optional
            The final diameter of the crater in meters.
        final_diameter : float, optional
            The final diameter of the crater in meters. This is an alias of "diameter"
        radius : float, optional
            The final radius of the crater in meters.
        final_radius : float, optional
            The final radius of the crater in meters. This is an alias of "radius"
        transient_diameter : float, optional
            The transient diameter of the crater in meters.
        transient_radius : float, optional
            The transient radius of the crater in meters.
        semimajor_axis : float, optional
            The semimajor axis of the crater in meters. Same as radius if the crater is circular.
        semiminor_axis : float, optional
            The semiminor axis of the crater in meters. Same as radius if the crater is circular.
        orientation : float, optional
            The orientation of the crater in degrees if it is elliptical.
        projectile_diameter : float, optional
            The diameter of the projectile in meters.
        projectile_radius : float, optional
            The radius of the projectile in meters.
        projectile_mass : float, optional
            The mass of the projectile in kilograms.
        projectile_density : float, optional
            The density of the projectile in kg/m^3. If not provided, it will be defined through the projectile population model provided.
        projectile_velocity : float, optional
            The total impact velocity of the projectile in m/s.
        projectile_mean_velocity : float, optional
            The mean velocity from which to sample a projectile velocity.
        projectile_vertical_velocity : float, optional
            The vertical component of the velocity in m/s.
        projectile_angle : float, optional
            The impact angle in degrees (0-90).
        projectile_direction : float, optional
            The direction of the impact in degrees (0-360).
        projectile_location : tuple of float, optional
            The (longitude, latitude) location of the projectile impact. This is equivalent to `location`, which takes precedence
        location : tuple of float, optional
            The (longitude, latitude) location of the crater.
        time : float, optional
            The time of the crater impact in Myr before present.
        measured_semimajor_axis : float, optional
            The measured semimajor axis of the crater in meters.
        measured_semiminor_axis : float, optional
            The measured semiminor axis of the crater in meters.
        measured_orientation : float, optional
            The measured orientation of the crater in degrees.
        measured_diameter : float, optional
            The measured diameter of the crater in meters.
        measured_radius : float, optional
            The measured radius of the crater in meters.
        measured_location : tuple of float, optional
            The measured (longitude, latitude) location of the crater.
        measured_rim_height : float, optional
            The measured rim height of the crater in meters.
        measured_floor_depth : float, optional
            The measured floor depth of the crater in meters.
        degradation_state : float, optional
            The current degradation state of the crater in meters squared.
        simdir : str | Path
            |simdir|
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        check_redundant_inputs : bool, optional
            If True, check for redundant inputs such as providing both diameter and radius. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        Crater
            A frozen Crater dataclass with derived attributes.

        Notes
        -----
        - Unless `check_redundant_inputs` is False, the following rules apply:
        - Exactly one of the following must be provided: `diameter`, `radius`, both of `semimajor_axis` and `semiminor_axis`, `transient_diameter`, `transient_radius`, `projectile_diameter`, `projectile_radius`, or `projectile_mass`.
        - Velocity may be specified in one of these ways:
        - `projectile_mean_velocity` alone (samples a velocity)
        - Any two of (`projectile_velocity`, `projectile_vertical_velocity`, `projectile_angle`). the third is inferred.
        - `projectile` is mutually exclusive with velocity-related inputs; if provided, it overrides velocity, angle, direction, and density unless explicitly set.
        - The `scaling`, and `rng` models are required for scaling and density inference, but are not stored in the returned Crater object.
        - If providing measured properties, either both `measured_semimajor_axis` and `measured_semiminor_axis`, or one of `measured_diameter` or `measured_radius` must be provided.
        - if `check_redundant_inputs` is True and any of the rules are violated, then it is assumed that the user wants a copy. In that case, some parameters will not be computed and only the user-input arguments will be used
        """
        from cratermaker.components.projectile import Projectile
        from cratermaker.components.scaling import Scaling
        from cratermaker.components.target import Target

        make_copy = crater is not None

        def _set_id(**kwargs: Any) -> np.uint32:
            """
            Sets the hash id of the crater based on input parameters.

            This is used as a unique identifier for the crater. To reduce storage requirements, the id is a 32-bit unsigned integer derived from the hash of the input parameters. There is a very small chance of hash collisions, but this is acceptable for the purposes of crater identification.

            Parameters
            ----------
            **kwargs : Any
                |kwargs|

            """
            id_args = [
                "semimajor_axis",
                "semiminor_axis",
                "orientation",
                "transient_diameter",
                "projectile_diameter",
                "projectile_velocity",
                "projectile_angle",
                "projectile_mass",
                "location",
            ]
            combined_args = [k for k in kwargs if k in id_args]
            combined_args.sort()  # Sort to ensure consistent ordering
            combined = "::".join(f"{k}:{kwargs[k]}" for k in combined_args)
            hexid = hashlib.shake_128(combined.encode()).hexdigest(4)
            return np.uint32(int(f"0x{hexid}", 16))

        # Convert from the old API "final_diameter/final_radius" to "diameter/radius"
        if final_diameter is not None:
            diameter = final_diameter
        if final_radius is not None:
            radius = final_radius
        if radius is not None and diameter is not None:
            if check_redundant_inputs:
                raise ValueError("Only one of diameter or radius may be set.")
            else:
                diameter = None  # Give priority to radius if both are set and we aren't checking for redundant inputs

        if location is None:
            location = projectile_location

        # Validate ellipticity parameters
        if semimajor_axis is not None or semiminor_axis is not None:
            if semimajor_axis is None or semiminor_axis is None:  # We will assume circular if only one is set
                radius = semimajor_axis if semiminor_axis is None else semiminor_axis
                semimajor_axis = radius
                semiminor_axis = radius
            if check_redundant_inputs and (diameter is not None or radius is not None):
                raise ValueError(
                    "diameter or radius cannot be used for elliptical craters; use semimajor_axis and semiminor_axis instead."
                )

        # Turn any circular crater arguments into elliptial ones
        if semimajor_axis is None and semiminor_axis is None:
            if radius is not None:
                semimajor_axis = radius
                semiminor_axis = radius
            elif diameter is not None:
                semimajor_axis = diameter / 2.0
                semiminor_axis = diameter / 2.0

        # For now, crater orientation and projectile direction are linked
        if check_redundant_inputs and orientation is not None and projectile_direction is not None:
            raise ValueError("Only one of orientation or projectile_direction may be set for elliptical craters.")
        if orientation is not None and projectile_direction is None:
            projectile_direction = orientation
        if projectile_direction is not None and orientation is None:
            orientation = projectile_direction

        velocity_inputs = {
            "projectile_velocity": projectile_velocity,
            "projectile_vertical_velocity": projectile_vertical_velocity,
            "projectile_angle": projectile_angle,
        }

        n_velocity_inputs = sum(x is not None for x in velocity_inputs.values())
        if n_velocity_inputs > 2:
            if check_redundant_inputs:
                raise ValueError(f"Only two of {', '.join(k for k, v in velocity_inputs.items() if v is not None)} may be set.")
            else:
                make_copy = True

        if projectile_mean_velocity is not None and (projectile_velocity is not None or projectile_vertical_velocity is not None):
            if check_redundant_inputs:
                raise ValueError("projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity")
            else:
                make_copy = True

        size_inputs = {
            "semimajor_axis": semimajor_axis,
            "transient_diameter": transient_diameter,
            "transient_radius": transient_radius,
            "projectile_diameter": projectile_diameter,
            "projectile_radius": projectile_radius,
            "projectile_mass": projectile_mass,
        }

        # Assemble initial arguments
        args = {
            "semimajor_axis": semimajor_axis,
            "semiminor_axis": semiminor_axis,
            "orientation": orientation,
            "transient_diameter": transient_diameter,
            "projectile_diameter": projectile_diameter,
            "projectile_mass": projectile_mass,
            "projectile_velocity": projectile_velocity,
            "projectile_angle": projectile_angle,
            "morphology_type": "Not Set",
            "location": location,
            "time": time,
        }

        if measured_diameter is not None:
            args["measured_diameter"] = measured_diameter
        if measured_radius is not None:
            args["measured_radius"] = measured_radius
        if measured_semimajor_axis is not None:
            args["measured_semimajor_axis"] = measured_semimajor_axis
        if measured_semiminor_axis is not None:
            args["measured_semiminor_axis"] = measured_semiminor_axis
        if measured_orientation is not None:
            args["measured_orientation"] = measured_orientation
        if measured_location is not None:
            args["measured_location"] = measured_location
        if measured_rim_height is not None:
            args["measured_rim_height"] = measured_rim_height
        if measured_floor_depth is not None:
            args["measured_floor_depth"] = measured_floor_depth
        if degradation_state is not None:
            args["degradation_state"] = degradation_state

        n_size_inputs = sum(v is not None for v in size_inputs.values())

        # Process input crater object if provided
        if crater is not None:
            if not isinstance(crater, Crater):
                raise TypeError("crater must be a Crater object.")
            old_parameters = {}
            for field in crater.as_dict():
                if field in locals() and locals()[field] is None and getattr(crater, field) is not None:
                    old_parameters[field] = getattr(crater, field)
            if (
                n_size_inputs == 0
            ):  # The user has not passed any size parameters, so we will use the size parameters from the crater object
                n_size_inputs = 1  # Make sure we don't trigger the error below
                make_copy = True
            else:  # The user is passing a size parameter, so we cannot use any of the values from the crater object
                make_copy = False
                for field in size_inputs:
                    old_parameters.pop(field, None)
            if projectile_mean_velocity is None:
                # Be sure to keep only two velocitity components, with a preference for angle over vertical velocity
                if n_velocity_inputs == 0:
                    old_parameters.pop("projectile_vertical_velocity", None)
                elif n_velocity_inputs == 1:
                    for k, v in velocity_inputs.items():
                        if v is not None:
                            old_parameters.pop(k, None)
                            if k == "projectile_velocity" or k == "projectile_angle":
                                old_parameters.pop("projectile_vertical_velocity", None)
                            elif k == "projectile_vertical_velocity":
                                old_parameters.pop("projectile_velocity", None)
                elif n_velocity_inputs == 2:
                    for field in velocity_inputs:
                        old_parameters.pop(field, None)
                # Make sure we don't override the projectile_velocity and projectile_angle
                n_velocity_inputs = 2

            # Now set the local arguments to be what's left from the old_parameters
            for field in old_parameters:
                if field in args and args[field] is None:
                    args[field] = old_parameters[field]

        if n_size_inputs != 1:
            if check_redundant_inputs:
                raise ValueError(f"Exactly one of {', '.join(k for k, v in size_inputs.items() if v is not None)} must be set.")
            else:
                make_copy = True

        if make_copy:
            # The following will change the id
            a = args.pop("semimajor_axis")
            b = args.pop("semiminor_axis")
            td = args.pop("transient_diameter")
            pdir = args.pop("orientation")
            pd = args.pop("projectile_diameter")
            pm = args.pop("projectile_mass")
            pv = args.pop("projectile_velocity")
            pang = args.pop("projectile_angle")
            location = args.pop("location")
            if "measured_diameter" in args:
                measured_diameter = args.pop("measured_diameter")
                measured_semimajor_axis = measured_diameter / 2.0
                measured_semiminor_axis = measured_diameter / 2.0
                measured_diameter = None
                measured_radius = None
            elif "measured_radius" in args:
                measured_radius = args.pop("measured_radius")
                measured_semimajor_axis = measured_radius
                measured_semiminor_axis = measured_radius
                measured_diameter = None
                measured_radius = None
            elif "measured_semimajor_axis" in args or "measured_semiminor_axis" in args:
                measured_semimajor_axis = args.pop("measured_semimajor_axis", None)
                measured_semiminor_axis = args.pop("measured_semiminor_axis", None)

            measured_orientation = args.pop("measured_orientation", None)
            measured_location = args.pop("measured_location", None)

            if crater is not None:
                mt = crater.morphology_type
            else:
                mt = kwargs.pop("morphology_type", "UNKNOWN")

        else:
            # --- Merge the CommonArgs with any arguments passed here using CratermakerBase ---
            argproc = CratermakerBase(simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
            if scaling is not None and isinstance(scaling, Scaling):
                if projectile is None:
                    projectile = scaling.projectile
                if target is None:
                    target = scaling.target

            kwargs = {**kwargs, **vars(argproc.common_args)}
            target = Target.maker(target, **kwargs)

            projectile_args = {
                "velocity": args["projectile_velocity"],
                "mass": args["projectile_mass"],
                "angle": args["projectile_angle"],
                "direction": args["orientation"],
                "location": args["location"],
            }
            projectile = Projectile.maker(
                projectile,
                target=target,
                **kwargs,
            ).new_projectile(
                **projectile_args,
            )

            scaling = Scaling.maker(
                scaling,
                target=target,
                projectile=projectile,
                **kwargs,
            )
            projectile = scaling.projectile
            scaling.recompute()
            target = scaling.target
            pv = projectile.velocity
            pvv = projectile.vertical_velocity
            pang = projectile.angle
            pdir = projectile.direction
            prho = projectile.density
            location = projectile.location

            # --- Ensure velocity/angle are all set ---
            n_set = sum(x is not None for x in [pv, pvv, pang])
            if n_set != 3:
                raise ValueError("Not enough information to infer a projectile velocity.")

            # --- Compute derived quantities ---
            fr = radius
            fd = diameter
            a = args["semimajor_axis"]
            b = args["semiminor_axis"]
            if a is not None and b is not None:
                fr = math.sqrt(a * b)

            if fr is not None:
                fd = 2 * fr
            elif fd is not None:
                fr = fd / 2.0

            td = args["transient_diameter"]
            tr = transient_radius
            pd = args["projectile_diameter"]
            pr = projectile_radius
            pm = projectile_mass
            mt = None

            if pr is not None:
                pd = 2 * pr
            if fr is not None:
                fd = 2 * fr
            if tr is not None:
                td = 2 * tr

            if fd is not None:
                td, mt = scaling.final_to_transient(fd)
                pd = scaling.transient_to_projectile(td)
            elif td is not None:
                fd, mt = scaling.transient_to_final(td)
                pd = scaling.transient_to_projectile(td)

            elif pd is not None:
                td = scaling.projectile_to_transient(pd)
                fd, mt = scaling.transient_to_final(td)

            elif pm is not None:
                pr = ((3.0 * pm) / (4.0 * math.pi * prho)) ** (1.0 / 3.0)
                pd = 2.0 * pr
                td = scaling.projectile_to_transient(pd)
                fd, mt = scaling.transient_to_final(td)

            pr = pd / 2
            tr = td / 2
            fr = fd / 2
            pm = (4.0 / 3.0) * math.pi * pr**3 * prho

            if a is None:
                a = fr
            if b is None:
                b = fr

            if a < 0 or b < 0 or td < 0 or pd < 0:
                raise ValueError("Crater and projectile sizes must be non-negative.")
            if pv < 0:
                raise ValueError("Projectile velocity must be non-negative.")
            if prho < 0:
                raise ValueError("Projectile density must be non-negative.")

        if check_redundant_inputs and measured_radius is not None and measured_diameter is not None:
            raise ValueError("Only one of measured_diameter or measured_radius may be set.")
        else:
            if measured_diameter is not None:
                measured_semimajor_axis = measured_diameter / 2.0
                measured_semiminor_axis = measured_diameter / 2.0
                measured_diameter = None
            if measured_radius is not None:
                measured_semimajor_axis = measured_radius
                measured_semiminor_axis = measured_radius
                measured_radius = None

        if measured_semimajor_axis is not None or measured_semiminor_axis is not None:
            if measured_semimajor_axis is None or measured_semiminor_axis is None:
                if check_redundant_inputs and (measured_diameter is not None or measured_radius is not None):
                    raise ValueError(
                        "If providing measured properties for an elliptical crater, either both measured_semimajor_axis and measured_semiminor_axis, or one of measured_diameter or measured_radius must be provided."
                    )
                measured_radius = measured_semimajor_axis if measured_semiminor_axis is None else measured_semiminor_axis
                measured_semimajor_axis = measured_radius
                measured_semiminor_axis = measured_radius
            if check_redundant_inputs and (measured_diameter is not None or measured_radius is not None):
                raise ValueError(
                    "measured_diameter or measured_radius cannot be used for elliptical craters; use measured_semimajor_axis and measured_semiminor_axis instead."
                )
        elif measured_radius is not None:
            measured_semimajor_axis = measured_radius
            measured_semiminor_axis = measured_radius
        else:
            measured_semimajor_axis = a if a is not None else None
            measured_semiminor_axis = b if b is not None else None

        if measured_orientation is None:
            measured_orientation = pdir if pdir is not None else None
        if measured_location is None:
            measured_location = (float(location[0]), float(location[1])) if location is not None else None

        # Assemble the final set of arguments to pass to the Crater constructor
        args["semimajor_axis"] = float(a) if a is not None else None
        args["semiminor_axis"] = float(b) if b is not None else None
        args["transient_diameter"] = float(td) if td is not None else None
        args["orientation"] = float(pdir) if pdir is not None else None
        args["projectile_diameter"] = float(pd) if pd is not None else None
        args["projectile_mass"] = float(pm) if pm is not None else None
        args["projectile_velocity"] = float(pv) if pv is not None else None
        args["projectile_angle"] = float(pang) if pang is not None else None
        args["location"] = (float(location[0]), float(location[1]))
        args["morphology_type"] = str(mt) if mt is not None else None
        if measured_location is not None:
            args["measured_location"] = (float(measured_location[0]), float(measured_location[1]))
        args["measured_semimajor_axis"] = float(measured_semimajor_axis) if measured_semimajor_axis is not None else None
        args["measured_semiminor_axis"] = float(measured_semiminor_axis) if measured_semiminor_axis is not None else None
        args["measured_orientation"] = float(measured_orientation) if measured_orientation is not None else None
        args["id"] = _set_id(**args)
        return cls(**args, **kwargs)

    def to_geoseries(
        self,
        n: int = 150,
        surface: Surface | None = None,
        split_antimeridian: bool = True,
        use_measured_properties: bool = True,
        **kwargs: Any,
    ) -> GeoSeries:
        """
        Geodesic ellipse on a sphere: for each bearing theta from the center, we shoot a geodesic with distance r(theta) = (a*b)/sqrt((b*ct)^2 + (a*st)^2), then rotate all bearings by the crater's orientation.

        Parameters
        ----------
        n : int, optional
            Number of points to use for the polygon, by default 150.
        surface : Surface | None, optional
            Surface object providing planetary radius and CRS.
        split_antimeridian : bool, optional
            If True, split the polygon into a GeometryCollection if it crosses the antimeridian, by default True.
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        A GeoSeries containing a Shapely Polygon in lon/lat degrees.
        """
        from cratermaker._cratermaker import counting_bindings
        from pyproj import Geod
        from shapely.geometry import GeometryCollection, LineString, Polygon
        from shapely.ops import split, transform

        from cratermaker.components.surface import Surface

        surface = Surface.maker(surface)
        geod = Geod(a=surface.target.radius, b=surface.target.radius)
        theta = np.linspace(0.0, 360.0, num=n, endpoint=False)
        if use_measured_properties:
            a = self.measured_semimajor_axis
            b = self.measured_semiminor_axis
            lon, lat = self.measured_location
            phi = self.measured_orientation
        else:
            a = self.semimajor_axis
            b = self.semiminor_axis
            lon, lat = self.location
            phi = self.orientation

        # Measure the rim height so that the polygon sits on to of the surface rather than underneath
        region = surface.extract_region(location=self.location, region_radius=self.measured_radius, at_least_one_face=True)
        rim_height = np.max(region.face_elevation)
        z = np.full_like(theta, rim_height)

        # Polar radius of an axis-aligned ellipse in a Euclidean tangent plane
        ct = np.cos(np.deg2rad(theta))
        st = np.sin(np.deg2rad(theta))
        r = (a * b) / np.sqrt((b * ct) ** 2 + (a * st) ** 2)

        # Bearings (from east, CCW) rotated by azimuth
        bearings = theta + phi % 360.0

        # Forward geodesic for each bearing/distance
        poly_lon, poly_lat, _ = geod.fwd(lon * np.ones_like(bearings), lat * np.ones_like(bearings), bearings, r)

        # Correct for potential antimeridian crossing
        if split_antimeridian and np.ptp(poly_lon) > 180.0:
            center_sign = np.sign(lon)
            poly_lon = np.where(np.sign(poly_lon) != center_sign, poly_lon + 360.0 * center_sign, poly_lon)
            poly = Polygon(zip(poly_lon, poly_lat, z, strict=True))
            merdian_lon = 180.0 * center_sign
            meridian = LineString([(merdian_lon, -90.0), (merdian_lon, 90.0)])
            poly = split(poly, meridian)

            def lon_flip(lon, lat, z=None):
                lon = np.where(np.abs(lon) >= 180.0, lon - 360.0 * np.sign(lon), lon)
                return lon, lat

            new_geoms = []
            for p in poly.geoms:
                if np.abs(p.centroid.x) > 180.0:
                    p = transform(lon_flip, p)
                new_geoms.append(p)
            poly = GeometryCollection(new_geoms)

        else:
            poly = Polygon(zip(poly_lon, poly_lat, z, strict=True))

        return GeoSeries([poly], crs=surface.crs)

    @property
    def final_diameter(self) -> float | None:
        """Final diameter of the crater in meters."""
        return self.diameter

    @property
    def final_radius(self) -> float | None:
        """Final radius of the crater in meters."""
        return self.radius

    def remove_complex_data(self, **kwargs) -> None:
        """
        Remove complex data types from the crater variable properties.

        This is useful for storing a lightweight representation, as it removes complex data types that can be recomputed from the fixed properties and the morphology model when needed. The base  Crater type does not have any complex data types, but is used by derived types that do, so this method is included here for consistency and to allow for future expansion of the base Crater type without breaking the API.
        """
        pass
