from __future__ import annotations

import hashlib
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from numpy.random import Generator

from ..utils.general_utils import format_large_units
from .base import CratermakerBase

if TYPE_CHECKING:
    from ..components.projectile import Projectile
    from ..components.scaling import Scaling
    from ..components.target import Target


@dataclass(frozen=True, slots=True)
class Crater:
    final_diameter: float | None = None
    transient_diameter: float | None = None
    projectile_diameter: float | None = None
    projectile_velocity: float | None = None
    projectile_angle: float | None = None
    projectile_direction: float | None = None
    projectile_density: float | None = None
    location: tuple[float, float] | None = None
    morphology_type: str | None = None
    age: float | None = None

    def __str__(self):
        if self.age is None:
            agetext = "Not set"
        else:
            agetext = f"{format_large_units(self.age, quantity='time')}"
        return (
            f"final_diameter: {format_large_units(self.final_diameter, quantity='length')}\n"
            f"transient_diameter: {format_large_units(self.transient_diameter, quantity='length')}\n"
            f"projectile_diameter: {format_large_units(self.projectile_diameter, quantity='length')}\n"
            f"projectile_mass: {self.projectile_mass:.4e} kg\n"
            f"projectile_density: {self.projectile_density:.0f} kg/m^3\n"
            f"projectile_velocity: {format_large_units(self.projectile_velocity, quantity='velocity')}\n"
            f"projectile_angle: {self.projectile_angle:.1f}°\n"
            f"projectile_direction: {self.projectile_direction:.1f}°\n"
            f"location (lon,lat): ({self.location[0]:.4f}°, {self.location[1]:.4f}°)\n"
            f"morphology_type: {self.morphology_type}\n"
            f"age: {agetext}"
        )

    @property
    def _id(self):
        """
        The hash id of the crater. This is used as a unique identifier for the crater.
        """
        # Ensure deterministic field order by sorting keys (optional but safer)
        data = asdict(self)
        combined = ":".join(str(data[k]) for k in sorted(data))
        return hashlib.sha256(combined.encode()).hexdigest()

    @property
    def final_radius(self) -> float | None:
        """Final radius of the crater in meters."""
        return self.final_diameter / 2.0 if self.final_diameter is not None else None

    @property
    def transient_radius(self) -> float | None:
        """Transient radius of the crater in meters."""
        return (
            self.transient_diameter / 2.0
            if self.transient_diameter is not None
            else None
        )

    @property
    def projectile_radius(self) -> float | None:
        """Projectile radius in meters."""
        return (
            self.projectile_diameter / 2.0
            if self.projectile_diameter is not None
            else None
        )

    @property
    def projectile_mass(self) -> float | None:
        """Projectile mass in kilograms."""
        if self.projectile_density is not None and self.projectile_radius is not None:
            return (
                (4.0 / 3.0)
                * math.pi
                * self.projectile_radius**3
                * self.projectile_density
            )
        return None

    @property
    def projectile_vertical_velocity(self) -> float | None:
        """Projectile vertical velocity in m/s."""
        if self.projectile_velocity is not None and self.projectile_angle is not None:
            return self.projectile_velocity * math.sin(
                math.radians(self.projectile_angle)
            )
        return None

    @classmethod
    def maker(
        cls: type[Crater],
        crater: Crater | None = None,
        scaling: str | Scaling | None = None,
        target: str | Target | None = None,
        projectile: str | Projectile | None = None,
        final_diameter: float | None = None,
        final_radius: float | None = None,
        transient_diameter: float | None = None,
        transient_radius: float | None = None,
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
        age: float | None = None,
        simdir: str | Path | None = None,
        rng: Generator = None,
        rng_seed: str | int | None = None,
        rng_state: dict | None = None,
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
        final_diameter : float, optional
            The final diameter of the crater in meters.
        final_radius : float, optional
            The final radius of the crater in meters.
        transient_diameter : float, optional
            The transient diameter of the crater in meters.
        transient_radius : float, optional
            The transient radius of the crater in meters.
        projectile_diameter : float, optional
            The diameter of the projectile in meters.
        projectile_radius : float, optional
            The radius of the projectile in meters.
        projectile_mass : float, optional
            The mass of the projectile in kilograms.
        projectile_density : float, optional
            The density of the projectile in kg/m^3. If not provided, the target's density is used.
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
        age : float, optional
            The age of the crater in Myr.
        simdir : str | Path
            The main project simulation directory. Defaults to the current working directory if None.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments for subclasses.

        Returns
        -------
        Crater
            A frozen Crater dataclass with derived attributes.

        Notes
        -----
        - Exactly one of the following must be provided: `final_diameter`, `final_radius`, `transient_diameter`, `transient_radius`, `projectile_diameter`, `projectile_radius`, or `projectile_mass`.
        - Velocity may be specified in one of these ways:
        - `projectile_mean_velocity` alone (samples a velocity)
        - Any two of (`projectile_velocity`, `projectile_vertical_velocity`, `projectile_angle`). the third is inferred.
        - `projectile` is mutually exclusive with velocity-related inputs; if provided, it overrides velocity, angle, direction, and density unless explicitly set.
        - The `scaling`, and `rng` models are required for scaling and density inference, but are not stored in the returned Crater object.
        """
        from ..components.projectile import Projectile
        from ..components.scaling import Scaling
        from ..components.target import Target

        # Validate that mutually exclusive arguments hve not been passed
        size_inputs = {
            "final_diameter": final_diameter,
            "final_radius": final_radius,
            "transient_diameter": transient_diameter,
            "transient_radius": transient_radius,
            "projectile_diameter": projectile_diameter,
            "projectile_radius": projectile_radius,
            "projectile_mass": projectile_mass,
        }

        velocity_inputs = {
            "projectile_velocity": projectile_velocity,
            "projectile_vertical_velocity": projectile_vertical_velocity,
            "projectile_angle": projectile_angle,
        }

        if location is None:
            location = projectile_location

        args = {
            "final_diameter": final_diameter,
            "transient_diameter": transient_diameter,
            "projectile_diameter": projectile_diameter,
            "projectile_density": projectile_density,
            "projectile_velocity": projectile_velocity,
            "projectile_angle": projectile_angle,
            "projectile_direction": projectile_direction,
            "morphology_type": "Not Set",
            "location": location,
            "age": age,
        }

        n_velocity_inputs = sum(x is not None for x in velocity_inputs.values())
        if n_velocity_inputs > 2:
            raise ValueError(
                f"Only two of {', '.join(k for k, v in velocity_inputs.items() if v is not None)} may be set."
            )

        if projectile_mean_velocity is not None:
            if (
                projectile_velocity is not None
                or projectile_vertical_velocity is not None
            ):
                raise ValueError(
                    "projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity"
                )

        n_size_inputs = sum(v is not None for v in size_inputs.values())

        # Process input crater object if provided
        if crater is not None:
            if not isinstance(crater, Crater):
                raise TypeError("crater must be a Crater object.")
            old_parameters = {}
            for field in cls.__dataclass_fields__:
                if field in locals() and locals()[field] is None:
                    old_parameters[field] = getattr(crater, field)
            if (
                n_size_inputs == 0
            ):  # The user has not passed any size parameters, so we will use the final diameter from the crater object
                for field in size_inputs:
                    if field != "final_diameter":
                        old_parameters.pop(field, None)
                n_size_inputs = 1  # Make sure we don't trigger the error below
            else:  # The user is passing a size parameter, so we cannot use any of the values from the crater object
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
            raise ValueError(
                f"Exactly one of {', '.join(k for k, v in size_inputs.items() if v is not None)} must be set."
            )

        # --- Normalize RNG, rng_seed, simdir using CratermakerBase ---
        argproc = CratermakerBase(
            simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state
        )
        if scaling is not None and isinstance(scaling, Scaling):
            if projectile is None:
                projectile = scaling.projectile
            if target is None:
                target = scaling.target

        target = Target.maker(target, **vars(argproc.common_args), **kwargs)

        projectile_args = {
            "velocity": args["projectile_velocity"],
            "density": args["projectile_density"],
            "angle": args["projectile_angle"],
            "direction": args["projectile_direction"],
            "location": args["location"],
        }
        projectile = Projectile.maker(
            projectile,
            target=target,
            **vars(argproc.common_args),
            **kwargs,
        ).new_projectile(
            **projectile_args,
        )

        scaling = Scaling.maker(
            scaling,
            target=target,
            projectile=projectile,
            **vars(argproc.common_args),
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
        fd = args["final_diameter"]
        fr = final_radius
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

        # Assemble final arguments
        args = {
            "final_diameter": float(fd) if fd is not None else None,
            "transient_diameter": float(td) if td is not None else None,
            "projectile_diameter": float(pd) if pd is not None else None,
            "projectile_density": float(prho) if prho is not None else None,
            "projectile_velocity": float(pv) if pv is not None else None,
            "projectile_angle": float(pang) if pang is not None else None,
            "projectile_direction": float(pdir) if pdir is not None else None,
            "location": (float(location[0]), float(location[1])),
            "morphology_type": str(mt) if mt is not None else None,
            "age": float(age) if age is not None else None,
        }
        return cls(**args)
