from __future__ import annotations
import numpy as np
from numpy.random import Generator
from pathlib import Path
from dataclasses import dataclass
from ..components.target import Target
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..components.scaling import Scaling
from ..components.projectile import Projectile
from typing import Any
from .base import CratermakerBase

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

    def __repr__(self):
        return (f"final_diameter={self.final_diameter} m, "
                f"transient_diameter={self.transient_diameter} m, "
                f"morphology_type={self.morphology_type} "
                f"projectile_diameter={self.projectile_diameter} m "
                f"projectile_mass={self.projectile_mass} kg " 
                f"projectile_density={self.projectile_density} kg/m^3 "
                f"projectile_velocity={self.projectile_velocity} m/s "
                f"projectile_angle={self.projectile_angle} deg, "
                f"projectile_direction={self.projectile_direction} deg, "
                f"lon: {self.location[0]}, lat {self.location[1]} "
                f"age={self.age} My")
    
    @property
    def final_radius(self) -> float | None:
        """Final radius of the crater in meters."""
        return self.final_diameter / 2.0 if self.final_diameter is not None else None
    @property
    def transient_radius(self) -> float | None:
        """Transient radius of the crater in meters."""
        return self.transient_diameter / 2.0 if self.transient_diameter is not None else None
    @property
    def projectile_radius(self) -> float | None:
        """Projectile radius in meters."""
        return self.projectile_diameter / 2.0 if self.projectile_diameter is not None else None
    @property
    def projectile_mass(self) -> float | None:
        """Projectile mass in kilograms."""
        if self.projectile_density is not None and self.projectile_radius is not None:
            return (4.0 / 3.0) * np.pi * self.projectile_radius**3 * self.projectile_density
        return None
    @property
    def projectile_vertical_velocity(self) -> float | None:
        """Projectile vertical velocity in m/s."""
        if self.projectile_velocity is not None and self.projectile_angle is not None:
            return self.projectile_velocity * np.sin(np.deg2rad(self.projectile_angle))
        return None

    @classmethod
    def maker(cls : type[Crater],
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
            location: tuple[float, float] | None = None,
            age: float | None = None,
            scaling: str | Scaling = "richardson2009",
            projectile: str | Projectile = "asteroids",
            target: str | Target = "Moon",
            simdir: str | Path | None = None,
            rng: Generator = None,
            rng_seed: str | int | None = None,
            rng_state: dict | None = None,
            **kwargs: Any): 
        """
        Create a Crater object with the given parameters.

        Parameters
        ----------
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
            The impact angle in degrees (0–90).
        projectile_direction : float, optional
            The direction of the impact in degrees (0–360).
        location : tuple of float, optional
            The (longitude, latitude) location of the impact.
        age : float, optional
            The age of the crater in Myr.
        scaling : str or Scaling, optional
            A string key or instance of a scaling model.
        projectile : str or Projectile, optional
            A string key or instance of an projectile model.
        target : str or Target, optional
            The target body name or object. Used internally, not stored on Crater.
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
        - The `target`, `scaling`, and `rng` models are required for scaling and density inference, but are not stored in the returned Crater object.

        """
        # --- Normalize RNG, rng_seed, simdir using CratermakerBase ---
        argproc = CratermakerBase(simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state)
        rng = argproc.rng

        target = Target.maker(target, **vars(argproc.common_args), **kwargs)
        projectile = Projectile.maker(projectile, target=target, **vars(argproc.common_args), **kwargs)

        # --- Normalize location and age ---
        if location is None:
            location = mc.get_random_location(rng=rng)
        else:
            location = validate_and_convert_location(location)
        if age is not None:
            age = float(age)

        # --- Handle projectile vs. raw velocity input ---
        pmv = projectile_mean_velocity
        pv = projectile_velocity
        pvv = projectile_vertical_velocity
        pang = projectile_angle
        pdir = projectile_direction
        prho = projectile_density


        # --- Resolve velocity input combinations ---
        if pmv is not None:
            if pv is not None or pvv is not None:
                raise ValueError("projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity")
            pmv = float(pmv)
            if pmv <= 0.0:
                raise ValueError("projectile_mean_velocity must be positive.")
            if pmv > target.escape_velocity:
                vencounter_mean = np.sqrt(pmv ** 2 - target.escape_velocity ** 2)
                vencounter = mc.get_random_velocity(vencounter_mean, rng=rng)
                pv = float(np.sqrt(vencounter ** 2 + target.escape_velocity ** 2))
            else:
                while True:
                    pv = float(mc.get_random_velocity(pmv, rng=rng))
                    if pv < target.escape_velocity:
                        break
        n_set = sum(x is not None for x in [pv, pvv, pang])
        if n_set == 0:
            projectile.new_projectile()
            pv = projectile.velocity
            pvv = projectile.vertical_velocity
            pang = projectile.angle
            pdir = projectile.direction
            prho = prho or projectile.density
        elif n_set > 2:
            raise ValueError("Only two of projectile_velocity, projectile_vertical_velocity, projectile_angle may be set")
        else:
            if pv is not None:
                pv = float(pv)
                if pv <= 0.0:
                    raise ValueError("projectile_velocity must be positive.")
            if pvv is not None:
                pvv = float(pvv)
                if pvv <= 0.0:
                    raise ValueError("projectile_vertical_velocity must be positive.")
            if pang is not None:
                pang = float(pang)
                if not (0.0 <= pang <= 90.0):
                    raise ValueError("projectile_angle must be between 0 and 90 degrees")
            if pv is not None and pang is not None:
                pvv = pv * np.sin(np.deg2rad(pang))
            elif pvv is not None and pang is not None:
                pv = pvv / np.sin(np.deg2rad(pang))
            elif pv is not None and pvv is not None:
                pang = np.rad2deg(np.arcsin(pvv / pv))
            elif pv is not None and pang is None:
                pang = mc.get_random_impact_angle(rng=rng)
                pvv = pv * np.sin(np.deg2rad(pang))
            elif pvv is not None and pang is None:
                pang = mc.get_random_impact_angle(rng=rng)
                pv = pvv / np.sin(np.deg2rad(pang))
            # Direction
            if pdir is None:
                pdir = float(rng.uniform(0.0, 360.0))
            else:
                pdir = float(pdir) % 360.0
            # Get or infer projectile density
            prho = prho or target.density
            projectile = Projectile(velocity=pv, angle=pang, density=prho, direction=pdir, sample_velocities=False, sample_angles=False, sample_directions=False, sample_direction=False)

        scaling = Scaling.maker(scaling, target=target, projectile=projectile, **vars(argproc.common_args), **kwargs)
        prho = scaling.projectile.density

        # --- Ensure velocity/angle are all set ---
        n_set = sum(x is not None for x in [pv, pvv, pang])
        if n_set != 3:
            raise ValueError("Not enough information to infer a projectile velocity.")

        # --- Resolve projectile size/mass inputs ---
        size_inputs = {
            "final_diameter": final_diameter,
            "final_radius": final_radius,
            "transient_diameter": transient_diameter,
            "transient_radius": transient_radius,
            "projectile_diameter": projectile_diameter,
            "projectile_radius": projectile_radius,
            "projectile_mass": projectile_mass
        }
        n_set = sum(v is not None for v in size_inputs.values())
        if n_set != 1:
            raise ValueError("Exactly one of final_diameter, final_radius, transient_diameter, transient_radius, projectile_diameter, projectile_radius, or projectile_mass must be set.")

        # --- Compute derived quantities ---
        fd = final_diameter
        fr = final_radius
        td = transient_diameter
        tr = transient_radius
        pd = projectile_diameter
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
            pr = ((3.0 * pm) / (4.0 * np.pi * prho)) ** (1.0 / 3.0)
            pd = 2.0 * pr
            td = scaling.projectile_to_transient(pd)
            fd, mt = scaling.transient_to_final(td)

        pr = pd / 2
        tr = td / 2
        fr = fd / 2
        pm = (4.0 / 3.0) * np.pi * pr**3 * prho

        return cls(final_diameter=fd,
                  transient_diameter=td,
                  projectile_diameter=pd,
                  projectile_density=prho,
                  projectile_velocity=pv,
                  projectile_angle=pang,
                  projectile_direction=pdir,
                  morphology_type=mt,
                  location=location,
                  age=age)
