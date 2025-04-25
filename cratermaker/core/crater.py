import numpy as np
from numpy.random import Generator
from pathlib import Path
from dataclasses import dataclass
from .target import Target
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..components.scaling import ScalingModel, get_scaling_model
from ..components.impactor import ImpactorModel, get_impactor_model
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


def make_crater(final_diameter: float | None = None,
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
                scale: str | ScalingModel = "richardson2009",
                impactor: str | ImpactorModel = "asteroids",
                target: str | Target = "Moon",
                seed: str | int | None = None,
                simdir: str | Path = Path.cwd(),
                rng: Generator = None,
                **kwargs: Any 
                ) -> Crater:
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
    scale : str or ScalingModel, optional
        A string key or instance of a scaling model.
    impactor : str or ImpactorModel, optional
        A string key or instance of an impactor model.
    target : str or Target, optional
        The target body name or object. Used internally, not stored on Crater.
    simdir : str | Path
        The main project simulation directory.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the seed if it is provided.
    seed : int | None
        The random seed for the simulation if rng is not provided. If None, a random seed is used.
    kwargs : Any
        Additional keyword arguments for subclasses.

    Returns
    -------
    Crater
        A frozen Crater dataclass with derived attributes.

    Notes
    -----
    - Exactly one of the following must be provided: 
    `final_diameter`, `final_radius`, `transient_diameter`, `transient_radius`, 
    `projectile_diameter`, `projectile_radius`, or `projectile_mass`.

    - Velocity may be specified in one of these ways:
    - `projectile_mean_velocity` alone (samples a velocity)
    - Any two of (`projectile_velocity`, `projectile_vertical_velocity`, `projectile_angle`)
        — the third is inferred.

    - `impactor` is mutually exclusive with velocity-related inputs; if provided, 
    it overrides velocity, angle, direction, and density unless explicitly set.

    - The `target`, `scale`, and `rng` models are required for scaling and density inference, but are not stored
    in the returned Crater object.
    """
    # --- Normalize RNG, seed, simdir using CratermakerBase ---
    argproc = CratermakerBase(simdir=simdir, rng=rng, seed=seed, **kwargs)
    rng = argproc.rng
    simdir = argproc.simdir
    seed = argproc.seed

    # --- Normalize target ---
    if isinstance(target, str):
        try:
            target = Target(target, **vars(argproc.common_args), **kwargs)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a valid name of a target body")

    # --- Normalize scale ---
    if isinstance(scale, str):
        try:
            scale = get_scaling_model(scale)(target_name=target.name, **vars(argproc.common_args), **kwargs)
        except:
            raise ValueError(f"Invalid scale name {scale}")
    elif not isinstance(scale, ScalingModel):
        raise TypeError("scale must be an instance of ScalingModel or a valid name of a scaling model")
    

    # --- Normalize impactor ---
    if isinstance(impactor, str):
        try:
            impactor = get_impactor_model(impactor)(target_name=target.name, **vars(argproc.common_args), **kwargs)
        except:
            raise ValueError(f"Invalid impactor model name {impactor}")
    elif not isinstance(impactor, ImpactorModel):
        raise TypeError("impactor must be an instance of ImpactorModel or a valid name of an impactor model")
    

    # --- Normalize location and age ---
    if location is None:
        location = mc.get_random_location(rng=rng)
    else:
        location = validate_and_convert_location(location)
    if age is not None:
        age = float(age)

    # --- Handle impactor vs. raw velocity input ---
    pmv = projectile_mean_velocity
    pv = projectile_velocity
    pvv = projectile_vertical_velocity
    pang = projectile_angle
    pdir = projectile_direction
    prho = projectile_density

    if impactor is not None:
        vargs = sum(x is not None for x in [pv, pvv, pmv])
        if vargs > 0:
            raise ValueError("projectile_velocity, projectile_vertical_velocity, and projectile_mean_velocity cannot be used with an impactor model")
        if isinstance(impactor, str):
            try:
                impactor = get_impactor_model(impactor)(target_name=target.name, rng=rng)
            except ValueError as e:
                raise ValueError(f"Invalid impactor name {impactor}") from e
        elif not isinstance(impactor, ImpactorModel):
            raise TypeError("impactor must be an instance of ImpactorModel or a valid name of an impactor model")
        impactor.new_projectile()
        pv = impactor.velocity
        pvv = impactor.vertical_velocity
        if pang is None:
            pang = impactor.angle
        if prho is None:
            prho = impactor.density
        if pdir is None:
            pdir = impactor.direction
    else:
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
        else:
            n_set = sum(x is not None for x in [pv, pvv, pang])
            if n_set == 0:
                impactor = get_impactor_model("asteroids")(target_name=target.name, rng=rng)
                impactor.new_projectile()
                pv = impactor.velocity
                pvv = impactor.vertical_velocity
                pang = impactor.angle
                pdir = impactor.direction
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
        if prho is None:
            prho = target.density
        impactor = ImpactorModel(velocity=pv, angle=pang, density=prho, direction=pdir, sample_velocities=False, sample_angles=False, sample_directions=False, sample_direction=False)

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
        td, mt = scale.final_to_transient(fd)
        pd = scale.transient_to_projectile(td)
    elif td is not None:
        fd, mt = scale.transient_to_final(td)
        pd = scale.transient_to_projectile(td)
    elif pd is not None:
        td = scale.projectile_to_transient(pd)
        fd, mt = scale.transient_to_final(td)
    elif pm is not None:
        pr = ((3.0 * pm) / (4.0 * np.pi * prho)) ** (1.0 / 3.0)
        pd = 2.0 * pr
        td = scale.projectile_to_transient(pd)
        fd, mt = scale.transient_to_final(td)

    pr = pd / 2
    tr = td / 2
    fr = fd / 2
    pm = (4.0 / 3.0) * np.pi * pr**3 * prho

    return Crater(final_diameter=fd,
                  transient_diameter=td,
                  projectile_diameter=pd,
                  projectile_density=prho,
                  projectile_velocity=pv,
                  projectile_angle=pang,
                  projectile_direction=pdir,
                  morphology_type=mt,
                  location=location,
                  age=age)
