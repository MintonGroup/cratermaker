import numpy as np
from numpy.random import Generator
from dataclasses import dataclass, field
from .target import Target
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..components.scaling import ScalingModel, get_scaling_model
from ..components.impactor import ImpactorModel, get_impactor_model


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
        return (f"Crater(final_diameter={self.final_diameter} m, "
                f"transient_diameter={self.transient_diameter} m, "
                f"morphology_type={self.morphology_type} "
                f"projectile_diameter={self.projectile_diameter} m, "
                f"projectile_mass={self.projectile_mass} kg, projectile_density={self.projectile_density} kg/m^3, "
                f"projectile_velocity={self.projectile_velocity} m/s, projectile_angle={self.projectile_angle} deg, "
                f"projectile_direction={self.projectile_direction} deg, "
                f"lon: {self.location[0]}, lat {self.location[1]} age={self.age} My)")
    
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
                scale: str | ScalingModel | None = None,
                impactor: str | ImpactorModel | None = None,
                target: str | Target = "Moon",
                rng: Generator = np.random.default_rng(),
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
    rng : numpy.random.Generator, optional
        Random number generator.

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

    if isinstance(target, str):
        try:
            target = Target(target)
        except:
            raise ValueError(f"Invalid target name {target}")
    elif not isinstance(target, Target):
        raise TypeError("target must be an instance of Target or a valid name of a target body")

    # Validate location
    if location is None:
        location = mc.get_random_location(rng=rng)
    else:
        location = validate_and_convert_location(location)

    # Validate age
    if age is not None:
        age = float(age)


    # Handle velocities and angles
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

        # Only override with values from the impactor model if they are not already set
        pv = impactor.velocity
        pvv = impactor.vertical_velocity
        if pang is None:
            pang = impactor.angle
        if prho is None:
            prho = impactor.density
        if pdir is None:
            pdir = impactor.direction
    else:
        if pmv is not None:
            if pv is not None or pvv is not None:
                raise ValueError("projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity")
            pmv = float(pmv)
            if pmv <= 0.0:
                raise ValueError("projectile_mean_velocity must be positive.")
            # Sample velocity from mean
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
            # Only two of pv, pvv, pang can be set
            n_set = sum(x is not None for x in [pv, pvv, pang])
            if n_set > 2:
                raise ValueError("Only two of projectile_velocity, projectile_vertical_velocity, projectile_angle may be set")
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

            # Infer missing velocity/angle
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

        # Create a basic impactor model
        impactor = ImpactorModel(velocity=pv, angle=pang, density=prho, direction=pdir, sample_velocities=False, sample_angles=False, sample_directions=False, sample_direction=False)

    # Ensure that none of the values of pv, pvv, and pang are still None
    n_set = sum(x is not None for x in [pv, pvv, pang])
    if n_set != 3:
        raise ValueError("Not enough information to infer a projectile velocity.")

    # Validate and resolve which size/mass is given
    fd = final_diameter
    fr = final_radius
    td = transient_diameter
    tr = transient_radius
    pd = projectile_diameter
    pr = projectile_radius
    pm = projectile_mass
    mt = None
    # Only one of fd, td, pd, pm can be set
    n_set = sum(x is not None for x in [fd, fr, td, tr, pd, pr, pm])
    if n_set != 1:
        raise ValueError("Exactly one of final_diameter, final_radius, transient_diameter, transient_radius, projectile_diameter, projectile_radius, or projectile_mass must be set.")

    # Compute all related sizes
    if pr is not None:
        pd = 2 * pr
    if fr is not None:
        fd = 2 * fr
    if tr is not None:
        td = 2 * tr

    if fd is not None:
        # Final diameter is given
        fd = float(fd)
        if fd <= 0.0:
            raise ValueError("final_diameter must be positive.")
        fr = fd / 2.0
        td, mt = scale.final_to_transient(fd)
        tr = td / 2.0
        pd = scale.transient_to_projectile(td)
        pr = pd / 2.0
        pm = 4.0 / 3.0 * np.pi * pr**3 * prho
    elif td is not None:
        # Transient diameter is given
        td = float(td)
        if td <= 0.0:
            raise ValueError("transient_diameter must be positive.")
        tr = td / 2.0
        fd, mt = scale.transient_to_final(td)
        fr = fd / 2.0
        pd = scale.transient_to_projectile(td)
        pr = pd / 2.0
        pm = 4.0 / 3.0 * np.pi * pr**3 * prho
    elif pd is not None:
        # Projectile diameter is given
        pd = float(pd)
        if pd <= 0.0:
            raise ValueError("projectile_diameter must be positive.")
        pr = pd / 2.0
        pm = 4.0 / 3.0 * np.pi * pr**3 * prho
        td = scale.projectile_to_transient(pd)
        tr = td / 2.0
        fd, mt = scale.transient_to_final(td)
        fr = fd / 2.0
    elif pm is not None:
        # Projectile mass is given
        pm = float(pm)
        if pm <= 0.0:
            raise ValueError("projectile_mass must be positive.")
        pr = ((3.0 * pm) / (4.0 * np.pi * prho)) ** (1.0 / 3.0)
        pd = 2.0 * pr
        td = scale.projectile_to_transient(pd)
        tr = td / 2.0
        fd, mt = scale.transient_to_final(td)
        fr = fd / 2.0
    else:
        raise RuntimeError("Failed to infer crater/projectile properties.")

    return Crater(final_diameter=fd, 
                  transient_diameter=td,
                  projectile_diameter=pd,
                  projectile_density=prho,
                  projectile_velocity=pv,
                  projectile_angle=pang,
                  morphology_type=mt,
                  location=location,
                    age=age)
