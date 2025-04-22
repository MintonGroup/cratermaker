import numpy as np
from numpy.random import Generator
from dataclasses import dataclass, field
from .target import Target
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..components.scaling import ScalingModel, get_scaling_model

@dataclass(frozen=True, slots=True)
class Crater:
    final_diameter: float | None = None
    transient_diameter: float | None = None
    projectile_diameter: float | None = None
    projectile_mass: float | None = None
    projectile_density: float | None = None
    projectile_velocity: float | None = None
    projectile_vertical_velocity: float | None = None
    projectile_mean_velocity: float | None = None
    projectile_angle: float | None = None
    projectile_direction: float | None = None
    location: tuple[float, float] | None = None
    age: float | None = None
    target: Target = field(default_factory=lambda: Target(name="Moon"))
    scale: str | ScalingModel | None = None
    rng: Generator = field(default_factory=np.random.default_rng)
    # computed fields
    final_radius: float = field(init=False)
    transient_radius: float = field(init=False)
    projectile_radius: float = field(init=False)
    morphology_type: str = field(init=False)

    def __post_init__(self):
        # Validate location
        loc = self.location
        if loc is None:
            loc = mc.get_random_location(rng=self.rng)
        else:
            loc = validate_and_convert_location(loc)
        object.__setattr__(self, 'location', loc)

        # Validate age
        age = self.age
        if age is not None:
            age = float(age)
        object.__setattr__(self, 'age', age)

        # Handle velocities and angles
        pmv = self.projectile_mean_velocity
        pv = self.projectile_velocity
        pvv = self.projectile_vertical_velocity
        pang = self.projectile_angle
        pdir = self.projectile_direction
        # Validate velocity/angle input
        if pmv is not None:
            if pv is not None or pvv is not None:
                raise ValueError("projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity")
            pmv = float(pmv)
            if pmv <= 0.0:
                raise ValueError("projectile_mean_velocity must be positive.")
            # Sample velocity from mean
            if pmv > self.target.escape_velocity:
                vencounter_mean = np.sqrt(pmv ** 2 - self.target.escape_velocity ** 2)
                vencounter = mc.get_random_velocity(vencounter_mean, rng=self.rng)
                pv = float(np.sqrt(vencounter ** 2 + self.target.escape_velocity ** 2))
            else:
                while True:
                    pv = float(mc.get_random_velocity(pmv, rng=self.rng))
                    if pv < self.target.escape_velocity:
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
                pang = mc.get_random_impact_angle(rng=self.rng)
                pvv = pv * np.sin(np.deg2rad(pang))
            elif pvv is not None and pang is None:
                pang = mc.get_random_impact_angle(rng=self.rng)
                pv = pvv / np.sin(np.deg2rad(pang))

        # Direction
        if pdir is None:
            pdir = float(self.rng.uniform(0.0, 360.0))
        else:
            pdir = float(pdir) % 360.0

        # Get or infer projectile density
        prho = self.projectile_density

        # Set up the scaling model
        match self.scale:
            case None:
                scale_obj = get_scaling_model("richardson2009")(target=self.target, rng=self.rng, projectile_density=prho, projectile_vertical_velocity=pvv)
            case str() as name:
                scale_obj = get_scaling_model(name)(target=self.target, rng=self.rng, projectile_density=prho, projectile_vertical_velocity=pvv)
            case ScalingModel():
                scale_obj = self.scale
            case _:
                raise TypeError("scale must be None, a ScalingModel instance, or a string model name")
            
        if prho is None:
            if hasattr(scale_obj, 'projectile_density') and scale_obj.projectile_density is not None:
                prho = scale_obj.projectile_density
            else:
                raise ValueError("Not enough information to infer a projectile density.")

        if pvv is None: 
            if hasattr(scale_obj, 'projectile_vertical_velocity') and scale_obj.projectile_vertical_velocity is not None:
                pvv = scale_obj.projectile_vertical_velocity
            else:
                raise ValueError("Not enough information to infer a projectile velocity.")

        if pang is None:
            if pv is not None:
                pang = np.rad2deg(np.arcsin(pvv / pv))
            else:
                pang = mc.get_random_impact_angle(rng=self.rng)
                pv = pvv / np.sin(np.deg2rad(pang))

        # Ensure that none of the values of pv, pvv, and pang are still None
        n_set = sum(x is not None for x in [pv, pvv, pang])
        if n_set != 3:
            raise ValueError("Not enough information to infer a projectile velocity.")

        object.__setattr__(self, 'scale', scale_obj)

        # Validate and resolve which size/mass is given
        fd = self.final_diameter
        td = self.transient_diameter
        pd = self.projectile_diameter
        pm = self.projectile_mass
        pr = None
        tr = None
        fr = None
        mt = None
        # Only one of fd, td, pd, pm can be set
        n_set = sum(x is not None for x in [fd, td, pd, pm])
        if n_set != 1:
            raise ValueError("Exactly one of final_diameter, transient_diameter, projectile_diameter, or projectile_mass must be set.")

        # Compute all related sizes
        if fd is not None:
            # Final diameter is given
            fd = float(fd)
            if fd <= 0.0:
                raise ValueError("final_diameter must be positive.")
            fr = fd / 2.0
            td, mt = self.scale.final_to_transient(fd)
            tr = td / 2.0
            pd = self.scale.transient_to_projectile(td)
            pr = pd / 2.0
            pm = 4.0 / 3.0 * np.pi * pr**3 * prho
        elif td is not None:
            # Transient diameter is given
            td = float(td)
            if td <= 0.0:
                raise ValueError("transient_diameter must be positive.")
            tr = td / 2.0
            fd, mt = self.scale.transient_to_final(td)
            fr = fd / 2.0
            pd = self.scale.transient_to_projectile(td)
            pr = pd / 2.0
            pm = 4.0 / 3.0 * np.pi * pr**3 * prho
        elif pd is not None:
            # Projectile diameter is given
            pd = float(pd)
            if pd <= 0.0:
                raise ValueError("projectile_diameter must be positive.")
            pr = pd / 2.0
            pm = 4.0 / 3.0 * np.pi * pr**3 * prho
            td = self.scale.projectile_to_transient(pd)
            tr = td / 2.0
            fd, mt = self.scale.transient_to_final(td)
            fr = fd / 2.0
        elif pm is not None:
            # Projectile mass is given
            pm = float(pm)
            if pm <= 0.0:
                raise ValueError("projectile_mass must be positive.")
            pr = ((3.0 * pm) / (4.0 * np.pi * prho)) ** (1.0 / 3.0)
            pd = 2.0 * pr
            td = self.scale.projectile_to_transient(pd)
            tr = td / 2.0
            fd, mt = self.scale.transient_to_final(td)
            fr = fd / 2.0
        else:
            raise RuntimeError("Failed to infer crater/projectile properties.")

        object.__setattr__(self, 'final_diameter', float(fd))
        object.__setattr__(self, 'final_radius', float(fr))
        object.__setattr__(self, 'transient_diameter', float(td))
        object.__setattr__(self, 'transient_radius', float(tr))
        object.__setattr__(self, 'projectile_diameter', float(pd))
        object.__setattr__(self, 'projectile_radius', float(pr))
        object.__setattr__(self, 'projectile_mass', float(pm))
        object.__setattr__(self, 'projectile_density', float(prho))
        object.__setattr__(self, 'morphology_type', str(mt))
        object.__setattr__(self, 'projectile_velocity', float(pv))
        object.__setattr__(self, 'projectile_vertical_velocity', float(pvv))
        object.__setattr__(self, 'projectile_angle', float(pang))
        object.__setattr__(self, 'projectile_direction', float(pdir))

    def __repr__(self):
        return (f"Crater(final_diameter={self.final_diameter} m, "
                f"transient_diameter={self.transient_diameter} m, "
                f"morphology_type={self.morphology_type} "
                f"projectile_diameter={self.projectile_diameter} m, "
                f"projectile_mass={self.projectile_mass} kg, projectile_density={self.projectile_density} kg/m^3, "
                f"projectile_velocity={self.projectile_velocity} m/s, projectile_angle={self.projectile_angle} deg, "
                f"projectile_vertical_velocity={self.projectile_vertical_velocity} m/s, projectile_direction={self.projectile_direction} deg, "
                f"lon: {self.location[0]}, lat {self.location[1]} age={self.age} My)")