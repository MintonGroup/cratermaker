import numpy as np
from numpy.random import Generator
from dataclasses import dataclass, field
from .target import Target
from .crater import Crater
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc

@dataclass    
class Projectile:
    """
    Represents the self.in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    Attributes
    ----------
    diameter : float
        The diameter of the projectile in m.
    radius : float
        The radius of the projectile in m.
    velocity : float
        The velocity of the projectile upon impact, in m/s.
    vertical_velocity : float
        The vertical component of the projectile velocity upon impact, in m/s.
    angle : float
        The angle of impact, in degrees.
    location : (2,) float
        The lat. and lon. of the impact point    
    """
    
    _diameter: np.float64 = field(default=None, repr=False)
    _radius: np.float64 = field(default=None, repr=False)
    _density: np.float64 = field(default=None, repr=False)
    _mass: np.float64 = field(default=None, repr=False)
    _velocity: np.float64 = field(default=None, repr=False)
    _vertical_velocity: np.float64 = field(default=None, repr=False)
    _location: np.ndarray = field(default=None, repr=False)
    
    def __post_init__(self, target: Target, rng: Generator = None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        
        # Evaluate and check diameter/radius inputs
        values_set = sum(x is not None for x in [self.diameter, self.radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius may be set")
        
        # Evaluate and check mass/density/radius inputs 
        values_set = sum(x is not None for x in [self.mass, self.density, self.radius])
        if values_set > 2:
            raise ValueError("Only two of mass, density, radius/diameter may be set")
        if self.density is None:
            if self.mass is None or self.radius is  None: # Default to target density if we are given no way to figure it out
                self.density = target.material.density 
                
        # Evaluate and check velocity/impact angle inputs
        values_set = sum(x is not None for x in [self.velocity, self.vertical_velocity, self.angle])
        if values_set > 2:
            raise ValueError("Only two of velocity, vertical_velocity, angle may be set")
            
        self._initialize_velocities(target,rng)
        self._initialize_location(rng)
        
        return
    
    
    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of projectile must be finite and positive!")
            self._radius = value / 2
            self._update_mass()
        return

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of projectile must be finite and positive!")            
            self._diameter = value * 2
            self._update_mass()
        return

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        if value is not None and value <= 0.0:
            raise ValueError("Mass of projectile must be finite and positive!")            
        self._update_volume_based_properties()
        return

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value
        if value is not None and value <= 0.0:
            raise ValueError("Density of projectile must be finite and positive!")                 
        self._update_mass()

    def _update_mass(self):
        if self._density is not None and self._radius is not None:
            self._mass = 4.0/3.0 * np.pi * self._radius**3 * self._density

    def _update_volume_based_properties(self):
        if self._mass is not None and self._density is not None:
            volume = self._mass / self._density
            self._radius = ((3.0 * volume) / (4.0 * np.pi))**(1.0/3.0)
            self._diameter = self._radius * 2

    @property
    def velocity(self):
        if self._velocity is None and self._vertical_velocity is not None and self._angle is not None:
            self._velocity = self._vertical_velocity / np.sin(self._angle)
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Velocity of projectile must be finite and positive!")        
        self._velocity = value
            
        if value is not None and self._angle is not None:
            self._vertical_velocity = value * np.sin(self._angle)

    @property
    def vertical_velocity(self):
        if self._vertical_velocity is None and self._velocity is not None and self._angle is not None:
            self._vertical_velocity = self._velocity * np.sin(self._angle)
        return self._vertical_velocity

    @vertical_velocity.setter
    def vertical_velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Vertical velocity of projectile must be finite and positive!")        
        self._vertical_velocity = value
        
        # Update velocity only if angle is already set
        if self._angle is not None and value is not None:
            try:
                self._velocity = value / np.sin(self._angle)
            except ValueError:
                raise ValueError(f"Invalid vertical velocity value {value} for a given angle value {np.rad2deg(self._angle)}!")

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        if value is not None:
            if value < 0.0 or value > 90.0:
                raise ValueError("Impact angle of projectile must be between 0 and 90 degrees")
            self._angle = np.deg2rad(value)
            # Update vertical_velocity only if velocity is already set
            if self._velocity is not None:
                self._vertical_velocity = self._velocity * np.sin(self._angle)

    def _initialize_velocities(self, target, rng):
        if self._velocity is None:
            vencounter_mean = np.sqrt(target.mean_impact_velocity**2 - target.escape_velocity**2)
            vencounter = mc.get_random_velocity(vencounter_mean)
            self.velocity = np.sqrt(vencounter**2 + target.escape_velocity**2)

        if self._angle is None:
            self.angle = mc.get_random_impact_angle(rng)
        return


    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, value):
        self._location = value
    
    
    def _initialize_location(self, rng):
        if self._location is None:
            self.location = mc.get_random_location(rng=rng)
        else:    
            self.location = validate_and_convert_location(self.location)


class ProjectileScaling:
    """
    A class for handling the scaling relationships between impactors and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties, 
    as well as determining crater morphology based on size and target properties.

    Parameters
    ----------
    target : Target
        The target body for the impact simulation.
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    """

    def __init__(self, target: Target, crater: Crater, rng: Generator = None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self.target = target
        self.rng = rng if rng else np.random.default_rng()
        
        return
        
    def projectile_to_crater(self, projectile: Projectile, target: Target) -> Crater:
        """
        Convert a projectile to its corresponding crater.

        Parameters
        ----------
        projectile : Projectile
            The projectile to be converted.
        target : Target
            The target body being impacted
        Returns
        -------
        Crater
            The crater resulting from the impact of the projectile.
        """
        transient_diameter = projectile_to_transient(profjectile, target)
        crater = Crater(transient_diameter=transient_diameter)

        return crater


    def crater_to_projectile(self, crater: Crater) -> Projectile:
        """
        Convert a crater back to its corresponding projectile.
        This operation is more hypothetical and approximates the possible projectile that created the crater.

        Parameters
        ----------
        crater : Crater
            The crater to be converted.

        Returns
        -------
        Projectile
            The estimated projectile that could have caused the crater.
        """
         
        
        return projectile

    @staticmethod
    def projectile_to_transient(projectile: Projectile, target: Target) -> np.float64:
        if not isinstance(projectile, Projectile):
            raise TypeError("target must be an instance of Projectile")
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
            
        # Compute some auxiliary quantites
        projectile.mass = 4.0/3.0 * np.pi * projectile.density * (projectile.radius)**3
        mu = target.material.mu
        kv = target.material.K1
        c1 = 1.0 + 0.5 * mu
        c2 = (-3 * mu)/(2.0 + mu)

        # Find dimensionless quantities
        pitwo = (target.gravity * projectile.radius)/(projectile.vertical_velocity**2)
        pithree = target.material.Ybar / (target.material.density * (projectile.vertical_velocity**2))
        pifour = target.material.density / projectile.density
        pivol = kv * ((pitwo * (pifour**(-1.0/3.0))) + (pithree**c1))**c2
        pivolg = kv * (pitwo * (pifour**(-1.0/3.0)))**c2
        
        # find transient crater volume and radii (depth = 1/3 diameter)
        cvol = pivol * (projectile.mass / target.material.density)
        cvolg = pivolg * (projectile.mass / target.material.density)
        transient_radius = (3 * cvol / np.pi)**(1.0/3.0)
        transient_radius_gravscale = (3 * cvolg / np.pi)**(1.0/3.0)
        
        # TODO: Compute whether or not to use gravity scaling or strength scaling
        transient_diameter = transient_radius * 2
        
        return transient_diameter


    @staticmethod
    def transient_to_projectile(crater: Crater, target: Target, rng: Generator = None) -> Projectile:
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
                
        
        # We'll create a Projectile object that will allow us to set velocity
        projectile = Projectile(diameter=crater.transient_diameter, target=target, location=crater.location, rng=rng)
        
        def root_func(projectile_diameter: float_like, 
                      projectile: Projectile, 
                      crater: Crater,
                      target: Target) -> np.float64:
            
            projectile.diameter = projectile_diameter
            transient_diameter = self.projectile_to_transient(projectile, target)
            return transient_diameter - crater.transient_diameter 
        
        sol = root_scalar(lambda x, *args: root_func(x, *args),bracket=(1e-5*crater.transient_diameter,crater.transient_diameter), args=(projectile, crater, target))
        
        # Regenerate the projectile with the new diameter value
        projectile = Projectile(diameter=sol.root, target=target, location=projectile.location, velocity=projectile.velocity, angle=projectile.angle, rng=rng)
        
        return projectile