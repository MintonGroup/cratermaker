from dataclasses import dataclass, field
import numpy as np
from ..models.production_function import ProductionFunction
from ..utils.general_utils import validate_and_convert_location
@dataclass    
class Projectile:
    """
    Represents the projectile in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    Attributes
    ----------
    diameter : float
        The diameter of the projectile, in m.
    radius : float
        The radius of the projectile, in m.
    velocity : float
        The velocity of the projectile upon impact, in m/s.
    vertical_velocity : float
        The vertical component of the projectile velocity upon impact, in m/s.
    angle : float
        The angle of impact, in degrees.
    location : (2,) float
        The lat. and lon. of the impact point    
    """
    diameter: np.float64 = field(default=None)
    radius: np.float64 = field(default=None)
    density: np.float64 = field(default=None)
    mass: np.float64 = field(default=None)
    velocity: np.float64 = field(default=None)
    vertical_velocity: np.float64 = field(default=None)
    angle: np.float64 = field(default=None)
    location: np.ndarray = field(default=None)    
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.diameter, self.radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius may be set")
        if self.diameter is not None:
            self.radius = self.diameter / 2
        elif self.radius is not None:
            self.diameter = self.radius * 2
            
        values_set = sum(x is not None for x in [self.mass, self.density, self.radius])
        if values_set > 2:
            raise ValueError("Only two of mass, density, radius/diameter may be set")

        values_set = sum(x is not None for x in [self.velocity, self.vertical_velocity, self.angle])
        if values_set > 2:
            raise ValueError("Only two of velocity, vertical_velocity, angle may be set")
        if self.velocity is not None and self.velocity <= 0.0:
            raise ValueError("velocity must be positive!")
        if self.vertical_velocity is not None and self.vertical_velocity <= 0.0:
            raise ValueError("vertical_velocity must be positive!")
        if self.vertical_velocity is not None and self.velocity is not None:
            if self.vertical_velocity > self.velocity:
                raise ValueError("Vertical component of velocity must be less than or equal to velocity!")

        if self.location is not None:
            self.location = validate_and_convert_location(self.location)
        
        
        return
