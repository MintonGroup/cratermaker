@dataclass    
class Projectile:
    """
    Represents the projectile in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    Attributes
    ----------
    velocity : float
        The velocity of the projectile upon impact, in m/s.
    angle : float
        The angle of impact, in degrees.
    material : Material
        The material composition of the projectile. 
    """
    production: ProductionFunction = field(default=None)
    diameter: np.float64 = field(default=None)
    radius: np.float64 = field(default=None)
    velocity: np.float64 = field(default=None)
    density: np.float64 = field(default=None)
    location: np.ndarray = field(default=None)    
    angle: np.float64 = field(default=None)
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.production, self.diameter, self.radius])

        if values_set > 1:
            raise ValueError("Only one of production, diameter, radius may be set")

        if self.diameter is not None:
            self.radius = self.diameter / 2
        elif self.radius is not None:
            self.diameter = self.radius * 2

        if self.location is not None:
            if not isinstance(self.location, np.ndarray):
                self.location = np.array(self.location, dtype=np.float64)
            if self.location.shape != (2,):
                raise ValueError("location must be a 2-element array")
        
        return    