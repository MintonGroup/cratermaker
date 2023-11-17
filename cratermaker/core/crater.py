@dataclass
class Crater:
    """
    Represents a crater formed by an impact in the simulation.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features.

    Attributes
    ----------
    TBD
    """
    production: ProductionFunction = field(default=None)
    diameter: np.float64 = field(default=None)
    radius: np.float64 = field(default=None)
    transient_diameter: np.float64 = field(default=None)
    transient_radius: np.float64 = field(default=None)
    location: np.ndarray = field(default=None)
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.production, self.diameter, self.radius, 
                                                 self.transient_diameter, self.transient_radius])

        if values_set > 1:
            raise ValueError("Only one of production, diameter, radius, transient_diameter, transient_radius may be set")

        if self.diameter is not None:
            self.radius = self.diameter / 2
        elif self.radius is not None:
            self.diameter = self.radius * 2
        elif self.transient_diameter is not None:
            self.transient_radius = self.transient_diameter / 2
        elif self.transient_radius is not None:
            self.transient_diameter = self.transient_radius * 2

        if self.location is not None:
            if not isinstance(self.location, np.ndarray):
                self.location = np.array(self.location, dtype=np.float64)
            if self.location.shape != (2,):
                raise ValueError("location must be a 2-element array")
        
        return
   
