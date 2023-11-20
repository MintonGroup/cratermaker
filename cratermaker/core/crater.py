import numpy as np
from numpy.random import Generator
from dataclasses import dataclass, field
from ..models.production_function import ProductionFunction
from ..models import crater_scaling 
from ..core.target import Target
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
    morphology_type: str = field(default=None)
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.production, self.diameter, self.radius, 
                                                 self.transient_diameter, self.transient_radius])
        if values_set > 1:
            raise ValueError("Only one of production, diameter, radius, transient_diameter, transient_radius may be set")
                
        if self.location is not None:
            if not isinstance(self.location, np.ndarray):
                self.location = np.array(self.location, dtype=np.float64)
            if self.location.shape != (2,):
                raise ValueError("location must be a 2-element array")
        
        return
    

    def final_to_transient(self, target: Target, rng: Generator=None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        if rng is None:
            rng = np.random.default_rng() 
        self.transient_diameter = crater_scaling.final_to_transient(self.diameter,target,rng)
        self.transient_radius = self.transient_diameter / 2
        return 


    def transient_to_final(self, target: Target, rng: Generator=None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        if rng is None:
            rng = np.random.default_rng()
        self.diameter, self.morphology_type = crater_scaling.transient_to_final(self.transient_diameter,target,rng)    
        self.radius = self.diameter / 2
        return 


     
   
