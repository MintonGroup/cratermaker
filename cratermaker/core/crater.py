import numpy as np
from numpy.random import Generator
from dataclasses import dataclass, field
from ..models.production_function import ProductionFunction
from ..models import craterscaling 
from ..core.target import Target
from ..utils.general_utils import validate_and_convert_location
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
    diameter: np.float64 = field(default=None)
    radius: np.float64 = field(default=None)
    transient_diameter: np.float64 = field(default=None)
    transient_radius: np.float64 = field(default=None)
    location: np.ndarray = field(default=None)
    _morphology_type: str = field(default=None) # This will be computed and is not meant to be set by the user
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.diameter, self.radius, 
                                                 self.transient_diameter, self.transient_radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius, transient_diameter, transient_radius may be set")
               
        if self.location is not None:
            self.location = validate_and_convert_location(self.location)
        
        return
    
    @property
    def morphology_type(self):
        return self._morphology_type    
    
    
    def set_morphology_type(self, target: Target, rng: Generator=None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        if rng is None:
            rng = np.random.default_rng()    
        
        _, transition_nominal, *_ = craterscaling.get_simple_to_complex_transition_factors(target,rng)
        # Use the 1/2x to 2x the nominal value of the simple->complex transition diameter to get the range of the "transitional" morphology type. This is supported by: Schenk et al. (2004) and Pike (1980) in particular  
        transition_range = (0.5*transition_nominal,2*transition_nominal)
        
        if self.diameter < transition_range[0]:
            self._morphology_type = "simple" 
        elif self.diameter > transition_range[1]:
            self._morphology_type = "complex"
        else:
            # We'll uses the distance from the nominal transition diameter to set a probability of being either simple, complex, or transitional.
            if self.diameter < transition_nominal:
                p = (transition_nominal - self.diameter)/(transition_nominal - transition_range[0])
                categories = ["simple","transitional"]
                prob = [p, 1.0-p] 
                self._morphology_type = rng.choice(categories,p=prob)
            else:
                p = (self.diameter - transition_nominal)/(transition_range[1] - transition_nominal)
                categories = ["complex","transitional"]
                prob = [p, 1.0-p] 
                self._morphology_type = rng.choice(categories,p=prob)                
        return