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
            raise ValueError("Only one of production, diameter, radius, transient_diameter, transient_radius may be set")
                
        if self.location is not None:
            if not isinstance(self.location, np.ndarray):
                self.location = np.array(self.location, dtype=np.float64)
            if self.location.shape != (2,):
                raise ValueError("location must be a 2-element array")
        
        return
    
    @property
    def morphology_type(self):
        return self._morphology_type    
    
    
    def set_morphology_type(self, target: Target, rng: Generator=None):
        _, transition_nominal, *_ = crater_scaling.get_simple_to_complex_transition_factors(target,rng)
        # Use the 1/2x to 2x the nominal value of the simple->complex transition diameter to get the range of the "transitional" morphology type. This is supported by: Schenk et al. (2004) and Pike (1980) in particular  
        transition_range = (0.5*transition_nominal,2*transition_nominal)
        print(f"Transition range: {transition_range}")
        if self.diameter < transition_range[0]:
            self._morphology_type = "simple" 
        elif self.diameter > transition_range[1]:
            self._morphology_type = "complex"
        else:
            self._morphology_type = "transitional"
        print(f"Morphology type: {self._morphology_type}")
        return
    
     
    def final_to_transient(self, target: Target, rng: Generator=None):
        self.transient_diameter = crater_scaling.final_to_transient(self.diameter,target,rng)
        self.transient_radius = self.transient_diameter / 2
        self.set_morphology_type(target,rng)
        return 


    def transient_to_final(self, target: Target, rng: Generator=None):
        self.diameter = crater_scaling.transient_to_final(self.transient_diameter,target,rng)    
        self.radius = self.diameter / 2
        self.set_morphology_type(target,rng)
        return 
