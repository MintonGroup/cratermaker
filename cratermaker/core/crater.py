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
    diameter: np.float64 = None
    radius: np.float64 = None
    transient_diameter: np.float64 = None
    transient_radius: np.float64 = None
    location: np.ndarray = None
    _morphology_type: str = None # This will be computed and is not meant to be set by the user
    
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
    
    @morphology_type.setter
    def morphology_type(self, value):
        self._morphology_type = value
        return 
