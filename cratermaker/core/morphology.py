import numpy as np
from numpy.random import Generator
from typing import Tuple
from .target import Target
from ..utils.general_utils import float_like
from ..utils import montecarlo as mc
from .surface import Surface
from .crater import Crater
from .target import Target

class Morphology:
    """
    An operations class for computing the morphology of a crater based on its size and target properties.

    This class encapsulates the logic for altering the topography of the surface based on the crater properties.

    Parameters
    ----------
    crater : Crater
        The crater to be created.
    target : Target
        The target body for the impact simulation.
    surf : Surface
        The surface to be altered.
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    """
    
    def __init__(self, crater: Crater, target: Target, surf: Surface, rng=None) -> None:
        self.crater = crater
        self.target = target
        self.surf = surf
        self.rng = rng if rng is not None else np.random.default_rng()
       
       
    def form_crater_interior(self):
        """
        Form the interior of the crater.

        This method forms the interior of the crater by removing the material from the surface mesh. It also updates the surface elevation to account for the removed material.
        """
         
        return  