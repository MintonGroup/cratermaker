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
       
        # Initialize morphology properties
        self.rimheight = None
        self.rimwidth = None
        self.floordepth = None
        self.floordiam = None
        self.peakheight = None

        # Set the morphology based on crater type
        self.set_morphology()        
        
    def set_morphology_parameters(self):
        """
        Sets the morphology parameters based on the crater type.
        """
        fcrat_scaled = self.crater.fcrat * 1e-3  # Convert to km for t hese models

        if self.crater.morphology_type in ["simple", "transition"]:
            # A hybrid model between Pike (1977) and Fassett & Thomson (2014)
            self.rimheight = 0.043 * fcrat_scaled**1.014 * 1e3  # Closer to Fassett & Thomson
            self.rimwidth = 0.257 * fcrat_scaled**1.011 * 1e3   # Pike model
            self.floordepth = 0.224 * fcrat_scaled**1.010 * 1e3 # Closer to Fassett & Thomson
            self.floordiam = 0.200 * fcrat_scaled**1.143 * 1e3  # Fassett & Thomson for D~1km, Pike for D~20km

        elif self.crater.morphology_type in ["complex", "peakring", "multiring"]:
            # Following Pike (1977)
            self.rimheight = 0.236 * fcrat_scaled**0.399 * 1e3  # Pike model
            self.rimwidth = 0.467 * fcrat_scaled**0.836 * 1e3   # Pike model
            self.floordepth = 1.044 * fcrat_scaled**0.301 * 1e3 # Pike model
            # Fassett & Thomson for D~1km, Pike for D~20km, but limited to 90% of fcrat
            self.floordiam = min(0.187 * fcrat_scaled**1.249 * 1e3, 0.9 * self.crater.fcrat)
            self.peakheight = 0.032 * fcrat_scaled**0.900 * 1e3  # Pike model


   
    @staticmethod 
    def craterform(r: float_like, crater: Crater) -> np.float64:
        """
        Calculate the elevation of a crater as a function of distance from the center.

        Parameters:
        ----------
        r : float-like
            Radial distance from the crater center in meters.
        crater: Crater
            The crater to be created.

        Returns:
        ----------
        elevation: np.float64
            Elevation of the crater relative to a reference surface.
        """

        # Empirical crater shape parameters from Fassett et al. (2014)
        r_floor = 0.2
        simple_depth_diam = 0.181
        r_rim = 0.98
        inner_c0 = -0.229
        inner_c1 = 0.228
        inner_c2 = 0.083
        inner_c3 = -0.039

        outer_c0 = 0.188
        outer_c1 = -0.187
        outer_c2 = 0.018
        outer_c3 = 0.015

          
    def form_crater_interior(self):
        """
        Form the interior of the crater.

        This method forms the interior of the crater by removing the material from the surface mesh. It also updates the surface elevation to account for the removed material.
        """
         
        return  