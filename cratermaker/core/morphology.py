import numpy as np
from numpy.random import Generator
from typing import Tuple, Any
from .target import Target
from ..utils.custom_types import FloatLike
from ..utils import montecarlo as mc
from .surface import Surface
from .target import Target

RIMDROP = 4.20

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
    
    def __init__(self, 
                 crater,  
                 target: Target | None=None, 
                 rng: Generator | None=None
                 ):

        self.crater = crater 
        self.target = target
        self.rng = rng
        self.rimheight = None
        self.rimwidth = None
        self.floordepth = None
        self.floordiam = None

        # Set the morphology based on crater type
        self.set_morphology_parameters()        
       
    def __repr__(self):
        return (f"Morphology(morphology_type={self.morphology_type}, diameter={self.diameter}, "
                f"rimheight: {self.rimheight}, rimwidth: {self.rimwidth}, floordepth: {self.floordepth}, floordiam: {self.floordiam})") 
               
    def set_morphology_parameters(self):
        """
        Sets the morphology parameters based on the crater type.
        """
        diameter_km = self.diameter * 1e-3  # Convert to km for these models

        if self.morphology_type in ["simple", "transitional"]:
            # A hybrid model between Pike (1977) and Fassett & Thomson (2014)
            self.rimheight = 0.043 * diameter_km**1.014 * 1e3  # Closer to Fassett & Thomson
            self.rimwidth = 0.257 * diameter_km**1.011 * 1e3   # Pike model
            self.floordepth = 0.224 * diameter_km**1.010 * 1e3 # Closer to Fassett & Thomson
            self.floordiam = 0.200 * diameter_km**1.143 * 1e3  # Fassett & Thomson for D~1km, Pike for D~20km

        elif self.morphology_type in ["complex", "peakring", "multiring"]:
            # Following Pike (1977)
            self.rimheight = 0.236 * diameter_km**0.399 * 1e3  # Pike model
            self.rimwidth = 0.467 * diameter_km**0.836 * 1e3   # Pike model
            self.floordepth = 1.044 * diameter_km**0.301 * 1e3 # Pike model
            # Fassett & Thomson for D~1km, Pike for D~20km, but limited to 90% of diameter
            self.floordiam = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * self.diameter)
            self.peakheight = 0.032 * diameter_km**0.900 * 1e3  # Pike model
        else:
            raise ValueError(f"Unknown morphology type: {self.morphology_type}")
            
        self.ejrim = 0.14 * (self.diameter * 0.5)**(0.74) # McGetchin et al. (1973) Thickness of ejecta at rim

   
    def crater_profile(self, r: FloatLike) -> np.float64:
        """
        Calculate the elevation of a crater as a function of distance from the center.

        Parameters
        -----------
        r : float-like
            Radial distance from the crater center in meters.
        crater: Crater
            The crater to be created.

        Returns
        --------
        np.float64
            Elevation of the crater relative to a reference surface.
        """

        # Constants
        A = 4.0 / 11.0
        B = -32.0 / 187.0

        # Calculate the floor radius relative to the final crater radius
        flrad = self.floordiam / self.diameter

        # Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
        c1 = (-self.floordepth - self.rimheight) / (flrad - 1.0 + A * (flrad**2 - 1.0) + B * (flrad**3 - 1.0))
        c0 = self.rimheight - c1 * (1.0 + A + B)
        c2 = A * c1
        c3 = B * c1
        
        r = np.abs(r) / self.radius

        # Compute the height based on the relative radial distance
        if r < flrad:
            h = -self.floordepth
        elif r >= 1.0:
            h = (self.rimheight - self.ejrim) * (r**(-RIMDROP))
        else:
            h = c0 + c1 * r + c2 * r**2 + c3 * r**3

        return h        

    
    def ejecta_profile(self, r: FloatLike) -> np.float64:
        """
        Calculate the thickness of ejecat as a function of distance from the center of the crater.

        Parameters
        ----------
        r : float-like
            Radial distance from the crater center in meters.
        crater: Crater
            The crater to be created.

        Returns
        -------
        np.float64
            Elevation of the crater relative to a reference surface.
        """
        
        ejprofile = 3.0
        r = np.abs(r) / self.radius
        thick = self.ejrim * (r)**(-ejprofile)
        return thick
          
    def form_crater(self, surf: Surface) -> None:
        """
        Form the interior of the crater.

        This method forms the interior of the crater by altering the elevation variable of the surface mesh.
        """
       
        def _crater_profile(r):
            h = self.crater_profile(r) 
            if r > self.crater.radius:
                h += self.ejecta_profile(r)
            if h > self.floordepth:
                h = self.floordepth
            return h
    
        try: 
            surf['elevation'] += np.vectorize(_crater_profile)(surf['crater_distance']) 
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
         
        return  
    
    @property
    def diameter(self) -> float:
        """
        Return the diameter of the crater in meters.
        """
        return self.crater.diameter

    @property
    def radius(self) -> float:
        """
        Return the radius of the crater in meters.
        """
        return self.crater.radius
   
    @property 
    def morphology_type(self) -> str:
        """
        Return the morphology type of the crater.
        """
        return self.crater.morphology_type 
    
    @property
    def rimheight(self) -> float:
        """
        Return the height of the crater rim in meters.
        """
        return self.crater.rimheight
    
    @rimheight.setter
    def rimheight(self, value: FloatLike) -> None:
        """
        Set the height of the crater rim in meters.
        """
        self.crater.rimheight = value
        
    @property
    def rimwidth(self) -> float:
        """
        Return the width of the crater rim in meters.
        """
        return self.crater.rimwidth
    
    @rimwidth.setter
    def rimwidth(self, value: FloatLike) -> None:
        """
        Set the width of the crater rim in meters.
        """
        self.crater.rimwidth = value
        
    @property
    def peakheight(self) -> float:
        """
        Return the height of the central peak in meters.
        """
        return self.crater.peakheight
    
    @peakheight.setter
    def peakheight(self, value: FloatLike) -> None:
        """
        Set the height of the central peak in meters.
        """
        self.crater.peakheight = value

    @property
    def floordiam(self) -> float:
        """
        Return the diameter of the crater floor in meters.
        """
        return self.crater.floordiam
    
    @floordiam.setter
    def floordiam(self, value: FloatLike) -> None:
        self.crater.floordiam = np.float64(value)
        
    @property
    def floordepth(self) -> np.float64:
        """
        Return the depth of the crater floor in m
        
        Returns
        -------
        np.float64
        """
        return self.crater.floordepth
    
    @floordepth.setter
    def floordepth(self, value: FloatLike) -> None:
        self.crater.floordepth = np.float64(value)
        
    @property
    def ejrim(self) -> np.float64:
        """
        The thickness of ejecta at the rim in m.
        
        Returns
        -------
        np.float64
        """
        return self.crater.ejrim
    
    @ejrim.setter
    def ejrim(self, value: FloatLike) -> None:
        self.crater.ejrim = np.float64(value)
        
        
    @property
    def crater(self):
        """
        The crater to be created.
        
        Returns
        -------
        Crater
        """ 
        return self._crater
    
    @crater.setter
    def crater(self, value):
        from .impact import Crater
        if value is not None and not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value
        return 
        
    @property
    def target(self):
        """
        The target body for the impact.
        
        Returns
        -------
        Target
        """ 
        return self._target
    
    @target.setter
    def target(self, value):
        from .target import Target
        if value is None:
            self._target = Target(name="Moon")
            return
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value
        return        
        
        
    @property
    def rng(self):
        """
        A random number generator instance.
        
        Returns
        -------
        Generator
        """ 
        return self._rng
    
    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()   
        