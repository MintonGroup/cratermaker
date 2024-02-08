import numpy as np
from numpy.random import Generator
import xarray as xr
import uxarray as uxr
from scipy.optimize import fsolve
from typing import Tuple, Any
from numpy.typing import ArrayLike
from .target import Target
from ..utils.custom_types import FloatLike
from ..utils import montecarlo as mc
from .surface import Surface
from .target import Target
from ..cython import morphology

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
    ejecta_truncation : float, optional
        The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a 
        truncation distance based on where the ejecta thickness reaches a small value.         
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    """
    
    def __init__(self, 
                 crater,  
                 target: Target | None = None, 
                 rng: Generator | None = None,
                 ejecta_truncation: FloatLike | None = None
                 ):

        self.crater = crater 
        self.target = target
        self.rng = rng
        self.ejecta_truncation = ejecta_truncation

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

    def profile(self, r: ArrayLike, r_ref: ArrayLike) -> np.float64:
        elevation = morphology.profile(r,
                                       r_ref, 
                                       self.diameter, 
                                       self.floordepth, 
                                       self.floordiam, 
                                       self.rimheight, 
                                       self.ejrim, 
                                       RIMDROP)
         
        return np.array(elevation, dtype=np.float64)
    
    def form_crater(self, 
                    surf: Surface,
                    **kwargs) -> None:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.
        
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here).
        """
      
        rmax = self.compute_rmax(minimum_thickness=surf.smallest_length)
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return

        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        #region_surf['node_crater_bearing'], region_surf['face_crater_bearing']  = surf.get_initial_bearing(self.crater.location)
        region_surf.get_reference_surface(self.crater.location, self.crater.radius)
        
        try:
            node_elevation = self.profile(region_surf['node_crater_distance'].values, region_surf['reference_node_elevation'].values)
            surf['node_elevation'].loc[{'n_node': region_surf.uxgrid._ds["subgrid_node_indices"]}] = node_elevation
            
            face_elevation = self.profile(region_surf['face_crater_distance'].values, region_surf['reference_face_elevation'].values)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] = face_elevation
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return  
   
    def compute_rmax(self,
                     minimum_thickness: np.float64,
                     ) -> np.float64:
        """
        Compute the maximum extent of the crater based on the minimum thickness of the ejecta blanket or the ejecta_truncation factor,
        whichever is smaller.
        
        Parameters
        ----------
        minimum_thickness : np.float64
            The minimum thickness of the ejecta blanket in meters.
        """ 
        
        # Compute the reference surface for the crater 
        def _profile_invert(r):
            return self.profile(r, np.zeros(1)) - minimum_thickness
       
        # Get the maximum extent 
        rmax = fsolve(_profile_invert, x0=self.radius)[0]
        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.radius)        
        # Get the maximum extent 
        rmax = fsolve(_profile_invert, x0=self.radius)[0]
        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.radius)
            
        return rmax

    @property
    def diameter(self) -> np.float64:
        """
        The diameter of the crater in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.diameter

    @property
    def radius(self) -> np.float64:
        """
        The radius of the crater in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.radius
   
    @property 
    def morphology_type(self) -> str:
        """
        The morphology type of the crater.
        
        Returns
        -------
        str 
        """
        return self.crater.morphology_type 
    
    @property
    def rimheight(self) -> np.float64:
        """
        The height of the crater rim in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.rimheight
    
    @rimheight.setter
    def rimheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimheight must be of type FloatLike")
        self.crater.rimheight = np.float64(value)
        
    @property
    def rimwidth(self) -> np.float64:
        """
        The width of the crater rim in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.rimwidth
    
    @rimwidth.setter
    def rimwidth(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimwidth must be of type FloatLike") 
        self.crater.rimwidth = np.float64(value)
        
    @property
    def peakheight(self) -> np.float64:
        """
        The height of the central peak in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.peakheight
    
    @peakheight.setter
    def peakheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("peakheight must be of type FloatLike") 
        self.crater.peakheight = np.float64(value)

    @property
    def floordiam(self) -> np.float64:
        """
        The diameter of the crater floor in meters.
        
        Returns
        -------
        np.float64
        """
        return self.crater.floordiam
    
    @floordiam.setter
    def floordiam(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("floordiam must be of type FloatLike")
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
        if not isinstance(value, FloatLike):
            raise TypeError("floordepth must be of type FloatLike")
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
        if not isinstance(value, FloatLike):
            raise TypeError("ejrim must be of type FloatLike")
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
       
    @property
    def ejecta_truncation(self) -> np.float64:
        """
        The radius at which the crater is truncated relative to the crater radius.
        
        Returns
        -------
        np.float64 or None
        """
        return self.crater.ejecta_truncation 
    
    @ejecta_truncation.setter
    def ejecta_truncation(self, value: FloatLike | None): 
        if value:
            if not isinstance(value, FloatLike):
                raise TypeError("truction_radius must be of type FloatLike")
            self.crater.ejecta_truncation = np.float64(value)
        else:
            self.crater.ejecta_truncation = None