import numpy as np
from numpy.random import Generator
from scipy.optimize import fsolve
from numpy.typing import ArrayLike
from typing import Any
from .target import Target
from ..utils.custom_types import FloatLike
from .surface import Surface
from .target import Target
from .. import crater, ejecta

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
    ejecta_truncation : float, optional
        The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a 
        truncation distance based on where the ejecta thickness reaches a small value.         
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    dorays : bool, optional
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket, default is True.
    **kwargs : Any
        Additional keyword arguments to be passed to internal functions.
    """
    
    def __init__(self, 
                 crater,  
                 target: Target | None = None, 
                 ejecta_truncation: FloatLike | None = None,
                 rng: Generator | None = None,
                 dorays: bool = True,
                 **kwargs: Any 
                 ):

        self.crater = crater 
        self.target = target
        self.rng = rng
        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays

        # Set the morphology based on crater type
        self._set_morphology_parameters()
       
    def __repr__(self):
        return (f"Morphology(morphology_type={self.morphology_type}, diameter={self.diameter}, "
                f"rimheight: {self.rimheight}, rimwidth: {self.rimwidth}, floordepth: {self.floordepth}, floordiam: {self.floordiam})") 
    
    def _set_morphology_parameters(self):
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


    def crater_profile(self, r: ArrayLike, r_ref: ArrayLike) -> np.float64:
        elevation = crater.profile(r,
                                   r_ref, 
                                   self.diameter, 
                                   self.floordepth, 
                                   self.floordiam, 
                                   self.rimheight, 
                                   self.ejrim
                                )
        
        return np.array(elevation, dtype=np.float64)
    

    def ejecta_profile(self, r: ArrayLike) -> np.float64:
        elevation = ejecta.profile(r,
                                   self.diameter, 
                                   self.ejrim
                                )
        elevation = np.array(elevation, dtype=np.float64)
        return elevation
   
    
    def ejecta_distribution(self, r: ArrayLike, theta: ArrayLike) -> np.float64:
        thickness = ejecta.distribution(r, theta,
                                       self.diameter, 
                                       self.ejrim, 
                                       self.ejecta_truncation,
                                       self.dorays
                                    )
        thickness = np.array(thickness, dtype=np.float64)
        return thickness
  
    def ray_intensity(self, r: ArrayLike, theta: ArrayLike) -> np.float64:
        intensity = ejecta.ray_intensity(r, theta,
                                       self.diameter, 
                                       self.ejrim, 
                                       self.ejecta_truncation,
                                    )
        intensity = np.array(intensity, dtype=np.float64)
        return intensity
           
    def compute_rmax(self, 
                     minimum_thickness: np.float64,
                     feature: str = "ejecta") -> np.float64:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor,
        whichever is smaller.
        
        Parameters
        ----------
        minimum_thickness : np.float64
            The minimum thickness of the feature blanket in meters.
        feature : str, optional, default = "ejecta"
            The feature to compute the maximum extent. Either "crater" or "ejecta". If "crater" is chosen, the rmax is based
            on where the raised rim is smaller than minimum thickness. 
        """ 
        
        # Compute the reference surface for the crater 

        if feature == "ejecta":
            def _profile_invert(r):
                return self.ejecta_profile(r) - minimum_thickness
        elif feature == "crater":
            def _profile_invert(r):
                return self.crater_profile(r, np.zeros(1)) - minimum_thickness
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")
    
        # Get the maximum extent 
        rmax = fsolve(_profile_invert, x0=self.radius*1.01)[0]
        
        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.radius)
            
        if feature == "crater":
            rmax = max(rmax, self.radius)

        return rmax     
     
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
        
        # Test if the crater is big enough to modify the surface
        rmax = self.compute_rmax(minimum_thickness=surf.smallest_length)
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        crater_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the crater area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > crater_area:
            return
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf.get_reference_surface(self.crater.location, self.crater.radius)
        
        try:
            node_elevation = self.crater_profile(region_surf['node_crater_distance'].values, 
                                                 region_surf['reference_node_elevation'].values)
            surf['node_elevation'].loc[{'n_node': region_surf.uxgrid._ds["subgrid_node_indices"]}] = node_elevation
            
            face_elevation = self.crater_profile(region_surf['face_crater_distance'].values, 
                                                 region_surf['reference_face_elevation'].values)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] = face_elevation
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return  
    
     
    def form_ejecta(self,
                    surf: Surface,
                    **kwargs) -> None:
        """
        This method forms the ejecta blanket around the crater by altering the elevation variable of the surface mesh.
       
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here). 
        """
                 
        # Test if the ejecta is big enough to modify the surface
        rmax = self.compute_rmax(minimum_thickness=surf.smallest_length) 
        if not self.ejecta_truncation:
            self.ejecta_truncation = rmax / self.radius
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        ejecta_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the ejecta blanket area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > ejecta_area:
            return                  
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf['node_crater_bearing'], region_surf['face_crater_bearing']  = region_surf.get_initial_bearing(self.crater.location)
        
        try:
            node_thickness = self.ejecta_distribution(region_surf['node_crater_distance'].values, 
                                                     region_surf['node_crater_bearing'].values)
            surf['node_elevation'].loc[{'n_node': region_surf.uxgrid._ds["subgrid_node_indices"]}] += node_thickness
            
            face_thickness = self.ejecta_distribution(region_surf['face_crater_distance'].values, 
                                                     region_surf['face_crater_bearing'].values)
            surf['face_elevation'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_thickness
            surf['ejecta_thickness'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_thickness
            
            if self.dorays: 
                face_intensity = self.ray_intensity(region_surf['face_crater_distance'].values, 
                                                     region_surf['face_crater_bearing'].values)
                surf['ray_intensity'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_intensity
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return  
    
    def form_secondaries(self,
                         surf: Surface,
                         **kwargs) -> None:
        """
        This method forms secondary craters around the primary crater. Currently it only generates the ray intensity function.
        Maximum ray length formula is from [1]_.
       
        Parameters
        ----------
        surf : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here). 
            
        References
        ----------
        .. [1] Elliott, J.R., Huang, Y.-H., Minton, D.A., Freed, A.M., 2018. The length of lunar crater rays explained using secondary crater scaling. Icarus 312, 231–246. https://doi.org/10.1016/j.icarus.2018.04.015

        """
        
        # Elliott et al. (2018) eq. 2
        A = 6.59 + self.rng.normal(loc=0.0,scale=1.45) 
        p = 1.27 + self.rng.normal(loc=0.0,scale=0.06)
        rmax = A * (self.radius / 1000) **p * 1000
        if not self.ejecta_truncation:
            self.ejecta_truncation = rmax / self.radius
        region_surf = surf.extract_region(self.crater.location, rmax)
        if not region_surf: # The crater is too small to change the surface
            return
        ray_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the ray effectt area
        if surf['face_areas'].isel(n_face=self.crater.face_index) > ray_area:
            return                  
        
        region_surf['node_crater_distance'], region_surf['face_crater_distance'] = region_surf.get_distance(self.crater.location)
        region_surf['node_crater_bearing'], region_surf['face_crater_bearing']  = region_surf.get_initial_bearing(self.crater.location)
        
        try:
            face_intensity = self.ray_intensity(region_surf['face_crater_distance'].values, 
                                                region_surf['face_crater_bearing'].values)
            surf['ray_intensity'].loc[{'n_face': region_surf.uxgrid._ds["subgrid_face_indices"]}] += face_intensity
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")
                 
        return      
            
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
            
    @property
    def dorays(self) -> bool:
        """
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket.
        
        Returns
        -------
        bool
        """
        return self._dorays
    
    @dorays.setter
    def dorays(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("dorays must be of type bool")
        self._dorays = value
        return