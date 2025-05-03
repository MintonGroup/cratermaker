import numpy as np
from numpy.random import Generator
from scipy.optimize import fsolve
from numpy.typing import NDArray, ArrayLike
from typing import Any
from cratermaker.core.crater import Crater
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.surface import Surface
from cratermaker.components.morphology import Morphology
from cratermaker.utils.general_utils import parameter
from cratermaker._cratermaker import crater_functions, ejecta_functions, morphology_functions

@Morphology.register("simplemoon")
class SimpleMoon(Morphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh.

    Parameters
    ----------
    ejecta_truncation : float, optional
        The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a 
        truncation distance based on where the ejecta thickness reaches a small value.         
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    dorays : bool, optional
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket, default is True.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """
    
    def __init__(self, 
                 ejecta_truncation: FloatLike | None = None,
                 dorays: bool = True,
                 rng: Generator | None = None,
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs: Any 
                 ):
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_ejecta_truncation" , None)
        object.__setattr__(self, "_rimheight" , None)
        object.__setattr__(self, "_rimwidth" , None)
        object.__setattr__(self, "_peakheight" , None)
        object.__setattr__(self, "_floor_diameter" , None)
        object.__setattr__(self, "_floordepth" , None)
        object.__setattr__(self, "_ejrim" , None)
        object.__setattr__(self, "_node" , None)

        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays

    def __repr__(self) -> str:
        base = super().__repr__()
        if self.ejecta_truncation is not None:
            base += f"\nEjecta Trunction: {self.ejecta_truncation:.2f} * crater.final_radius"
        else:
            base += f"\nEjecta Truncation: Off"
        return (
            f"{base}\n"
            f"Ejecta Rays: {self.dorays}"
        )

    def _set_morphology_config(self) -> None:
        """
        This method adds a crater to the morphology model and sets the parameters for the morphology based on the crater type. 
        """
        # Set the morphology based on crater type
        diameter_m = self.crater.final_diameter
        diameter_km = diameter_m * 1e-3  # Convert to km for some models

        if self.crater.morphology_type in ["simple", "transitional"]:
            # A hybrid model between Pike (1977) and Fassett & Thomson (2014)
            self.rimheight = 0.043 * diameter_km**1.014 * 1e3  # Closer to Fassett & Thomson
            self.rimwidth = 0.257 * diameter_km**1.011 * 1e3   # Pike model
            self.floordepth = 0.224 * diameter_km**1.010 * 1e3 # Closer to Fassett & Thomson
            self.floor_diameter = 0.200 * diameter_km**1.143 * 1e3  # Fassett & Thomson for D~1km, Pike for D~20km

        elif self.crater.morphology_type in ["complex", "peakring", "multiring"]:
            # Following Pike (1977)
            self.rimheight = 0.236 * diameter_km**0.399 * 1e3  # Pike model
            self.rimwidth = 0.467 * diameter_km**0.836 * 1e3   # Pike model
            self.floordepth = 1.044 * diameter_km**0.301 * 1e3 # Pike model
            # Fassett & Thomson for D~1km, Pike for D~20km, but limited to 90% of diameter
            self.floor_diameter = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m)
            self.peakheight = 0.032 * diameter_km**0.900 * 1e3  # Pike model
        else:
            raise ValueError(f"Unknown morphology type: {self.crater.morphology_type}")
            
        self.ejrim = 0.14 * (diameter_m * 0.5)**(0.74) # McGetchin et al. (1973) Thickness of ejecta at rim
        return


    def form_crater(self, 
                    surface: Surface,
                    crater: Crater | None = None, 
                    **kwargs) -> None:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.
        
        Parameters
        ----------
        surface : Surface
            The surface to be altered.
        crater : Crater
            The crater object to be formed. This is optional if it has already been added
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here).
        """
        self.crater = Crater.maker(crater, **kwargs)

        if not isinstance(surface, Surface):
            raise TypeError("surface must be an instance of Surface")
        self.node_index, self.face_index = surface.find_nearest_index(self.crater.location)

        # Test if the crater is big enough to modify the surface
        rmax = self.rmax(minimum_thickness=surface.smallest_length)
        region_view = surface.extract_region(self.crater.location, rmax)
        if region_view is None: # The crater is too small to change the surface
            return
        crater_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the crater area
        if surface.face_areas[self.face_index] > crater_area:
            return
        
        node_crater_distance, face_crater_distance = surface.get_distance(region_view, self.crater.location)
        reference_face_elevation, reference_node_elevation = surface.get_reference_surface(region_view, face_crater_distance, node_crater_distance, self.crater.location, self.crater.final_radius)
        
        try:
            node_elevation = self.crater_profile(node_crater_distance, 
                                                 reference_node_elevation)
            surface.node_elevation[region_view.node_indices] = node_elevation
            
            face_elevation = self.crater_profile(face_crater_distance, 
                                                 reference_face_elevation)
            surface.face_elevation[region_view.face_indices] = face_elevation
        except:
            print(self)
            raise ValueError("Something went wrong with this crater!")

        self.form_ejecta(surface, crater=self.crater, **kwargs) 
        return  


    def form_ejecta(self,
                    surface: Surface,
                    crater: Crater | None = None,
                    **kwargs) -> None:
        """
        This method forms the ejecta blanket around the crater by altering the elevation variable of the surface mesh.
       
        Parameters
        ----------
        surface : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here). 
        """
        if crater:
            self.crater = crater

        if not isinstance(surface, Surface):
            raise TypeError("surface must be an instance of Surface")
        self.node_index, self.face_index = surface.find_nearest_index(self.crater.location) 

        # Test if the ejecta is big enough to modify the surface
        rmax = self.rmax(minimum_thickness=surface.smallest_length) 
        if not self.ejecta_truncation:
            ejecta_truncation = rmax / self.crater.final_radius
        else:
            ejecta_truncation = self.ejecta_truncation
        region_view = surface.extract_region(self.crater.location, rmax)
        if region_view is None: # The crater is too small to change the surface
            return
        ejecta_area = np.pi * rmax**2
        
        # Check to make sure that the face at the crater location is not smaller than the ejecta blanket area
        if surface.face_areas[self.face_index] > ejecta_area:
            return
        
        node_crater_distance, face_crater_distance = surface.get_distance(region_view, self.crater.location)
        node_crater_bearing, face_crater_bearing  = surface.get_initial_bearing(region_view, self.crater.location)
        morphology_functions.form_ejecta(self, 
                                         region_view, 
                                         node_crater_distance, 
                                         face_crater_distance, 
                                         node_crater_bearing,
                                         face_crater_bearing,
                                         ejecta_truncation,
                                         surface.node_elevation,
                                         surface.face_elevation,
                                         surface.ejecta_thickness,
                                         surface.ray_intensity)
        return


    def crater_profile(self, r: ArrayLike, r_ref: ArrayLike) -> NDArray[np.float64]:
        elevation = crater_functions.profile(r,
                                   r_ref, 
                                   self.crater.final_diameter, 
                                   self.floordepth, 
                                   self.floor_diameter, 
                                   self.rimheight, 
                                   self.ejrim
                                )
        
        return np.array(elevation, dtype=np.float64)
    

    def ejecta_profile(self, r: ArrayLike) -> NDArray[np.float64]:
        elevation = ejecta_functions.profile(r,
                                   self.crater.final_diameter, 
                                   self.ejrim
                                )
        elevation = np.array(elevation, dtype=np.float64)
        return elevation
   
    
    def ejecta_distribution(self, r: ArrayLike, theta: ArrayLike) -> NDArray[np.float64]:
        thickness = ejecta_functions.distribution(r, theta,
                                       self.crater.final_diameter, 
                                       self.ejrim, 
                                       self.ejecta_truncation,
                                       self.dorays
                                    )
        thickness = np.array(thickness, dtype=np.float64)
        return thickness


    def ray_intensity(self, r: ArrayLike, theta: ArrayLike) -> NDArray[np.float64]:
        intensity = ejecta_functions.ray_intensity(r, theta,
                                       self.crater.final_diameter, 
                                       self.ejecta_truncation,
                                    )
        intensity = np.array(intensity, dtype=np.float64)
        return intensity


    def rmax(self, 
            minimum_thickness: FloatLike,
            feature: str = "ejecta",
            crater: Crater | None = None) -> float:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor,
        whichever is smaller.
        
        Parameters
        ----------
        minimum_thickness : FloatLike
            The minimum thickness of the feature blanket in meters.
        feature : str, optional, default = "ejecta"
            The feature to compute the maximum extent. Either "crater" or "ejecta". If "crater" is chosen, the rmax is based
            on where the raised rim is smaller than minimum thickness. 
        crater : Crater, optional
            The crater object to be used. If None, the current crater object is used. If passed, the current crater object is replaced.
        Returns
        -------
        float
            The maximum extent of the crater or ejecta blanket in meters.
        """ 
        
        # Compute the reference surface for the crater 
        if crater is not None:
            self.crater = crater

        if feature == "ejecta":
            def _profile_invert(r):
                return self.ejecta_profile(r) - minimum_thickness
        elif feature == "crater":
            def _profile_invert(r):
                return self.crater_profile(r, np.zeros(1)) - minimum_thickness
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")
    
        # Get the maximum extent 
        rmax = fsolve(_profile_invert, x0=self.crater.final_radius*1.01)[0]
        
        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.crater.final_radius)
            
        if feature == "crater":
            rmax = max(rmax, self.crater.final_radius)

        return float(rmax)

    @property
    def rimheight(self) -> float:
        """
        The height of the crater rim in meters.
        
        Returns
        -------
        float
        """
        return self._rimheight
    
    @rimheight.setter
    def rimheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimheight must be of type FloatLike")
        self._rimheight = float(value)
        
    @property
    def rimwidth(self) -> float:
        """
        The width of the crater rim in meters.
        
        Returns
        -------
        float
        """
        return self._rimwidth
    
    @rimwidth.setter
    def rimwidth(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("rimwidth must be of type FloatLike") 
        self._rimwidth = float(value)
        
    @property
    def peakheight(self) -> float:
        """
        The height of the central peak in meters.
        
        Returns
        -------
        float
        """
        return self._peakheight
    
    @peakheight.setter
    def peakheight(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("peakheight must be of type FloatLike") 
        self._peakheight = float(value)

    @property
    def floor_diameter(self) -> float:
        """
        The diameter of the crater floor in meters.
        
        Returns
        -------
        float
        """
        return self._floor_diameter
    
    @floor_diameter.setter
    def floor_diameter(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("floor_diameter must be of type FloatLike")
        self._floor_diameter = float(value)
        
    @property
    def floordepth(self) -> float:
        """
        Return the depth of the crater floor in m
        
        Returns
        -------
        float
        """
        return self._floordepth
    
    @floordepth.setter
    def floordepth(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("floordepth must be of type FloatLike")
        self._floordepth = float(value)
        
    @property
    def ejrim(self) -> float:
        """
        The thickness of ejecta at the rim in m.
        
        Returns
        -------
        float
        """
        return self._ejrim
    
    @ejrim.setter
    def ejrim(self, value: FloatLike) -> None:
        if not isinstance(value, FloatLike):
            raise TypeError("ejrim must be of type FloatLike")
        self._ejrim = float(value)
        
    @property
    def crater(self):
        return super().crater

    @crater.setter
    def crater(self, value):
        Morphology.crater.fset(self, value)
        self._set_morphology_config()

    @property
    def node_index(self):
        """
        The index of the node closest to the crater location.
        
        Returns
        -------
        int
        """
        return self._node_index
    
    @node_index.setter
    def node_index(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError("node_index must be of type int")
        self._node_index = value

    @property
    def face_index(self):
        """
        The index of the face closest to the crater location.
        
        Returns
        -------
        int
        """
        return self._face_index
    
    @face_index.setter
    def face_index(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError("face_index must be of type int")
        self._face_index = value

    @parameter
    def ejecta_truncation(self) -> float:
        """
        The radius at which the crater is truncated relative to the crater radius.
        
        Returns
        -------
        float or None
        """
        return self._ejecta_truncation 
    
    @ejecta_truncation.setter
    def ejecta_truncation(self, value: FloatLike | None): 
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("truction_radius must be of type FloatLike")
            self._ejecta_truncation = float(value)
        else:
            self._ejecta_truncation = None
            
    @parameter
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
