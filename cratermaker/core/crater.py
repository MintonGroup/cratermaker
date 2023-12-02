import numpy as np
from numpy.random import Generator
from scipy.optimize import root_scalar
from typing import Tuple
from .target import Target
from ..utils.general_utils import validate_and_convert_location, float_like
from ..utils import montecarlo as mc

class Crater:
    """
    Represents a crater formed by an impact in the simulation.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features.

    Attributes
    ----------
    diameter : float
        The diameter of the crater rim in m.
    radius : float
        The radius of the crater rim in m.
    transient_diameter : float
        The diameter of the transient crater in m.
    transient_radius : float
        The radius of the transient crater in m.    
    location : (2,) float
        The lat. and lon. of the impact point    
    """

    def __init__(self, 
                diameter: float_like = None,
                radius: float_like = None,
                transient_diameter: float_like = None,
                transient_radius: float_like = None,
                location: np.ndarray = None,
                target: Target = None, 
                rng: Generator = None):
       
        if target is None:
            target = Target(name="Moon")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        
        #  Evaluate and check diameter/radius values 
        values_set = sum(x is not None for x in [diameter, radius, transient_diameter, transient_radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius, transient_diameter, transient_radius may be set")
        elif values_set == 0:
            raise ValueError("A crater must include one of diameter, radius, transient_diameter, or transient_radius!")
        
        # Initialize the crater scaling operations class
        self.scale = CraterScaling(target, rng)
        
        # Now call the setters to ensure that all related calculations and checks are performed
        self._diameter = diameter
        self._radius = radius
        self._transient_diameter = transient_diameter
        self._transient_radius = transient_radius    
        self._location = location
        
        if diameter is not None:
            self.diameter = diameter
        elif radius is not None:
            self.radius = radius
        elif transient_diameter is not None:
            self.transient_diameter = transient_diameter
        elif transient_radius is not None:
            self.transient_radius = transient_radius

        # Set location last since it does not depend on other properties
        if location is not None:
            self.location = location              
        self._initialize_location(rng) 
            
        return
    
    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of crater rim must be finite and positive!")
            self._diameter = value
            self._radius = value / 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2
        return
    
    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of crater rim must be finite and positive!")            
            self._radius = value
            self._diameter = value * 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2
            
    @property
    def transient_diameter(self):
        return self._transient_diameter

    @transient_diameter.setter
    def transient_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of transient crater must be finite and positive!")
            self._transient_diameter = value
            self._transient_radius = value / 2
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
        return
    
    @property
    def transient_radius(self):
        return self._transient_radius

    @transient_radius.setter
    def transient_radius(self, value):
        self._transient_radius = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of transient crater must be finite and positive!")            
            self._transient_diameter = value * 2            
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
    
    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, value):
        self._location = value
    
    def _initialize_location(self, rng):
        if self._location is None:
            self.location = mc.get_random_location(rng=rng)
        else:    
            self.location = validate_and_convert_location(self.location)


class CraterScaling:
    """
    A class for handling the scaling relationships between impactors and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties, 
    as well as determining crater morphology based on size and target properties.

    Parameters
    ----------
    target : Target
        The target body for the impact simulation.
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    """

    def __init__(self, target: Target, rng: Generator = None):
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self.target = target
        self.rng = rng if rng else np.random.default_rng()
        
        # Initialize additional attributes for simple->complex transition scale factors. These are set to None here just for clarity
        self.transition_diameter = None
        self.transition_nominal = None
        self.simple_enlargement_factor = None
        self.complex_enlargement_factor = None
        self.final_exp = None

        # Initialize transition factors
        self._compute_simple_to_complex_transition_factors() 
        return
        
        
    def _compute_simple_to_complex_transition_factors(self):
        """
        Computes and sets the internal attributes for transition factors between simple and complex craters.
        """    
        # These terms are used to compute the ratio of the transient crater to simple crater size       
        simple_enlargement_mean = 0.84 # See Melosh (1989) pg. 129 just following eq. 8.2.1
        simple_enlargement_std = 0.04 # Just using Pike (1980) fig. 9 the crater depth varies by about the same amount on either side of the transition so this is a reasonable assumption
        
        # These terms are used in the exponent in the final rim radius/ simple crater radius vs  final radius / transition radius relationship
        # See Holsapple (1993) eq. 28
        final_exp_mean = 0.079    
        final_exp_std = 0.0001 # We add noise because this is nature and nature messy
        complex_enlargement_factor = 1.02
    
        # These terms are used to compute the transition diameter as a function of gravity
        # The transition values come from CTEM and are a synthesis of Pike (1980), Croft (1985), Schenk et al. (2004).
        if self.target.transition_scale_type == "silicate":
            simple_complex_exp = -1.0303 
            simple_complex_mean = 2*16533.8 
            simple_complex_std = 0.04
        elif self.target.transition_scale_type == "ice":
            simple_complex_exp = -1.22486
            simple_complex_mean = 2*3081.39
            simple_complex_std = 0.04
        
        # The nominal value will be used for determining the range of the "transitional" morphology type
        transition_nominal= simple_complex_mean * self.target.gravity**simple_complex_exp
        
        # Draw from a truncated normal distribution for each component of the model
        simple_enlargement_factor = 1.0 / mc.bounded_norm(simple_enlargement_mean, simple_enlargement_std)
        final_exp = mc.bounded_norm(final_exp_mean, final_exp_std)
        simple_complex_fac = simple_complex_mean * np.exp(self.rng.normal(loc=0.0,scale=simple_complex_std))
        transition_diameter = simple_complex_fac * self.target.gravity**simple_complex_exp
        self.transition_diameter = transition_diameter
        self.transition_nominal=transition_nominal
        self.simple_enlargement_factor = simple_enlargement_factor
        self.complex_enlargement_factor = complex_enlargement_factor
        self.final_exp = final_exp
        return 


    def get_morphology_type(self, final_diameter: float_like) -> str:
        """
        Computes and the morphology type of a crater and returns a string corresponding to its type.

        Parameters
        ----------
        final_diameter : float
            The diameter of the crater to compute
        
        Returns
        ----------
        str
            The type of crater "simple", "complex", or "transitional" 
        """
        
        # Use the 1/2x to 2x the nominal value of the simple->complex transition diameter to get the range of the "transitional" morphology type. This is supported by: Schenk et al. (2004) and Pike (1980) in particular  
        transition_range = (0.5*self.transition_nominal,2*self.transition_nominal)
        
        if final_diameter < transition_range[0]:
            morphology_type = "simple" 
        elif final_diameter > transition_range[1]:
            morphology_type = "complex"
        else:
            # We'll uses the distance from the nominal transition diameter to set a probability of being either simple, complex, or transitional.
            if final_diameter < self.transition_nominal:
                p = (self.transition_nominal - diameter)/(self.transition_nominal - transition_range[0])
                categories = ["simple","transitional"]
                prob = [p, 1.0-p] 
                morphology_type = self.rng.choice(categories,p=prob)
            else:
                p = (final_diameter - self.transition_nominal)/(transition_range[1] - self.transition_nominal)
                categories = ["complex","transitional"]
                prob = [p, 1.0-p] 
                morphology_type = self.rng.choice(categories,p=prob)                
        
        return morphology_type


    def f2t_simple(self, Df):
        return Df / self.simple_enlargement_factor
    
    
    def f2t_complex(self, Df):
        return Df / (self.simple_enlargement_factor * self.complex_enlargement_factor) * (Df / self.transition_diameter)**-self.final_exp
    
    
    def final_to_transient(self, final_diameter: float_like, morphology_type: str | None = None) -> np.float64:
        """
        Computes the transient diameter of a crater based on its final diameter and morphology type.

        This method first ensures that the morphology type of the crater is computed. It then calculates
        the transient crater diameter based on the final diameter using scaling factors for simple or complex
        crater morphologies.

        Parameters
        ----------
        final_diameter : float-like
            The final crater diameter in meters for which to compute the transient diameter.
        morphology_type : str, optional
            The morphology type of the crater ("simple", "complex", "transitional")

        Returns
        -------
        np.float64
            Returns the crater transient diameter in meters
        """        
        if not morphology_type:
            morphology_type = self.get_morphology_type(final_diameter) 
        
        if morphology_type == "simple": 
            transient_diameter = self.f2t_simple(final_diameter)
        else:
            transient_diameter = self.f2t_complex(final_diameter)

        transient_diameter = np.float64(transient_diameter)
        return transient_diameter, morphology_type


    def transient_to_final(self, transient_diameter: float_like) -> Tuple[np.float64, str]:
        """
        Computes the final diameter of a crater based on its transient diameter and morphology type.

        This method first ensures that the morphology type of the crater is computed. It then calculates
        the final crater diameter based on the transient diameter using scaling factors for simple or complex
        crater morphologies. This is a bit more complicated than the final->transient calculation  because In 
        the transition region, a particular transient crater diameter could be associate with simple, complex, 
        or transitional crater morphologies. Therefore we need to monte carlo our way into a solution to avoid 
        biasing in favor of one or another in the transient->final computation

        Parameters
        ----------
        transient_diameter : float-like
            The transient diameter in meters of the crater to convert to final

        Returns
        -------
        np.float64
            The final crater diameter
        str
            The morphology type of the crater
        """ 
        
        # Invert the final -> transient functions for  each crater type
        final_diameter_simple = transient_diameter * self.simple_enlargement_factor
        def root_func(final_diameter,Dt,scale):
            return scale.f2t_complex(final_diameter) - Dt
            
        sol = root_scalar(lambda x, *args: root_func(x, *args),bracket=(0.1*final_diameter_simple,10*final_diameter_simple), args=(transient_diameter, self))
        final_diameter_complex = sol.root
        
        # Evaluate the potential morphology that this transient crater could be consistent with. If both potential diameter values are unambigusously is unambiguosuly simple or complex, go with that.
        # If there is disagreement, then we'll draw the answer from a hat and just check to make sure that final_diameter > transient_diameter 
        morphology_options = [self.get_morphology_type(final_diameter_simple),self.get_morphology_type(final_diameter_complex)]
        
        if len(set(morphology_options)) == 1: # We have agreement!
            morphology_type = morphology_options[0]
            if morphology_type == "simple":
                final_diameter = final_diameter_simple
            else:
                final_diameter = final_diameter_complex # this includes transitional types as well
        else: 
            if "simple" in morphology_options: # The disagreement is between simple/complex or simple/transitional
                if morphology_options[0] == "simple":
                    sind = 0
                    cind = 1 
                else:
                    sind = 1
                    cind = 0
                    
                # Randomly draw a morphology based on weighting by whichever option is closest to the transition 
                is_simple = self.rng.random() < np.abs(final_diameter_complex - self.transition_diameter) / np.abs(final_diameter_simple - final_diameter_complex)
                if is_simple:
                    final_diameter = final_diameter_simple
                    morphology_type = morphology_options[sind] 
                else:
                    final_diameter = final_diameter_complex
                    morphology_type = morphology_options[cind]
            else:
                final_diameter = final_diameter_complex
                morphology_type = self.rng.choice(morphology_options)
        
        final_diameter = np.float64(final_diameter)
        morphology_type = morphology_type
        return final_diameter, morphology_type

