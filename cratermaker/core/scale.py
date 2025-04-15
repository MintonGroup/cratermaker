import numpy as np
from numpy.random import Generator
from typing import Tuple, Any, Dict
from scipy.optimize import root_scalar
from .target import Target
from ..utils.custom_types import FloatLike
from ..utils import montecarlo as mc
from ..utils.general_utils import set_properties, create_catalogue, check_properties

class Material:
    """
    Represents the material properties relevant to the crater simulation.

    This class defines various physical properties of the material involved in the cratering process.

    """
    config_ignore = ['catalogue']  # Instance variables to ignore when saving to file
    def __init__(self,
                 name: str | None = None,
                 K1: FloatLike | None = None,
                 mu: FloatLike | None = None,
                 Ybar: FloatLike | None = None,
                 density: FloatLike | None = None,
                 catalogue: Dict[str, Dict[str, FloatLike]] | None = None,
                 **kwargs: Any,
                 ):
        """
        Initialize the target object, setting properties from the provided arguments,
        and creating a catalogue of known solar system targets if not provided.
        
        Parameters
        ----------
        name : str
            The name of the material. If the material is matched to one that is present in the catalogue, the rest of the properties 
            will be retrieved for it unless specified. If the name is not known from the catalogue, then all other properties must 
            be supplied and in order to build a custom material.
        K1 : FloatLike
            Variable used in crater scaling (see _[1])
        mu : FloatLike
            Variable used in crater scaling (see _[1])
        Ybar : FloatLike
            The strength of the material, (Pa)
        density : FloatLike
            Volumentric density of material, (kg/m^3)
        catalogue : Dict[str, Dict[str, FloatLike]]
            An optional dictionary containing a catalogue of known materials to use for the simulation. The catalogue should be 
            constructed using a nested dictionary, where the first level of keys are the names of the materials, and the second level
            are the corresponding property names (K1, mu, Ybar, density). 
            If not provided, a default catalogue will be used.
        **kwargs : Any
            Additional keyword argumments that could be set by the user.
            
        Notes
        -----
        The material properties defined here include crater scaling relationship values that were used in CTEM from Richardson (2009) [1]_.  These values used are from Holsapple (1993) [2]_ and Kraus et al. (2011) [3]_.
        
        References
        ----------
        .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029
        .. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333–373. https://doi.org/10.1146/annurev.ea.21.050193.002001
        .. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724–738. https://doi.org/10.1016/j.icarus.2011.05.016        

        """ 
        # Set the attributes for this class
        self._name = name
        self._K1 = None
        self._mu = None
        self._Ybar = None
        self._density = None
        self._catalogue = None
       
        self.catalogue = catalogue
        
        # Set properties for the Material object based on the arguments passed to the function
        self.set_properties(name=name,
                            K1=K1,
                            mu=mu,
                            Ybar=Ybar,
                            density=density,
                            catalogue=self.catalogue,
                            **kwargs) 
        
        # Check to make sure all required properties are set 
        check_properties(self)
        
        return    

    @property
    def catalogue(self):
        """
        A nested dictionary containing a catalogue of known materials to use for the simulation. 
        
        Returns
        -------
        Dict[str, Dict[str, FloatLike]] or None
            A catalogue of known materials to use for the simulation.
        """
        return self._catalogue
    
    @catalogue.setter
    def catalogue(self, value):
                      
        if not isinstance(value, dict) and value is not None:
            raise TypeError("catalogue must be a dict or None") 
             
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "density" 
        ]
        material_values = [
            ("Water",     2.30,     0.55,   0.0,        1000.0),
            ("Sand",      0.24,     0.41,   0.0,        1750.0),
            ("Dry Soil",  0.24,     0.41,   0.18e6,     1500.0),
            ("Wet Soil",  0.20,     0.55,   1.14e6,     2000.0),
            ("Soft Rock", 0.20,     0.55,   7.60e6,     2250.0),
            ("Hard Rock", 0.20,     0.55,   18.0e6,     2500.0),
            ("Ice",       15.625,   0.48,   0.0,        900.0), 
        ]        
       
        if value is None: 
            self._catalogue = create_catalogue(material_properties, material_values)
        else:
            self._catalogue = value    

    @property
    def name(self):
        """
        The name of the material.
        
        Returns
        -------
        str 
            Name of the material.
        """
        return self._name
    
    @name.setter
    def name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("name must be a string or None")
        self._name = value
        
    @property
    def K1(self):
        """
        K1 crater scaling relationship term. 
        
        Returns
        -------
        np.float64 
        
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._K1
    
    @K1.setter
    def K1(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("K1 must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("K1 must be a positive number")
        self._K1 = np.float64(value)
        
    @property
    def mu(self):
        """
        mu crater scaling relationship term.
        
        Returns
        -------
        np.float64 
        
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._mu
    
    @mu.setter
    def mu(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("mu must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("mu must be a positive number")
        self._mu = np.float64(value)
        
    @property
    def Ybar(self):
        """
        The strength of the material in Pa.
        
        Returns
        -------
        np.float64 
            
                    
        References
        ----------
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697–715. https://doi.org/10.1016/j.icarus.2009.07.029 
        """
        return self._Ybar
    
    @Ybar.setter
    def Ybar(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("Ybar must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("Ybar must be a positive number")
        self._Ybar = np.float64(value)
        
    @property
    def density(self):
        """
            Volumentric density of material in kg/m^3.
        
        Returns
        -------
        np.float64 
        """
        return self._density
    
    @density.setter
    def density(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("density must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("density must be a positive number")
        self._density = np.float64(value)
        
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """         
        set_properties(self,**kwargs)
        return
    
      
class Scale():
    """
    An operations class for computing the scaling relationships between impactors and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties, 
    as well as determining crater morphology based on size and target properties.

    """

    def __init__(self, 
                 target: Target,
                 material: Material | None = None,
                 material_name: str | None = None,                  
                 rng: Generator | None = None, 
                 **kwargs):
        """
        Create an operations class for computing the scaling relationships between impactors and craters.
        
        Parameters
        ----------
        target : Target
            The target body for the impact simulation.
        material : Material or None
            Material composition of the target body. This will override the material supplied by name given by the target object. 
            Only one of material or material_name can be provided.
        material_name : str or None
            Name of the material composition of the target body. This will override the material supplied by name given by the
            target object. Only one of material or material_name can be provided.
        rng : Generator, optional
            A random number generator instance. If not provided, the default numpy RNG will be used. 
        """
        
        self.rng = rng
        self.target = target
        self._material = None
        self._material_name = None
        
        if material is not None and material_name is not None:
            raise ValueError("Only one of material or material_name can be provided")
        
        if self.target.material_name is None:
            if material_name is None or material is None:
                raise ValueError("The target does not provide a material name. Either material or material_name must be provided")
       
        if material is not None:
            self.material = material
        elif material_name is not None: 
            self.material_name = material_name
            self.material = Material(name=self.material_name)        
        else:
            self.material_name = self.target.material_name
            self.material = Material(name=self.material_name)
        
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


    def get_morphology_type(self, final_diameter: FloatLike) -> str:
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
                p = (self.transition_nominal - final_diameter)/(self.transition_nominal - transition_range[0])
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
    
    
    def final_to_transient(self, final_diameter: FloatLike, morphology_type: str | None = None, **kwargs) -> np.float64:
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


    def transient_to_final(self, transient_diameter: FloatLike) -> Tuple[np.float64, str]:
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
 
        
    def projectile_to_crater(self, projectile, **kwargs):
        """
        Convert a projectile to its corresponding crater.

        Parameters
        ----------
        projectile : Projectile
            The projectile to be converted.
        target : Target
            The target body being impacted
        Returns
        -------
        Crater
            The crater resulting from the impact of the projectile.
        """
        from .impact import Crater
        transient_diameter = self.projectile_to_transient(projectile, target=self.target, material=self.material, rng=self.rng)
        crater = Crater(transient_diameter=transient_diameter, target=self.target, rng=self.rng, **kwargs, location=projectile.location, age=projectile.age)

        return crater


    def crater_to_projectile(self, crater, **kwargs):
        """
        Convert a crater back to its corresponding projectile.
        This operation is more hypothetical and approximates the possible projectile that created the crater.

        Parameters
        ----------
        crater : Crater
            The crater to be converted.

        Returns
        -------
        Projectile
            The estimated projectile that could have caused the crater.
        """
        projectile = self.transient_to_projectile(crater, rng=self.rng, **kwargs)
        
        return projectile


    @staticmethod
    def projectile_to_transient(projectile, 
                                target: Target,
                                material: Material,
                                **kwargs: Any) -> np.float64:
        """
        Calculate the transient diameter of a crater based on the properties of the projectile and target.

        Parameters
        ----------
        projectile : Projectile
            The projectile responsible for the impact.
        target : Target
            The target body for the impact simulation.
        material : Material
            The material composition of the target body.
        **kwargs : Any
            Additional keyword arguments that might influence the calculation.

        Returns
        -------
        np.float64
            The calculated transient diameter of the crater resulting from the impact.
        """
        from .impact import Projectile
        if not isinstance(projectile, Projectile):
            raise TypeError("projectile must be an instance of Projectile")
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        if not isinstance(material, Material):
            raise TypeError("material must be an instance of Material")
            
        # Compute some auxiliary quantites
        projectile.mass = 4.0/3.0 * np.pi * projectile.density * (projectile.radius)**3
        mu = material.mu
        kv = material.K1
        c1 = 1.0 + 0.5 * mu
        c2 = (-3 * mu)/(2.0 + mu)

        # Find dimensionless quantities
        pitwo = (target.gravity * projectile.radius)/(projectile.vertical_velocity**2)
        pithree = material.Ybar / (material.density * (projectile.vertical_velocity**2))
        pifour = material.density / projectile.density
        pivol = kv * ((pitwo * (pifour**(-1.0/3.0))) + (pithree**c1))**c2
        pivolg = kv * (pitwo * (pifour**(-1.0/3.0)))**c2
        
        # find transient crater volume and radii (depth = 1/3 diameter)
        cvol = pivol * (projectile.mass / material.density)
        cvolg = pivolg * (projectile.mass / material.density)
        transient_radius = (3 * cvol / np.pi)**(1.0/3.0)
        transient_radius_gravscale = (3 * cvolg / np.pi)**(1.0/3.0)
        
        transient_diameter = transient_radius * 2
        
        if transient_diameter < projectile.diameter:
            transient_diameter = projectile.diameter
        
        return transient_diameter


    def transient_to_projectile(self, 
                                crater, 
                                rng: Generator = None, 
                                **kwargs: Any): 
        """
        Estimate the characteristics of the projectile that could have created a given crater.

        This method approximates the properties of a hypothetical projectile based on the characteristics
        of a known crater.

        Parameters
        ----------
        crater : Crater
            The crater for which to estimate the projectile.
        rng : Generator, optional
            Random number generator instance used for any probabilistic calculations.
        **kwargs : Any
            Additional keyword arguments that might influence the calculation.

        Returns
        -------
        Projectile
            The computed projectile for the crater.
        """
        from .impact import Crater, Projectile
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
                
        # We'll create a Projectile object that will allow us to set velocity
        # First pop off any pre-computed values so that we use the transient_diameter alone
        kwargs.pop("diameter",None)
        kwargs.pop("location",None)
        kwargs.pop("age",None)
        projectile = Projectile(diameter=crater.transient_diameter, target=self.target, location=crater.location, rng=rng, age=crater.age, **kwargs)
        
        def root_func(projectile_diameter: FloatLike, 
                      projectile: Projectile, 
                      target: Target,
                      material: Material
                      ) -> np.float64:
            
            projectile.diameter = projectile_diameter
            transient_diameter = self.projectile_to_transient(projectile, target, material)
            return transient_diameter - crater.transient_diameter 
        
        sol = root_scalar(lambda x, *args: root_func(x, *args),bracket=(1e-5*crater.transient_diameter,1.2*crater.transient_diameter), args=(projectile, self.target, self.material))
        
        # Regenerate the projectile with the new diameter value
        projectile = Projectile(diameter=sol.root, target=self.target, location=projectile.location, velocity=projectile.velocity, angle=projectile.angle, direction=projectile.direction, rng=rng, age=projectile.age)
        
        return projectile

    @property
    def transition_diameter(self) -> np.float64:
        """
        The transition diameter between simple and complex craters in m.
        
        Returns
        -------
        np.float64
        """
        return self._transition_diameter
    
    @transition_diameter.setter
    def transition_diameter(self, value: FloatLike) -> None:
        self._transition_diameter = np.float64(value)
        
    @property
    def transition_nominal(self) -> np.float64:
        """
        The nominal transition diameter for crater morphology in m.
        
        Returns
        -------
        np.float64
        """
        return self._transition_nominal
    
    @transition_nominal.setter
    def transition_nominal(self, value: FloatLike) -> None:
        self._transition_nominal = np.float64(value)

    @property
    def simple_enlargement_factor(self) -> np.float64:
        """
        The enlargement factor for simple craters.
        
        Returns
        -------
        np.float64
        """
        return self._simple_enlargement_factor
    
    @simple_enlargement_factor.setter
    def simple_enlargement_factor(self, value: FloatLike) -> None:
        self._simple_enlargement_factor = np.float64(value)

    @property
    def complex_enlargement_factor(self) -> np.float64:
        """
        The enlargement factor for complex craters.
        
        Returns
        -------
        np.float64
        """
        return self._complex_enlargement_factor
    
    @complex_enlargement_factor.setter
    def complex_enlargement_factor(self, value: FloatLike) -> None:
        self._complex_enlargement_factor = np.float64(value)

    @property
    def final_exp(self) -> np.float64:
        """
        The exponent used in the final rim radius to simple crater radius relationship.
        
        Returns
        -------
        np.float64
        """
        return self._final_exp
    
    @final_exp.setter
    def final_exp(self, value: FloatLike) -> None:
        self._final_exp = np.float64(value)

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
        if value is None:
            self._target = Target(name="Moon")
            return
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value
        return 
    

    @property
    def material_name(self):
        """
        The name of the material composition of the target body.
        
        Returns
        -------
        str 
        """
        return self._material_name

    @material_name.setter
    def material_name(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError("material_name must be a string or None")
        self._material_name = value
        
        
    @property
    def material(self):
        """
        The material properties associated with the target body.
        
        Returns
        -------
        Material
        """
        return self._material
    
    @material.setter
    def material(self, value):
        if not isinstance(value, Material) and value is not None:
            raise TypeError("material must be a Material object or None")
        self._material = value    
    
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
      
