import numpy as np
from numpy.random import Generator
from typing import Any
from scipy.optimize import root_scalar
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils import montecarlo as mc
from cratermaker.utils.general_utils import _set_properties
from cratermaker.components.scaling import register_scaling_model, ScalingModel
from cratermaker.utils.general_utils import _create_catalogue
from cratermaker.core.target import Target

@register_scaling_model("richardson2009")
class Richardson2009(ScalingModel):
    """
    This is an operations class for computing the scaling relationships between impactors and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties, 
    as well as determining crater morphology based on size and target propertiesImplements the scaling laws described in Richardson (2009) [1]_ that were implemented in CTEM.
        
    Parameters
    ----------
    target : Target or str, required
        The target body for the impact simulation, or the name of a target. If the name is provided, you can also provide additional keyword arguments that will be passed to the Target class.
    material_name : str or None
        Name of the target material composition of the target body to look up from the built-in catalogue. Options include "water", "sand", "dry soil", "wet soil", "soft rock", "hard rock", and "ice".
    K1 : FloatLike, optional
        Variable used in crater scaling (see _[1])
    mu : FloatLike, optional
        Variable used in crater scaling (see _[1])
    Ybar : FloatLike, optional
        The strength of the target material, (Pa)
    target_density : FloatLike, optional
        Volumentric density of target material, (kg/m^3)
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used. 

    Notes
    -----
    - The `target` parameter is required and must be an instance of the `Target` class.
    - The `material_name` parameter is optional. If not provided, it will be retrieved from `target`. Setting it explicitly will override the value in `target`.
    - The `K1`, `mu`, `Ybar`, and `target_density` parameters are optional. If not provided, they will be retrieved from the material catalogue based on the `material_name`. Setting them explicitly will override the values in the catalogue.
    - The built-in material property values are from Holsapple (1993) [2]_ and Kraus et al. (2011) [3]_. 
    
    References
    ----------
    .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
    .. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
    .. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016        
    """  

    def __init__(self, 
                 target: str | Target = None,
                 material_name: str | None = None,
                 K1: FloatLike | None = None,
                 mu: FloatLike | None = None,
                 Ybar: FloatLike | None = None,
                 target_density: FloatLike | None = None,
                 rng: Generator | None = None, 
                 **kwargs):
        """

        References
        ----------
        .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
        """
        super().__init__()
        if isinstance(target, str):
            try:
                target = Target(target,material_name=material_name,**kwargs)
            except:
                raise ValueError(f"Invalid target name {target}")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target or a valid name of a target body")

        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_material_name", None)
        object.__setattr__(self, "_K1", None)
        object.__setattr__(self, "_mu", None)
        object.__setattr__(self, "_Ybar", None)
        object.__setattr__(self, "_target_density", None)
        object.__setattr__(self, "_rng", None)

        if material_name is None and target.material_name is not None: 
            material_name = target.material_name

        material_catalogue = self._create_material_catalogue() 

        _set_properties(self,
                        target=target,
                        name=material_name,
                        K1=K1,
                        mu=mu,
                        Ybar=Ybar,
                        target_density=target_density,
                        catalogue=material_catalogue,
                        rng=rng,
                        **kwargs
                    ) 
        
        # Initialize additional attributes for simple->complex transition scale factors. These are set to None here just for clarity
        self.transition_diameter = None
        self.transition_nominal = None
        self.simple_enlargement_factor = None
        self.complex_enlargement_factor = None
        self.final_exp = None

        # Initialize transition factors
        self._compute_simple_to_complex_transition_factors() 
        return
    
    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        include_list=("material_name", "K1", "mu", "Ybar", "target_density")
        # Add it to the set of user-defined parameters if it is in the list of parameters
        public_name = name.lstrip("_")
        if public_name in include_list:
            self._user_defined.add(public_name)


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
        from cratermaker.core.impact import Crater
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
            transient_diameter = self._f2t_simple(final_diameter)
        else:
            transient_diameter = self._f2t_complex(final_diameter)

        transient_diameter = np.float64(transient_diameter)
        return transient_diameter, morphology_type


    def transient_to_final(self, transient_diameter: FloatLike) -> tuple[np.float64, str]:
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
            return scale._f2t_complex(final_diameter) - Dt
            
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

    @staticmethod
    def projectile_to_transient(projectile, 
                                target: Target,
                                K1: FloatLike,
                                mu: FloatLike,
                                Ybar: FloatLike,
                                target_density: FloatLike,
                                **kwargs: Any) -> np.float64:
        """
        Calculate the transient diameter of a crater based on the properties of the projectile and target.

        Parameters
        ----------
        projectile : Projectile
            The projectile responsible for the impact.
        target : Target
            The target body for the impact simulation.
        K1 : FloatLike, optional
            Variable used in crater scaling (see _[1])
        mu : FloatLike, optional
            Variable used in crater scaling (see _[1])
        Ybar : FloatLike, optional
            The strength of the target material, (Pa)
        target_density : FloatLike, optional
            Volumentric density of target material, (kg/m^3)

        **kwargs : Any
            Additional keyword arguments that might influence the calculation.

        Returns
        -------
        np.float64
            The calculated transient diameter of the crater resulting from the impact.

        References
        ----------
            .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
        """
        from cratermaker.core.impact import Projectile
        if not isinstance(projectile, Projectile):
            raise TypeError("projectile must be an instance of Projectile")
        if not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
            
        # Compute some auxiliary quantites
        projectile.mass = 4.0/3.0 * np.pi * projectile.density * (projectile.radius)**3
        c1 = 1.0 + 0.5 * mu
        c2 = (-3 * mu)/(2.0 + mu)

        # Find dimensionless quantities
        pitwo = (target.gravity * projectile.radius)/(projectile.vertical_velocity**2)
        pithree = Ybar / (target_density * (projectile.vertical_velocity**2))
        pifour = target_density / projectile.density
        pivol = K1 * ((pitwo * (pifour**(-1.0/3.0))) + (pithree**c1))**c2
        pivolg = K1 * (pitwo * (pifour**(-1.0/3.0)))**c2
        
        # find transient crater volume and radii (depth = 1/3 diameter)
        cvol = pivol * (projectile.mass / target_density)
        cvolg = pivolg * (projectile.mass / target_density)
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
        from cratermaker.core.impact import Crater, Projectile
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
        
        def root_func(
                        projectile_diameter: FloatLike, 
                        projectile: Projectile, 
                        target: Target,
                        K1: FloatLike,
                        mu: FloatLike,
                        Ybar: FloatLike,
                        target_density: FloatLike,
                      ) -> np.float64:
            
            projectile.diameter = projectile_diameter
            transient_diameter = self.projectile_to_transient(projectile, target, K1, mu, Ybar, target_density)
            return transient_diameter - crater.transient_diameter 
        
        sol = root_scalar(lambda x, *args: root_func(x, *args),bracket=(1e-5*crater.transient_diameter,1.2*crater.transient_diameter), args=(projectile, self.target, self.K1, self.mu, self.Ybar, self.target_density))
        
        # Regenerate the projectile with the new diameter value
        projectile = Projectile(diameter=sol.root, target=self.target, location=projectile.location, velocity=projectile.velocity, angle=projectile.angle, direction=projectile.direction, rng=rng, age=projectile.age)
        
        return projectile

    def _create_material_catalogue(self):
        
        # Define some built-in catalogue values for known solar system materials of interest
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "target_density" 
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
        return _create_catalogue(material_properties, material_values)
    

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

    def _f2t_simple(self, Df):
        return Df / self.simple_enlargement_factor
    
    def _f2t_complex(self, Df):
        return Df / (self.simple_enlargement_factor * self.complex_enlargement_factor) * (Df / self.transition_diameter)**-self.final_exp

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
    def name(self):
        """
        The name of the material composition of the target body.
        
        Returns
        -------
        str 
        """
        return self._material_name

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
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029 
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
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029 
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
        Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029 
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
    def target_density(self):
        """
        Volumentric density of material in kg/m^3.
        
        Returns
        -------
        np.float64 
        """
        return self._target_density
    
    @target_density.setter
    def target_density(self, value):
        if not isinstance(value, FloatLike) and value is not None:
            raise TypeError("target_density must be a numeric value or None")
        if value is not None and value < 0:
            raise ValueError("target_density must be a positive number")
        self._target_density = np.float64(value)
        
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
