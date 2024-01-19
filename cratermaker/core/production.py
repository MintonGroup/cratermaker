import numpy as np
from numpy.random import Generator
from scipy.optimize import root_scalar
from cratermaker.utils.montecarlo import get_random_size
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import R_to_CSFD
from numpy.typing import ArrayLike
from typing import Union, Sequence, Tuple, Any
import warnings

class Production():
    """
    An operations class for computing the production function for craters and impactors. This implements a very simple power law 
    production function that can be used as either a crater or projectile production function. The production function is defined as
    the cumulative number of craters greater than a given diameter per unit m^2 surface area.

    """
    def __init__(self, 
                rng: Generator | None = None,
                **kwargs: Any):
          
        self.valid_generator_types = ["crater", "projectile"]
        self.rng = rng 
        self.set_model_parameters(**kwargs)
       
        return
       
    
    def set_model_parameters(self, **kwargs: Any) -> None:
        """
        Set the parameters for the power law production function.
        
        Parameters
        ----------
        **kwargs : Any
            This function accepts the following keyword arguments:
          
        model : str
            The specific model to use for the production function. Defaults to "Powerlaw". 
        generator_type : str
            The type of generator to use. This can be either "crater" or "projectile". Defaults to "crater". 
        N1_coef : float
            The coefficient for the power law production function at 1 m diameter per 1 My. 
            Defaults to 7.9.e-3 (lunar craters) or 2.2e-8 (lunar impactors) based on fits to the NPF on the Moon.
        slope : float
            The slope of the power law production function. 
            Defaults to -3.33 (lunar craters) or -2.26 (lunar impactors) based on fits to the NPF on the Moon.
        """
        if not hasattr(self, "valid_models"):
            self.valid_models = ["Powerlaw"] 
        model = kwargs.get("model", "Powerlaw")
        self.model = model
        
        # Set the generator type. For the default generator, it can be either "crater" or "projectile" 
        generator_type = kwargs.get("generator_type", "crater")
        self.generator_type = generator_type
        self.valid_time = (0,None)  # Range over which the production function is valid       
        
        # Default values that are approximately equal to the NPF for the Moon
        default_N1_coef = {
            "crater" : 7.883e-3, 
            "projectile" : 7.989e-7
            }
       
        default_slope = {
            "crater" : -3.328, 
            "projectile" : -2.634
            } 
        # Set the power law parameters for the production function along with defaults 
        N1_coef = kwargs.get("N1_coef",default_N1_coef[self.generator_type] )
        
        if not isinstance(N1_coef, FloatLike):
            raise ValueError("N1_coef must be a float")
        if N1_coef < 0.0:
            raise ValueError("N1_coef must be positive")
        self.N1_coef = N1_coef
       
        # Set the power law exponent for the production function along with defaults 
        slope = kwargs.get("slope", default_slope[self.generator_type])
        if not isinstance(slope, FloatLike):
            raise ValueError("slope must be a float")   
        elif slope > 0.0: # Slope must be negative, but convention in the field is mixed. So we flip the sign if it is positive.
            slope *= -1
        self.slope = slope 
        
        self.impact_velocity_model = kwargs.get("impact_velocity_model", "Moon")
        
    
    def set_mean_impact_velocity(self, 
                                 **kwargs: Any
                                ) -> None:
        """
        Sets an appropriate mean impact velocity for a target body.
        
        Parameters
        ----------
        **kwargs : Any
            This function accepts the following keyword arguments:        
        impact_velocity_model : str
            The name of the mean impact velocity model to use for the impact simulation.  Valid options are "Mercury_MBA", "Venus_MBA", "Earth_MBA", "Moon_MBA", "Mars_MBA", and "MBA_MBA". If None, the mean_impact_velocity is not set and a warning will be raised and Projectile <-> Crater scalings will not be performed.
            
        Returns
        -------
        np.float64 or None
            The mean impact velocity in m/s. If the target_name is None, returns None.
        
        
        Notes
        ----- 
        Mean velocities for terrestrial planets and the Moon are based on analysis of simulations of main-belt derived asteroids from Minton & Malhotra (2010) [1]_  and Yue et al. (2013) [2]_. Mean velocities for the asteroids are from Bottke et al. (1994) [3]_.
        
        References
        ----------
        .. [1] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744–757. https://doi.org/10.1016/j.icarus.2009.12.008
        .. [2] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
        .. [3] Bottke, W.F., Nolan, M.C., Greenberg, R., Kolvoord, R.A., 1994. Velocity distributions among colliding asteroids. Icarus 107, 255–268. https://doi.org/10.1006/icar.1994.1021
        """
     
        predefined_models = ['Mercury_MBA', 'Venus_MBA', 'Earth_MBA', 'Moon_MBA', 'Mars_MBA', 'MBA_MBA']
        predefined_velocities = [41100.0, 29100.0, 24600.0, 22100.0, 10700.0, 5300.0]
        predefined = dict(zip(predefined_models, predefined_velocities))
       
        impact_velocity_model = kwargs.get("impact_velocity_model", None) 
        if impact_velocity_model is None:
            warnings.warn("No impact_velocity_model provided. mean_impact_velocity not set.")
            return
        
        if impact_velocity_model not in predefined_models:
            warnings.warn(f"impact_velocity_model {impact_velocity_model} is not one of {predefined_models}. mean_impact_velocity not set.")

            return
        return np.float64(predefined[impact_velocity_model])
            
       
    def function(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
        """
        Return the cumulative size-frequency distribution of craters over a given age range and crater diameter for a simple power
        law model.

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
        age : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike or ArrayLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        **kwargs : Any
            Any additional keywords. These are not used in this base class, but included here so that any extended class can share
            the same function signature.

        Returns
        -------
        FloatLike or numpy array of FloatLike
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given age range.
        """         
        diameter, _ = self._validate_csfd(diameter=diameter)   
        age, age_end = self._validate_age(age, age_end)
        return self.N1_coef * diameter**self.slope * (age - age_end) 
    

    def function_inverse(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike,
             cumulative_number: FloatLike | Sequence[FloatLike] | ArrayLike,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
      
        """
        Return the age in My for a given number density of craters and diameter 

        Parameters
        ----------
        diameter : float-lik or  array-like
            diameter of the crater in m
        cumulative_number : float-like or array-like
            number density of craters per m^2 surface area greater than the input diameter
        **kwargs: Any
            Any additional keywords that are passed to the function method.

        Returns
        -------
        float_like or numpy array
            The age in My for the given relative number density of craters. 
        """
       
        diameter, cumulative_number = self._validate_csfd(diameter=diameter, cumulative_number=cumulative_number) 
        
        def _root_func(t,D,N):
            retval = self.function(diameter=D,age=t,check_valid_time=False,**kwargs) - N
            return retval
             
        xtol = 1e-10
        x0 = 4400.0
        retval = []
        darr = np.array(diameter)
        converged = []
        flag = []
        for i,d in np.ndenumerate(darr):
            sol = root_scalar(lambda x: _root_func(x,d,cumulative_number[i]), x0=x0, xtol=xtol, method='brentq', bracket=[0,10000]) 
            retval.append(sol.root)
            converged.append(sol.converged)
            flag.append(sol.flag)
        retval = np.array(retval)
        converged = np.array(converged)
        flag = np.array(flag)
        if np.all(converged):
            return retval.item() if np.isscalar(diameter) else retval
        else:
            raise ValueError(f"The root finding algorithm did not converge for all values of diameter and cumulative_number. Flag {flag}")


    def sample(self,
               age: FloatLike | None = None,
               age_end: FloatLike | None = None,
               cumulative_number_at_diameter: PairOfFloats | None = None,
               reference_cumulative_number_at_diameter: PairOfFloats | None = None,
               diameter_range: PairOfFloats | None = None,
               area: FloatLike | None = None, 
               ) -> np.ndarray:
        
        """
        Sample crater diameters and ages from the production function. This function can either sample from a given age range or
        from a given cumulative number/diameter pair (but not both). 
       
        Parameters
        ----------
        age : FloatLike or ArrayLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike or ArrayLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        cumulative_number_at_diameter : PairOfFloats, optional
            A pair of cumulative number and diameter values, in the form of a (N, D). If provided, the function convert this value
            to a corresponding age and use the production function for a given age.
        reference_cumulative_number_at_diameter : PairOfFloats, optional
            A pair of cumulative number and diameter values, in the form of a (N, D). If provided, the function will convert this
            value to a corresponding reference age and use the production function for a given age.
        diameter_range : PairOfFloats
            The minimum and maximum crater diameter to sample from in meters. 
        area : FloatLike, optional
            The area in m^2 over which the production function is evaluated to generate the expected number, which is the production
            function over the input age/cumulative number range at the minimum diameter.
            
        Returns
        -------
        numpy array
            The sampled diameter values
        numpy array
            The sampled age values
            
        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - Both the age and cumulative_number_at_diameter arguments are provided
            - The diameter_range is not provided.
            - If diameter_range is not a pair of values.
            - If the minimum diameter is less than or equal to 0.
            - If the maximum diameter is less than or equal the minimum.
            - If the cumulative_number_at_diameter argument is not a pair of values, or any of them are less than 0
            - If the reference_cumulative_number_at_diameter argument is not a pair of values, or any of them are less than 0
            - If if the area argument is less than 1. 

        """
        # Validate all the input arguments
        if age is None and cumulative_number_at_diameter is None:
            raise ValueError("Either the 'age' or 'cumulative_number_at_diameter' must be provided")
        elif age is not None and cumulative_number_at_diameter is not None:
            raise ValueError("Only one of the 'age' or 'cumulative_number_at_diameter' arguments can be provided")
        if age_end is not None and reference_cumulative_number_at_diameter is not None: 
            raise ValueError("Only one of the 'age_end' or 'reference_cumulative_number_at_diameter' arguments can be provided") 
        
        if age is not None and not np.isscalar(age):
            raise ValueError("The 'age' must be a scalar")
        
        if diameter_range is None:
            raise ValueError("The 'diameter_range' must be provided")
        if len(diameter_range) != 2:
            raise ValueError("The 'diameter_range' must be a pair of values")
        if diameter_range[0] <= 0:
            raise ValueError(f"Diameter range invalid: {diameter_range}. The minimum diameter must be greater than 0")
        if diameter_range[1] < diameter_range[0]:
            raise ValueError(f"Diameter range invalid: {diameter_range}. The maximum diameter must be greater than or equal to the minimum diameter")
        
        if area is not None:
            if not np.isscalar(area):
                raise ValueError("The 'area' must be a scalar")
            if area < 0.0:
                raise ValueError("The 'area' must be greater than 0")
       
        if reference_cumulative_number_at_diameter is not None:
            if len(reference_cumulative_number_at_diameter) != 2:
                raise ValueError("The 'reference_cumulative_number_at_diameter' must be a pair of values")
            reference_cumulative_number_at_diameter = self._validate_csfd(*reference_cumulative_number_at_diameter)
            age_end = self.function_inverse(diameter=reference_cumulative_number_at_diameter[1], cumulative_number=reference_cumulative_number_at_diameter[0])
            # Check to be sure that the diameter in the reference_cumulative_number_at_diameter is the same as cumulative_number_at_diameter. 
            # If not, we need to adjust it so that they match 
            if cumulative_number_at_diameter is not None:
                if len(cumulative_number_at_diameter) != 2:
                    raise ValueError("The 'cumulative_number_at_diameter' must be a pair of values")
                if reference_cumulative_number_at_diameter[1] != cumulative_number_at_diameter[1]:
                    reference_cumulative_number_at_diameter[0] = self.function(diameter=cumulative_number_at_diameter[1], age=age_end)
                    reference_cumulative_number_at_diameter[1] = cumulative_number_at_diameter[1]
                # Now adjust the N value so that when we convert it to age, we get the age value relative to the present
                cumulative_number_at_diameter[0] += reference_cumulative_number_at_diameter[0] 
            
            cumulative_number_at_diameter = self._validate_csfd(*cumulative_number_at_diameter)
            age = self.function_inverse(diameter=cumulative_number_at_diameter[1], cumulative_number=cumulative_number_at_diameter[0])
        
        age, age_end = self._validate_age(age, age_end)
           
        # Build the cumulative distribution function from which we will sample 
        input_diameters = np.logspace(np.log10(diameter_range[0]), np.log10(diameter_range[1]))
        cdf = self.function(diameter=input_diameters, age=age, age_end=age_end)
        expected_num = cdf[0] * area if area is not None else None
        diameters = get_random_size(diameters=input_diameters, cdf=cdf, mu=expected_num, rng=self.rng)
        
        # Get the corresponding age values for the sampled diameters
        ages = np.empty_like(diameters)
        
        # First divide up the interval into smaller subintervals to capture the curvature of the production function
        age_subinterval = np.linspace(age_end, age, num=1000) 
        for i, d in enumerate(diameters):
            N_vs_age = self.function(diameter=d, age=age_subinterval) # Compute the chronology curve for this size crater over the interval
            weights = N_vs_age - np.min(N_vs_age) # Shift the curve so that the minimum value is 0
            weights = weights[:-1] # Remove the last value
            weights /= np.sum(weights) # Normalize the weights so that they sum to 1
            chosen_subinterval_index = self.rng.choice(range(len(age_subinterval) - 1), p=weights) # Choose a subinterval based on the weights
            age_lo = age_subinterval[chosen_subinterval_index]
            age_hi = age_subinterval[chosen_subinterval_index + 1]
            # Now sample uniformly from the subinterval
            ages[i] = self.rng.uniform(low=age_lo, high=age_hi)
            
        # Sort the ages and diameters so that they are in order of increasing age
        sort_indices = np.argsort(ages)
        diameters = diameters[sort_indices]
        ages = ages[sort_indices]
        
        return diameters, ages
    

    @property
    def N1_coef(self):
        """Get the N1 coefficient of the power law production function."""
        return self._N1_coef

    @N1_coef.setter
    def N1_coef(self, value):
        """Set the N1 coefficient of the power law production function."""
        if not isinstance(value, FloatLike):
            raise TypeError("N1_coef must be a numeric value (float or int)")
        if value < 0:
            raise ValueError("N1_coef must be positive")
        self._N1_coef = value

    @property
    def slope(self):
        """Get the slope of the power law production function."""
        return self._slope

    @slope.setter
    def slope(self, value):
        """Set the slope of the power law production function."""
        if not isinstance(value, FloatLike):
            raise TypeError("slope must be a numeric value (float or int)")
        self._slope = value

    @property
    def valid_generator_types(self):
        """Get the list of valid generator types."""
        return self._valid_generator_types

    @valid_generator_types.setter
    def valid_generator_types(self, value):
        """Set the list of valid generator types."""
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            raise TypeError("valid_generator_types must be a list of strings")
        self._valid_generator_types = value

    @property
    def impact_velocity_model(self):
        """Get the impact velocity model name."""
        return self._impact_velocity_model

    @impact_velocity_model.setter
    def impact_velocity_model(self, value):
        """Set the impact velocity model name."""
        if not isinstance(value, str):
            raise TypeError("impact_velocity_model must be a string")
        self._impact_velocity_model = value

    @property
    def valid_time(self):
        """Get the valid time range for the production function."""
        return self._valid_time

    @valid_time.setter
    def valid_time(self, value):
        """Set the valid time range for the production function."""
        if not (isinstance(value, tuple) and len(value) == 2):
            raise ValueError("valid_time must be a tuple of two values")
        self._valid_time = value
        
    

    def _validate_age(self, 
                       age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
                       age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
                       ) -> Union[FloatLike, ArrayLike]:
        """
        Processes the age argument and age_end arguments. Checks that they are valid and returns a tuple of age and age_end.

        Parameters
        ----------        
        age : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike or ArrayLike, optional
            The reference used when computing age in My. If none is passed, it will be set either 0 or an array of zeros, depending on 
            the size of age.  

        Returns
        ------- 
        Tuple of np.float64
            The start and end ages in units of My.
            
        Raises
        ------
        ValueError
            If the the start age is greater than the end age or the age variable is not a scalar or a sequence of 2 values.
        """
          
        if np.isscalar(age):
            if not isinstance(age, FloatLike):
                raise TypeError("age must be a numeric value (float or int)")
            age = np.float64(age)
            
            if age_end is None:
                age_end = np.float64(0.0)
            else:
                if not np.isscalar(age_end):
                    raise ValueError("If age is a scalar, age_end must be a scalar")
                elif not isinstance(age_end, FloatLike):
                    raise TypeError("age_end must be a numeric value (float or int)")
                age_end = np.float64(age_end)
            if age < age_end:
                raise ValueError("age must be greater than or equal to the age_end")
        elif isinstance(age, (list, tuple, np.ndarray)):
            age = np.array(age, dtype=np.float64)
            if age_end is None:  
                age_end = np.zeros_like(age)
            elif isinstance(age_end, (list, tuple, np.ndarray)):
                age_end = np.array(age_end, dtype=np.float64)
            else:
                raise ValueError("If age is a sequence, age_end must be a sequence")
            if age.size != age_end.size:
                raise ValueError("If age is a sequence, age_end must be a sequence of the same size")
            if np.any(age < age_end):
                raise ValueError("age must be greater than the age_end")
        else:
            raise ValueError("age must be a scalar or a sequence")
       
        return age, age_end 


    def _validate_csfd(self,
                        diameter: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
                        cumulative_number: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
                       ) -> Tuple[Union[FloatLike, ArrayLike], Union[FloatLike, ArrayLike]]:
        """
        Validates the diameter and cumulative_number arguments. Both arguments can be either
        scalar or array-like, but they must be both scalars or both arrays of the same length.
        Values must be non-negative.

        Parameters
        ----------
        diameter : float-like or array-like, optional
            Diameter of the crater in meters.
        cumulative_number : float-like or array-like, optional
            Number density of craters per square meter surface area greater than the input diameter.

        Raises
        ------
        ValueError
            If any of the conditions on the CSFD are not met.
        """
        
        if diameter is None and cumulative_number is None:
            raise ValueError("Either the 'diameter' or 'cumulative_number' must be provided")
        
        # Convert inputs to numpy arrays for uniform processing
        if diameter is not None:
            diameter_array = np.atleast_1d(diameter)
        else:
            diameter_array = None
        if cumulative_number is not None:
            cumulative_number_array = np.atleast_1d(cumulative_number)   
        else:
            cumulative_number_array = None

        # Check if both are provided, they should be of the same length
        if diameter_array is not None and cumulative_number_array is not None:
            # Check for length consistency
            if (len(diameter_array) != len(cumulative_number_array)):
                raise ValueError("The 'diameter' and 'cumulative_number' must have the same length when both are provided")
            
        # Validate non-negative values
        if diameter_array is not None and np.any(diameter_array < 0):
            raise ValueError("All values in 'diameter' must be non-negative")
        if cumulative_number_array is not None and np.any(cumulative_number_array < 0):
            raise ValueError("All values in 'cumulative_number' must be non-negative")

        if diameter is not None and not np.isscalar(diameter):
            diameter = diameter_array
        if cumulative_number is not None and not np.isscalar(cumulative_number):
            cumulative_number = cumulative_number_array
        return diameter, cumulative_number
            
    @property
    def rng(self):
        """
        A random number generator instance. If not provided, the default numpy RNG will be used.
        """
        return self._rng

    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()

    @property
    def model(self):
        """
        The specific model to use for the production function. Defaults to 'Powerlaw'.
        """
        return self._model

    @model.setter
    def model(self, value):
        """
        Validates the given model string against the list of valid models.

        Parameters
        ----------
        value : str | None
            The model name to validate. If None, the first model in the valid models list is returned.

        Raises
        ------
        ValueError
            If the model is not a string or if the model is not in the list of valid models.
        """       
        if not value:
            self._model = self.valid_models[0]
            return
        if not isinstance(value, str):
            raise ValueError("model must be a string")
        value = value.capitalize()
        if value not in self.valid_models:
            raise ValueError(f"Invalid model {value}. Must be one of {self.valid_models}")
        self._model = value
    
    @property
    def valid_models(self):
        """
        A list of valid models for the production function. 
        These models define the parameters used for calculating the size-frequency distribution 
        of craters and impactors
        """
        return self._valid_models

    @valid_models.setter
    def valid_models(self, value):
        if not isinstance(value, list):
            raise TypeError("valid_models must be a list of strings")
        if any(not isinstance(model, str) for model in value):
            raise ValueError("All items in valid_models must be strings")
        self._valid_models = value    
       
    @property
    def generator_type(self):
        """
        The type of generator to use. This can be either "crater" or "projectile". 
        This determines the nature of the production function, differentiating between
        'crater' and 'projectile' types, which affect the calculations and interpretations
        of the production function.
        """
        return self._generator_type

    @generator_type.setter
    def generator_type(self, value):
        if not value:
            self._generator_type = self.valid_generator_types[0]
            return
        if not isinstance(value, str):
            raise ValueError("generator_type must be a string")
        if value not in self.valid_generator_types:
            raise ValueError(f"Invalid generator_type {value}. Must be one of {self.valid_generator_types}")
        self._generator_type = value
        return 
        
class NeukumProduction(Production):
    """
    An operations class for computing the the Neukum production function for the Moon and Mars.

    Parameters
    ----------
    model : {"Moon", "Mars", "Projectile"}, optional
        The specific model to use for the production function. "Moon" and "Mars" are both crater production functions, and
        "Projectile" is a projectile function. Defaults to "Moon".
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    **kwargs : Any
        Includes arguments that were called from the parent class. These are not used in this class.
        
    Notes
    ----- 
    The CSFD is computed using the model of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 for the Moon and Mars, with 
    minor changes. Notably, there is a typo in the chronology function (Eq. 5) of the original paper. The linear term in the paper
    is given as 8.38e-4. The value should be 10^(a0), and therefore the number given in the paper is based on the "Old"
    coefficients from Neukum (1983). The correct value is 10^(-3.0876) = 8.17e-4. We compute the value from the coefficients 
    in our implementation of the chronology function.       
   
    References
    ---------- 
    Lunar PF from: Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to 
        the Lunar Reference System. Space Science Reviews 96, 55–86. https://doi.org/10.1023/A:1011989004263
    Mars PF from: Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. Space Science Reviews 96, 87–104.
        https://doi.org/10.1023/A:1011941121102
    Projectile PF from: Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters 
        and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1–34. 
        https://doi.org/10.1007/978-94-010-0712-2_1
    """
        
    def set_model_parameters(self, **kwargs: Any) -> None:
        """
        Set the parameters for Neukum production. This will set the following attributes based on the value of the keyword argument
        `model`, which is either "Moon", "Mars", or "Projectile".
        
        - sfd_coef : the coefficients for the size-frequency distribution function (See Table 1 of Neukum et al. 2001)
        - sfd_range : the range of diameters over which the size-frequency distribution function is valid 
        - valid_time : the range of ages over which the chronology function is valid
        - tau : the time constant for the chronology function (See Eq. 5 of Neukum et al. 2001)
        - Cexp : the coefficient for the exponential componenbt of the chronology function (See Eq. 5 of Neukum et al. 2001)
        - Clin : the coefficient for the linear component of the chronology function (See Eq. 5 of Neukum et al. 2001, but our implementation corrects the typo in that expression)
        
        Parameters
        ----------
        **kwargs : Any
            This function accepts the following keyword arguments:
          
            model : str, {"Moon", "Mars", "Projectile"}
                The specific model to use for the production function. Defaults to "Moon" 
        """
        # Set the generator type. For the default generator, it can be either "crater" or "projectile" 
        if not hasattr(self, "valid_models"):
            self.valid_models = ["Moon", "Mars", "Projectile"]
        model = kwargs.get("model", "Moon")
        self.model = model
        
        if self.model == "Projectile":
            self.generator_type = "projectile"
        else:
            self.generator_type = "crater"

        sfd_coef = {
                "Moon" : np.array(
                    [
                        -3.0876,
                        -3.557528,
                        +0.781027,
                        +1.021521,
                        -0.156012,
                        -0.444058,
                        +0.019977,
                        +0.086850,
                        -0.005874,
                        -0.006809,
                        +8.25e-4, 
                        +5.54e-5
                    ]),
                "Mars" : np.array(
                    [
                        -3.384, 
                        -3.197,
                        +1.257,
                        +0.7915,
                        -0.4861,
                        -0.3630,
                        +0.1016,
                        +6.756e-2,
                        -1.181e-2,
                        -4.753e-3,
                        +6.233e-4,
                        +5.805e-5
                    ]),
                "Projectile" : np.array(
                    [
                        0,
                        +1.375458,
                        +1.272521e-1,
                        -1.282166,
                        -3.074558e-1,
                        +4.149280e-1,
                        +1.910668e-1,
                        -4.260980e-2,
                        -3.976305e-2,
                        -3.180179e-3,
                        +2.799369e-3,
                        +6.892223e-4,
                        +2.614385e-6,
                        -1.416178e-5,
                        -1.191124e-6
                    ]
                )
            }
        self.sfd_coef = sfd_coef[self.model]
        sfd_range = {
                "Moon" : np.array([0.01,1000]),
                "Mars" : np.array([0.015,362]),
                "Projectile" : np.array([0.0001, 200.0]) # Estimated based on Fig. 16 of Ivanov et al. (2001)
            }
        self.sfd_range = sfd_range[self.model]
        
        # Chronology function parameters
        self.valid_time = (0,4500)  # Range over which the production function is valid
        self.tau = 1.0 / 6.93
        Cexp_moon = 5.44e-14
        Clin = {
                "Moon" : 10**(sfd_coef.get("Moon")[0]),
                "Mars" : 10**(sfd_coef.get("Mars")[0]),
                "Projectile": 10**(sfd_coef.get("Projectile")[0]),
        }
        Cexp = {
                "Moon" : Cexp_moon,
                "Mars" : Cexp_moon * Clin["Mars"] / Clin["Moon"],
                "Projectile": Cexp_moon * Clin["Projectile"] / Clin["Moon"],
            }   
        self.Cexp = Cexp[self.model]
        self.Clin = Clin[self.model]
        

    def function(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
             check_valid_time: bool=True,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
        """
        Return the cumulative size-frequency distribution of craters over a given age range and crater diameter.

        Parameters
        ----------
        diameter : FloatLike or numpy array
            Crater diameter(s) in units of meters to compute corresponding cumulative number density value.
        age : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike or ArrayLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        check_valid_time : bool, optional (default=True)
            If True, return NaN for age values outside the valid age range

        Returns
        -------
        FloatLike or numpy array of FloatLike
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given age range.
        """
        age, age_end = self._validate_age(age, age_end) 
        diameter, _ = self._validate_csfd(diameter=diameter)

        return self._size_frequency_distribution(diameter) * (self._chronology(age,check_valid_time) - self._chronology(age_end,check_valid_time))
    
     
    def _chronology(self,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             check_valid_time: bool=True
             ) -> Union[FloatLike, ArrayLike]:
        """
        Returns the relative number of craters produced over a given age range. This implements the chronology function given in
        Eq. 5 of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86, but takes in the age argument in the Cratermaker unit 
        system of My instead of Gy. The returned value is normalized to the number of craters greater than 1 km in diameter at the
        reference time of 1 Gy. 

        Parameters
        ----------
        age : FloatLike or ArrayLike, default=1.0
            Age in the past relative to the present day to compute cumulative SFD in units of My. 
        check_valid_time : bool, optional (default=True)
            If True, return NaN for age values outside the valid age range

        Returns
        -------
        FloatLike or numpy array of FloatLike
            The cumulative number of craters per square meter greater than the input diameter that would be expected to form on a 
            surface over the given age range.
            
        """     
        time_Gy = np.array(age) * 1e-3  # Convert age range from My to Gy ago for internal functions
        
        def _N1(age: FloatLike | Sequence[FloatLike] | ArrayLike,
                check_valid_time:bool=True
                ) -> Union[FloatLike, ArrayLike]:
            """
            Return the cumulative number of 1 km craters as a function of age in Gy. This is a direct implementation of Eq. 5 in
            Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 (with corrected coefficient for the linear term).  

            Parameters
            ----------
            age : FloatLike or numpy array
                Time ago in units of Gy
            check_valid_time : bool, optional (default=True)
                If True, return NaN for age values outside the valid age range        

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than 1 km in diameter
            """
            N1 = self.Cexp * (np.exp(age/self.tau) - 1.0) + self.Clin * age
            if check_valid_time:
                if self.valid_time[0] is not None:
                    min_time = self.valid_time[0] * 1e-3
                    N1 = np.where(age >= min_time, N1, np.nan)
                if self.valid_time[1] is not None:
                    max_time = self.valid_time[1] * 1e-3
                    N1 = np.where(age <= max_time, N1, np.nan) 
            return N1.item() if np.isscalar(age) else N1
        
        N1_reference = _N1(1.0) 
        N1_values = _N1(time_Gy,check_valid_time)
        N1_values /= N1_reference
        
        return  N1_values 
   
    
    def _size_frequency_distribution(self,diameter: FloatLike | ArrayLike,) -> Union[FloatLike, ArrayLike]:
        """
        Return the cumulative size frequency distribution of craters at a given age relative to age = 1 Gy ago per m^2.

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            Time in units of meter 
            
        Returns
        -------
        FloatLike or numpy array
           Cumulative number density of craters per square meter greater than the input diameter.
        """        
        def _extrapolate_sfd(side: str = "lo") -> Union[FloatLike, ArrayLike]:
            """
            Return the exponent, p, and and proportionality constant, A, for  the extrapolated 
            CSFD in the form N(D) = A * D**-p. 

            Parameters
            ----------
            side : str
                The side of the range to extrapolate. Valid values are "lo" and "hi"

            Returns
            -------
            A, p
            """    
            if side == "lo":
                idx = 0
            elif side == "hi":
                idx = 1
            else:
                raise ValueError("side must be 'lo' or 'hi'")
            p = _dNdD(self.sfd_range[idx])
            A = _CSFD(self.sfd_range[idx])
            return A, p


        def _dNdD(Dkm: FloatLike | Sequence[FloatLike] | ArrayLike) -> Union[FloatLike, ArrayLike]:
            """
            Return the derivative of the cumulative size-frequency distribution as a function of diameter. For diameter values outside 
            the range of the NPF, the derivative is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at age = 1 Gy ago.
            """        
            def _dNdD_scalar(Dkm): 
                dcoef = self.sfd_coef[1:]
                if Dkm < self.sfd_range[0]:
                    _extrapolate_sfd(side="lo")
                    return A * (p / Dkm) * (Dkm / self.sfd_range[0]) ** p 
                elif Dkm > self.sfd_range[1]:
                    _extrapolate_sfd(side="hi")
                    return A * (p / Dkm) * (Dkm / self.sfd_range[0]) ** p 
                else:
                    return sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))
            
            return _dNdD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_dNdD_scalar)(Dkm)

    
        def _CSFD(Dkm: FloatLike | Sequence[FloatLike] | ArrayLike) -> Union[FloatLike, ArrayLike]:
            """
            Return the cumulative size-frequency distribution at the reference age of 1 Gy ago. For diameter values outside 
            the range of the NPF, the CSFD is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than Dkm in diameter at age=1 Gy ago.
            """
            def _CSFD_scalar(Dkm):
                if Dkm < self.sfd_range[0]:
                    A, p = _extrapolate_sfd(side="lo")
                    return A * (Dkm / self.sfd_range[0]) ** p
                elif Dkm > self.sfd_range[1]:
                    A, p = _extrapolate_sfd(side="hi")
                    p -= 2.0 # Steepen the upper branch of the SFD to prevent anomolously large craters from forming
                    return A * (Dkm / self.sfd_range[1]) ** p
                else:
                    logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(self.sfd_coef))
                    return 10 ** logCSFD
        
            return _CSFD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_CSFD_scalar)(Dkm)


        if np.any(diameter < 0.0):
            raise ValueError("diameter must be greater than or equal to 0.0")
              
        Dkm = diameter * 1e-3 # Convert m to km for internal functions

        if self.model == "Projectile":
            Ncumulative = R_to_CSFD(R=_CSFD, D=Dkm) * 2.94e-5 # This is a multiplication factor that gets the projectile CSFD to approximately match the lunar crater CSFD 
        else:
            Ncumulative = _CSFD(Dkm) 
            
        return Ncumulative * 1e-6 # convert from km^-2 to m^-2    

    @property
    def sfd_coef(self):
        """
        Coefficients for the size-frequency distribution function used in the Neukum production model.
        """
        return self._sfd_coef

    @sfd_coef.setter
    def sfd_coef(self, value):
        """
        Set the coefficients for the size-frequency distribution function used in the Neukum production model.
        """
        if not isinstance(value, np.ndarray):
            raise TypeError("sfd_coef must be a numpy array")
        self._sfd_coef = value

    @property
    def sfd_range(self):
        """
        Range of diameters over which the size-frequency distribution function is valid in the Neukum production model.
        """
        return self._sfd_range

    @sfd_range.setter
    def sfd_range(self, value):
        """
        Set the range of diameters over which the size-frequency distribution function is valid in the Neukum production model.
        """
        if not (isinstance(value, (list, tuple, np.ndarray)) and len(value) == 2):
            raise ValueError("sfd_range must be a list, tuple, or numpy array of length 2")
        self._sfd_range = value

    @property
    def valid_time(self):
        """
        Range of ages over which the chronology function is valid in the Neukum production model.
        """
        return self._valid_time

    @valid_time.setter
    def valid_time(self, value):
        """
        Set the range of ages over which the chronology function is valid in the Neukum production model.
        """
        if not (isinstance(value, tuple) and len(value) == 2):
            raise ValueError("valid_time must be a tuple of two values")
        self._valid_time = value

    @property
    def tau(self):
        """Get the time constant for the chronology function."""
        return self._tau

    @tau.setter
    def tau(self, value):
        """Set the time constant for the chronology function."""
        if not isinstance(value, (float, int)):
            raise TypeError("tau must be a numeric value")
        self._tau = value

    @property
    def Cexp(self):
        """Get the coefficient for the exponential component of the chronology function."""
        return self._Cexp

    @Cexp.setter
    def Cexp(self, value):
        """Set the coefficient for the exponential component of the chronology function."""
        if not isinstance(value, (float, int)):
            raise TypeError("Cexp must be a numeric value")
        self._Cexp = value

    @property
    def Clin(self):
        """Get the coefficient for the linear component of the chronology function."""
        return self._Clin

    @Clin.setter
    def Clin(self, value):
        """Set the coefficient for the linear component of the chronology function."""
        if not isinstance(value, (float, int)):
            raise TypeError("Clin must be a numeric value")
        self._Clin = value
        
if __name__ == "__main__":
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from cratermaker.core.target import Target
    
    def plot_npf_csfd():
        fig = plt.figure(1, figsize=(8, 7))
        ax = {'Moon': fig.add_subplot(121),
            'Mars': fig.add_subplot(122)}

        tvals = [0.01,1.0,4.0]
        x_min = 1e-3
        x_max = 1e4
        y_min = 1e-12
        y_max = 1e6
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        for key in ax:
            production = NeukumProduction(model=key)
            ax[key].title.set_text(key)
            ax[key].set_xscale('log')
            ax[key].set_yscale('log')
            ax[key].set_ylabel('$\\mathregular{N_{>D} (km^{-2})}$')
            ax[key].set_xlabel('Diameter (km)')
            ax[key].set_xlim(x_min, x_max)
            ax[key].set_ylim(y_min, y_max)
            ax[key].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].grid(True,which="minor",ls="-",lw=0.5,zorder=5)
            ax[key].grid(True,which="major",ls="-",lw=1,zorder=10)
            inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
            lo = Dvals < production.sfd_range[0]
            hi = Dvals > production.sfd_range[1]
            for t in tvals:
                Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
                Nvals *= 1e6 # convert from m^-2 to km^-2
                ax[key].plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
                ax[key].plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
                ax[key].plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)
                labeli = int(0.25*nD)
                ax[key].text(Dvals[labeli],3*Nvals[labeli],f"{t:.2f} ", ha="left", va="top",rotation=-72)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()

        
    def plot_npf_N1_vs_T():
        fig = plt.figure(1, figsize=(8, 4))
        ax = fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_ylabel('$\\mathregular{N(1) (km^{-2})}$')
        ax.set_xlabel('Time (Gy ago)')
        ax.set_xlim(4.5, 0)
        moon = NeukumProduction(model="Moon")
        mars = NeukumProduction(model="Mars")
        tvals = np.linspace(4.5, 0.0, num=1000)
        N1_moon = moon.function(diameter=1000.0, age=tvals*1e3)*1e6
        N1_mars = mars.function(diameter=1000.0, age=tvals*1e3)*1e6
        ax.plot(tvals, N1_moon, '-', color='dimgrey', linewidth=2.0, zorder=50, label="Moon")
        ax.plot(tvals, N1_mars, '-', color='orange', linewidth=2.0, zorder=50, label="Mars")
        ax.legend()
        plt.tight_layout()
        plt.show() 


    def plot_npf_fit():
        # Define the power-law function
        def power_law(D, N1_coef, slope):
            return N1_coef * D ** slope
       
        # Create the lunar crater and projectile production functions 
        crater_production = NeukumProduction(model="Moon")
        projectile_production = NeukumProduction(model="Projectile")

        # Make a population of sizes that spans the size of interest
        Dc = np.logspace(-1,2)
        Di = np.logspace(-2,1)

        # Compute the reference N values (N>D at 1 My)
        Nc = crater_production.function(diameter=Dc,age=1.0)
        Ni = projectile_production.function(diameter=Di,age=1.0)
       
        D1proj = 22.3 # Approximate diameter of crater made by a 1 m projectile using the default Cratermaker scaling relationships
         
        N1_c = crater_production.function(diameter = D1proj, age = 1.0)
        N1_p = projectile_production.function(diameter = 1.0, age = 1.0)        

        # Fit the power-law function to the crater data
        params_crater, _ = curve_fit(power_law, Dc, Nc)

        # Fit the power-law function to the projectile data
        params_projectile, _ = curve_fit(power_law, Di, Ni * N1_c / N1_p)

        # Extracting the parameters
        N1_coef_crater, slope_crater = params_crater
        N1_coef_projectile, slope_projectile = params_projectile

        # Output the results
        print("Crater Production Fit Parameters N1_coef =", N1_coef_crater, ", slope =", slope_crater)
        print("Projectile Production Fit Parameters N1_coef =", N1_coef_projectile, ", slope =", slope_projectile)
        
        # Create the fitted curve data
        fitted_curve_crater = power_law(Dc, *params_crater)
        fitted_curve_projectile = power_law(Di, *params_projectile)

        # Plotting the crater data
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.loglog(Dc, Nc, 'o', label='Original Data (Crater)')
        plt.loglog(Dc, fitted_curve_crater, '-', label='Fitted Curve (Crater)')
        plt.xlabel('Crater Diameter (m)')
        plt.ylabel('N>D')
        plt.title('Crater Production')
        plt.legend()

        # Plotting the projectile data
        plt.subplot(1, 2, 2)
        plt.loglog(Di, Ni * N1_c/N1_p, 'o', label='Original Data (Projectile)')
        plt.loglog(Di, fitted_curve_projectile, '-', label='Fitted Curve (Projectile)')
        plt.xlabel('Projectile Diameter (m)')
        plt.ylabel('N>D')
        plt.title('Projectile Production')
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()        


    def plot_npf_proj_csfd():
        fig = plt.figure(1, figsize=(4, 7))
        ax = fig.add_subplot(111)

        x_min = 1e-5
        x_max = 1e4
        y_min = 1e-7
        y_max = 1e13
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        production = NeukumProduction(model="Projectile")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('$\\mathregular{N_{>D}}$')
        ax.set_xlabel('Projectile Diameter (km)')
        ax.set_xlim(x_min, x_max)
        #ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
        ax.grid(True,which="minor",ls="-",lw=0.5,zorder=5)
        ax.grid(True,which="major",ls="-",lw=1,zorder=10)
        inrange = (Dvals >= production.sfd_range[0]) & (Dvals <= production.sfd_range[1])
        lo = Dvals < production.sfd_range[0]
        hi = Dvals > production.sfd_range[1]
        t = 1.0
        Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
        Nvals *= 1e6 # convert from m^-2 to km^-2
        ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
        ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()    
   
   
    def plot_sampled_csfd():
        fig = plt.figure(1, figsize=(8, 7))
        ax = {'Power Law': fig.add_subplot(121),
            'NPF (Moon)': fig.add_subplot(122)}

        production = {
                    'Power Law': Production(),
                    'NPF (Moon)': NeukumProduction(model="Moon")
                    }
        
        target = Target("Moon")
        area = 4*np.pi*target.radius**2 
        age = 4100.0
        x_min = 1e0
        x_max = 1e5
        y_min = 1e0
        y_max = 1e4
        diameter_range = (2e3,10000e3) # Range of diameters to generate in m
        nD = 1000
        Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
        Nevaluations = 100
        for key in ax:
            ax[key].title.set_text(key)
            ax[key].set_xscale('log')
            ax[key].set_yscale('log')
            ax[key].set_ylabel('$\\mathregular{N_{>D}}$')
            ax[key].set_xlabel('Diameter (km)')
            ax[key].set_xlim(x_min, x_max)
            ax[key].set_ylim(y_min, y_max)
            ax[key].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            ax[key].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
            ax[key].xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
            # Plot the sampled values
            for i in range(Nevaluations):
                Dsampled = production[key].sample(age=age, diameter_range=diameter_range, area=area)
                Dsampled = np.sort(Dsampled)[::-1]
                Nsampled = range(1,len(Dsampled)+1) 
                ax[key].plot(Dsampled*1e-3, Nsampled, '-', color='cornflowerblue', linewidth=2.0, zorder=50, alpha=0.2)
               
            # Plot the production function 
            Nvals = production[key].function(diameter=Dvals*1e3,age=age)
            Nvals *= area # convert from per unit area to total number
            ax[key].plot(Dvals, Nvals, '-', color='black', linewidth=3.0, zorder=50)

        plt.tick_params(axis='y', which='minor')
        plt.tight_layout()
        plt.show()
    
            
    plot_npf_csfd()
    plot_npf_N1_vs_T()
    plot_npf_fit()    
    plot_npf_proj_csfd()
    plot_sampled_csfd()
