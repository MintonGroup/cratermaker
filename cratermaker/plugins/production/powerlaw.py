import numpy as np
from numpy.random import Generator
from scipy.optimize import root_scalar
from cratermaker.plugins.production import register_production_model, ProductionModel
from cratermaker.utils.montecarlo import get_random_size
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import R_to_CSFD
from numpy.typing import ArrayLike
from collections.abc import Sequence
from typing import Any, Union

@register_production_model("powerlaw")
class PowerLawProduction(ProductionModel):
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
        self._set_model_parameters(**kwargs)
       
        return
       
    def _set_model_parameters(self, **kwargs: Any) -> None:
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
        mean_velocity : float
            The mean impact velocity to use for the impact simulation. Either mean_velocity or impact_velocity_model must be provided.
        impact_velocity_model : str
            The name of the mean impact velocity model to use for the impact simulation.  Valid options are "Mercury_MBA", "Venus_MBA", "Earth_MBA", "Moon_MBA", "Mars_MBA", and "MBA_MBA". 
            Either mean_velocity or impact_velocity_model must be provided.
        """
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
       
        if "mean_velocity" in kwargs and "impact_velocity_model" in kwargs:
            raise ValueError("Only one of 'mean_velocity' or 'impact_velocity_model' can be provided")
         
        if "mean_velocity" in kwargs:
            self.mean_velocity = kwargs["mean_velocity"]
        elif "impact_velocity_model" in kwargs:
            self.impact_velocity_model = kwargs.get("impact_velocity_model")
        else:
            raise ValueError("Either 'mean_velocity' or 'impact_velocity_model' must be provided")
      
       
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
        
        n_array = np.asarray(self.N1_coef * diameter**self.slope)
        age_difference = np.asarray(age - age_end)
        
        if n_array.ndim > 0 and age_difference.ndim > 0:
            return n_array[:, None] * age_difference
        else: 
            return n_array * age_difference
    

    def function_inverse(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike,
             cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]:
      
        """
        Return the age in My for a given number density of craters and diameter 

        Parameters
        ----------
        diameter : float-lik or  array-like
            diameter of the crater in m
        cumulative_number_density : float-like or array-like
            number density of craters per m^2 surface area greater than the input diameter
        **kwargs: Any
            Any additional keywords that are passed to the function method.

        Returns
        -------
        float_like or numpy array
            The age in My for the given relative number density of craters. 
        """
       
        diameter, cumulative_number_density = self._validate_csfd(diameter=diameter, cumulative_number_density=cumulative_number_density) 
        
        def _root_func(t,D,N):
            retval = self.function(diameter=D,age=t,check_valid_time=False,**kwargs) - N
            return retval
             
        xtol = 1e-10
        x0 = 4400.0
        retval = []
        darr = np.array(diameter)
        narr = np.array(cumulative_number_density)
        converged = []
        flag = []
        for i,d in np.ndenumerate(darr):
            sol = root_scalar(lambda x: _root_func(x,d,narr[i]), x0=x0, xtol=xtol, method='brentq', bracket=[0,10000]) 
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
               diameter_number: PairOfFloats | None = None,
               diameter_number_end: PairOfFloats | None = None,
               diameter_range: PairOfFloats | None = None,
               area: FloatLike | None = None, 
               return_age: bool = True
               ) -> np.ndarray:
        
        """
        Sample diameters and ages from the production function. This function can either sample from a given age range or
        from a given cumulative number/diameter pair (but not both). 
       
        Parameters
        ----------
        age : FloatLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        diameter_number : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function convert this value
            to a corresponding age and use the production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function will convert this
            value to a corresponding reference age and use the production function for a given age.
        diameter_range : PairOfFloats
            The minimum and maximum crater diameter to sample from in meters. 
        area : FloatLike, optional
            The area in m^2 over which the production function is evaluated to generate the expected number, which is the production
            function over the input age/cumulative number range at the minimum diameter.
        return_age : bool, optional
            If True, the function will return the sampled ages in addition to the diameters. The default is True.
            
        Returns
        -------
        numpy array
            The sampled diameter values
        numpy array or None
            The sampled age values if return_age is True, otherwise None.
            
        """
        arguments = locals().copy()
        arguments.pop("self")
        arguments = self._validate_sample_args(**arguments)
        age = arguments["age"]
        age_end = arguments["age_end"]
        diameter_range = arguments["diameter_range"]
        area = arguments["area"]
        return_age = arguments["return_age"]
           
        # Build the cumulative distribution function from which we will sample 
        input_diameters = np.logspace(np.log10(diameter_range[0]), np.log10(diameter_range[1]))
        cdf = self.function(diameter=input_diameters, age=age, age_end=age_end)
        expected_num = cdf[0] * area if area is not None else None
        diameters = np.asarray(get_random_size(diameters=input_diameters, cdf=cdf, mu=expected_num, rng=self.rng))
        if diameters.size == 0:
            return np.empty(0), np.empty(0)
        elif diameters.size == 1:
            diameters = np.array([diameters])
       
        if return_age: 
            age_subinterval = np.linspace(age_end, age, num=1000) 
            N_vs_age = np.asarray(self.function(diameter=diameters, age=age_subinterval))
            
            # Normalize the weights for each diameter
            if N_vs_age.ndim > 1:
                N_vs_age -= np.min(N_vs_age, axis=1)[:, None]  # Subtract the min along age dimension
                N_vs_age = N_vs_age[:, :-1]  # Remove the last value to avoid the cumulative sum reaching 1
                weights_sum = np.sum(N_vs_age, axis=1)[:, None]  # Sum of weights for each diameter
                weights = N_vs_age / weights_sum  # Normalize weights for each diameter
                
                # Compute the CDF for each diameter
                cdf = np.cumsum(weights, axis=1)

            else:
                N_vs_age -= np.min(N_vs_age)  # Subtract the min along age dimension
                N_vs_age = N_vs_age[:-1]  # Remove the last value to avoid the cumulative sum reaching 1
                weights_sum = np.sum(N_vs_age) # Sum of weights for each diameter
                weights = N_vs_age / weights_sum  # Normalize weights for each diameter
                
                # Compute the CDF for each diameter
                cdf = np.cumsum(weights)    

            # Generate uniform random numbers for each diameter
            random_values = self.rng.uniform(0, 1, size=diameters.size)
            if N_vs_age.ndim > 1:
                chosen_subinterval_indices = np.zeros(diameters.size, dtype=int)
                # For each diameter, find the corresponding index in the CDF
                for i, cdf_row in enumerate(cdf):
                    # Find the corresponding index in the CDF for each diameter
                    chosen_subinterval_indices[i] = np.searchsorted(cdf_row, random_values[i], side='right') - 1
            else:
                chosen_subinterval_indices = np.searchsorted(cdf, random_values[:], side='right') - 1 
                
            # Ensure indices are within valid range
            chosen_subinterval_indices = np.clip(chosen_subinterval_indices, 0, len(age_subinterval) - 2)
            # Sample a random age within the selected subinterval for each diameter
            age_lo = age_subinterval[chosen_subinterval_indices]
            age_hi = age_subinterval[chosen_subinterval_indices + 1]
            ages = self.rng.uniform(low=age_lo, high=age_hi) 

            # Sort the ages and diameters so that they are in order of decreasing age
            if ages.size > 1:
                sort_indices = np.argsort(ages)[::-1]
                diameters = diameters[sort_indices]
                ages = ages[sort_indices]
        else:
            ages = np.empty(0)
            
        return diameters, ages
    
    
    def _validate_sample_args(self,**kwargs: dict) -> dict:
        """
        Validate all the input arguments to the sample method. This function will raise a ValueError if any of the arguments are invalid.
        It will also convert age arguments to diameter_number and vice versa.
        
        Parameters
        ----------
        age : FloatLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. 
            The default is 0 (present day).
        diameter_number : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N), which gives the total cumulative number of 
            impactors, N, larger than diameter, D. If provided, the function convert this value to a corresponding age and use the 
            production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N), which gives the total cumulative number of 
            impactors, N, larger than diameter, D.. If provided, the function will convert this value to a corresponding age_end
            and use the production function for a given age. The default is (1000.0, 0) (present day). 
        diameter_range : PairOfFloats
            The minimum and maximum crater diameter to sample from in meters.
        area : FloatLike, optional
            The area in m^2 over which the production function is evaluated to generate the expected number, which is the production
            function over the input age/cumulative number range at the minimum diameter.
        return_age : bool, optional
            If True, the function will return the sampled ages in addition to the diameters. The default is True.
                
        Returns
        -------
        FloatLike
            The start age in units of My relative to the present.
        FloatLike
            The end age in units of My relative to the present.
        PairOfFloats
            A pair of diameter and cumulative number values, in the form of a (D, N) representing the cumulative number and diameter values for the start age.
        PairOfFloats
            A pair of diameter and cumulative number values, in the form of a (D, N) representing the cumulative number and diameter values for the end age.
        PairOfFloats
            The minimum and maximum diameter values to sample from in meters.
        FloatLike
            The area in m^2 over which the production function is evaluated to generate the expected number, which is the production
            function over the input age/cumulative number range at the minimum diameter.
            
        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - Neither the age nore the diameter_number argument is provided.
            - Both the age and diameter_number arguments are provided.
            - Both the age_end and diameter_number_end arguments are provided.
            - The age argument is provided but is not a scalar.
            - The age_end argument is provided but is not a scalar. 
            - The diameter_number argument is not a pair of values, or any of them are less than 0
            - The diameter_number_end argument is not a pair of values, or any of them are less than 0
            - The diameter_range is not provided.
            - The diameter_range is not a pair of values.
            - The minimum diameter is less than or equal to 0.
            - The maximum diameter is less than or equal the minimum.
            - The area argument is not a scalar or is less than 0.             
        """ 
        _REF_DIAM = 1000.0 # The diameter value used if age arguments are provided
        age = kwargs.get("age", None)
        age_end = kwargs.get("age_end", None)
        diameter_number = kwargs.get("diameter_number", None)
        diameter_number_end = kwargs.get("diameter_number_end", None)
        diameter_range = kwargs.get("diameter_range", None)
        area = kwargs.get("area", None)
        return_age = kwargs.get("return_age", True)
        
        if age is None and diameter_number is None:
            raise ValueError("Either the 'age' or 'diameter_number' must be provided")
        elif age is not None and diameter_number is not None:
            raise ValueError("Only one of the 'age' or 'diameter_number' arguments can be provided")
        if age_end is not None and diameter_number_end is not None: 
            raise ValueError("Only one of the 'age_end' or 'diameter_number_end' arguments can be provided")
        if age is not None and not np.isscalar(age):
            raise ValueError("The 'age' must be a scalar")
        if age_end is not None and not np.isscalar(age_end):
            raise ValueError("The 'age_end' must be a scalar")
            
        if diameter_range is None:
            raise ValueError("The 'diameter_range' must be provided")
        if len(diameter_range) != 2:
            raise ValueError("The 'diameter_range' must be a pair of values")
        if diameter_range[0] <= 0:
            raise ValueError(f"Diameter range invalid: {diameter_range}. The minimum diameter must be greater than 0")
        if diameter_range[1] < diameter_range[0]:
            raise ValueError(f"Diameter range invalid: {diameter_range}. The maximum diameter must be greater than or equal to the minimum diameter")
      
        if area is None:
            raise ValueError("The 'area' must be provided") 
        else:
            if not np.isscalar(area):
                raise ValueError("The 'area' must be a scalar")
            if area < 0.0:
                raise ValueError("The 'area' must be greater than 0")
           
        if diameter_number is not None:
            if len(diameter_number) != 2:
                raise ValueError("The 'diameter_number' must be a pair of values in the form (D, N)") 
            diameter_number = self._validate_csfd(*diameter_number)
            diameter_number_density = (diameter_number[0], diameter_number[1] / area)
            age = self.function_inverse(*diameter_number_density)
        else:
            diameter_number_density = (_REF_DIAM, self.function(diameter=_REF_DIAM, age=age))
            diameter_number = (diameter_number_density[0], diameter_number_density[1] * area)
            
        if age_end is None and diameter_number_end is None:
            diameter_number_end = (_REF_DIAM, 0.0)
            diameter_number_density_end = diameter_number_end
            age_end = 0.0
            
        if diameter_number_end is not None:
            if len(diameter_number_end) != 2:
                raise ValueError("The 'diameter_number_end' must be a pair of values in the form (D,N)")
            diameter_number_end = self._validate_csfd(*diameter_number_end)
            diameter_number_density_end = (diameter_number_end[0], diameter_number_end[1] / area)
            # Check to be sure that the diameter in the diameter_number_end is the same as diameter_number. 
            # If not, we need to adjust the end diameter value it so that they match 
            if diameter_number is not None:
                if diameter_number_density_end[0] != diameter_number_density[0]:
                    diameter_number_density_end = (diameter_number_density[0], self.function(diameter=diameter_number_density[0], age=age_end))
                    diameter_number_end = (diameter_number[0], diameter_number_density_end[1] * area)
        
        if age_end is None:
            diameter_number_end = self._validate_csfd(*diameter_number_end)
            diameter_number_density_end = (diameter_number_end[0], diameter_number_end[1] / area)
            age_end = self.function_inverse(*diameter_number_density_end)            
        
        age, age_end = self._validate_age(age, age_end)
        
        if diameter_number_end is None:
            diameter_number_end = (diameter_number[0], self.function(diameter=diameter_number[0], age=age_end) * area)
        
        if not isinstance(return_age, bool):
            raise ValueError("The 'return_age' argument must be a boolean")
           
        kwargs["age"] = age
        kwargs["age_end"] = age_end
        kwargs["diameter_number"] = diameter_number
        kwargs["diameter_number_end"] = diameter_number_end
        kwargs["diameter_range"] = diameter_range
        kwargs["area"] = area
        kwargs["return_age"] = return_age
        
        return kwargs
    
    
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
        tuple of np.float64
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
                        cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
                       ) -> tuple[Union[FloatLike, ArrayLike], Union[FloatLike, ArrayLike]]:
        """
        Validates the diameter and cumulative_number arguments. Both arguments can be either
        scalar or array-like, but they must be both scalars or both arrays of the same length.
        Values must be non-negative.

        Parameters
        ----------
        diameter : float-like or array-like, optional
            Diameter of the impacts 
        cumulative_number_density: float-like or array-like, optional
            Cumulative number of impacts larger than `diameter`

        Raises
        ------
        ValueError
            If any of the conditions on the CSFD are not met.
        """
        
        if diameter is None and cumulative_number_density is None:
            raise ValueError("Either the 'diameter' or 'cumulative_number_density' must be provided")
        
        # Convert inputs to numpy arrays for uniform processing
        if diameter is not None:
            diameter_array = np.atleast_1d(diameter)
        else:
            diameter_array = None
        if cumulative_number_density is not None:
            cumulative_number_density_array = np.atleast_1d(cumulative_number_density)   
        else:
            cumulative_number_density_array = None

        # Check if both are provided, they should be of the same length
        if diameter_array is not None and cumulative_number_density_array is not None:
            # Check for length consistency
            if (len(diameter_array) != len(cumulative_number_density_array)):
                raise ValueError("The 'diameter' and 'cumulative_number_density' must have the same length when both are provided")
            
        # Validate non-negative values
        if diameter_array is not None and np.any(diameter_array < 0):
            raise ValueError("All values in 'diameter' must be non-negative")
        if cumulative_number_density_array is not None and np.any(cumulative_number_density_array < 0):
            raise ValueError("All values in 'cumulative_number_density' must be non-negative")

        if diameter is not None and not np.isscalar(diameter):
            diameter = diameter_array
        if cumulative_number_density is not None and not np.isscalar(cumulative_number_density):
            cumulative_number_density = cumulative_number_density_array
        return diameter, cumulative_number_density
    
    
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
        """"
        Set the name of the mean impact velocity model to use for the impact simulation.  Valid options are "Mercury_MBA", "Venus_MBA", "Earth_MBA", "Moon_MBA", "Mars_MBA", and "MBA_MBA". 
            
        Notes
        ----- 
        Mean velocities for terrestrial planets and the Moon are based on analysis of simulations of main-belt derived asteroids from Minton & Malhotra (2010) [1]_  and Yue et al. (2013) [2]_. Mean velocities for the asteroids are from Bottke et al. (1994) [3]_.
        
        References
        ----------
        .. [1] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744-757. https://doi.org/10.1016/j.icarus.2009.12.008
        .. [2] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
        .. [3] Bottke, W.F., Nolan, M.C., Greenberg, R., Kolvoord, R.A., 1994. Velocity distributions among colliding asteroids. Icarus 107, 255-268. https://doi.org/10.1006/icar.1994.1021
        """
     
        predefined_models = ['Mercury_MBA', 'Venus_MBA', 'Earth_MBA', 'Moon_MBA', 'Mars_MBA', 'MBA_MBA']
        predefined_velocities = [41100.0, 29100.0, 24600.0, 22100.0, 10700.0, 5300.0]
        predefined = dict(zip(predefined_models, predefined_velocities))
        if value is None:
            if self._mean_velocity is None:
                raise ValueError("impact_velocity_model must be set if mean_velocity is not set")
        elif not isinstance(value, str):
            raise TypeError("impact_velocity_model must be a string")
        elif value not in predefined_models:
            raise ValueError(f"impact_velocity_model {value} is not one of {predefined_models}")
        self._impact_velocity_model = value
        self._mean_velocity = np.float64(predefined[value])
       
    @property 
    def mean_velocity(self):
        """The mean impact velocity for the production function."""
        return self._mean_velocity
    
    @mean_velocity.setter
    def mean_velocity(self, value):
        if value is None:
            if self._impact_velocity_model is None:
                raise ValueError("mean_velocity must be set if impact_velocity_model is not set")
        if not isinstance(value, FloatLike):
            raise TypeError("mean_velocity must be a numeric value (float or int)")
        if value <= 0.0:
            raise ValueError("mean_velocity must be finite and positive")
        self._mean_velocity = np.float64(value)

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