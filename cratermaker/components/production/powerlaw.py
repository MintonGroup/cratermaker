import numpy as np
from numpy.random import Generator
from cratermaker.components.production import register_production_model, ProductionModel
from cratermaker.utils.custom_types import FloatLike
from numpy.typing import ArrayLike
from collections.abc import Sequence
from typing import Any, Union

@register_production_model("powerlaw")
class PowerLawProduction(ProductionModel):
    """
    An operations class for computing the production function for craters and impactors. This impliments a very simple power law 
    production function that can be used as either a crater or projectile production function. The production function is defined as
    the cumulative number of craters greater than a given diameter per unit m^2 surface area.
        
    Parameters
    ----------
    generator_type : str, optional
        The type of generator to use. This can be either "crater" or "projectile". Default is "crater". 

    mean_velocity : float, optional
        The mean impact velocity to use for the impact simulation. Only one of either mean_velocity or impact_velocity_model can be provided.
    impact_velocity_model : str, optional
        The name of the mean impact velocity model to use for the impact simulation.  Valid options are "Mercury_MBA", "Venus_MBA", "Earth_MBA", "Moon_MBA", "Mars_MBA", and "MBA_MBA". 
        Only one of either mean_velocity or impact_velocity_model can be provided. Default is "Moon_MBA"
    N1_coef : float, optional
        The coefficient for the power law production function at 1 m diameter per 1 My. 
        Defaults to 7.9.e-3 (lunar craters) or 2.2e-8 (lunar impactors) based on fits to the NPF on the Moon.
    slope : float, optional
        The slope of the power law production function. 
        Defaults to -3.33 (lunar craters) or -2.26 (lunar impactors) based on fits to the NPF on the Moon.
    rng : numpy.random.Generator, optional
        A random number generator to use for sampling. If None, a default generator will be used.
    """
    def __init__(self, 
                rng: Generator | None = None,
                generator_type: str = "crater",
                mean_velocity: FloatLike | None = None,
                impact_velocity_model: str | None = None,
                N1_coef: FloatLike | None = None,
                slope: FloatLike | None = None,
                **kwargs: Any):

        super().__init__(rng=rng, 
                         mean_velocity=mean_velocity, 
                         impact_velocity_model=impact_velocity_model, 
                         **kwargs) 
        self.generator_type = generator_type
        
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
        if N1_coef is None:
            N1_coef = default_N1_coef[self.generator_type] 
        
        if not isinstance(N1_coef, FloatLike):
            raise ValueError("N1_coef must be a float")
        if N1_coef < 0.0:
            raise ValueError("N1_coef must be positive")
        self.N1_coef = N1_coef
       
        # Set the power law exponent for the production function along with defaults 
        if slope is None:
            slope = default_slope[self.generator_type]

        if not isinstance(slope, FloatLike):
            raise ValueError("slope must be a float")   
        elif slope > 0.0: # Slope must be negative, but convention in the field is mixed. So we flip the sign if it is positive.
            slope *= -1
        self.slope = slope 


      
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
    
    def chronology(self,
                   age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0, 
                   **kwargs: Any) -> Union[FloatLike, ArrayLike]:
        """
        Returns the age in My. Because the powerlaw model assumes constant impact rate, the returned age is the same as the input age.
        
        Parameters
        ----------
        age : FloatLike or ArrayLike, default=1.0
            Age in the past relative to the present day to compute cumulative SFD in units of My. 
        **kwargs: Any
            Any additional keywords that are passed to the function method.
            
        Returns
        -------
        FloatLike or numpy array of FloatLike
            The age in My for the given relative number density of craters. 
        """
        return age

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


