from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike

from cratermaker.components.production import Production
from cratermaker.constants import FloatLike


@Production.register("powerlaw")
class PowerLawProduction(Production):
    """
    An operations class for computing the production function for craters and projectiles. This impliments a very simple power law
    production function that can be used as either a crater or projectile production function. The production function is defined as
    the cumulative number of craters greater than a given diameter per unit m^2 surface area.

    Parameters
    ----------
    generator_type : str, optional
        The type of generator to use. This can be either "crater" or "projectile". Default is "crater".
    N1_coef : float, optional
        The coefficient for the power law production function at 1 m diameter per 1 My.
        Defaults to 7.9.e-3 (lunar craters) or 2.2e-8 (lunar projectiles) based on fits to the NPF on the Moon.
    slope : float, optional
        The slope of the power law production function.
        Defaults to -3.33 (lunar craters) or -2.26 (lunar projectiles) based on fits to the NPF on the Moon.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        generator_type: str = "crater",
        N1_coef: FloatLike | None = None,
        slope: FloatLike | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        self.generator_type = generator_type

        # Default values that are approximately equal to the NPF for the Moon
        default_N1_coef = {"crater": 7.883e-3, "projectile": 7.989e-7}

        default_slope = {"crater": -3.328, "projectile": -2.634}

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
        elif (
            slope > 0.0
        ):  # Slope must be negative, but convention in the field is mixed. So we flip the sign if it is positive.
            slope *= -1
        self.slope = slope

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nN1 Coefficient: {self.N1_coef:.2e}\nSlope: {self.slope:.3f}"

    def function(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike:
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

        n_array = np.asarray(self.csfd(diameter))
        age_difference = np.asarray(age - age_end)

        if n_array.ndim > 0 and age_difference.ndim > 0:
            return n_array[:, None] * age_difference
        else:
            return n_array * age_difference

    def chronology(
        self, age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0, **kwargs: Any
    ) -> FloatLike | ArrayLike:
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

    def csfd(
        self, diameter: FloatLike | ArrayLike, **kwargs: Any
    ) -> FloatLike | ArrayLike:
        """
        Return the cumulative size frequency distribution of craters at a given age relative to age = 1 My ago per m^2.

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            units of meter

        Returns
        -------
        FloatLike or numpy array
           Cumulative number density of per square meter greater than the input diameter.
        """
        return self.N1_coef * diameter**self.slope

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
