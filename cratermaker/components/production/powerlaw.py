from typing import Any

from numpy.random import Generator
from numpy.typing import ArrayLike

from cratermaker.components.production import Production
from cratermaker.constants import FloatLike


@Production.register("powerlaw")
class PowerLawProduction(Production):
    """
    An operations class for computing the production function for craters and projectiles.

    This impliments a very simple power law production function that can be used as either a crater or projectile production function. The production function is defined as the cumulative number of craters greater than a given diameter per unit m² surface area.

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
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        generator_type : str, optional
            The type of generator to use. This can be either "crater" or "projectile". Default is "crater".
        N1_coef : float, optional
            The coefficient for the power law production function (value of N>1 m diameter per 1 My per m²).
            Default is 7.9.e-3 (lunar craters) or 2.2e-8 (lunar projectiles) based on fits to the NPF on the Moon.
        slope : float, optional
            The slope of the power law production function.
            Default is -3.33 (lunar craters) or -2.26 (lunar projectiles) based on fits to the NPF on the Moon.
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        self.generator_type = generator_type

        # Default values that are approximately equal to the NPF for the Moon
        default_N1_coef = {"crater": 7.883e-3, "projectile": 7.989e-7}

        default_slope = {"crater": -3.328, "projectile": -2.634}

        # Set the power law parameters for the production function along with defaults
        if N1_coef is None:
            N1_coef = default_N1_coef[self.generator_type]

        if not isinstance(N1_coef, FloatLike):
            raise TypeError("N1_coef must be a float")
        if N1_coef < 0.0:
            raise ValueError("N1_coef must be positive")
        self.N1_coef = N1_coef

        # Set the power law exponent for the production function along with defaults
        if slope is None:
            slope = default_slope[self.generator_type]

        if not isinstance(slope, FloatLike):
            raise TypeError("slope must be a float")
        elif slope > 0.0:  # Slope must be negative, but convention in the field is mixed. So we flip the sign if it is positive.
            slope *= -1
        self.slope = slope

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += f"N1 Coefficient: {self.N1_coef:.2e}\n"
        str_repr += f"Slope: {self.slope:.3f}\n"
        return str_repr

    def _chronology(
        self,
        time_start: FloatLike | ArrayLike,
        time_end: FloatLike | ArrayLike,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike:
        """
        The chronology function, which returns the number of craters relative to the amount produced in the last 1 My.

        Parameters
        ----------
        time_start : FloatLike or ArrayLike
            Age in the past relative to the present day to compute cumulative SFD in units of My.
        time_end : FloatLike or ArrayLike
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD.
        **kwargs: Any
            |kwargs|

        Returns
        -------
        FloatLike or ArrayLike
            The number of craters relative to the amount produced in the last 1 My.

        Notes
        -----
        For a simple power law production function, the chronology function assumes a constant impact rate and so crater accumulation is linear with time. Therefore, the chronology function is simply the ratio of the cumulative SFD at the given age to the cumulative SFD at 1 My.
        """
        return time_start - time_end

    def _csfd(self, diameter: FloatLike | ArrayLike, **kwargs: Any) -> FloatLike | ArrayLike:
        """
        The cumulative size frequency distribution of craters at a given age relative to age = 1 My ago per m².

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            units of meter
        **kwargs: Any
            |kwargs|

        Returns
        -------
        FloatLike or numpy array
           Cumulative number density of per square meter greater than the input diameter.
        """
        return self.N1_coef * diameter**self.slope

    @property
    def N1_coef(self):
        """The N1 coefficient of the power law production function (N>1 m per My per m²)."""
        return self._N1_coef

    @N1_coef.setter
    def N1_coef(self, value):
        if not isinstance(value, FloatLike):
            raise TypeError("N1_coef must be a numeric value (float or int)")
        if value < 0:
            raise ValueError("N1_coef must be positive")
        self._N1_coef = value

    @property
    def slope(self):
        """The slope of the power law production function."""
        return self._slope

    @slope.setter
    def slope(self, value):
        if not isinstance(value, FloatLike):
            raise TypeError("slope must be a numeric value (float or int)")
        self._slope = value
