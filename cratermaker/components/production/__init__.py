from __future__ import annotations

from abc import abstractmethod
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike
from scipy.optimize import root_scalar

from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import parameter


class Production(ComponentBase):
    _registry: dict[str, Production] = {}

    def __init__(
        self,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        An abstract operations class that forms the base of classes that compute the production function for craters and projectiles.  The production function is defined as
        the cumulative number of craters greater than a given diameter per unit m^2 surface area per unit My time.

        Parameters
        ----------
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_valid_generator_types", ["crater", "projectile"])

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nGenerator type: {self.generator_type}"

    @classmethod
    def maker(
        cls,
        production: str | Production | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ) -> Production:
        """
        This is a helper function that can be used to validate and initialize the production model.

        Parameters
        ----------
        production : str | Production | None, optional
            The production model to use. This can be either a string or a Production instance.
            If None, the default production model is "neukum" and the version is based on the target (if provided), either Moon, Mars, or Projectile for all other bodies. Default is "Moon"
        target : Target | str | None, optional
            The target body for the impact. Can be a Target object or a string representing the target name.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        Returns
        -------
        Production
            An instance of the specified production model.
        Raises
        ------
        KeyError
            If the specified production model name is not found in the registry.
        TypeError
            If the specified production model is not a string or a subclass of Production.
        ValueError
            If there is an error initializing the production model.

        """

        version = kwargs.pop("version", None)
        if production is None or production == "neukum":
            if target is None and version in ["Moon", "Mars"]:
                target = version
            target = Target.maker(target, **kwargs)
            if target.name in ["Mercury", "Venus", "Earth", "Moon", "Mars"]:
                production = "neukum"
                if version is None:
                    if target.name in ["Moon", "Mars"]:
                        version = target.name
                    else:
                        version = "projectile"
            else:
                production = "powerlaw"
        return super().maker(
            component=production,
            version=version,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )

    def sample(
        self,
        age: FloatLike | None = None,
        age_end: FloatLike | None = None,
        diameter_number: PairOfFloats | None = None,
        diameter_number_end: PairOfFloats | None = None,
        diameter_range: PairOfFloats | None = None,
        area: FloatLike | None = None,
        return_age: bool = True,
        **kwargs: Any,
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
        **kwargs: Any
            Any additional keywords that are passed to the function method.

        Returns
        -------
        numpy array
            The sampled diameter values
        numpy array or None
            The sampled age values if return_age is True, otherwise None.

        """
        arguments = {
            "age": age,
            "age_end": age_end,
            "diameter_number": diameter_number,
            "diameter_number_end": diameter_number_end,
            "diameter_range": diameter_range,
            "area": area,
            "return_age": return_age,
            **kwargs,
        }
        arguments = self._validate_sample_args(**arguments)
        age = arguments["age"]
        age_end = arguments["age_end"]
        diameter_range = arguments["diameter_range"]
        area = arguments["area"]
        return_age = arguments["return_age"]

        # Build the cumulative distribution function from which we will sample
        input_diameters = np.logspace(
            np.log10(diameter_range[0]), np.log10(diameter_range[1])
        )
        cdf = self.function(
            diameter=input_diameters, age=age, age_end=age_end, **kwargs
        )
        expected_num = cdf[0] * area if area is not None else None
        diameters = mc.get_random_size(
            diameters=input_diameters,
            cdf=cdf,
            mu=expected_num,
            **vars(self.common_args),
        )
        if diameters.size == 0:
            return np.empty(0), np.empty(0)

        if return_age:
            age_subinterval = np.linspace(age_end, age, num=1000)
            N_vs_age = self.function(diameter=diameters, age=age_subinterval, **kwargs)

            # Normalize the weights for each diameter
            if N_vs_age.ndim > 1:
                N_vs_age -= np.min(N_vs_age, axis=1)[
                    :, None
                ]  # Subtract the min along age dimension
                N_vs_age = N_vs_age[
                    :, :-1
                ]  # Remove the last value to avoid the cumulative sum reaching 1
                weights_sum = np.sum(N_vs_age, axis=1)[
                    :, None
                ]  # Sum of weights for each diameter
                weights = N_vs_age / weights_sum  # Normalize weights for each diameter

                # Compute the CDF for each diameter
                cdf = np.cumsum(weights, axis=1)

            else:
                N_vs_age -= np.min(N_vs_age)  # Subtract the min along age dimension
                N_vs_age = N_vs_age[
                    :-1
                ]  # Remove the last value to avoid the cumulative sum reaching 1
                weights_sum = np.sum(N_vs_age)  # Sum of weights for each diameter
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
                    chosen_subinterval_indices[i] = (
                        np.searchsorted(cdf_row, random_values[i], side="right") - 1
                    )
            else:
                chosen_subinterval_indices = (
                    np.searchsorted(cdf, random_values[:], side="right") - 1
                )

            # Ensure indices are within valid range
            chosen_subinterval_indices = np.clip(
                chosen_subinterval_indices, 0, len(age_subinterval) - 2
            )
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

    def function_inverse(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike,
        cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike:
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

        diameter, cumulative_number_density = self._validate_csfd(
            diameter=diameter, cumulative_number_density=cumulative_number_density
        )

        def _root_func(t, D, N):
            kwargs.pop("check_valid_age", None)
            retval = (
                self.function(diameter=D, age=t, check_valid_age=False, **kwargs) - N
            )
            return retval

        xtol = 1e-8
        x0 = 4000.0
        if self.valid_age[0] is not None:
            xlo = self.valid_age[0]
        else:
            xlo = 0.0
        if self.valid_age[1] is not None:
            xhi = self.valid_age[1]
        else:
            xhi = 1e6
        retval = []
        darr = np.array(diameter)
        narr = np.array(cumulative_number_density)
        converged = []
        flag = []
        for i, d in np.ndenumerate(darr):
            sol = root_scalar(
                lambda x: _root_func(x, d, narr[i]),
                x0=x0,
                xtol=xtol,
                method="brentq",
                bracket=[xlo, xhi],
            )
            retval.append(sol.root)
            converged.append(sol.converged)
            flag.append(sol.flag)
        retval = np.array(retval)
        converged = np.array(converged)
        flag = np.array(flag)
        if np.all(converged):
            return retval.item() if np.isscalar(diameter) else retval
        else:
            raise ValueError(
                f"The root finding algorithm did not converge for all values of diameter and cumulative_number. Flag {flag}"
            )

    def get_projectile_properties(
        self, projectile_mean_velocity: FloatLike | None = None, **kwargs
    ) -> dict:
        """
        Generates basic projectile properties including density, velocity, and angle.

        Returns
        -------
        dict
            A dictionary containing the projectile properties.
        """

        return {
            "projectile_density": 0.0,
            "projectile_velocity": 0.0,
            "projectile_angle": 0.0,
        }

    @abstractmethod
    def chronology(
        self,
        age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        check_valid_age: bool = True,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike: ...

    @abstractmethod
    def function(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        check_valid_age: bool = True,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike: ...

    @abstractmethod
    def csfd(
        self, diameter: FloatLike | ArrayLike, **kwargs: Any
    ) -> FloatLike | ArrayLike: ...

    def _validate_csfd(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        cumulative_number_density: FloatLike
        | Sequence[FloatLike]
        | ArrayLike
        | None = None,
    ) -> tuple[FloatLike | ArrayLike, FloatLike | ArrayLike]:
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
            raise ValueError(
                "Either the 'diameter' or 'cumulative_number_density' must be provided"
            )

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
            if len(diameter_array) != len(cumulative_number_density_array):
                raise ValueError(
                    "The 'diameter' and 'cumulative_number_density' must have the same length when both are provided"
                )

        # Validate non-negative values
        if diameter_array is not None and np.any(diameter_array < 0):
            raise ValueError("All values in 'diameter' must be non-negative")
        if cumulative_number_density_array is not None and np.any(
            cumulative_number_density_array < 0
        ):
            raise ValueError(
                "All values in 'cumulative_number_density' must be non-negative"
            )

        if diameter is not None and not np.isscalar(diameter):
            diameter = diameter_array
        if cumulative_number_density is not None and not np.isscalar(
            cumulative_number_density
        ):
            cumulative_number_density = cumulative_number_density_array
        return diameter, cumulative_number_density

    def _validate_sample_args(self, **kwargs: dict) -> dict:
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
            projectiles, N, larger than diameter, D. If provided, the function convert this value to a corresponding age and use the
            production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N), which gives the total cumulative number of
            projectiles, N, larger than diameter, D.. If provided, the function will convert this value to a corresponding age_end
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
        _REF_DIAM = 1000.0  # The diameter value used if age arguments are provided
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
            raise ValueError(
                "Only one of the 'age' or 'diameter_number' arguments can be provided"
            )
        if age_end is not None and diameter_number_end is not None:
            raise ValueError(
                "Only one of the 'age_end' or 'diameter_number_end' arguments can be provided"
            )
        if age is not None and not np.isscalar(age):
            raise ValueError("The 'age' must be a scalar")
        if age_end is not None and not np.isscalar(age_end):
            raise ValueError("The 'age_end' must be a scalar")

        if diameter_range is None:
            raise ValueError("The 'diameter_range' must be provided")
        if len(diameter_range) != 2:
            raise ValueError("The 'diameter_range' must be a pair of values")
        if diameter_range[0] <= 0:
            raise ValueError(
                f"Diameter range invalid: {diameter_range}. The minimum diameter must be greater than 0"
            )
        if diameter_range[1] < diameter_range[0]:
            raise ValueError(
                f"Diameter range invalid: {diameter_range}. The maximum diameter must be greater than or equal to the minimum diameter"
            )

        if area is None:
            raise ValueError("The 'area' must be provided")
        else:
            if not np.isscalar(area):
                raise ValueError("The 'area' must be a scalar")
            if area < 0.0:
                raise ValueError("The 'area' must be greater than 0")

        if diameter_number is not None:
            if len(diameter_number) != 2:
                raise ValueError(
                    "The 'diameter_number' must be a pair of values in the form (D, N)"
                )
            diameter_number = self._validate_csfd(*diameter_number)
            diameter_number_density = (diameter_number[0], diameter_number[1] / area)
            age = self.function_inverse(*diameter_number_density)
        else:
            diameter_number_density = (
                _REF_DIAM,
                self.function(diameter=_REF_DIAM, age=age),
            )
            diameter_number = (
                diameter_number_density[0],
                diameter_number_density[1] * area,
            )

        if age_end is None and diameter_number_end is None:
            diameter_number_end = (_REF_DIAM, 0.0)
            diameter_number_density_end = diameter_number_end
            age_end = 0.0

        if diameter_number_end is not None:
            if len(diameter_number_end) != 2:
                raise ValueError(
                    "The 'diameter_number_end' must be a pair of values in the form (D,N)"
                )
            diameter_number_end = self._validate_csfd(*diameter_number_end)
            diameter_number_density_end = (
                diameter_number_end[0],
                diameter_number_end[1] / area,
            )
            # Check to be sure that the diameter in the diameter_number_end is the same as diameter_number.
            # If not, we need to adjust the end diameter value it so that they match
            if diameter_number is not None:
                if diameter_number_density_end[0] != diameter_number_density[0]:
                    diameter_number_density_end = (
                        diameter_number_density[0],
                        self.function(diameter=diameter_number_density[0], age=age_end),
                    )
                    diameter_number_end = (
                        diameter_number[0],
                        diameter_number_density_end[1] * area,
                    )

        if age_end is None:
            diameter_number_end = self._validate_csfd(*diameter_number_end)
            diameter_number_density_end = (
                diameter_number_end[0],
                diameter_number_end[1] / area,
            )
            age_end = self.function_inverse(*diameter_number_density_end)

        age, age_end = self._validate_age(age, age_end)

        if diameter_number_end is None:
            diameter_number_end = (
                diameter_number[0],
                self.function(diameter=diameter_number[0], age=age_end) * area,
            )

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

    def _validate_age(
        self,
        age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        check_valid_age: bool = True,
    ) -> FloatLike | ArrayLike:
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
        tuple of float
            The start and end ages in units of My.

        Raises
        ------
        ValueError
            If the the start age is greater than the end age or the age variable is not a scalar or a sequence of 2 values.
        """

        if np.isscalar(age):
            if not isinstance(age, FloatLike):
                raise TypeError("age must be a numeric value (float or int)")
            age = float(age)

            if age_end is None:
                age_end = float(0.0)
            else:
                if not np.isscalar(age_end):
                    raise ValueError("If age is a scalar, age_end must be a scalar")
                elif not isinstance(age_end, FloatLike):
                    raise TypeError("age_end must be a numeric value (float or int)")
                age_end = float(age_end)
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
                raise ValueError(
                    "If age is a sequence, age_end must be a sequence of the same size"
                )
            if np.any(age < age_end):
                raise ValueError("age must be greater than the age_end")
        else:
            raise ValueError("age must be a scalar or a sequence")

        if check_valid_age:
            if self.valid_age[0] is not None:
                if np.any(age < self.valid_age[0]):
                    raise ValueError(
                        f"age must be greater than the minimum valid age {self.valid_age[0]}"
                    )
            if self.valid_age[1] is not None:
                if np.any(age > self.valid_age[1]):
                    raise ValueError(
                        f"age must be less than the maximum valid age {self.valid_age[1]}"
                    )

        return age, age_end

    @parameter
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
            self._generator_type = self._valid_generator_types[0]
            return
        if not isinstance(value, str):
            raise ValueError("generator_type must be a string")
        if value.lower() not in self._valid_generator_types:
            raise ValueError(
                f"Invalid generator_type {value}. Must be one of {self._valid_generator_types}"
            )
        self._generator_type = value.lower()
        return

    @property
    def valid_age(self) -> tuple[float, float]:
        """
        The range of ages over which the production function is valid. The range is given in My.

        Returns
        -------
        tuple
            The lower and upper bounds of the valid time range in My, or None if not applicable.
        """
        return (None, None)


import_components(__name__, __path__, ignore_private=True)
