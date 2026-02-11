from __future__ import annotations

from abc import abstractmethod
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import root_scalar

from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.utils.general_utils import parameter


class Production(ComponentBase):
    """
    The base class for computing the production function for craters and projectiles.
    """

    _registry: dict[str, Production] = {}

    def __init__(
        self,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
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
        Initialize a Production model with the given name or instance.

        Parameters
        ----------
        production : str | Production | None, optional
            The production model to use. This can be either a string or a Production instance.
            If None, the default production model is "neukum" and the version is based on the target (if provided), either Moon, Mars, or Projectile for all other bodies. Default is "Moon"
        target : Target | str | None, optional
            The target body for the impact. Can be a Target object or a string representing the target name.
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|

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
        time_start: FloatLike | None = None,
        time_end: FloatLike | None = None,
        diameter_number: PairOfFloats | None = None,
        diameter_number_end: PairOfFloats | None = None,
        diameter_range: PairOfFloats | None = None,
        area: FloatLike | None = None,
        return_age: bool = True,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        Sample diameters and ages from the production function.

        This function can either sample from a given age range or
        from a given cumulative number/diameter pair (but not both).

        Parameters
        ----------
        time_start : FloatLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        time_end, FloatLike, optional
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
        if "age" in kwargs:
            time_start = kwargs.pop("age")
        arguments = {
            "time_start": time_start,
            "time_end": time_end,
            "diameter_number": diameter_number,
            "diameter_number_end": diameter_number_end,
            "diameter_range": diameter_range,
            "area": area,
            "return_age": return_age,
            **kwargs,
        }
        arguments = self._validate_sample_args(**arguments)
        time_start = arguments["time_start"]
        time_end = arguments["time_end"]
        diameter_range = arguments["diameter_range"]
        area = arguments["area"]
        return_age = arguments["return_age"]
        validate_inputs = kwargs.pop("validate_inputs", False)

        # Build the cumulative distribution function from which we will sample
        input_diameters = np.logspace(np.log10(diameter_range[0]), np.log10(diameter_range[1]))
        cdf = self.function(
            diameter=input_diameters,
            time_start=time_start,
            time_end=time_end,
            validate_inputs=validate_inputs,
            **kwargs,
        )
        if cdf.ndim > 1:
            cdf = cdf[:, 0]
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
            time_subinterval = np.linspace(time_end, time_start, num=1000)
            N_vs_age = self.function(
                diameter=diameters,
                time=time_subinterval,
                validate_inputs=validate_inputs,
                **kwargs,
            )

            # Normalize the weights for each diameter
            if N_vs_age.ndim > 1:
                N_vs_age -= np.min(N_vs_age, axis=1)[:, None]  # Subtract the min along age dimension
                N_vs_age = N_vs_age[:, :-1]  # Remove the last value to avoid the cumulative sum reaching 1
                weights_sum = np.sum(N_vs_age, axis=1)[:, None]  # Sum of weights for each diameter
                weights = N_vs_age / weights_sum  # Normalize weights for each diameter

                # Compute the CDF for each diameter
                cdf = np.cumsum(weights, axis=1)

            else:
                if len(N_vs_age) > 0:
                    N_vs_age -= np.min(N_vs_age)  # Subtract the min along age dimension
                    N_vs_age = N_vs_age[:-1]  # Remove the last value to avoid the cumulative sum reaching 1
                    weights_sum = np.sum(N_vs_age)  # Sum of weights for each diameter
                    weights = N_vs_age / weights_sum  # Normalize weights for each diameter
                    # Compute the CDF for each diameter
                    cdf = np.cumsum(weights)
                else:
                    cdf = np.empty(0)

            # Generate uniform random numbers for each diameter
            random_values = self.rng.uniform(0, 1, size=diameters.size)
            if N_vs_age.ndim > 1:
                chosen_subinterval_indices = np.zeros(diameters.size, dtype=int)
                # For each diameter, find the corresponding index in the CDF
                for i, cdf_row in enumerate(cdf):
                    # Find the corresponding index in the CDF for each diameter
                    chosen_subinterval_indices[i] = np.searchsorted(cdf_row, random_values[i], side="right") - 1
            else:
                chosen_subinterval_indices = np.searchsorted(cdf, random_values[:], side="right") - 1

            # Ensure indices are within valid range
            chosen_subinterval_indices = np.clip(chosen_subinterval_indices, 0, len(time_subinterval) - 2)
            # Sample a random age within the selected subinterval for each diameter
            time_lo = time_subinterval[chosen_subinterval_indices]
            time_hi = time_subinterval[chosen_subinterval_indices + 1]
            timevals = self.rng.uniform(low=time_lo, high=time_hi)

            # Sort the ages and diameters so that they are in order of decreasing age
            if timevals.size > 1:
                sort_indices = np.argsort(timevals)[::-1]
                diameters = diameters[sort_indices]
                timevals = timevals[sort_indices]
        else:
            timevals = np.empty(0)

        return diameters, timevals

    @abstractmethod
    def function(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        time_start: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        time_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        **kwargs: Any,
    ) -> NDArray[np.float64]: ...

    def age_from_D_N(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike,
        cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike,
        validate_inputs: bool = True,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike:
        """
        Return the age in My for a given number density of craters and diameter.

        Parameters
        ----------
        diameter : float-like or array-like
            diameter of the crater in m
        cumulative_number_density : float-like or array-like
            number density of craters per m^2 surface area greater than the input diameter
        validate_inputs : bool, optional
            If True, the function will validate the inputs. The default is True.
        **kwargs: Any
            Any additional keywords that are passed to the function method.

        Returns
        -------
        float_like or numpy array
            The age in My for the given relative number density of craters.
        """
        if validate_inputs:
            diameter, cumulative_number_density = self._validate_csfd(
                diameter=diameter, cumulative_number_density=cumulative_number_density
            )

        def _root_func(t, D, N):
            kwargs.pop("validate_inputs", None)
            retval = self.function(diameter=D, time_start=t, validate_inputs=False, **kwargs) - N
            return retval

        xtol = 1e-8
        x0 = 4000.0
        if self.valid_time[0] is not None:
            xlo = self.valid_time[0]
        else:
            xlo = 0.0
        if self.valid_time[1] is not None:
            xhi = self.valid_time[1]
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

    def D_from_N_age(
        self,
        cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike,
        time_start: FloatLike | Sequence[FloatLike] | ArrayLike,
        time_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        validate_inputs: bool = True,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike:
        """
        Return the diameter for a given number density of craters and age.

        Parameters
        ----------
        cumulative_number_density : float-like or array-like
            number density of craters per m^2 surface area greater than the input diameter
        time_start : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        time_end, FloatLike or ArrayLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        validate_inputs : bool, optional
            If True, the function will validate the inputs. The default is True.
        **kwargs: Any
            Any additional keywords that are passed to the function method.

        Returns
        -------
        float_like or numpy array
            The diameter for the given number density and age.
        """
        if "age" in kwargs:
            time_start = kwargs.pop("age")
        if validate_inputs:
            time_start, time_end = self._validate_age(time_start, time_end)

        def _root_func(D, N, time_start, time_end):
            kwargs.pop("validate_inputs", None)
            retval = (
                self.function(
                    diameter=D,
                    time_start=time_start,
                    time_end=time_end,
                    validate_inputs=False,
                    **kwargs,
                )
                - N
            )
            return retval

        xtol = 1e-8
        x0 = 1.0  # Initial guess for the diameter
        xlo = 1e-6
        xhi = 1e6
        retval = []
        narr = np.array(cumulative_number_density)
        converged = []
        flag = []
        for _, n in np.ndenumerate(narr):
            sol = root_scalar(
                lambda x: _root_func(x, n, time_start, time_end),
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
            return retval.item() if np.isscalar(cumulative_number_density) else retval
        else:
            raise ValueError(
                f"The root finding algorithm did not converge for all values of cumulative_number and age. Flag {flag}"
            )

    @abstractmethod
    def chronology(
        self,
        time_start: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        time_end: FloatLike | Sequence[FloatLike] | ArrayLike = 0.0,
        validate_inputs: bool = True,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike: ...

    @abstractmethod
    def csfd(self, diameter: FloatLike | ArrayLike, **kwargs: Any) -> FloatLike | ArrayLike: ...

    def _validate_csfd(
        self,
        diameter: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
        cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
    ) -> tuple[FloatLike | ArrayLike, FloatLike | ArrayLike]:
        """
        Validates the diameter and cumulative_number arguments.

        Both arguments can be either scalar or array-like, but they must be both scalars or both arrays of the same length.  Values must be non-negative.

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
            if len(diameter_array) != len(cumulative_number_density_array):
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

    def _validate_sample_args(self, **kwargs: dict) -> dict:
        """
        Validate all the input arguments to the sample method.

        This function will raise a ValueError if any of the arguments are invalid.  It will also convert age arguments to diameter_number and vice versa.

        Parameters
        ----------
        time_start : FloatLike, optional
            The starting time in units of My relative to the present, which is used to compute the
        time_end, FloatLike, optional
            The ending time in units of My relative to the present, which is used to compute the cumulative SFD.
            The default is 0 (present day).
        diameter_number : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N), which gives the total cumulative number of
            projectiles, N, larger than diameter, D. If provided, the function convert this value to a corresponding age and use the
            production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N), which gives the total cumulative number of
            projectiles, N, larger than diameter, D.. If provided, the function will convert this value to a corresponding time_end
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
            The start time in units of My relative to the present.
        FloatLike
            The end time in units of My relative to the present.
        PairOfFloats
            A pair of diameter and cumulative number values, in the form of a (D, N) representing the cumulative number and diameter values for a surface at the age given by the start time
        PairOfFloats
            A pair of diameter and cumulative number values, in the form of a (D, N) representing the cumulative number and diameter values for a surface at the age given by the end time.
        PairOfFloats
            The minimum and maximum diameter values to sample from in meters.
        FloatLike
            The area in m^2 over which the production function is evaluated to generate the expected number, which is the production
            function over the input time/cumulative number range at the minimum diameter.

        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - Neither the age nore the diameter_number argument is provided.
            - Both the age and diameter_number arguments are provided.
            - Both the time_end and diameter_number_end arguments are provided.
            - The time_start or age argument is provided but is not a scalar.
            - The time_end argument is provided but is not a scalar.
            - The diameter_number argument is not a pair of values, or any of them are less than 0
            - The diameter_number_end argument is not a pair of values, or any of them are less than 0
            - The diameter_range is not provided.
            - The diameter_range is not a pair of values.
            - The minimum diameter is less than or equal to 0.
            - The maximum diameter is less than or equal the minimum.
            - The area argument is not a scalar or is less than 0.
        """
        _REF_DIAM = 1000.0  # The diameter value used if age arguments are provided
        age = kwargs.pop("age", None)
        if age is not None:
            kwargs["time_start"] = age
        time_start = kwargs.get("time_start")
        time_end = kwargs.get("time_end")
        diameter_number = kwargs.get("diameter_number")
        diameter_number_end = kwargs.get("diameter_number_end")
        diameter_range = kwargs.get("diameter_range")
        area = kwargs.get("area")
        return_age = kwargs.get("return_age", True)

        if time_start is None and diameter_number is None:
            raise ValueError("Either the 'time_start' or 'diameter_number' must be provided")
        elif time_start is not None and diameter_number is not None:
            raise ValueError("Only one of the 'time_start' or 'diameter_number' arguments can be provided")
        if time_end is not None and diameter_number_end is not None:
            raise ValueError("Only one of the 'time_end' or 'diameter_number_end' arguments can be provided")
        if time_start is not None and not np.isscalar(time_start):
            raise ValueError("The 'time_start' must be a scalar")
        if time_end is not None and not np.isscalar(time_end):
            raise ValueError("The 'time_end' must be a scalar")

        if diameter_range is None:
            raise ValueError("The 'diameter_range' must be provided")
        if len(diameter_range) != 2:
            raise ValueError("The 'diameter_range' must be a pair of values")
        if diameter_range[0] <= 0:
            raise ValueError(f"Diameter range invalid: {diameter_range}. The minimum diameter must be greater than 0")
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
                raise ValueError("The 'diameter_number' must be a pair of values in the form (D, N)")
            diameter_number = self._validate_csfd(*diameter_number)
            diameter_number_density = (diameter_number[0], diameter_number[1] / area)
            time_start = self.age_from_D_N(*diameter_number_density)
        else:
            diameter_number_density = (
                _REF_DIAM,
                self.function(diameter=_REF_DIAM, age=time_start),
            )
            diameter_number = (
                diameter_number_density[0],
                diameter_number_density[1] * area,
            )

        if time_end is None and diameter_number_end is None:
            diameter_number_end = (_REF_DIAM, 0.0)
            diameter_number_density_end = diameter_number_end
            time_end = 0.0

        if diameter_number_end is not None:
            if len(diameter_number_end) != 2:
                raise ValueError("The 'diameter_number_end' must be a pair of values in the form (D,N)")
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
                        self.function(diameter=diameter_number_density[0], time_start=time_end),
                    )
                    diameter_number_end = (
                        diameter_number[0],
                        diameter_number_density_end[1] * area,
                    )

        if time_end is None:
            diameter_number_end = self._validate_csfd(*diameter_number_end)
            diameter_number_density_end = (
                diameter_number_end[0],
                diameter_number_end[1] / area,
            )
            time_end = self.age_from_D_N(*diameter_number_density_end)

        time_start, time_end = self._validate_age(time_start, time_end)

        if diameter_number_end is None:
            diameter_number_end = (
                diameter_number[0],
                self.function(diameter=diameter_number[0], time_start=time_end) * area,
            )

        if not isinstance(return_age, bool):
            raise ValueError("The 'return_age' argument must be a boolean")

        kwargs["time_start"] = time_start
        kwargs["time_end"] = time_end
        kwargs["diameter_number"] = diameter_number
        kwargs["diameter_number_end"] = diameter_number_end
        kwargs["diameter_range"] = diameter_range
        kwargs["area"] = area
        kwargs["return_age"] = return_age

        return kwargs

    def _validate_age(
        self,
        time_start: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
        time_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
    ) -> FloatLike | ArrayLike:
        """
        Processes the time_start argument and time_end arguments. Checks that they are valid and returns a tuple of age and time_end.

        Parameters
        ----------
        time_start : FloatLike or ArrayLike, default=1.0
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        time_end, FloatLike or ArrayLike, optional
            The reference used when computing age in My. If none is passed, it will be set either 0 or an array of zeros, depending on
            the size of time_start.

        Returns
        -------
        tuple of float
            The start and end ages in units of My.

        Raises
        ------
        ValueError
            If the the start age is greater than the end age or the age variable is not a scalar or a sequence of 2 values.
        """
        if np.isscalar(time_start):
            if not isinstance(time_start, FloatLike):
                raise TypeError("time_start must be a numeric value (float or int)")
            if time_end is None:
                time_end = 0.0
            elif not np.isscalar(time_end) or not isinstance(time_end, FloatLike):
                raise TypeError("time_end must be a numeric value (float or int)")
            if time_start < time_end:
                raise ValueError("time_start must be greater than the time_end")
            if self.valid_time[0] is not None and time_start < self.valid_time[0]:
                raise ValueError(
                    f"time_start must be greater than the minimum valid time {self.valid_time[0]} (it is units of My before present)"
                )
            if self.valid_time[1] is not None and time_start > self.valid_time[1]:
                raise ValueError(
                    f"time_start must be less than the maximum valid time {self.valid_time[1]} (it is in units of My before present)"
                )
        else:
            if isinstance(time_start, (list, tuple)):
                time_start = np.array(time_start, dtype=np.float64)
            elif not isinstance(time_start, np.ndarray):
                raise TypeError("time_start must be a numeric value (float or int) or an array")

            if time_end is None:
                time_end = np.zeros_like(time_start, dtype=np.float64)
            elif np.isscalar(time_end):
                time_end = np.full_like(time_start, time_end, dtype=np.float64)
            elif isinstance(time_end, (list, tuple)):
                time_end = np.array(time_end, dtype=np.float64)
            elif not isinstance(time_end, np.ndarray):
                raise TypeError("time_end must be a numeric value (float or int) or an array")

            if time_start.size != time_end.size:
                raise ValueError("The 'time_start' and 'time_end' arguments must be the same size if both are provided")

            if np.any(time_start < time_end):
                raise ValueError("time_start must be greater than the time_end")

            if self.valid_time[0] is not None:
                if np.any(time_start < self.valid_time[0]) or np.any(time_end < self.valid_time[0]):
                    raise ValueError(
                        f"time_start must be greater than the minimum valid age {self.valid_time[0]} (it is in units of My before present)"
                    )
            if self.valid_time[1] is not None:
                if np.any(time_start > self.valid_time[1]) or np.any(time_end > self.valid_time[1]):
                    raise ValueError(
                        f"time_start must be less than the maximum valid age {self.valid_time[1]} (it is in units of My before present)"
                    )

        return time_start, time_end

    @parameter
    def generator_type(self):
        """
        The type of generator to use.

        This can be either "crater" or "projectile".  This determines the nature of the production function, differentiating between 'crater' and 'projectile' types, which affect the calculations and interpretations of the production function.
        """
        return self._generator_type

    @generator_type.setter
    def generator_type(self, value):
        if not value:
            self._generator_type = self._valid_generator_types[0]
            return
        if not isinstance(value, str):
            raise TypeError("generator_type must be a string")
        if value.lower() not in self._valid_generator_types:
            raise ValueError(f"Invalid generator_type {value}. Must be one of {self._valid_generator_types}")
        self._generator_type = value.lower()
        return

    @property
    def valid_time(self) -> tuple[float, float]:
        """
        The range of ages over which the production function is valid. The range is given in My.

        Returns
        -------
        tuple
            The lower and upper bounds of the valid time range in My, or None if not applicable.
        """
        return (None, None)


import_components(__name__, __path__)
