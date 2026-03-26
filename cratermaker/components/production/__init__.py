from __future__ import annotations

from abc import abstractmethod
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import root_scalar

from cratermaker.components.crater import Crater
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.utils.general_utils import format_large_units, parameter
from cratermaker.utils.montecarlo_utils import bounded_norm


class Production(ComponentBase):
    """
    The base class for computing the production function for craters and projectiles.
    """

    _registry: dict[str, Production] = {}

    def __init__(
        self,
        quasimc_file: str | Path | None = None,
        quasimc_craters: list[Crater] | None = None,
        diameter_range: PairOfFloats | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        quasimc_file : str | Path, optional
            Path to a file (CSV or NetCDF) containing the parameters used for craters emplaced using the quasi-Monte Carlo method. This file should contain at a minimum the diameter (or radius) of each crater, its location (lon, lat), and one of either production_time or production_D and production_N (D in km, N in units given by the ND_conversion_factor attribute).
        quasimc_craters : list[Crater], optional
            A list of Crater objects that are emplaced using the quasi-Monte Carlo method. This is an alternative to providing a quasimc_file. Only one of either quasimc_file or quasimc_craters should be provided.
        diameter_range : PairOfFloats, optional
            The minimum and maximum crater diameter to sample from in meters. If not provided, the default is (0, inf), unless quasimc_file or quasimc_craters is provided, in which case the upper range will be set based on the smallest diameter in the quasi-Monte Carlo data.
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
        object.__setattr__(self, "_quasimc_file", None)
        object.__setattr__(self, "_quasimc_craters", None)
        object.__setattr__(self, "_ND_conversion_factor", None)
        object.__setattr__(self, "_ND_unit", None)
        object.__setattr__(self, "_diameter_range", None)
        self.diameter_range = diameter_range

        if quasimc_file is not None and quasimc_craters is not None:
            raise ValueError("Cannot provide both quasimc_file and quasimc_craters")

        if quasimc_file is not None:
            if Path(quasimc_file).exists():
                self.read_quasimc_file(filename=quasimc_file, **kwargs)
            else:
                raise FileNotFoundError(f"quasimc_file {quasimc_file} not found")
        elif quasimc_craters is not None:
            self.process_quasimc_craters(craters=quasimc_craters, **kwargs)

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nGenerator type: {self.generator_type}"

    @classmethod
    def maker(
        cls,
        production: str | Production | None = None,
        target: Target | str | None = None,
        quasimc_file: str | Path | None = None,
        quasimc_craters: list[Crater] | None = None,
        diameter_range: PairOfFloats | None = None,
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
        quasimc_file : str | Path, optional
            Path to a file (CSV or NetCDF) containing the parameters used for craters emplaced using the quasi-Monte Carlo method. This file should contain at a minimum the diameter (or radius) of each crater, its location (lon, lat), and one of either production_time or production_D and production_N (D in km, N in units given by the ND_conversion_factor attribute).
        quasimc_craters : list[Crater], optional
            A list of Crater objects that are emplaced using the quasi-Monte Carlo method. This is an alternative to providing a quasimc_file. Only one of either quasimc_file or quasimc_craters should be provided.
        diameter_range : PairOfFloats, optional
            The minimum and maximum crater diameter to sample from in meters. If not provided, the default is (0, inf), unless quasimc_file or quasimc_craters is provided, in which case the upper range will be set based on the smallest diameter in the quasi-Monte Carlo data.
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
        set_version = False
        version = None
        if production is None:
            set_version = True
        elif isinstance(production, str):
            production = production.lower()
            if production == "neukum":
                set_version = True
            if set_version:
                version = kwargs.pop("version", None)
                target = Target.maker(target, **kwargs)
                if target.name in ["Mercury", "Venus", "Earth", "Moon", "Mars"]:
                    production = "neukum"
                    if version is None or version != "projectile":
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
            quasimc_file=quasimc_file,
            quasimc_craters=quasimc_craters,
            diameter_range=diameter_range,
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
        diameter_range : PairOfFloats, optional
            The minimum and maximum crater diameter to sample from in meters. If not provided, the :py:attr:`~cratermaker.components.production.Production.diameter_range` value will be used.
        area : FloatLike, optional
            The area in m² over which the production function is evaluated to generate the expected number, which is the production
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
                time_start=time_subinterval,
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
            number density of craters per m² surface area greater than the input diameter
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

        def _root_func(t, d, n):
            kwargs.pop("validate_inputs", None)
            retval = self.function(diameter=d, time_start=t, validate_inputs=False, **kwargs) - n
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
                lambda t, d=d, n_i=narr[i]: _root_func(t, d, n_i),
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
            number density of craters per m² surface area greater than the input diameter
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

        def _root_func(d, n, time_start, time_end):
            kwargs.pop("validate_inputs", None)
            retval = (
                self.function(
                    diameter=d,
                    time_start=time_start,
                    time_end=time_end,
                    validate_inputs=False,
                    **kwargs,
                )
                - n
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
                lambda d, n=n: _root_func(d, n, time_start, time_end),
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
            The minimum and maximum crater diameter to sample from in meters. If not provided, the :py:attr:`~cratermaker.components.production.Production.diameter_range`
        area : FloatLike, optional
            The area in m² over which the production function is evaluated to generate the expected number, which is the production
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
            The area in m² over which the production function is evaluated to generate the expected number, which is the production
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
        diameter_range = kwargs.pop("diameter_range", self.diameter_range)
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
                self.function(diameter=_REF_DIAM, time_start=time_start),
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
            raise TypeError("The 'return_age' argument must be a boolean")

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
            if isinstance(time_start, (list, tuple, ArrayLike)):
                time_start = np.array(time_start, dtype=np.float64)
            elif not isinstance(time_start, np.ndarray):
                raise TypeError("time_start must be a numeric value (float or int) or an array")

            if time_end is None:
                time_end = np.zeros_like(time_start, dtype=np.float64)
            elif np.isscalar(time_end):
                time_end = np.full_like(time_start, time_end, dtype=np.float64)
            elif isinstance(time_end, (list, tuple, ArrayLike)):
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

    @parameter
    def quasimc_file(self) -> Path:
        """
        File containing the quasi-Monte Carlo craters.
        """
        return self._quasimc_file

    @quasimc_file.setter
    def quasimc_file(self, value) -> None:
        if value is not None:
            if isinstance(value, str):
                value = Path(value)
            if isinstance(value, Path):
                if not value.exists():
                    raise FileNotFoundError(f"quasimc_file {str(value)} not found")
                self._quasimc_file = value
                return

    @property
    def quasimc_craters(self) -> list[Crater]:
        """
        The list of processed quasi-Monte Carlo craters.
        """
        return self._quasimc_craters

    @quasimc_craters.setter
    def quasimc_craters(self, value: list[Crater]):
        if not isinstance(value, list) or not all(isinstance(c, Crater) for c in value):
            raise TypeError("quasimc_craters must be a list of Crater objects")
        self._quasimc_craters = value

    def process_quasimc_craters(
        self, craters: list[Crater], set_max_diameter_from_quasimc: bool = True, **kwargs: Any
    ) -> list[Crater]:
        """
        Process the quasi-Monte Carlo craters by computing production_time if present, production_ND if present, and will drop craters that have neither.

        Parameters
        ----------
        craters : list[Crater]
            The list of Crater objects to process.
        set_max_diameter_from_quasimc : bool, optional
            If True, the largest_crater attribute of the Simulation will be set to the smallest crater in the quasi-Monte Carlo file. Default is True.

        Returns
        -------
        list[Crater]
            The list of processed Crater objects with production_time and production_ND computed as needed, and any craters that have neither dropped.
        """
        _DSTD = 1e3  # Standard diameter to use for conversion between time and N(D) values when not provided

        def _draw_quasimc_time(
            crater: Crater,
            sequence_extents: dict[int, float] | None = None,
        ) -> float | None:
            """Draw one time sample for a crater, or None if no production metadata."""
            if sequence_extents is not None:
                s = crater.production_sequence
                seq_nlo, seq_nhi = sequence_extents[s]
                seq_nlo = max(seq_nlo, 0.0)
                seq_nmean = (seq_nlo + seq_nhi) / 2
                seq_nsig = (seq_nhi - seq_nlo) / 2
                if seq_nsig / seq_nmean < 1e-12:
                    seq_nsig = 0.0
                seq_tlo = self.age_from_D_N(diameter=_DSTD, cumulative_number_density=seq_nlo)
                seq_thi = self.age_from_D_N(diameter=_DSTD, cumulative_number_density=seq_nhi)
                seq_tmean = (seq_tlo + seq_thi) / 2
            else:
                s = None

            has_sequence = s is not None and sequence_extents is not None

            if crater.production_ND is None and crater.production_time is None:
                if not has_sequence:
                    raise ValueError(
                        "Craters without production_time or production_ND need a production_sequence and sequence_extents"
                    )

                if seq_nsig > 0.0:
                    nval = bounded_norm(loc=seq_nmean, scale=seq_nsig, lower_bound=seq_nlo, upper_bound=seq_nhi, rng=self.rng)
                else:
                    nval = seq_nmean
                nval = max(nval, 0.0)
                return float(self.age_from_D_N(diameter=_DSTD, cumulative_number_density=nval))
            elif crater.production_ND is not None and crater.production_ND != [None, None, None]:
                # N(D) values take precedence over time values
                prod_diam, nmean, nsig = crater.production_ND
                prod_diam *= 1e3
                nmean /= self.ND_conversion_factor
                nsig /= self.ND_conversion_factor
                # Convert from N(prod_diam) to N(1) for alignment
                fac = self.csfd(_DSTD) / self.csfd(prod_diam)
                nmean *= fac
                nsig *= fac
                if has_sequence:
                    if seq_nsig > 0.0:
                        nval = bounded_norm(loc=nmean, scale=nsig, lower_bound=seq_nlo, upper_bound=seq_nhi, rng=self.rng)
                    else:
                        nval = seq_nmean
                else:
                    nval = self.rng.normal(loc=nmean, scale=nsig) if nsig > 0 else nmean
                nval = max(0.0, nval)
                return float(self.age_from_D_N(diameter=_DSTD, cumulative_number_density=nval))
            elif crater.production_time is not None and crater.production_time != [None, None]:
                tmean, tsig = crater.production_time
                if has_sequence:
                    if seq_nsig > 0.0:
                        if tsig > 0.0:
                            t = bounded_norm(loc=tmean, scale=tsig, lower_bound=seq_tlo, upper_bound=seq_thi, rng=self.rng)
                        else:
                            t = tmean
                    else:
                        t = seq_tmean
                else:
                    t = self.rng.normal(loc=tmean, scale=tsig) if tsig > 0 else tmean
                max(t, 0.0)
                return float(t)
            else:
                return None

        def _sequence_is_consistent(
            times_by_idx: dict[int, float],
            sequence_groups: dict[int, list[int]],
        ) -> bool:
            """
            For increasing tag value, times must strictly decrease: sequence(a) < sequence(b)  ==>  time(a) > time(b).
            """
            if any(t is None for t in times_by_idx.values()):
                return False
            ordered_sequence = sorted(sequence_groups.keys())[:-1]
            for t_lo, t_hi in zip(ordered_sequence[:-1], ordered_sequence[1:], strict=True):
                lo_min = min(
                    times_by_idx[i] for i in sequence_groups[t_lo]
                )  # all lower-sequence times must exceed higher-sequence times
                hi_max = max(times_by_idx[i] for i in sequence_groups[t_hi])
                if not (lo_min >= hi_max):
                    return False
            return True

        def _get_sequence_extents(valid_craters: list[Crater], sequence_groups: dict[int, list[int]]) -> dict[int, float]:
            """
            Determine the N(1) extents of each sequence number from the provided production metadata.

            Parameters
            ----------
            valid_craters : list[Crater]
                The list of Crater objects to process
            sequence_groups : dict[int, list[int]]
                A dictionary with each key being a sequence number and the values are lists of indices to the entries in the valid_craters list. It is assumed that the sequence_groups are ordered by key so that the first entry is the lowest sequence number

            Returns
            -------
            sequence_extents : dict[int, float]
                A dictionary with each key being a sequence number and values are tuples with the upper and lower N(1) extent of the associated sequence
            """
            if len(sequence_groups) == 0:
                return {}
            # First go through the sequence groups and get their production time/N(D) metadata, converting time values to N(1) values
            sequence_extents: dict[int, float] = {}
            nvals = {}
            nsigvals = {}
            sequence_numbers = list(sequence_groups.keys())
            lowest_sequence_number = sequence_numbers[0]
            highest_sequence_number = sequence_numbers[-1]

            for s, idxs in sequence_groups.items():
                n = []
                sig = []
                for i in idxs:
                    crater = valid_craters[i]
                    if crater.production_ND is not None and crater.production_ND != [None, None, None]:
                        prod_diam, nmean, nstdev = crater.production_ND
                        prod_diam *= 1e3
                        nmean /= self.ND_conversion_factor
                        nstdev /= self.ND_conversion_factor
                        # Convert from N(prod_diam) to N(1) for alignment
                        fac = self.csfd(_DSTD) / self.csfd(prod_diam)
                        nmean *= fac
                        nstdev *= fac
                        n.append(nmean)
                        sig.append(nstdev)
                    elif crater.production_time is not None and crater.production_time != [None, None]:
                        # Time doesn't map linearly with N(1), but to keep things consistent here for computing sequence boundaries we will convert the upper/lower 1 sigma bounds of production_time into upper/lower 1 sigma bounds on N(1). Later, craters with production_time+/-production_time_stdev values will be drawn using time-based values
                        tmean, tstdev = crater.production_time
                        tlo = max(tmean - tstdev, 0.0)
                        thi = tmean + tstdev
                        nlo = self.function(diameter=_DSTD, time_start=tlo)
                        nhi = self.function(diameter=_DSTD, time_start=thi)
                        nmean = (nhi + nlo) / 2
                        nstdev = (nhi - nlo) / 2
                        n.append(nmean)
                        sig.append(nstdev)
                if len(n) > 0:
                    sig = np.sqrt(np.sum([u**2 for u in sig])) / len(n)
                    n = np.mean(n)
                    nvals[s] = n
                    nsigvals[s] = sig
                else:
                    # Verify that the first entry has valid production metadata, otherwise we have no way to anchor the sequences to values of N(1)
                    if s == lowest_sequence_number:
                        raise ValueError(
                            "At least one crater with the lowest production_sequence number must have either a production_time or production_ND value"
                        )
                    elif s == highest_sequence_number:
                        nvals[s] = 0.0
                        nsigvals[s] = 0.0
                    else:
                        # Hold off on calculating the empty slots until all of the sequence numbers with production metadata have been computed
                        nvals[s] = None
                        nsigvals[s] = 0.0

            # Linearlly interpolate using the N(1) vs sequence_number value to compute the mean N(1) value for any that doen't have one
            empties = [s for s in sequence_numbers if nvals[s] is None]
            filled = [s for s in sequence_numbers if nvals[s] is not None]
            for s in empties:
                nvals[s] = np.interp(s, xp=filled, fp=[nvals[s] for s in filled])

            # Iteratively compute upper and lower extents until there is no overlap between any sequences
            for _ in range(10):
                for i, s in enumerate(sequence_numbers):
                    n = nvals[s]
                    sig = nsigvals[s]
                    nlo = max(n - sig, 0.0)
                    nhi = n + sig
                    if s != lowest_sequence_number:
                        sprev = sequence_numbers[i - 1]
                        nloprev = max(nvals[sprev] - nsigvals[sprev], 0.0)
                        nhi += (nloprev - nhi) / 2
                    if s != highest_sequence_number:
                        snext = sequence_numbers[i + 1]
                        nhinext = nvals[snext] + nsigvals[snext]
                        nlo += (nhinext - nlo) / 2
                        nlo = max(nlo, 0.0)
                    if nlo > nhi:
                        tmp = nlo
                        nlo = nhi
                        nhi = tmp
                    sequence_extents[s] = (float(nlo), float(nhi))

                nvals = {k: (nlo + nhi) / 2 for k, (nlo, nhi) in sequence_extents.items()}
                nsigvals = {k: (nhi - nlo) / 2 for k, (nlo, nhi) in sequence_extents.items()}

                # Check that both nlo and nhi are monotonically decreasing with increasing sequence number and break out if they are
                nlo_values = [v[0] for v in sequence_extents.values()]
                nhi_values = [v[1] for v in sequence_extents.values()]
                if all(x >= y for x, y in zip(nlo_values[:-1], nlo_values[1:], strict=True)) and all(
                    x >= y for x, y in zip(nhi_values[:-1], nhi_values[1:], strict=True)
                ):
                    break

            # Ensure that all boundaries are consistent such that nlo of one sequences is the same value as nhi of the next
            for i, s in enumerate(sequence_numbers):
                nlo, nhi = sequence_extents[s]
                if s != lowest_sequence_number:
                    sprev = sequence_numbers[i - 1]
                    nlo_prev, nhi_prev = sequence_extents[sprev]
                    if nhi != nlo_prev:
                        nhi = nlo_prev
                if s != highest_sequence_number:
                    snext = sequence_numbers[i + 1]
                    nlo_next, nhi_next = sequence_extents[snext]
                    if nlo != nhi_next:
                        nlo = nhi_next
                sequence_extents[s] = (nlo, nhi)

            return sequence_extents

        # Keep only craters that can produce a time
        valid_craters: list[Crater] = []
        for c in craters:
            if (
                (c.production_ND is not None and c.production_ND != [None, None, None])
                or (c.production_time is not None and c.production_time != [None, None])
                or (c.production_sequence is not None)
            ):
                valid_craters.append(c)
            else:
                print(
                    f"Skipping crater {c.name if c.name is not None else format_large_units(c.diameter, quantity='length')} because it has no production metadata"
                )

        # Partition craters between those tagged with a sequence number vs untagged
        tagged_idx: list[int] = []
        untagged_idx: list[int] = []
        sequence_groups: dict[int, list[int]] = {}

        # Group tagged craters
        for i, c in enumerate(valid_craters):
            sequence = c.production_sequence
            if sequence is None:
                untagged_idx.append(i)
            else:
                sequence = int(sequence)
                tagged_idx.append(i)
                sequence_groups.setdefault(sequence, []).append(i)

        sequence_groups = dict(sorted(sequence_groups.items()))
        sequence_extents = _get_sequence_extents(valid_craters, sequence_groups)

        # Draw untagged once (unconstrained)
        times_by_idx: dict[int, float] = {}
        for i in untagged_idx:
            times_by_idx[i] = _draw_quasimc_time(valid_craters[i])
        for i in tagged_idx:
            times_by_idx[i] = None

        # Draw tagged iteratively until constraints are satisfied
        if tagged_idx:
            for i in tagged_idx:
                times_by_idx[i] = _draw_quasimc_time(valid_craters[i], sequence_extents)

            # Check for consistency, then repeat if necessary
            if not _sequence_is_consistent(times_by_idx, sequence_groups):
                raise RuntimeError("Could not satisfy production_sequence ordering constraints ")

        # Build output crater list
        processed: list[Crater] = []
        smallest_diameter = self.diameter_range[1]
        for i, c in enumerate(valid_craters):
            t = times_by_idx.get(i)
            if t is None:
                continue
            processed.append(Crater.maker(crater=c, time=t))
            if self.generator_type == "crater":
                smallest_diameter = min(smallest_diameter, c.diameter)
            else:
                smallest_diameter = min(smallest_diameter, c.projectile_diameter)

        if set_max_diameter_from_quasimc and processed:
            self.diameter_range = (self.diameter_range[0], smallest_diameter)

        processed.sort(key=lambda c: c.time, reverse=True)
        self._quasimc_craters = processed
        return self._quasimc_craters

    def read_quasimc_file(
        self,
        filename: Path | str | None = None,
        **kwargs: Any,
    ) -> list[Crater]:
        """
        Reads in the quasi-Monte Carlo crater from file, processes the ages and sorts the crater list by age in order of decreasing age (increasing time) and computes production_time if present, production_ND if present, and will drop craters that have neither.

        Parameters
        ----------
        filename : Path | str | None, optional
            The path to the quasi-Monte Carlo file. If not provided, the Simulation must have a quasimc_file attribute set. This function will replace the stored quasimc_file attribute.
        set_max_diameter_from_quasimc : bool, optional
            If True, the largest_crater attribute of the Simulation will be set to the smallest crater in the quasi-Monte Carlo file. Default is True.
        """
        if filename is not None:
            self._quasimc_file = filename
        input_craters = Crater.from_file(self.quasimc_file)
        return self.process_quasimc_craters(input_craters, **kwargs)

    @property
    def diameter_range(self) -> tuple[float, float]:
        """
        The lower and upper range of diameters that will be returned from :py:meth:`~cratermaker.components.production.Production.sample`.
        """
        return self._diameter_range

    @diameter_range.setter
    def diameter_range(self, value: PairOfFloats | None):
        # Reset to default
        if value is None:
            self._diameter_range = (0, np.inf)
        elif not isinstance(value, (list, tuple, ArrayLike)) or len(value) != 2:
            raise TypeError("diameter_range must be a list, tuple, or array of length 2")
        else:
            self._diameter_range = (float(value[0]), float(value[1]))
        return

    @parameter
    def ND_conversion_factor(self) -> float:
        """
        The conversion factor used to convert N(D) format into number of craters per m².

        Only three options are allowed, 1.0 for units of '# per m²', 1e6 for units of '# per km²', or 1e12 for units of per 10⁶ km². The default value is 1e12."
        """
        if self._ND_conversion_factor is None:
            return 1e12

    @ND_conversion_factor.setter
    def ND_conversion_factor(self, value: float | None):
        if value is None or value == 1e12:
            self._ND_conversion_factor = None
            self._ND_unit = None
        elif value == 1e6:
            self._ND_conversion_factor = value
            self._ND_unit = "# per km²"
        elif value == 1.0:
            self._ND_conversion_factor = value
            self._ND_unit = "# per m²"
        else:
            raise ValueError("ND_conversion factor can only be 1.0, 1e6, or 1e12")

    @property
    def ND_unit(self) -> float:
        if self._ND_unit is None:
            return "# per 10⁶ km²"
        else:
            return self._ND_unit


import_components(__name__, __path__)
