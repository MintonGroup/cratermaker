from contextlib import suppress
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from matplotlib.axes import Axes
from numpy.random import Generator
from numpy.typing import ArrayLike
from tqdm import tqdm

from cratermaker.components.counting import Counting
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology
from cratermaker.components.production import Production
from cratermaker.components.projectile import Projectile
from cratermaker.components.scaling import Scaling
from cratermaker.components.surface import Surface
from cratermaker.components.surface.hireslocal import HiResLocalSurface
from cratermaker.components.target import Target
from cratermaker.constants import (
    _COMPONENT_NAMES,
    _CONFIG_FILE_NAME,
    FloatLike,
    PairOfFloats,
)
from cratermaker.core.base import CratermakerBase, _convert_for_yaml
from cratermaker.utils.general_utils import _set_properties, format_large_units, parameter


class Simulation(CratermakerBase):
    """
    Creates a simulation of a crater population on a target body. It allows for the generation of craters based on a variety of parameters, including the target body, scaling laws, production functions, and morphology models.

    Parameters
    ----------
    target: Target or str, optional, default "Moon"
        Name target body for the simulation, default is "Moon".
    scaling : Scaling or str, optional
        The projectile->crater size scaling model to use from the components library. The default is "montecarlo".
    production: Production or str, optional
        The production function model to use from the components library that defines the production function used to populate the surface with craters. If none provided,
        then the default will be based on the target body, with the NeukumProduction crater-based scaling law used if the target
        body is the Moon or Mars, the NeukumProduction projectile-based scaling law if the target body is Mercury, Venus, or
        Earth, and a simple power law model otherwise.
    morphology : str, optional
        The model used to generate the morphology of the crater. If none provided, then the default will "simplemoon", which is similar to the one used by CTEM.
    projectile : str, optional
        The projectile model to use from the components library, which is used to generate the projectile properties for the simulation, such as velocity and density. The default is "asteroids" when target is Mercury, Venus, Earth, Moon, Mars, Ceres, or Vesta, and "comets" otherwise.
    surface : str, optional
        The name of the surface used for the surface. Default is "icosphere".
    counting : Counting or str, optional
        The crater counting model to use from the components library. Default is "minton2019".
    simdir : str | Path
        |simdir|
    rng : numpy.random.Generator | None
        |rng|
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        |rng_seed|
    rng_state : dict, optional
        |rng_state|
    reset : bool, optional
        Flag to indicate whether to reset the simulation or resume from an old simulation. If False, the simulation will attempt to load the previous state from the config file. Default is False if `ask_overwrite=False` and a config file is detected, otherwise default is True.
    do_counting : bool, optional
        If True, the counting component will keep track of emplaced craters during the simulation. Default is True.
    ask_overwrite : bool, optional
        If True, the user will be prompted before overwriting any existing files. Default is True.
    **kwargs : Any
        |kwargs|, including those for component function constructors. Refer to the documentation of each component module for details.
    """

    def __init__(
        self,
        *,  # Enforce keyword-only arguments
        target: Target | str | None = None,
        scaling: Scaling | str | None = None,
        production: Production | str | None = None,
        morphology: Morphology | str | None = None,
        projectile: Projectile | str | None = None,
        surface: Surface | str | None = None,
        counting: Counting | str | None = None,
        simdir: str | Path | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        reset: bool = None,
        ask_overwrite: bool = True,
        do_counting: bool = True,
        **kwargs: Any,
    ):
        super().__init__(simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_scaling", scaling)
        object.__setattr__(self, "_production", production)
        object.__setattr__(self, "_morphology", morphology)
        object.__setattr__(self, "_projectile", projectile)
        object.__setattr__(self, "_surface", surface)
        object.__setattr__(self, "_counting", counting)
        object.__setattr__(self, "_interval", None)
        object.__setattr__(self, "_elapsed_time", None)
        object.__setattr__(self, "_time", None)
        object.__setattr__(self, "_elapsed_n1", None)
        object.__setattr__(self, "_smallest_crater", None)
        object.__setattr__(self, "_smallest_projectile", None)
        object.__setattr__(self, "_largest_crater", None)
        object.__setattr__(self, "_largest_projectile", None)
        object.__setattr__(self, "_ask_overwrite", None)
        object.__setattr__(self, "_config_readonly", False)
        if reset is None:
            if ask_overwrite and self.config_file.exists():
                response = input(
                    "Old run detected. Enter y to reset the simulation and n to resume from the previous state. To disable this message, pass `ask_overwrite=False` to this function. (y/[N]): "
                )
                if response.lower() == "y":
                    print("Resetting simulation.")
                    reset = True
                else:
                    reset = False
            else:
                reset = False

        if not reset and self.config_file.exists():
            config_file = self.config_file
            object.__setattr__(self, "_config_readonly", True)
        else:
            config_file = None

        config_override = {}
        for component in _COMPONENT_NAMES:
            # Set to true if a local variable from the argument list with the component name is set to something other than None, otherwise false
            config_override[component] = getattr(self, f"_{component}") is not None

        _, unmatched = _set_properties(
            self,
            target=target,
            scaling=scaling,
            production=production,
            morphology=morphology,
            projectile=projectile,
            surface=surface,
            counting=counting,
            config_file=config_file,
            ask_overwrite=ask_overwrite,
            **vars(self.common_args),
        )

        for component in _COMPONENT_NAMES:
            if config_override[component]:
                # If the component is set to something other than None, then remove it from the unmatched dictionary
                unmatched.pop(f"{component}_config", None)

        production_config = unmatched.pop("production_config", {})
        scaling_config = unmatched.pop("scaling_config", {})
        surface_config = unmatched.pop("surface_config", {})
        morphology_config = unmatched.pop("morphology_config", {})
        target_config = unmatched.pop("target_config", {})
        projectile_config = unmatched.pop("projectile_config", {})
        if do_counting:
            counting_config = unmatched.pop("counting_config", {})
        kwargs.update(unmatched)
        kwargs = {**kwargs, **vars(self.common_args)}

        target_config = {**target_config, **kwargs}
        self.target = Target.maker(self.target, **target_config)

        production_config = {**production_config, **kwargs}
        self.production = Production.maker(self.production, target=self.target, **production_config)

        projectile_config = {**projectile_config, **kwargs}
        self.projectile = Projectile.maker(self.projectile, target=self.target, **projectile_config)

        scaling_config = {**scaling_config, **kwargs}
        self.scaling = Scaling.maker(
            self.scaling,
            target=self.target,
            projectile=self.projectile,
            ask_overwrite=self.ask_overwrite,
            **scaling_config,
        )

        surface_config = {
            **surface_config,
            **kwargs,
            "ask_overwrite": self.ask_overwrite,
        }
        if "superdomain_scale_factor" not in surface_config:
            surface_config["superdomain_scale_factor"] = (
                None  # This will trigger setting of the superdomain after the Morphology and Scaling models are set
            )
        self.surface = Surface.maker(
            self.surface,
            target=self.target,
            reset=reset,
            **surface_config,
        )

        if do_counting:
            counting_config = {**counting_config, **kwargs}
            self.counting = Counting.maker(
                self.counting,
                surface=self.surface,
                reset=reset,
                **counting_config,
            )

        morphology_config = {**morphology_config, **kwargs}
        self.morphology = Morphology.maker(
            self.morphology,
            surface=self.surface,
            production=self.production,
            counting=self.counting,
            **morphology_config,
        )

        # If this is a variant of the HiResLocalSurface we need to check to see if it has a grid yet.
        # This is because when creating a new Surface object of this type, the grid generation is deferred until the Scaling and Morphology objects are initialized in order to set the superdomain properly.
        if issubclass(self.surface.__class__, HiResLocalSurface) and self.surface.uxgrid is None:
            self.surface.set_superdomain(
                scaling=self.scaling,
                morphology=self.morphology,
                reset=reset,
                **surface_config,
            )

        if reset:
            # The Surface has already had its reset method called.
            skip_components = ["surface"]
            if not do_counting:
                skip_components.append("counting")
            self.reset(ask_overwrite=ask_overwrite, skip_component=skip_components)
        else:
            # Now that all components are initialized, we turn off the config_readonly flag so that any changes to the simulation parameters will be saved to the config file.
            self._config_readonly = False

        self.to_config()

        return

    def __str__(self) -> str:
        """
        Returns a string representation of the Simulation object.
        """
        return (
            f"<Simulation>\n\n"
            f"{self.counting}\n\n"
            f"{self.morphology}\n\n"
            f"{self.production}\n\n"
            f"{self.projectile}\n\n"
            f"{self.scaling}\n\n"
            f"{self.surface}\n\n"
            f"{self.target}\n\n"
            f"<Current state>\n"
            f"Current time : {format_large_units(self.time, quantity='time')} before present\n"
            f"Elapsed time: {format_large_units(self.elapsed_time, quantity='time')}\n"
            f"Elapsed N_1 : {self.elapsed_n1} #/m^2\n"
            f"Interval    : {self.interval}\n"
            f"simdir      : {str(self.simdir)}\n"
        )

    def __repr__(self) -> str:
        config = self.to_config(save_to_file=False)
        txt = f"{self.__class__.__name__}("
        for k, v in config.items():
            if isinstance(v, str):
                v = f"'{v}'"
            txt += f"\n    {k}={v},"
        txt += "\n)"
        return txt

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        if name not in self._user_defined:
            return
        # Avoid recursive calls during initialization or early access
        if hasattr(self, "to_config") and callable(getattr(self, "to_config", None)) and _convert_for_yaml(value) is not None:
            with suppress(Exception):
                self.to_config()

    def run(
        self,
        age: FloatLike | None = None,
        time_start: FloatLike | None = None,
        time_end: FloatLike | None = None,
        diameter_number: PairOfFloats | None = None,
        time_interval: FloatLike | None = None,
        ninterval: int | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Run the simulation over a specified interval using the current production function.

        Parameters
        ----------
        age : FloatLike, optional
            Start age in My relative to the present for the simulation, used to compute the starting point of the production function.
            Default is None, which requires 'time_start' or `diameter_number` to be set.
        time_start : Floatlike, optional
            An alternative to `age` that specifies the starting time in My relative to the present for the simulation, used to compute the starting point of the production function. This is used in conjunction with `time_end` in order to allow for simulations that span a range of time rather than being of a specific age. Default is None, which requires either `age` or `diameter_number` to be set.
        time_end : FloatLike, optional
            Ending time in My relative to the present for the simulation, used to compute the ending point of the production function.
            Default is 0 (present day) if not provided but `time_start` is provided.
        diameter_number: PairOfFloats, optional
            Cumulative number and diameter pair (D, N) that defines the total crater accumulation for the given production function. Default is None, which requires `age` or `time_start` and `time_end` to be set.
        time_interval : FloatLike, optional
            Interval in My for outputting intermediate results. If not provided, calculated as `age` / `ninterval` or (`time_start - time_end`) / `ninterval` if `ninterval` is provided, otherwise set to the total simulation duration (e.g. `ninterval=1`).
        ninterval : int, optional
            Number of intervals for outputting results. This has a special use case where one can specify age-based inputs but output are computed in equal cumulative number intervals and vice versa.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        This function allows defining the simulation parameters either in terms of time or crater frequency (cumulative number). The
        arguments `age`, `time_start`, `time_end` are mutually exclusive with `diameter_number` and `time_interval` is mutually exclusive with `ninterval`.

        The initial state of the simulation (before any craters are emplaced) is always saved automatically.
        As a result, the total number of saved states will be `ninterval + 1`, where `ninterval` is the number
        of simulation intervals requested.

        Examples
        --------
        .. code-block:: python

            # Create a simulation object with default parameters (Moon, NeukumProduction, etc.)
            sim = cratermaker.Simulation()

            # Run the simulation for 3.8 billion years, saving the results every 100 million years
            sim.run(age=3.8e3, time_interval=100.0)

            # Run the simulation for 3.8 billion years, saving 100 intervals with equal cumulative number intervals
            sim.run(age=3.8e3, ninterval=100)

            # Run the simulation to create 80 craters larger than 300 km and output 100 equal cumulative number intervals
            sim.run(diameter_number=(300e3, 80), ninterval=100)

            # Run the simulation from 3.8 billion years to 3.0 billion years, saving the results every 100 million years
            sim.run(time_start=3.8e3, time_end=3.0e3, time_interval=100.0)
        """

        def _validate_run_args(**kwargs: Any) -> dict:
            """
            Validate all the input arguments to the sample method. This function will raise a ValueError if any of the arguments are invalid.

            Parameters
            ----------
            kwargs : Any
                A dictionary of all the arguments passed to the run method, including age, time_start, time_end, diameter_number, time_interval, ninterval, and any additional arguments passed through **kwargs.

            Returns
            -------
            A dict containing all arguments listed in Parameters above, as well as `is_time_interval`, which is a boolean flag indicating
            whether or not the simulation is being run in equal age intervals or equal number intervals.

            Raises
            ------
            ValueError
                If any of the following conditions are met:
                - Neither the age nore the diameter_number argument is provided.
                - Both the age and diameter_number arguments are provided.
                - Both the age and either time_start or time_end arguments are provided.
                - Both the time_start and diameter_number arguments are provided.
                - The age, time_start, or time_end arguments are provided but are not a scalar.
                - The time_interval is provided but is not a positive scalar.
                - The time_interval provided is greater than age or time_start - time_end
                - The diameter_number argument is not a pair of values, or any of them are less than 0
                - The time_interval and nintervaql arguments are both provided.
                - The ninterval is provided but is not an integer or is less than 1.
            """
            # Determine whether we are going to do equal time intervals or equal number intervals
            age = kwargs.pop("age", None)
            time_start = kwargs.pop("time_start", None)
            time_end = kwargs.pop("time_end", None)
            diameter_number = kwargs.pop("diameter_number", None)
            time_interval = kwargs.pop("time_interval", None)
            ninterval = kwargs.pop("ninterval", None)

            if age is not None:
                if time_start is not None or time_end is not None:
                    raise ValueError("Cannot specify both age and time_start or time_end")
                if diameter_number is not None:
                    raise ValueError("Cannot specify both age and diameter_number")
                if not np.isscalar(age):
                    raise ValueError("age must be a scalar value")
                # as age is just a convenience variable, we replace it with time_start = age and time_end = 0
                time_start = age
                time_end = 0.0
                del age
            if time_start is not None:
                if time_end is None:
                    time_end = 0.0
                if diameter_number is not None:
                    raise ValueError("Cannot specify both time_start and diameter_number")
                if not np.isscalar(time_start):
                    raise ValueError("time_start must be a scalar value")
                if not np.isscalar(time_end):
                    raise ValueError("time_end must be a scalar value")
            elif time_end is not None:
                if self._time is None:
                    raise ValueError("time_end cannot be used without time_start")
                time_start = self.time
            elif diameter_number is None:
                raise ValueError("Must provide one of age, time_start, or diameter_number")

            if ninterval is not None:
                if not isinstance(ninterval, int):
                    raise TypeError("ninterval must be an integer")
                if ninterval < 1:
                    raise ValueError("ninterval must be greater than zero")
                if time_interval is not None:
                    raise ValueError("Cannot specify both ninterval and time_interval")
            elif time_interval is None:
                ninterval = 1

            is_time_interval = time_interval is not None

            # Validate arguments using the production function validator first, which will convert time-based values to diameter_number-based ones
            kwargs["diameter_range"] = (
                self._get_smallest_diameter(),
                self._get_largest_diameter(),
            )
            kwargs["area"] = self.surface.area
            kwargs["time_start"] = time_start
            kwargs["time_end"] = time_end
            kwargs["diameter_number"] = diameter_number
            kwargs["ninterval"] = ninterval
            kwargs["time_interval"] = time_interval
            kwargs = self.production._validate_sample_args(**kwargs)

            if is_time_interval:
                if time_interval is None:
                    if ninterval is None:
                        ninterval = 1
                    time_interval = (time_start - time_end) / ninterval
                else:
                    if time_interval > time_start - time_end:
                        raise ValueError("time_interval must be less than age or time_start - time_end")
                    elif time_interval <= 0:
                        raise ValueError("time_interval must be greater than zero")
                    ninterval = int(np.ceil((time_start - time_end) / time_interval))

                kwargs["time_interval"] = time_interval
                kwargs["ninterval"] = ninterval
            else:
                diameter_number = kwargs.get("diameter_number", None)
                diameter_number_end = kwargs.get("diameter_number_end", None)
                if diameter_number is None:
                    raise ValueError("Something went wrong! diameter_number should be set by self.production_validate_sample_args")
                if diameter_number_end is None:
                    raise ValueError(
                        "Something went wrong! diameter_number_end should be set by self.production_validate_sample_args"
                    )
                if ninterval is None:
                    ninterval = 1
                diameter_number_interval = (
                    diameter_number[0],
                    (diameter_number[1] - diameter_number_end[1]) / ninterval,
                )
                kwargs["diameter_number_interval"] = diameter_number_interval
                kwargs["ninterval"] = ninterval

            kwargs["is_time_interval"] = is_time_interval

            # Remove unecessary arguments that came out of the production._validate_sample_args method
            kwargs.pop("diameter_range")
            kwargs.pop("area")
            kwargs.pop("return_age")

            return kwargs

        arguments = {
            "age": age,
            "time_start": time_start,
            "time_end": time_end,
            "time_interval": time_interval,
            "diameter_number": diameter_number,
            "ninterval": ninterval,
            **kwargs,
        }

        arguments = _validate_run_args(**arguments)
        time_start = arguments.pop("time_start", None)
        time_end = arguments.pop("time_end", None)
        time_interval = arguments.pop("time_interval", None)
        diameter_number = arguments.pop("diameter_number", None)
        diameter_number_interval = arguments.pop("diameter_number_interval", None)
        ninterval = arguments.pop("ninterval", None)
        is_time_interval = arguments.pop("is_time_interval", None)
        validate_inputs = kwargs.pop("validate_inputs", False)
        if ninterval == 1:
            is_time_interval = True
            time_interval = time_start - time_end

        if not is_time_interval:
            diameter_number_density_start = (
                diameter_number[0],
                diameter_number[1] / self.surface.area,
            )
            time_start = self.production.age_from_D_N(*diameter_number_density_start, validate_inputs=validate_inputs)

        # If this is a fresh run, we need to set the value of the current time based on the requested starting condition.
        if self._time is None:
            self.time = time_start
            self.interval = 0

        # If this is a restarted run, we need to distinguish between a true restart and a continuation with different parameters
        if time_start < self.time:
            raise RuntimeError(
                "Starting time cannot be later than the current time. Choose a starting time value equal to or larger than the current time, or reset this simulation."
            )
        if is_time_interval:
            initial_interval = int((time_start - self.time) / time_interval)
        else:
            delta_n1_start = self.production.function(
                diameter=1000.0,
                time_start=time_start,
                time_end=self.time,
                validate_inputs=validate_inputs,
            ).item()
            n1_interval = (
                self.production.function(
                    diameter=1000.0,
                    time_start=time_start,
                    time_end=time_end,
                    validate_inputs=validate_inputs,
                ).item()
                / ninterval
            )
            initial_interval = int(delta_n1_start / n1_interval)

        self.save(merge_with_existing=False, **kwargs)
        for i in tqdm(
            range(initial_interval, ninterval),
            total=ninterval,
            initial=initial_interval,
            desc="Simulation interval",
            unit="interval",
            position=3,
            leave=True,
        ):
            if self.do_counting:
                self.counting._emplaced = []
            if is_time_interval:
                time = time_start - i * time_interval
                current_time_end = time_start - (i + 1) * time_interval
                if current_time_end < time_end:
                    current_time_end = time_end
                self.populate(time_start=time, time_end=current_time_end)
            else:
                current_diameter_number = (
                    diameter_number[0],
                    diameter_number[1] - i * diameter_number_interval[1],
                )
                current_diameter_number_end = (
                    diameter_number[0],
                    diameter_number[1] - (i + 1) * diameter_number_interval[1],
                )
                self.populate(
                    diameter_number=current_diameter_number,
                    diameter_number_end=current_diameter_number_end,
                )
                current_diameter_number_density = (
                    current_diameter_number[0],
                    current_diameter_number[1] / self.surface.area,
                )
                time = self.production.age_from_D_N(*current_diameter_number_density, validate_inputs=validate_inputs)
                if current_diameter_number_end[1] > 0:
                    current_diameter_number_density_end = (
                        current_diameter_number_end[0],
                        current_diameter_number_end[1] / self.surface.area,
                    )

                    current_time_end = self.production.age_from_D_N(
                        *current_diameter_number_density_end,
                        validate_inputs=validate_inputs,
                    )
                else:
                    current_time_end = 0.0
                if current_time_end < 0.0:
                    current_time_end = 0.0
                time_interval = time - current_time_end
            self.elapsed_time += time_interval
            self.elapsed_n1 += self.production.function(
                diameter=1000.0,
                time_start=time,
                time_end=current_time_end,
                validate_inputs=validate_inputs,
            ).item()
            self.time = current_time_end
            self.interval += 1

            self.save(merge_with_existing=False, **kwargs)

        return

    def populate(
        self,
        age: FloatLike | None = None,
        time_start: FloatLike | None = None,
        time_end: FloatLike | None = None,
        diameter_number: PairOfFloats | None = None,
        diameter_number_end: PairOfFloats | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Populate the surface with craters over a specified interval using the current production function.

        Parameters
        ----------
        age : FloatLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        time_start : FloatLike, optional
            An alternative to `age` that specifies the starting time in My relative to the present for the simulation, used to compute the starting point of the production function. This is used in conjunction with `time_end` in order to allow for simulations that span a range of time rather than being of a specific age. Default is None, which requires either `age` or `diameter_number`
        time_end : FloatLike, optional
            The ending time in My relative to the present for the simulation, used to compute the ending point of the production function. Default is 0 (present day).
        diameter_number : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function will convert this value
            to a corresponding age and use the production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function will convert this
            value to a corresponding reference age and use the production function for a given age.
        """
        if not hasattr(self, "production"):
            raise RuntimeError("No production function defined for this simulation")
        elif not hasattr(self.production, "generator_type"):
            raise RuntimeError("The production function is not properly defined. Missing 'generator_type' attribute")
        elif self.production.generator_type not in ["crater", "projectile"]:
            raise RuntimeError(f"Invalid production function type {self.production.generator_type}")
        if age is not None:
            if time_start is not None:
                raise ValueError("Cannot specify both age and time_start")
            time_start = age
            del age

        from_projectile = self.production.generator_type == "projectile"
        diam_key = "projectile_diameter" if from_projectile else "diameter"
        diam_max = self._get_largest_diameter(from_projectile=from_projectile)
        diam_min = self._get_smallest_diameter(from_projectile=from_projectile)

        # Loop over each face in the mesh to build up a population of craters in this interval. This is done because faces may
        # not all have the same surface area, the range of crater sizes that can be formed on each face may be different.
        impact_diameters = []
        impact_times = []
        impact_locations = []

        # Process each bin
        for i, face_indices in enumerate(self.surface.face_bin_indices):
            total_bin_area = self.surface.face_bin_area[i]
            area_ratio = total_bin_area / self.surface.area

            diam_min = self._get_smallest_diameter(self.surface.face_bin_min_sizes[i], from_projectile=from_projectile)
            diameter_number_local = (diameter_number[0], diameter_number[1] * area_ratio) if diameter_number is not None else None

            if diameter_number_end is not None:
                diameter_number_end_local = (
                    diameter_number_end[0],
                    diameter_number_end[1] * area_ratio,
                )
            else:
                diameter_number_end_local = None

            diameters, times = self.production.sample(
                time_start=time_start,
                time_end=time_end,
                diameter_number=diameter_number_local,
                diameter_number_end=diameter_number_end_local,
                diameter_range=(diam_min, diam_max),
                area=total_bin_area,
                **kwargs,
            )
            if diameters.size > 0:
                impact_diameters.extend(diameters.tolist())
                impact_times.extend(times.tolist())

                # Get the relative probability of impact onto any particular face then get the locations of the impacts
                p = self.surface.face_area[face_indices] / total_bin_area
                face_indices = self.rng.choice(face_indices, size=diameters.shape, p=p)
                locations = self.surface.get_random_location_on_face(face_indices)
                impact_locations.extend(np.array(locations).T.tolist())

        if len(impact_diameters) > 0:
            craterlist = []
            # Sort the times, diameters, and locations so that they are in order of decreasing age
            sort_indices = np.argsort(impact_times)[::-1]
            impact_diameters = np.asarray(impact_diameters)[sort_indices]
            impact_times = np.asarray(impact_times)[sort_indices]
            impact_locations = np.array(impact_locations)[sort_indices]
            for diameter, location, time in tqdm(
                zip(impact_diameters, impact_locations, impact_times, strict=False),
                total=len(impact_diameters),
                desc="Generating crater population",
                unit="crater",
                position=0,
                leave=False,
            ):
                diam_arg = {diam_key: diameter}
                craterlist.append(
                    Crater.maker(
                        location=location,
                        age=time,
                        scaling=self.scaling,
                        **diam_arg,
                        **vars(self.common_args),
                        **kwargs,
                    )
                )
            self.emplace(craterlist)

        return

    def emplace(self, craters: list[Crater] | Crater | None = None, **kwargs: Any) -> list[Crater]:
        """
        Emplace one or more craters in the simulation.

        This method orchestrates the creation and placement of a crater in the
        simulation. It can create a crater directly or based on the characteristics
        of a projectile.

        Parameters
        ----------
        craters : Crater or list of Crater objects, optional
            The Crater object(s) to be emplaced. If provided, this will be used directly. Otherwise, a single crat er will be generated based on the keyword arguments.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        list[Crater]
            A list of the Crater objects that were emplaced in the simulation. Returns an empty list if no craters were emplaced.

        Notes
        -----
        The keyword arguments provided are passed down to :meth:`Crater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        Examples
        --------
        .. code-block:: python

            from cratermaker import Simulation, Crater

            sim = Simulation()

            # Create a crater with specific diameter
            sim.emplace(diameter=10.0e3)

            # Create a crater based on a projectile with given mass and projectile_velocity
            sim.emplace(projectile_mass=1e15, projectile_velocity=20e3)

            # Create a crater with a specific transient diameter and location
            sim.emplace(transient_diameter=50e3, location=(43.43, -86.92))

            # Create multiple craters
            craters = [Crater.maker(diameter=20.0e3), Crater.maker(diameter=20.0e3)]
            sim.emplace(craters)

        """
        crater_args = {**kwargs, **vars(self.common_args)}
        if craters is None and "scaling" not in crater_args:
            crater_args["scaling"] = self.scaling
        return self.morphology.emplace(craters=craters, **crater_args)

    def save(self, **kwargs: Any) -> None:
        """
        Save the current simulation state to a file.

        Parameters
        ----------
        **kwargs : Any
            Additional keyword argumments to pass to the component save methods.
        """
        self.surface.save(
            interval=self.interval,
            time_variables=self.time_variables,
            **kwargs,
        )

        if self.do_counting:
            self.counting.save(interval=self.interval, **kwargs)

        self.to_config(**kwargs)
        self.plot(show=False, save=True, **kwargs)

        return

    def export(
        self,
        driver: str = "OpenCraterTool",
        interval: int = -1,
        ask_overwrite: bool | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Export component output to a specified file format.

        Parameters
        ----------
        driver : str, optional
            The driver to use export the data to. Supported formats are 'OpenCraterTool', 'VTK' or a driver supported by GeoPandas ('GPKG', 'ESRI Shapefile', etc.). This is overridden if either the filename or file_extension parameters are provided. Default is 'OpenCraterTool'.
        interval : int, optional
            The interval number to export. Default is -1 (the most current interval saved in the simulation).
        ask_overwrite : bool, optional
            If True, the user will be prompted before overwriting any existing files. Default is set to the value provided when the Simulation object was created.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        The default driver is 'OpenCraterTool', which is designed to output data into a format that is relatively easy to import into QGIS with the OpenCraterTool plugin. This will create a GeoTIFF file representation of the surface, and a set of SCC files for the crater counting data if counting is enabled.
        """
        if ask_overwrite is None:
            ask_overwrite = self.ask_overwrite
        if interval < 0:
            interval = self.interval + 1 + interval
        self.save()
        if driver.lower() == "opencratertool":
            surface_driver = "GeoTIFF"
            counting_driver = "SCC"
        else:
            surface_driver = driver
            counting_driver = driver
        self.surface.export(
            driver=surface_driver,
            interval=interval,
            ask_overwrite=ask_overwrite,
            **kwargs,
        )

        if self.do_counting:
            self.counting.export(
                craters=self.counting.observed,
                interval=interval,
                driver=counting_driver,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )

        return

    def plot(self, include_counting: bool = False, **kwargs: Any) -> Axes:
        """
        Plot the current state of the surface.

        Parameters
        ----------
        include_counting : bool, optional
            If True, the counting data will be included in the plot if counting is enabled. Default is False
        **kwargs : Any
        |kwargs|

        Returns
        -------
        Axes
            The matplotlib Axes object created by the surface plot method.
        """
        label = kwargs.pop("label", f"Time: {self.time:.0f} My bp\nAge : {self.elapsed_time:.0f} My")
        plot_style = kwargs.pop("plot_style", "hillshade")
        if include_counting and self.do_counting:
            ax = self.counting.plot(interval=self.interval, plot_style=plot_style, label=label, **kwargs)
        else:
            ax = self.surface.plot(interval=self.interval, plot_style=plot_style, label=label, **kwargs)
        return ax

    def show(self, engine: str = "pyvista", **kwargs: Any) -> None:
        """
        Show the current state of the simulated surface.

        Parameters
        ----------
        engine : str, optional
            The engine to use for plotting. Currently, only "pyvista" is supported. Default is "pyvista".
        **kwargs : Any
        |kwargs|
        """
        if "interval" not in kwargs:
            kwargs["interval"] = self.interval
        if self.do_counting:
            self.counting.show(engine=engine, **kwargs)
        else:
            self.surface.show(engine=engine, **kwargs)

        return

    def to_config(self, save_to_file: bool = True, **kwargs: Any) -> dict:
        """
        Converts values to types that can be used in yaml.safe_dump.

        This will convert various types into a format that can be saved in a human-readable YAML file. This will consolidate all of the configuration
        parameters into a single dictionary that can be saved to a YAML file. This will also remove any common arguments from the individual configurations for each component model to avoid repeating them.

        Parameters
        ----------
        save_to_file : bool, optional
            If True, the configuration will be saved to a file. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        dict[str, Any]
            A dictionary of the object's attributes that can be serialized to YAML.

        Notes
        -----
        - The function will ignore any attributes that are not serializable to human-readable YAML. Therefore, it will ignore anything that cannot be converted into a str, int, float, or bool.
        - The function will convert Numpy types to their native Python types.
        """
        sim_config = super().to_config(remove_common_args=False)
        if self.target is not None:
            sim_config["target"] = self.target.name
            sim_config["target_config"] = self.target.to_config(remove_common_args=True)
        if self.scaling is not None:
            sim_config["scaling"] = self.scaling._component_name
            sim_config["scaling_config"] = self.scaling.to_config(remove_common_args=True)
        if self.production is not None:
            sim_config["production"] = self.production._component_name
            sim_config["production_config"] = self.production.to_config(remove_common_args=True)
        if self.surface is not None:
            sim_config["surface"] = self.surface._component_name
            sim_config["surface_config"] = self.surface.to_config(remove_common_args=True)
        if self.projectile is not None:
            sim_config["projectile"] = self.projectile._component_name
            sim_config["projectile_config"] = self.projectile.to_config(remove_common_args=True)
        if self.morphology is not None:
            sim_config["morphology"] = self.morphology._component_name
            sim_config["morphology_config"] = self.morphology.to_config(remove_common_args=True)
        if self.counting is not None:
            sim_config["counting"] = self.counting._component_name
            sim_config["counting_config"] = self.counting.to_config(remove_common_args=True)

        for config in [
            "target",
            "scaling",
            "production",
            "projectile",
            "morphology",
            "surface",
            "counting",
        ]:
            # drop any empty values or {} from either f"{config} or f"{config}_config" if when they are either None or empty
            if config in sim_config and (sim_config[config] is None or sim_config[config] == {}):
                sim_config.pop(config)
            if f"{config}_config" in sim_config and (
                sim_config[f"{config}_config"] is None or sim_config[f"{config}_config"] == {}
            ):
                sim_config.pop(f"{config}_config")

        # Write the combined configuration to a YAML file
        if save_to_file and not self._config_readonly:
            with Path.open(self.config_file, "w") as f:
                yaml.safe_dump(sim_config, f, indent=4)

        return sim_config

    def reset(
        self,
        ask_overwrite: bool | None = None,
        skip_component: str | list[str] | None = None,
    ) -> None:
        """
        Reset the simulation by clearing all data and files associated with it.

        Parameters
        ----------
        ask_overwrite : bool, optional
            If True, the user will be prompted before overwriting any existing files. Default is what is set during initialization, which is True unless specified otherwise.
        skip_component : str or list of str, optional
            List of component names to skip during the reset process. Default is an empty list, which means all components will be reset.
        """
        if skip_component is None:
            skip_component = []
        elif isinstance(skip_component, str):
            skip_component = [skip_component]
        elif not isinstance(skip_component, list) or not all(isinstance(c, str) for c in skip_component):
            raise TypeError("skip_component must be a string or a list of strings")
        if ask_overwrite is None:
            ask_overwrite = self.ask_overwrite

        if ask_overwrite:
            files_to_remove = []
            for component in _COMPONENT_NAMES:
                if component not in skip_component and hasattr(self, component):
                    files_to_remove += getattr(self, component).saved_output_files()
            if len(files_to_remove) > 0:
                print("The following files will be deleted:")
                for f in files_to_remove:
                    print(f"  {f}")
                print("To disable this message, pass `ask_overwrite=False` to this function.")
                response = input(f"Are you sure you want to delete {len(files_to_remove)} files? [y/N]: ")
                if response.lower() != "y":
                    raise RuntimeError("User aborted the reset operation.")
        for component in _COMPONENT_NAMES:
            if component not in skip_component and hasattr(self, component):
                getattr(self, component).reset(ask_overwrite=False)

        # Remove any old local surface output files
        files_to_remove = self.surface.output_dir.glob("local_*")
        for f in files_to_remove:
            f.unlink(missing_ok=True)

        self._interval = 0
        self._elapsed_time = None
        self._time = None
        self._elapsed_n1 = None
        self._smallest_crater = 0.0  # The smallest crater will be determined by the smallest face area
        self._smallest_projectile = 0.0  # The smallest crater will be determined by the smallest face area
        self._largest_crater = np.inf  # The largest crater will be determined by the target body radius
        self._largest_projectile = np.inf  # The largest projectile will be determined by the target body radius

        return

    def update_elevation(self, *args: Any, **kwargs: Any) -> None:
        """
        Set the elevation on the surface. Delegates to the Surface object.

        Parameters
        ----------
        *args: Variable length argument list to pass to self.surface.update_elevation.
        |kwargs|
        """
        return self.surface.update_elevation(*args, **kwargs)

    def _get_smallest_diameter(self, face_size: ArrayLike | None = None, from_projectile: bool = False) -> float:
        """
        Get the smallest possible crater or projectile be formed on a face.

        Parameters
        ----------
        face_size : FloatLike, optional
            The effective size of the face to determine the smallest crater size that can be formed on it. If None, the size of the smallest face on the surface will be used.
        from_projectile : bool, optional
            If True, the smallest projectile diameter will be returned instead of the smallest crater diameter. Default is False.

        Returns
        -------
        float
            The smallest possible crater or projectile diameter that can be formed on the surface.
        """
        if face_size is None:
            face_size = np.min(self.surface.face_size)
        if from_projectile:
            crater = Crater.maker(
                diameter=face_size,
                angle=90.0,
                projectile_velocity=self.scaling.projectile_mean_velocity * 10,
                scaling=self.scaling,
                **vars(self.common_args),
            )
            return crater.projectile_diameter
        else:
            return float(face_size)

    def _get_largest_diameter(self, from_projectile: bool = False) -> float:
        """
        Get the largest possible crater or projectile that can be formed on the surface.
        """
        largest_crater = self.target.radius * 2
        if from_projectile:
            crater = Crater.maker(
                diameter=largest_crater,
                angle=1.0,
                projectile_velocity=self.scaling.projectile_mean_velocity / 10.0,
                scaling=self.scaling,
                **vars(self.common_args),
            )
            return crater.projectile_diameter
        else:
            return largest_crater

    @property
    def target(self):
        """
        The target body for the impact simulation. Set during initialization.
        """
        return self._target

    @target.setter
    def target(self, value):
        if not isinstance(value, (Target | str)):
            raise TypeError("target must be an instance of Target or str")
        self._target = value

    @property
    def surface(self):
        """
        Surface mesh data for the simulation. Set during initialization.
        """
        return self._surface

    @surface.setter
    def surface(self, value):
        if not isinstance(value, (Surface | str | type)) or (isinstance(value, type) and not issubclass(value, Surface)):
            raise TypeError("surface must be an instance of Surface, a subclass of Surface, or str")
        self._surface = value

    @property
    def production(self):
        """
        The Production class instance used for crater production. Set during initialization.
        """
        return self._production

    @production.setter
    def production(self, value):
        if not isinstance(value, (Production | str)):
            raise TypeError("production must be a subclass of Production or str")
        self._production = value

    @property
    def scaling(self):
        """
        The Scaling object that defines the crater scaling relationships model. Set during initialization.
        """
        return self._scaling

    @scaling.setter
    def scaling(self, value):
        if not isinstance(value, (Scaling | str)):
            raise TypeError("scaling must be of Scaling type or str")
        self._scaling = value

    @property
    def morphology(self):
        """
        The crater morphology model. Set during initialization.
        """
        return self._morphology

    @morphology.setter
    def morphology(self, value):
        if not isinstance(value, (Morphology | str)):
            raise TypeError("morphology must be of Morphology type or str")
        self._morphology = value

    @property
    def projectile(self):
        """
        The crater projectile model. Set during initialization.
        """
        return self._projectile

    @projectile.setter
    def projectile(self, value):
        if not isinstance(value, (Projectile | str)):
            raise TypeError("projectile must be of Projectile type or str")
        self._projectile = value

    @property
    def counting(self):
        """
        The crater counting model. Set during initialization.
        """
        return self._counting

    @counting.setter
    def counting(self, value):
        if not isinstance(value, (Counting | str)):
            raise TypeError("counting must be of Counting type or str")
        self._counting = value

    @property
    def n_node(self):
        """
        Number of nodes in the simulation mesh. Dynamically set based on `surface` attribute.
        """
        return self.surface.uxgrid.n_node

    @property
    def n_face(self):
        """
        Number of faces in the simulation mesh. Dynamically set based on `surface` attribute.
        """
        return self.surface.uxgrid.n_face

    @parameter
    def interval(self):
        """
        The index of the current time step.
        """
        if self._interval is None:
            return 0
        return self._interval

    @interval.setter
    def interval(self, value):
        if not isinstance(value, int):
            raise TypeError("interval must be an integer")
        if value < 0:
            raise ValueError("interval must be greater than or equal to zero")

        self._interval = value

    @parameter
    def elapsed_time(self):
        """
        The elapsed time in My since the start of the simulation.
        """
        if self._elapsed_time is None:
            return 0.0
        return self._elapsed_time

    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = float(value)

    @parameter
    def time(self):
        """
        The age of the current time step in My relative to the present from the chronology of the production function.
        """
        if self._time is None:
            return -1.0
        return self._time

    @time.setter
    def time(self, value):
        self._time = float(value)

    @parameter
    def elapsed_n1(self):
        """
        The elapsed number of craters larger than 1 km in diameter.
        """
        if self._elapsed_n1 is None:
            return 0.0
        return self._elapsed_n1

    @elapsed_n1.setter
    def elapsed_n1(self, value):
        self._elapsed_n1 = float(value)

    @parameter
    def smallest_crater(self):
        """
        The smallest crater diameter in meters. Set during initialization.
        """
        return self._smallest_crater

    @smallest_crater.setter
    def smallest_crater(self, value):
        if value is None:
            self._smallest_crater = 0.0
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("smallest_crater must be a scalar value")
        elif value < 0:
            raise ValueError("smallest_crater must be greater than or equal to zero")
        elif self._largest_crater is not None and value > self._largest_crater:
            raise ValueError("smallest_crater must be less than or equal to largest_crater")
        self._smallest_crater = float(value)

    @parameter
    def largest_crater(self):
        """
        The largest crater diameter in meters. Set during initialization.
        """
        return self._largest_crater

    @largest_crater.setter
    def largest_crater(self, value):
        if value is None or np.isinf(value):
            self._largest_crater = np.inf
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("largest_crater must be a scalar value")
        elif value <= 0:
            raise ValueError("largest_crater must be greater than zero")
        elif self._smallest_crater is not None and value < self._smallest_crater:
            raise ValueError("largest_crater must be greater than or equal to smallest_crater")
        self._largest_crater = float(value)

    @parameter
    def smallest_projectile(self):
        """
        The smallest projectile diameter in meters. Set during initialization.
        """
        return self._smallest_projectile

    @smallest_projectile.setter
    def smallest_projectile(self, value):
        if value is None:
            self._smallest_projectile = 0.0
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("smallest_projectile must be a scalar value")
        elif value < 0:
            raise ValueError("smallest_projectile must be greater or equal to zero")
        elif self._largest_projectile is not None and value > self._largest_projectile:
            raise ValueError("smallest_projectile must be less than or equal to largest_projectile")
        self._smallest_projectile = float(value)

    @parameter
    def largest_projectile(self):
        """
        The largest projectile diameter in meters. Set during initialization.
        """
        return self._largest_projectile

    @largest_projectile.setter
    def largest_projectile(self, value):
        if value is None or np.isinf(value):
            self._largest_projectile = np.inf
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("largest_projectile must be a scalar value")
        elif value <= 0:
            raise ValueError("largest_projectile must be greater than zero")
        elif self._smallest_projectile is not None and value < self._smallest_projectile:
            raise ValueError("largest_projectile must be greater than or equal to smallest_projectile")
        self._largest_projectile = float(value)

    @property
    def name(self):
        """
        The name of the simulation.
        """
        return "Cratermaker Simulation object"

    @property
    def config_file(self):
        """
        The path to the configuration file for the simulation.
        """
        return self.simdir / _CONFIG_FILE_NAME

    @property
    def config_readonly(self) -> bool:
        """
        Flag indicating whether the configuration is read-only.
        """
        return self._config_readonly

    @parameter
    def ask_overwrite(self):
        """
        Flag to indicate whether the user should be prompted to overwrite any old files or not.
        """
        return self._ask_overwrite

    @ask_overwrite.setter
    def ask_overwrite(self, value):
        if not isinstance(value, bool):
            raise TypeError("ask_overwrite must be a bool")
        self._ask_overwrite = value

    @property
    def time_variables(self) -> dict[str, float]:
        """
        A dictionary of time-related variables for the simulation.

        Returns
        -------
        dict[str, float]
            A dictionary containing the following keys:
            - "time": The current age in My relative to the present.
            - "elapsed_time": The elapsed time in My since the start of the simulation.
            - "elapsed_n1": The elapsed number of craters larger than 1 km in diameter.
        """
        return {
            "time": self.time,
            "elapsed_time": self.elapsed_time,
            "elapsed_n1": self.elapsed_n1,
        }

    @property
    def do_counting(self) -> bool:
        """
        A boolean flag indicating whether or not counting is enabled for the simulation. This is determined by whether or not a counting model is present and has counting enabled.

        Returns
        -------
        bool
            True if counting is enabled, False otherwise.
        """
        return self.counting is not None and self.morphology.do_counting if self.morphology is not None else False

    @property
    def observed(self) -> dict[Crater] | None:
        """
        Pass-through to retrieve the current observed craters from the counting model, if it is enabled.
        """
        if self.do_counting:
            return self.counting.observed
        else:
            return None

    @property
    def emplaced(self) -> list[Crater] | None:
        """
        Pass-through to retrieve the current emplaced craters from the morphology model, if it is enabled.
        """
        if self.morphology is not None:
            return self.counting.emplaced
        else:
            return None

    @property
    def n_observed(self) -> int | None:
        """
        Pass-through to retrieve the current number of observed craters from the counting model, if it is enabled.
        """
        if self.do_counting:
            return self.counting.n_observed
        else:
            return None

    @property
    def n_emplaced(self) -> int | None:
        """
        Pass-through to retrieve the current number of emplaced craters from the counting model, if it is enabled.
        """
        if self.do_counting:
            return self.counting.n_emplaced
        else:
            return None
