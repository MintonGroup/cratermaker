import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from numpy.random import Generator
from numpy.typing import ArrayLike
from tqdm import tqdm

from ..components.morphology import Morphology
from ..components.production import Production
from ..components.projectile import Projectile
from ..components.scaling import Scaling
from ..components.surface import Surface
from ..components.target import Target
from ..constants import (
    _COMPONENT_NAMES,
    _CONFIG_FILE_NAME,
    _CRATER_DIR,
    FloatLike,
    PairOfFloats,
)
from ..utils import export
from ..utils.general_utils import _set_properties, format_large_units, parameter
from .base import CratermakerBase, _convert_for_yaml
from .crater import Crater


class Simulation(CratermakerBase):
    """
    This class is used to create a simulation of a crater population on a target body. It allows for the generation of craters based on a variety of parameters, including the target body, scaling laws, production functions, and morphology models.

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
    surface : str, optional
        The name of the surface used for the surface. Default is "icosphere".
    projectile : str, optional
        The projectile model to use from the components library, which is used to generate the projectile properties for the simulation, such as velocity and density. The default is "asteroids" when target is Mercury, Venus, Earth, Moon, Mars, Ceres, or Vesta, and "comets" otherwise.
    simdir : str | Path
        The main project simulation directory. Defaults to the current working directory if None.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    resume_old : bool, optional
        Flag to indicate whether to resume from an old simulation. If True, the simulation will attempt to load the previous state from the config file.
    **kwargs : Any
        Additional keyword arguments that can be passed to other cratermaker components, such as arguments to set the surface, scaling,
        morphology, or production function constructors. Refer to the documentation of each component module for details.
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
        simdir: str | Path | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        resume_old: bool = False,
        **kwargs: Any,
    ):
        super().__init__(
            simdir=simdir, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs
        )
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_scaling", scaling)
        object.__setattr__(self, "_production", production)
        object.__setattr__(self, "_morphology", morphology)
        object.__setattr__(self, "_projectile", projectile)
        object.__setattr__(self, "_surface", surface)
        object.__setattr__(self, "_craterlist", None)
        object.__setattr__(self, "_crater", None)
        object.__setattr__(self, "_interval_number", None)
        object.__setattr__(self, "_elapsed_time", None)
        object.__setattr__(self, "_current_age", None)
        object.__setattr__(self, "_elapsed_n1", None)
        object.__setattr__(self, "_smallest_crater", None)
        object.__setattr__(self, "_smallest_projectile", None)
        object.__setattr__(self, "_largest_crater", None)
        object.__setattr__(self, "_largest_projectile", None)
        object.__setattr__(self, "_true_crater_list", [])

        if self.config_file.exists():
            config_file = self.config_file
            if not resume_old:
                config_file = None

        else:
            config_file = None

        config_override = {}
        for component in _COMPONENT_NAMES:
            # Set to true if a local variable from the argument list with the component name is set to something other than None, otherwise false
            config_override[component] = getattr(self, f"_{component}") is not None

        _, unmatched = _set_properties(
            self,
            target=target,
            rng_seed=rng_seed,
            scaling=scaling,
            production=production,
            morphology=morphology,
            projectile=projectile,
            surface=surface,
            config_file=config_file,
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
        kwargs.update(unmatched)
        kwargs = {**kwargs, **vars(self.common_args)}
        if not resume_old:
            self.reset()

        target_config = {**target_config, **kwargs}
        self.target = Target.maker(self.target, **target_config)

        production_config = {**production_config, **kwargs}
        self.production = Production.maker(
            self.production, target=self.target, **production_config
        )

        projectile_config = {**projectile_config, **kwargs}
        self.projectile = Projectile.maker(
            self.projectile, target=self.target, **projectile_config
        )

        scaling_config = {**scaling_config, **kwargs}
        self.scaling = Scaling.maker(
            self.scaling,
            target=self.target,
            projectile=self.projectile,
            **scaling_config,
        )

        surface_config = {**surface_config, **kwargs}
        surface_config["reset"] = not resume_old
        self.surface = Surface.maker(
            self.surface,
            target=self.target,
            **surface_config,
        )

        morphology_config = {**morphology_config, **kwargs}
        self.morphology = Morphology.maker(
            self.morphology, surface=self.surface, **morphology_config
        )
        if self.surface.gridtype == "hireslocal":
            self.surface.set_superdomain(
                scaling=self.scaling, morphology=self.morphology, **surface_config
            )

        self._craterlist = []
        self._crater = None
        self._interval_number = 0
        self._elapsed_time = 0.0
        self._current_age = 0.0
        self._elapsed_n1 = 0.0
        self._smallest_crater = (
            0.0  # The smallest crater will be determined by the smallest face area
        )
        self._smallest_projectile = (
            0.0  # The smallest crater will be determined by the smallest face area
        )
        self._largest_crater = (
            np.inf
        )  # The largest crater will be determined by the target body radius
        self._largest_projectile = (
            np.inf
        )  # The largest projectile will be determined by the target body radius
        self.to_config()

        return

    def __str__(self) -> str:
        """
        Returns a string representation of the Simulation object.
        """
        return (
            f"<Simulation>\n\n"
            f"{self.target}\n\n"
            f"{self.projectile}\n\n"
            f"{self.scaling}\n\n"
            f"{self.production}\n\n"
            f"{self.morphology}\n\n"
            f"{self.surface}\n\n"
            f"<Current state>\n"
            f"Current age : {format_large_units(self.current_age, quantity='time')}\n"
            f"Elapsed time: {format_large_units(self.elapsed_time, quantity='time')}\n"
            f"Elapsed N_1 : {self.elapsed_n1} #/m^2\n"
            f"Interval    : {self.interval_number}\n"
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
        if hasattr(self, "to_config") and callable(getattr(self, "to_config", None)):
            if _convert_for_yaml(value) is not None:
                try:
                    self.to_config()
                except Exception:
                    pass

    def run(
        self,
        age: FloatLike | None = None,
        age_end: FloatLike | None = None,
        age_interval: FloatLike | None = None,
        diameter_number: PairOfFloats | None = None,
        diameter_number_end: PairOfFloats | None = None,
        diameter_number_interval: FloatLike | None = None,
        ninterval: int | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Run the simulation over a specified interval using the current production function.

        Parameters
        ----------
        age : FloatLike, optional
            Start age in My relative to the present for the simulation, used to compute the starting point of the production function.
            Default is None, which requires `diameter_number` to be set.
        age_end : FloatLike, optional
            End age in My relative to the present for the simulation, used to compute the ending point of the production function.
            Default is 0 (present day) if not provided.
        age_interval : FloatLike, optional
            Interval in My for outputting intermediate results. If not provided, calculated as `age - age_end` / `ninterval` if `ninterval` is provided, otherwise set to the total simulation duration.
        diameter_number : PairOfFloats, optional
            Starting cumulative number and diameter pair (D, N) to define the simulation start in terms of crater frequency.
            Default is None, which requires `age` to be set.
        diameter_number_end : PairOfFloats, optional
            Ending cumulative number and diameter pair (D, N) to define the simulation end in terms of crater frequency.
            Default is the present-day values if not provided.
        diameter_number_interval : PairOfFloats, optional
            Interval for outputting results in terms of cumulative number and diameter (D, N). Calculated based on `ninterval` if provided.
        ninterval : int, optional
            Number of intervals for outputting results. This has a special use case where one can specify age-based inputs but output
            in equal cumulative number intervals and vice versa.

        Notes
        -----
        This function allows defining the simulation parameters either in terms of age or crater frequency (cumulative number). The
        arguments `age`, `age_end`, `age_interval` are mutually exclusive with `diameter_number`, `diameter_number_end`, `diameter_number_interval`.

        The input diameter values used in diameter_number, diameter_number_end, and diameter_number_interval need not be the same.
        However, the production function will be used to convert all numbers to a common diameter value on output.

        The initial state of the simulation (before any craters are emplaced) is always saved automatically.
        As a result, the total number of saved states will be `ninterval + 1`, where `ninterval` is the number
        of simulation intervals requested.

        Examples
        --------
        .. code-block:: python

            # Create a simulation object with default parameters (Moon, NeukumProduction, etc.)
            sim = cratermaker.Simulation()

            # Run the simulation for 3.8 billion years, saving the results every 100 million years
            sim.run(age=3.8e3, age_interval=100.0)

            # Run the simulation to create 80 craters larger than 300 km in diameter and output the results after approximately 100
            # craters larger than 50km in diameter have formed
            sim.run(diameter_number=(300e3, 80), diameter_number_interval=(50e3, 100))

            # Run the simulation for 3.8 billion years, saving 100 intervals with equal cumulative number intervals
            sim.run(age=3.8e3, ninterval=100)

            # Run the simulation to create 80 craters larger than 300 km and output 100 intervals with equal age intervals
            sim.run(diameter_number=(300e3, 80), ninterval=100)

            # Run the simulation from 3.8 billion years to 3.0 billion years, saving the results every 100 million years
            sim.run(age=3.8e3, age_end=3.0e3, age_interval=100.0)

        """
        arguments = {
            "age": age,
            "age_end": age_end,
            "age_interval": age_interval,
            "diameter_number": diameter_number,
            "diameter_number_end": diameter_number_end,
            "diameter_number_interval": diameter_number_interval,
            "ninterval": ninterval,
            **kwargs,
        }
        arguments = self._validate_run_args(**arguments)
        age = arguments.pop("age", None)
        age_end = arguments.pop("age_end", None)
        age_interval = arguments.pop("age_interval", None)
        diameter_number = arguments.pop("diameter_number", None)
        diameter_number_end = arguments.pop("diameter_number_end", None)
        diameter_number_interval = arguments.pop("diameter_number_interval", None)
        ninterval = arguments.pop("ninterval", None)
        is_age_interval = arguments.pop("is_age_interval", None)
        if ninterval == 1:
            is_age_interval = True
            age_interval = age

        self.current_age = age
        self.elapsed_time = 0.0
        self.elapsed_n1 = 0.0
        self.save()
        for i in tqdm(
            range(ninterval),
            total=ninterval,
            desc="Simulation interval",
            unit="interval",
            position=1,
            leave=True,
        ):
            self.interval_number = i + 1
            if is_age_interval:
                current_age = age - i * age_interval
                current_age_end = age - (i + 1) * age_interval
                if current_age_end < 0.0:
                    current_age_end = 0.0
                self.populate(age=current_age, age_end=current_age_end)
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
                current_age = self.production.function_inverse(
                    *current_diameter_number_density
                )
                if current_diameter_number_end[1] > 0:
                    current_diameter_number_density_end = (
                        current_diameter_number_end[0],
                        current_diameter_number_end[1] / self.surface.area,
                    )

                    current_age_end = self.production.function_inverse(
                        *current_diameter_number_density_end
                    )
                else:
                    current_age_end = 0.0
                if current_age_end < 0.0:
                    current_age_end = 0.0
                age_interval = current_age - current_age_end
            self.elapsed_time += age_interval
            self.elapsed_n1 += self.production.function(
                diameter=1000.0, age=current_age
            ) - self.production.function(diameter=1000.0, age=current_age_end)
            self.current_age = current_age_end

            self.save()

        self.export("vtp")
        return

    def _validate_run_args(self, **kwargs) -> dict:
        """
        Validate all the input arguments to the sample method. This function will raise a ValueError if any of the arguments are invalid.
        It will also convert age arguments to diameter_number and vice versa.

        Parameters
        ----------
        age : FloatLike, optional
            Start age in My relative to the present for the simulation, used to compute the starting point of the production function.
            Default is None, which requires `diameter_number` to be set.
        age_end : FloatLike, optional
            End age in My relative to the present for the simulation, used to compute the ending point of the production function.
            Default is 0 (present day) if not provided.
        age_interval : FloatLike, optional
            Interval in My for outputting intermediate results. If not provided, calculated as `age - age_end` / `ninterval` if `ninterval` is provided, otherwise set to the total simulation duration.
        diameter_number : PairOfFloats, optional
            Starting cumulative number and diameter pair (D, N) to define the simulation start in terms of crater frequency.
            Default is None, which requires `age` to be set.
        diameter_number_end : PairOfFloats, optional
            Ending cumulative number and diameter pair (D, N) to define the simulation end in terms of crater frequency.
            Default is the present-day values if not provided.
        diameter_number_interval : PairOfFloats, optional
            Interval for outputting results in terms of cumulative number and diameter (D, N). Calculated based on `ninterval` if provided.
        ninterval : int, optional
            Number of intervals for outputting results. This has a special use case where one can specify age-based inputs but output
            in equal cumulative number intervals and vice versa.

        Returns
        -------
        A dict containing all arguments listed in Parameters above, as well as `is_age_interval`, which is a boolean flag indicating
        whether or not the simulation is being run in equal age intervals or equal number intervals.

        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - Neither the age nore the diameter_number argument is provided.
            - Both the age and diameter_number arguments are provided.
            - Both the age_end and diameter_number_end arguments are provided.
            - The age argument is provided but is not a scalar.
            - The age_end argument is provided but is not a scalar.
            - The age_interval is provided but is not a positive scalar
            - The age_interval provided is negative, or is greater than age - age_end
            - The diameter_number argument is not a pair of values, or any of them are less than 0
            - The diameter_number_end argument is not a pair of values, or any of them are less than 0
            - The diameter_number_interval argument is not a pair of values, or any of them are less than 0
            - The age_interval and diameter_number_interval arguments are both provided.
            - The diameter_number_interval provided is negative, or is greater than diameter_number - diameter_number_end
            - The ninterval is provided but is not an integer or is less than 1.
            - The ninterval is provided and either age_interval or diameter_number_interval is also provided
        """
        # Determine whether we are going to do equal time intervals or equal number intervals
        age = kwargs.get("age", None)
        age_interval = kwargs.get("age_interval", None)
        diameter_number_interval = kwargs.get("diameter_number_interval", None)
        ninterval = kwargs.pop("ninterval", None)
        if age_interval is not None and diameter_number_interval is not None:
            raise ValueError(
                "Cannot specify both ninterval and age_interval or diameter_number_interval"
            )
        if ninterval is not None:
            if not isinstance(ninterval, int):
                raise TypeError("ninterval must be an integer")
            if ninterval < 1:
                raise ValueError("ninterval must be greater than zero")
            if age_interval is not None or diameter_number_interval is not None:
                raise ValueError(
                    "Cannot specify both ninterval and age_interval or diameter_number_interval"
                )

        is_age_interval = age_interval is not None

        # Validate arguments using the production function validator first
        if "diameter_range" not in kwargs:
            kwargs["diameter_range"] = (
                self.get_smallest_diameter(),
                self.get_largest_diameter(),
            )
        if "area" not in kwargs:
            kwargs["area"] = self.surface.area.item()
        kwargs = self.production._validate_sample_args(**kwargs)

        if is_age_interval:
            age = kwargs.get("age", None)
            if age is None:
                raise ValueError(
                    "Something went wrong! age should be set by self.production_validate_sample_args"
                )
            age_end = kwargs.get("age_end", None)
            if age_end is None:
                raise ValueError(
                    "Something went wrong! age_end should be set by self.production_validate_sample_args"
                )
            if age_interval is None:
                if ninterval is None:
                    ninterval = 1
                age_interval = (age - kwargs["age_end"]) / ninterval
            else:
                if age_interval > age - age_end:
                    raise ValueError("age_interval must be less than age - age_end")
                elif age_interval <= 0:
                    raise ValueError("age_interval must be greater than zero")
                ninterval = int(np.ceil((age - age_end) / age_interval))

            kwargs["age_interval"] = age_interval
            kwargs["ninterval"] = ninterval
        else:
            diameter_number = kwargs.get("diameter_number", None)
            if diameter_number is None:
                raise ValueError(
                    "Something went wrong! diameter_number should be set by self.production_validate_sample_args"
                )
            diameter_number_end = kwargs.get("diameter_number_end", None)
            if diameter_number_end is None:
                raise ValueError(
                    "Something went wrong! diameter_number_end should be set by self.production_validate_sample_args"
                )
            if diameter_number_interval is None:
                if ninterval is None:
                    ninterval = 1
                diameter_number_interval = (
                    diameter_number[0],
                    (diameter_number[1] - diameter_number_end[1]) / ninterval,
                )
            else:
                if len(diameter_number_interval) != 2:
                    raise ValueError(
                        "The 'diameter_number_interval' must be a pair of values in the form (D,N)"
                    )
                # Check to be sure that the diameter in the diameter_number_interval is the same as diameter_number.
                # If not, we need to adjust the end diameter value it so that they match
                diameter_number_interval = self.production._validate_csfd(
                    *diameter_number_interval
                )
                if diameter_number_interval[0] != diameter_number[0]:
                    area = kwargs.get("area", None)
                    if area is None:
                        raise ValueError(
                            "Something went wrong! area should be set by self.production_validate_sample_args"
                        )
                    diameter_number_density_interval = (
                        diameter_number_interval[0],
                        diameter_number_interval[1] / area,
                    )
                    age_val = self.production.function_inverse(
                        *diameter_number_density_interval
                    )
                    diameter_number_density_interval = (
                        diameter_number[0],
                        self.production.function(
                            diameter=diameter_number[0], age=age_val
                        ),
                    )
                    diameter_number_interval = (
                        diameter_number[0],
                        diameter_number_density_interval[1] * area,
                    )

                if (
                    diameter_number_interval[1]
                    >= diameter_number[1] - diameter_number_end[1]
                ):
                    raise ValueError(
                        "diameter_number_interval must be less than diameter_number - diameter_number_end"
                    )
                if diameter_number_interval[1] <= 0:
                    raise ValueError(
                        "diameter_number_interval must be greater than zero"
                    )
                ninterval = int(
                    np.ceil(
                        (diameter_number[1] - diameter_number_end[1])
                        / diameter_number_interval[1]
                    )
                )

            kwargs["diameter_number_interval"] = diameter_number_interval
            kwargs["ninterval"] = ninterval

        kwargs["is_age_interval"] = is_age_interval

        # Remove unecessary arguments that came out of the production._validate_sample_args method
        kwargs.pop("diameter_range")
        kwargs.pop("area")
        kwargs.pop("return_age")

        return kwargs

    def populate(
        self,
        age: FloatLike | None = None,
        age_end: FloatLike | None = None,
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
        age_end, FloatLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
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
            raise RuntimeError(
                "The production function is not properly defined. Missing 'generator_type' attribute"
            )
        elif self.production.generator_type not in ["crater", "projectile"]:
            raise RuntimeError(
                f"Invalid production function type {self.production.generator_type}"
            )

        from_projectile = self.production.generator_type == "projectile"
        if from_projectile:
            diam_key = "projectile_diameter"
        else:
            diam_key = "final_diameter"
        Dmax = self.get_largest_diameter(from_projectile=from_projectile)
        Dmin = self.get_smallest_diameter(from_projectile=from_projectile)

        # Loop over each face in the mesh to build up a population of craters in this interval. This is done because faces may
        # not all have the same surface area, the range of crater sizes that can be formed on each face may be different.
        impact_diameters = []
        impact_ages = []
        impact_locations = []
        face_areas = self.surface.face_areas
        min_area = face_areas.min()
        surface_area = self.surface.area.item()

        # Group surfaces into bins based on their area. All bins within a factor of 2 in surface area are grouped together.
        max_bin_index = np.ceil(np.log2(face_areas.max() / min_area)).astype(int)
        bins = {i: [] for i in range(max_bin_index + 1)}

        for face_index, area in enumerate(face_areas):
            bin_index = np.floor(np.log2(area / min_area)).astype(int)
            bins[bin_index].append(face_index)

        # Process each bin
        for bin_index, face_indices in bins.items():
            if not face_indices:
                continue  # Skip empty bins
            face_indices = np.array(face_indices)

            bin_areas = face_areas[face_indices]
            total_bin_area = bin_areas.sum()
            area_ratio = total_bin_area / surface_area

            Dmin = self.get_smallest_diameter(
                bin_areas, from_projectile=from_projectile
            )
            if diameter_number is not None:
                diameter_number_local = (
                    diameter_number[0],
                    diameter_number[1] * area_ratio,
                )
            else:
                diameter_number_local = None
            if diameter_number_end is not None:
                diameter_number_end_local = (
                    diameter_number_end[0],
                    diameter_number_end[1] * area_ratio,
                )
            else:
                diameter_number_end_local = None

            diameters, ages = self.production.sample(
                age=age,
                age_end=age_end,
                diameter_number=diameter_number_local,
                diameter_number_end=diameter_number_end_local,
                diameter_range=(Dmin, Dmax),
                area=total_bin_area,
                **kwargs,
            )
            if diameters.size > 0:
                impact_diameters.extend(diameters.tolist())
                impact_ages.extend(ages.tolist())

                # Get the relative probability of impact onto any particular face then get the locations of the impacts
                p = bin_areas / total_bin_area
                face_indices = self.rng.choice(face_indices, size=diameters.shape, p=p)
                locations = self.surface.get_random_location_on_face(face_indices)
                impact_locations.extend(np.array(locations).T.tolist())
        if len(impact_diameters) > 0:
            craterlist = []
            # Sort the ages, diameters, and locations so that they are in order of decreasing age
            sort_indices = np.argsort(impact_ages)[::-1]
            impact_diameters = np.asarray(impact_diameters)[sort_indices]
            impact_ages = np.asarray(impact_ages)[sort_indices]
            impact_locations = np.array(impact_locations)[sort_indices]
            for diameter, location, age in zip(
                impact_diameters, impact_locations, impact_ages
            ):
                diam_arg = {diam_key: diameter}
                craterlist.append(
                    Crater.maker(
                        location=location,
                        age=age,
                        scaling=self.scaling,
                        **diam_arg,
                        **vars(self.common_args),
                        **kwargs,
                    )
                )
            self.emplace(craterlist)

        return

    def emplace(
        self, crater: Crater | list[Crater] | None = None, **kwargs: Any
    ) -> None:
        """
        Emplace one or more craters in the simulation.

        This method orchestrates the creation and placement of a crater in the
        simulation. It can create a crater directly or based on the characteristics
        of a projectile.

        Parameters
        ----------
        crater : Crater or list of Crater objects, optional
            The Crater object(s) to be emplaced. If provided, this will be used directly. Otherwise, a single will be generated based on the keyword arguments.
        **kwargs : Any
            Keyword arguments to pass to :class:`Crater.maker`.
            Refer to the documentation of this class for details on valid keyword arguments.

        Notes
        -----
        The keyword arguments provided are passed down to :meth:`Crater.maker`.
        Refer to its documentation for a detailed description of valid
        keyword arguments.

        Examples
        --------
        .. code-block:: python

            from cratermaker import Simulation, Crater
            sim = Simulation()

            # Create a crater with specific diameter
            sim.emplace(final_diameter=10.0e3)

            # Create a crater based on a projectile with given mass and projectile_velocity
            sim.emplace(projectile_mass=1e15, projectile_velocity=20e3)

            # Create a crater with a specific transient diameter and location
            sim.emplace(transient_diameter=50e3, location=(43.43, -86.92))

            # Create multiple craters
            craters = [Crater.maker(final_diameter=20.0e3), Crater.maker(final_diameter=20.0e3)]
            sim.emplace(craters)

        """
        if crater is None:
            crater_args = {**kwargs, **vars(self.common_args)}
            # Add scaling=self.scaling to the kwargs if it is not already present
            if "scaling" not in crater_args:
                crater_args["scaling"] = self.scaling
            crater = Crater.maker(**crater_args)
        elif isinstance(crater, list) and len(crater) > 0:
            self._true_crater_list.extend(crater)
            for c in crater:
                self._enqueue_crater(c)
            self._process_queue()
            return
        if isinstance(crater, Crater):
            self._true_crater_list.append(crater)
            self.morphology.emplace(crater, **kwargs)

        return

    def _enqueue_crater(self, crater: Crater | None = None, **kwargs) -> None:
        """
        Add a crater to the queue for later emplacement.

        Parameters
        ----------
        crater : Crater
            The crater object to enqueue.

        **kwargs : Any
            Additional keyword arguments for initializing the :class:`Crater`.

        Raises
        ------
        RuntimeError
            If the queue manager has not been initialized.
        """
        self.morphology._enqueue_crater(crater, **kwargs)
        return

    def _process_queue(self) -> None:
        """
        Process all queued craters in the order they were added, forming non-overlapping
        batches and applying each to the surface.


        Raises
        ------
        RuntimeError
            If the queue manager has not been initialized.
        """
        self.morphology._process_queue()
        return

    def save(self, **kwargs: Any) -> None:
        """
        Save the current simulation state to a file.
        """

        time_variables = {
            "current_age": self.current_age,
            "elapsed_time": self.elapsed_time,
            "elapsed_n1": self.elapsed_n1,
        }

        self.surface.save_to_files(
            interval_number=self.interval_number,
            time_variables=time_variables,
            **kwargs,
        )

        self.dump_crater_lists()

        self.to_config(**kwargs)

        return

    def dump_crater_lists(self) -> None:
        """
        Dump the crater lists to a file and reset the true crater list.
        """
        crater_dir = self.simdir / _CRATER_DIR
        crater_dir.mkdir(parents=True, exist_ok=True)
        truefilename = crater_dir / f"true_crater_list{self.interval_number:06d}.csv"

        # Convert current crater list to dicts, splitting location into longitude/latitude
        new_data = []
        for c in self._true_crater_list:
            d = asdict(c)
            if "location" in d:
                lon, lat = d.pop("location")
                d["longitude"] = lon
                d["latitude"] = lat
            new_data.append(d)

        # If the file already exists, read it and merge
        if truefilename.exists():
            import csv

            with truefilename.open("r", newline="") as f:
                reader = csv.DictReader(f)
                existing_data = list(reader)
            combined_data = existing_data + new_data
            # Sort by final_diameter descending
            combined_data = sorted(
                combined_data, key=lambda d: -float(d["final_diameter"])
            )
        else:
            combined_data = new_data

        # Write merged data back to file
        if combined_data:
            import csv

            with truefilename.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=combined_data[0].keys())
                writer.writeheader()
                writer.writerows(combined_data)

        self._true_crater_list = []
        return

    def to_config(self, save_to_file: bool = True, **kwargs: Any) -> dict:
        """
        Converts values to types that can be used in yaml.safe_dump. This will convert various types into a format that can be saved in a human-readable YAML file. This will consolidate all of the configuration
        parameters into a single dictionary that can be saved to a YAML file. This will also remove any common arguments from the individual configurations for each component model to avoid repeating them.

        Parameters
        ----------
        save_to_file : bool, optional
            If True, the configuration will be saved to a file. Default is True.
        **kwargs : Any
            Additional keyword arguments for subclasses.

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
        sim_config["target_config"] = self.target.to_config(remove_common_args=True)
        sim_config["scaling_config"] = self.scaling.to_config(remove_common_args=True)
        sim_config["production_config"] = self.production.to_config(
            remove_common_args=True
        )
        sim_config["surface_config"] = self.surface.to_config(remove_common_args=True)
        sim_config["projectile_config"] = self.projectile.to_config(
            remove_common_args=True
        )
        sim_config["morphology_config"] = self.morphology.to_config(
            remove_common_args=True
        )
        sim_config["target"] = self.target.name
        sim_config["scaling"] = self.scaling._component_name
        sim_config["production"] = self.production._component_name
        sim_config["projectile"] = self.projectile._component_name
        sim_config["morphology"] = self.morphology._component_name
        sim_config["surface"] = self.surface._component_name

        for config in [
            "target",
            "scaling",
            "production",
            "projectile",
            "morphology",
            "surface",
        ]:
            # drop any empty values or {} from either f"{config} or f"{config}_config" if when they are either None or empty
            if config in sim_config:
                if sim_config[config] is None or sim_config[config] == {}:
                    sim_config.pop(config)
            if f"{config}_config" in sim_config:
                if (
                    sim_config[f"{config}_config"] is None
                    or sim_config[f"{config}_config"] == {}
                ):
                    sim_config.pop(f"{config}_config")

        # Write the combined configuration to a YAML file
        if save_to_file:
            with open(self.config_file, "w") as f:
                yaml.safe_dump(sim_config, f, indent=4)

        return sim_config

    def export(self, format="vtp", *args, **kwargs) -> None:
        """
        Export the surface mesh to a file in the specified format. Currently only VTK is supported.
        """
        self.save()
        if format == "vtp" or format == "vtk":
            export.to_vtk(self.surface, *args, **kwargs)
        else:
            raise ValueError(f"Unsupported export format: {format}")

    def reset(self):
        crater_dir = self.simdir / _CRATER_DIR
        # delete crater_dir using Pathlib
        if crater_dir.exists():
            shutil.rmtree(crater_dir)

    def update_elevation(self, *args: Any, **kwargs: Any) -> None:
        """
        Set the elevation on the surface. Delegates to the Surface object.

        Parameters
        ----------
        *args: Variable length argument list to pass to self.surface.update_elevation.
        **kwargs: Arbitrary keyword arguments to pass to self.surface.update_elevation.
        """
        return self.surface.update_elevation(*args, **kwargs)

    def get_smallest_diameter(
        self, face_areas: ArrayLike | None = None, from_projectile: bool = False
    ) -> float:
        """
        Get the smallest possible crater or projectile be formed on the surface.
        """
        if face_areas is None:
            face_areas = self.surface.face_areas
        else:
            face_areas = np.asarray(face_areas)
        smallest_crater = np.sqrt(face_areas.min().item() / np.pi) * 2
        if from_projectile:
            crater = Crater.maker(
                final_diameter=smallest_crater,
                angle=90.0,
                projectile_velocity=self.scaling.projectile_mean_velocity * 10,
                scaling=self.scaling,
                **vars(self.common_args),
            )
            return crater.projectile_diameter
        else:
            return smallest_crater

    def get_largest_diameter(self, from_projectile: bool = False) -> float:
        """
        Get the largest possible crater or projectile that can be formed on the surface.
        """
        largest_crater = self.target.radius * 2
        if from_projectile:
            crater = Crater.maker(
                final_diameter=largest_crater,
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
        if not isinstance(value, (Target, str)):
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
        if not isinstance(value, (Surface, str)):
            raise TypeError("surface must be an instance of Surface or str")
        self._surface = value

    @property
    def production(self):
        """
        The Production class instance used for crater production. Set during initialization.
        """
        return self._production

    @production.setter
    def production(self, value):
        if not isinstance(value, (Production, str)):
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
        if not isinstance(value, (Scaling, str)):
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
        if not isinstance(value, (Morphology, str)):
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
        if not isinstance(value, (Projectile, str)):
            raise TypeError("projectile must be of Projectile type or str")
        self._projectile = value

    @property
    def data_dir(self):
        """
        Directory where the data files are stored. Dynamically set based on `surface` attribute.
        """
        return self.surface.data_dir

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

    @property
    def true_crater_list(self):
        """
        The list of craters that have been emplaced in the simulation.
        """
        return self._true_crater_list

    @parameter
    def interval_number(self):
        """
        The index of the current time step.
        """
        return self._interval_number

    @interval_number.setter
    def interval_number(self, value):
        if not isinstance(value, int):
            raise TypeError("interval_number must be an integer")
        if value < 0:
            raise ValueError("interval_number must be greater than or equal to zero")

        self._interval_number = value

    @parameter
    def elapsed_time(self):
        """
        The elapsed time in My since the start of the simulation.
        """
        return self._elapsed_time

    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = float(value)

    @parameter
    def current_age(self):
        """
        The age of the current time step in My relative to the present from the chronology of the production function.
        """
        return self._current_age

    @current_age.setter
    def current_age(self, value):
        self._current_age = float(value)

    @parameter
    def elapsed_n1(self):
        """
        The elapsed number of craters larger than 1 km in diameter.
        """
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
            raise ValueError(
                "smallest_crater must be less than or equal to largest_crater"
            )
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
            raise ValueError(
                "largest_crater must be greater than or equal to smallest_crater"
            )
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
            raise ValueError(
                "smallest_projectile must be less than or equal to largest_projectile"
            )
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
        elif (
            self._smallest_projectile is not None and value < self._smallest_projectile
        ):
            raise ValueError(
                "largest_projectile must be greater than or equal to smallest_projectile"
            )
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
