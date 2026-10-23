from __future__ import annotations

import sys
from typing import (
    TYPE_CHECKING,
    Any,
)

try:
    if sys.version_info >= (3, 11):
        from typing import Self, TypeAlias
    else:
        from typing import TypeAlias

        from typing_extensions import Self
except ImportError:
    if TYPE_CHECKING:
        raise
    else:
        Self: Any = None

import math
from typing import TYPE_CHECKING

import numpy as np
from numpy.random import Generator

from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.utils.general_utils import (
    format_large_units,
    parameter,
    validate_and_normalize_location,
)

if TYPE_CHECKING:
    from cratermaker.components.target import Target


class Projectile(ComponentBase):
    """
    An abstract base class for all projectile models. It defines the interface for generating projectile velocities, angles, and densities for a given target body.
    """

    _registry: dict[str, Projectile] = {}
    _catalogue = None

    def __init__(
        self,
        monte_carlo: bool = True,
        projectile_mean_velocity: FloatLike | None = None,
        projectile_velocity: FloatLike | None = None,
        projectile_density: FloatLike | None = None,
        projectile_angle: FloatLike | None = None,
        projectile_direction: FloatLike | None = None,
        projectile_location: tuple[float, float] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs,
    ):
        """

        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        monte_carlo : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to True, the `projectile_mean_velocity` argument is required. If set to False, the `projectile_velocity` argument is required.
        projectile_mean_velocity : float, optional
            The mean projectile_velocity of the projectile in m/s. Required if `monte_carlo` is True, ignored if `monte_carlo` is False.
        projectile_velocity : float | None
            The impact projectile_velocity in m/s. If `monte_carlo` is True, this value is ignored. If `monte_carlo` is False, this value is required.
        projectile_density : float, optional
            The projectile_density of the projectile in kg/m³.
        projectile_angle : float, optional
            The impact projectile_angle in degrees. Default is 90.0 degrees (vertical impact) if `monte_carlo` is False. If `monte_carlo` is True, this value is ignored.
        projectile_direction : float | None
            The impact projectile_direction in degrees. Default is 0.0 degrees (due North) if `monte_carlo` is False. If `monte_carlo` is True, this value is ignored.
        projectile_location : tuple[float, float] | None
            The projectile_location of the projectile on the target body in (lon, lat) coordinates. If None, the projectile_location will be sampled from a distribution.
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.components.target import Target

        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

        object.__setattr__(self, "_monte_carlo", monte_carlo)
        object.__setattr__(self, "_mean_velocity", projectile_mean_velocity)
        object.__setattr__(self, "_velocity", projectile_velocity)
        object.__setattr__(self, "_density", projectile_density)
        object.__setattr__(self, "_angle", projectile_angle)
        object.__setattr__(self, "_direction", projectile_direction)
        object.__setattr__(self, "_location", projectile_location)
        object.__setattr__(self, "_target", Target.maker(target, **kwargs))

        if self.monte_carlo:
            if self.projectile_mean_velocity is None:
                raise ValueError("projectile_mean_velocity must be provided when monte_carlo is True")
        else:
            if self.projectile_velocity is None:
                raise ValueError("projectile_velocity must be provided when monte_carlo is False")
            if self.projectile_angle is None:
                raise ValueError("projectile_angle must be provided when monte_carlo is False")
            if self.projectile_direction is None:
                raise ValueError("projectile_direction must be provided when monte_carlo is False")
            if self.projectile_location is None:
                raise ValueError("projectile_location must be provided when monte_carlo is False")

        if self.projectile_density is None:
            raise ValueError("projectile_density must be provided")
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        if self.monte_carlo:
            str_repr += f"Sample from distributions: {self.monte_carlo}\n"
            str_repr += f"Mean Velocity: {format_large_units(self.projectile_mean_velocity, quantity='velocity')}\n"
        else:
            str_repr += f"Velocity: {format_large_units(self.projectile_velocity, quantity='velocity')}\n"
            str_repr += f"Angle: {self.projectile_angle:.1f}°\n"
            str_repr += f"Direction: {self.projectile_direction:.1f}°\n"
        str_repr += f"Density: {self.projectile_density:.1f} kg/m³\n"
        return str_repr

    def _copy(self, deep: bool = True, memo: dict[int, Any] | None = None) -> Self:
        import copy
        import inspect

        copier = copy.deepcopy if deep else copy.copy
        memo = {} if memo is None else memo

        # Get all init parameters except 'self'
        cls = self.__class__
        sig = inspect.signature(cls.__init__)
        init_keys = sig.parameters.keys() - {"self"}

        # Extract values for those keys
        init_kwargs = {}
        for k in init_keys:
            if k == "target":
                init_kwargs[k] = self.target
            elif k == "rng":
                init_kwargs[k] = self.rng
            else:
                attr = getattr(self, k, None)
                init_kwargs[k] = copier(attr, memo) if deep else copier(attr)

        return cls(**init_kwargs)

    def __copy__(self) -> Self:
        return self._copy(deep=False)

    def __deepcopy__(self, memo: dict[int, Any] | None = None) -> Self:
        return self._copy(deep=True, memo=memo)

    @classmethod
    def maker(
        cls,
        projectile: Projectile | str | None = None,
        projectile_mean_velocity: FloatLike | None = None,
        projectile_density: FloatLike | None = None,
        monte_carlo: bool = True,
        projectile_angle: FloatLike | None = None,
        projectile_velocity: FloatLike | None = None,
        projectile_direction: FloatLike | None = None,
        projectile_location: tuple[float, float] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ) -> Projectile:
        """
        Initialize a Projectile model with the given name or instance.

        Parameters
        ----------
        projectile : Projectile or str
            The projectile model to initialize. Can be a class or a string representing the model name.
        projectile_mean_velocity : float
            The mean projectile_velocity of the projectile in m/s.
        projectile_density : float
            The projectile_density of the projectile in kg/m³.
        monte_carlo : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean projectile_velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.
        projectile_angle : float
            The impact projectile_angle in degrees.
        projectile_velocity : float | None
            The impact projectile_velocity in m/s. If None, the projectile_velocity will be sampled from a distribution.
        projectile_direction : float | None
            The impact projectile_direction in degrees. If None, the projectile_direction will be sampled from a distribution.
        projectile_location : tuple[float, float]
            The projectile_location of the projectile on the target body in (lon, lat) coordinates. If None, the projectile_location will be sampled from a distribution
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
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
        Projectile
            The initialized projectile model.

        Raises
        ------
        KeyError
            If the specified projectile model is not found in the registry.
        TypeError
            If the specified projectile component type is not a string or a subclass of Projectile.
        """
        from cratermaker.components.projectile.asteroids import AsteroidProjectiles
        from cratermaker.components.projectile.comets import CometProjectiles
        from cratermaker.components.target import Target

        target = Target.maker(target, **kwargs)
        if projectile is None:
            if target.name in AsteroidProjectiles._catalogue:
                projectile = "asteroids"
            elif target.name in CometProjectiles._catalogue:
                projectile = "comets"
            else:
                projectile = "generic"

        # if this is a brand new uninstantiated projectile, we need to flag it so that it its properties can be set propertly. Otherwise, it should just pass through as is.
        isfresh = (
            isinstance(projectile, str)
            or projectile_velocity is not None
            or projectile_angle is not None
            or projectile_direction is not None
            or projectile_location is not None
            or projectile_density is not None
        )

        projectile = super().maker(
            component=projectile,
            projectile_mean_velocity=projectile_mean_velocity,
            projectile_density=projectile_density,
            monte_carlo=monte_carlo,
            projectile_angle=projectile_angle,
            projectile_velocity=projectile_velocity,
            projectile_direction=projectile_direction,
            projectile_location=projectile_location,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )
        if isfresh:
            return projectile.new_projectile(
                projectile_velocity=projectile_velocity,
                projectile_angle=projectile_angle,
                projectile_direction=projectile_direction,
                projectile_location=projectile_location,
                projectile_density=projectile_density,
                **kwargs,
            )
        else:
            return projectile

    def new_projectile(
        self,
        projectile_velocity: FloatLike | None = None,
        projectile_angle: FloatLike | None = None,
        projectile_direction: FloatLike | None = None,
        projectile_location: tuple[float, float] | None = None,
        projectile_density: FloatLike | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Returns a new projectile instance with updated sampled or default values, based on the original instance.

        Parameters
        ----------
        projectile_velocity : float | None
            The impact projectile_velocity in m/s. If None, the projectile_velocity will be sampled from a distribution if `monte_carlo` is True, otherwise it will be set to the value of `projectile_mean_velocity`.
        projectile_angle : float, optional
            The impact projectile_angle in degrees. Default is 90.0 degrees (vertical impact) if `monte_carlo` is False.
        projectile_direction : float | None
            The impact projectile_direction in degrees. Default is 0.0 degrees (due North) if `monte_carlo` is False.
        projectile_location : tuple[float, float] | None
            The projectile_location of the projectile on the target body in (lon, lat) coordinates. If None, the projectile_location will be sampled from a distribution.
        projectile_density : float | None
            The projectile_density of the projectile in kg/m³. If None, the projectile_density will be the default for this Projectile type.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        Projectile
            A new Projectile instance with updated properties.
        """
        if projectile_velocity is None:
            if self.monte_carlo:
                projectile_velocity = float(
                    mc.get_random_velocity(
                        vmean=self.projectile_mean_velocity,
                        vescape=self.target.escape_velocity,
                        rng=self.rng,
                    )[0]
                )
            else:
                projectile_velocity = self.projectile_mean_velocity

        if projectile_angle is None:
            if self.monte_carlo:
                projectile_angle = float(mc.get_random_impact_angle(rng=self.rng)[0])
            elif self.projectile_angle is not None:
                projectile_angle = self.projectile_angle
            else:
                projectile_angle = 90.0

        if projectile_direction is None:
            if self.monte_carlo:
                projectile_direction = float(mc.get_random_impact_direction(rng=self.rng)[0])
            elif self.projectile_direction is not None:
                projectile_direction = self.projectile_direction
            else:
                projectile_direction = 0.0

        if projectile_location is None:
            if self.monte_carlo:
                projectile_location = mc.get_random_location(rng=self.rng)[0]
            elif self.projectile_location is not None:
                projectile_location = self.projectile_location
            else:
                projectile_location = (0, 0)

        if projectile_density is None:
            projectile_density = self.projectile_density

        return self.__class__(
            projectile_velocity=projectile_velocity,
            projectile_angle=projectile_angle,
            projectile_direction=projectile_direction,
            projectile_location=projectile_location,
            projectile_density=projectile_density,
            projectile_mean_velocity=self.projectile_mean_velocity,
            monte_carlo=self.monte_carlo,
            target=self.target,
            rng=self.rng,
            **kwargs,
        )

    @parameter
    def monte_carlo(self):
        """
        Flag that determines whether to sample velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean projectile_velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.

        Returns
        -------
        bool
        """
        return self._monte_carlo

    @monte_carlo.setter
    def monte_carlo(self, value):
        if not isinstance(value, bool):
            raise TypeError("monte_carlo must be a boolean value")
        self._monte_carlo = value
        return

    @parameter
    def projectile_mean_velocity(self):
        """
        The mean projectile_velocity of the projectile in m/s.

        Returns
        -------
        float
        """
        return self._mean_velocity

    @projectile_mean_velocity.setter
    def projectile_mean_velocity(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if isinstance(value, FloatLike):
            if value < 0:
                raise ValueError("projectile_mean_velocity must be a positive number")
            self._mean_velocity = float(value)
        else:
            raise TypeError("projectile_mean_velocity must be a numeric value")
        return

    @property
    def mean_velocity(self):
        """
        The mean projectile_velocity of the projectile in m/s.

        Returns
        -------
        float
        """
        return self._mean_velocity

    @mean_velocity.setter
    def mean_velocity(self, value):
        self.projectile_mean_velocity = value

    @property
    def projectile_angle(self):
        """The impact projectile_angle relative to horizontal in degrees."""
        return self._angle

    @projectile_angle.setter
    def projectile_angle(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("projectile_angle must be a numeric value")
        if value < 0:
            raise ValueError("projectile_angle must be a positive number")
        self._angle = float(value)

    @property
    def angle(self):
        """The impact projectile_angle relative to horizontal in degrees."""
        return self._angle

    @angle.setter
    def angle(self, value):
        self.projectile_angle = value

    @property
    def projectile_direction(self):
        """The impact projectile_direction measured clockwise from North in degrees."""
        return self._direction

    @projectile_direction.setter
    def projectile_direction(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("projectile_direction must be a numeric value")
        if value < 0:
            raise ValueError("projectile_direction must be a positive number")
        self._direction = float(value)

    @property
    def direction(self):
        """The impact projectile_direction measured clockwise from North in degrees."""
        return self._direction

    @direction.setter
    def direction(self, value):
        self.projectile_direction = value

    @property
    def projectile_location(self) -> tuple[float, float]:
        """The projectile_location of the projectile (longitude, latitude) in degrees on the target body."""
        return self._location

    @projectile_location.setter
    def projectile_location(self, value: tuple[float, float]):
        self._location = validate_and_normalize_location(value)

    @property
    def location(self) -> tuple[float, float]:
        """The projectile_location of the projectile (longitude, latitude) in degrees on the target body."""
        return self._location

    @location.setter
    def location(self, value):
        self.projectile_location = value

    @property
    def projectile_velocity(self):
        """The impact projectile_velocity in m/s."""
        return self._velocity

    @projectile_velocity.setter
    def projectile_velocity(self, value):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("projectile_velocity must be a numeric value")
            if value < 0:
                raise ValueError("projectile_velocity must be a positive number")
            self._velocity = float(value)
        else:
            self._velocity = None
        return

    @property
    def velocity(self):
        """The impact projectile_velocity in m/s."""
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        self.projectile_velocity = value

    @parameter
    def projectile_density(self):
        """The bulk projectile_density of the projectile in kg/m³."""
        return self._density

    @projectile_density.setter
    def projectile_density(self, value):
        if not isinstance(value, FloatLike):
            raise TypeError("projectile_density must be a numeric value")
        if value < 0:
            raise ValueError("projectile_density must be a positive number")
        self._density = float(value)

    @property
    def density(self):
        """The bulk projectile_density of the projectile in kg/m³."""
        return self._density

    @density.setter
    def density(self, value):
        self.projectile_density = value

    @property
    def vertical_velocity(self):
        """The vertical component of the impact projectile_velocity in m/s."""
        return self.projectile_velocity * math.sin(math.radians(self.projectile_angle))

    @property
    def projectile_vertical_velocity(self):
        """The vertical component of the impact projectile_velocity in m/s."""
        return self.vertical_velocity

    @property
    def population(self):
        """The name of the population of the projectile model."""
        return self._component_name

    @property
    def target(self):
        """The Target object associated with theprojectile model."""
        return self._target

    @target.setter
    def target(self, value):
        from cratermaker.components.target import Target

        self._target = Target.maker(value)

    @property
    def catalogue(self) -> str:
        """A string representation of the projectile catalogue for the target body."""
        from cratermaker.utils.general_utils import format_large_units

        if self.__class__._catalogue is None:
            return "This Projectile component does not have a catalogue."
        lines = []
        header1 = "Target"
        lines.append(f"\n{header1:<11}| Mean Velocity")
        lines.append(26 * "-")
        for name, projectile_velocity in self.__class__._catalogue.items():
            formatted_velocity = format_large_units(projectile_velocity, quantity="velocity")
            lines.append(f"{name:<11}|     {formatted_velocity}")
        return "\n".join(lines)


import_components(__name__, __path__)
