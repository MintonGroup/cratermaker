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
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import (
    format_large_units,
    parameter,
    validate_and_normalize_location,
)

if TYPE_CHECKING:
    from cratermaker.components.target import Target


class Projectile(ComponentBase):
    _registry: dict[str, Projectile] = {}
    _catalogue = None

    def __init__(
        self,
        sample: bool = True,
        mean_velocity: FloatLike | None = None,
        velocity: FloatLike | None = None,
        density: FloatLike | None = None,
        angle: FloatLike | None = None,
        direction: FloatLike | None = None,
        location: tuple[float, float] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs,
    ):
        """
        This is the abstract base class for all projectile models. It defines the interface for generating projectile velocities, angles, and densities for a given target body.

        Parameters
        ----------
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to True, the `mean_velocity` argument is required. If set to False, the `velocity` argument is required.
        mean_velocity : float, optional
            The mean velocity of the projectile in m/s. Required if `sample` is True, ignored if `sample` is False.
        velocity : float | None
            The impact velocity in m/s. If `sample` is True, this value is ignored. If `sample` is False, this value is required.
        density : float, optional
            The density of the projectile in kg/m^3.
        angle : float, optional
            The impact angle in degrees. Default is 90.0 degrees (vertical impact) if `sample` is False. If `sample` is True, this value is ignored.
        direction : float | None
            The impact direction in degrees. Default is 0.0 degrees (due North) if `sample` is False. If `sample` is True, this value is ignored.`
        location : tuple[float, float] | None
            The location of the projectile on the target body in (lon, lat) coordinates. If None, the location will be sampled from a distribution.
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """
        from cratermaker.components.target import Target

        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

        object.__setattr__(self, "_sample", sample)
        object.__setattr__(self, "_mean_velocity", mean_velocity)
        object.__setattr__(self, "_velocity", velocity)
        object.__setattr__(self, "_density", density)
        object.__setattr__(self, "_angle", angle)
        object.__setattr__(self, "_direction", direction)
        object.__setattr__(self, "_location", location)
        object.__setattr__(self, "_target", Target.maker(target, **kwargs))

        if self.sample:
            if self.mean_velocity is None:
                raise ValueError("mean_velocity must be provided when sample is True")
        else:
            if self.velocity is None:
                raise ValueError("velocity must be provided when sample is False")
            if self.angle is None:
                raise ValueError("angle must be provided when sample is False")
            if self.direction is None:
                raise ValueError("direction must be provided when sample is False")
            if self.location is None:
                raise ValueError("location must be provided when sample is False")

        if self.density is None:
            raise ValueError("density must be provided")
        return

    def __str__(self) -> str:
        base = super().__str__()
        if self.sample:
            params = f"\nMean Velocity: {format_large_units(self.mean_velocity, quantity='velocity')}"
        else:
            params = (
                f"\nVelocity: {format_large_units(self.velocity, quantity='velocity')}"
            )
            params += f"\nAngle: {self.angle:.1f} degrees"
            params += f"\nDirection: {self.direction:.1f} degrees"
        return (
            f"{base}\n"
            f"Sample from distributions: {self.sample}\n"
            f"{params}\n"
            f"Density: {self.density:.1f} kg/m³\n"
        )

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
        mean_velocity: FloatLike | None = None,
        density: FloatLike | None = None,
        sample: bool = True,
        angle: FloatLike | None = None,
        velocity: FloatLike | None = None,
        direction: FloatLike | None = None,
        location: tuple[float, float] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ) -> Projectile:
        """
        Initialize an projectile model based on the provided name or class.

        Parameters
        ----------
        projectile : Projectile or str
            The projectile model to initialize. Can be a class or a string representing the model name.
        mean_velocity : float
            The mean velocity of the projectile in m/s.
        density : float
            The density of the projectile in kg/m^3.
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.
        angle : float
            The impact angle in degrees.
        velocity : float | None
            The impact velocity in m/s. If None, the velocity will be sampled from a distribution.
        direction : float | None
            The impact direction in degrees. If None, the direction will be sampled from a distribution.
        location : tuple[float, float]
            The location of the projectile on the target body in (lon, lat) coordinates. If None, the location willb e sampled from a distribution
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
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
        Projectile
            The initialized projectile model.

        Raises
        ------
        KeyError
            If the specified projectile model is not found in the registry.
        TypeError
            If the specified projectile model is not a string or a subclass of Projectile.
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
                mean_velocity = 20.0e3 if mean_velocity is None else mean_velocity

        # if this is a brand new uninstantiated projectile, we need to flag it so that it its properties can be set propertly. Otherwise, it should just pass through as is.
        isfresh = isinstance(projectile, str)

        projectile = super().maker(
            component=projectile,
            mean_velocity=mean_velocity,
            density=density,
            sample=sample,
            angle=angle,
            velocity=velocity,
            direction=direction,
            location=location,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )
        if isfresh:
            return projectile.new_projectile(
                velocity=velocity,
                angle=angle,
                direction=direction,
                location=location,
                **kwargs,
            )
        else:
            return projectile

    def new_projectile(
        self,
        velocity: FloatLike | None = None,
        angle: FloatLike | None = None,
        direction: FloatLike | None = None,
        location: tuple[float, float] | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Returns a new projectile instance with updated sampled or default values,
        based on the original instance.

        Parameters
        ----------
        velocity : float | None
            The impact velocity in m/s. If None, the velocity will be sampled from a distribution if `sample` is True, otherwise it will be set to the value of `mean_velocity`.
        angle : float, optional
            The impact angle in degrees. Default is 90.0 degrees (vertical impact) if `sample` is False.
        direction : float | None
            The impact direction in degrees. Default is 0.0 degrees (due North) if `sample` is False.
        location : tuple[float, float] | None
            The location of the projectile on the target body in (lon, lat) coordinates. If None, the location will be sampled from a distribution.
        **kwargs : Any
            Additional keyword arguments to override attributes.

        Returns
        -------
        Projectile
            A new Projectile instance with updated properties.
        """
        import copy

        new_obj = copy.copy(self)

        if velocity is not None:
            new_obj.velocity = velocity
        elif new_obj.sample:
            new_obj._velocity = float(
                mc.get_random_velocity(
                    vmean=new_obj.mean_velocity,
                    vescape=new_obj.target.escape_velocity,
                    rng=new_obj.rng,
                )[0]
            )
        elif new_obj._velocity is None:
            new_obj._velocity = new_obj.mean_velocity

        if angle is not None:
            new_obj.angle = angle
        elif new_obj.sample:
            new_obj._angle = float(mc.get_random_impact_angle(rng=new_obj.rng)[0])
        elif new_obj._angle is None:
            new_obj._angle = 90.0

        if direction is not None:
            new_obj.direction = direction
        elif new_obj.sample:
            new_obj._direction = float(
                mc.get_random_impact_direction(rng=new_obj.rng)[0]
            )
        elif new_obj._direction is None:
            new_obj._direction = 0.0

        if location is not None:
            new_obj.location = location
        elif new_obj.sample:
            new_obj._location = mc.get_random_location(rng=new_obj.rng)[0]
        elif new_obj._location is None:
            new_obj._location = (0, 0)

        return new_obj

    @parameter
    def sample(self):
        """
        Flag that determines whether to sample velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.

        Returns
        -------
        bool
        """
        return self._sample

    @sample.setter
    def sample(self, value):
        if not isinstance(value, bool):
            raise TypeError("sample must be a boolean value")
        self._sample = value
        return

    @parameter
    def mean_velocity(self):
        """
        The mean velocity of the projectile in m/s.

        Returns
        -------
        float
        """
        return self._mean_velocity

    @mean_velocity.setter
    def mean_velocity(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if isinstance(value, FloatLike):
            if value < 0:
                raise ValueError("mean_velocity must be a positive number")
            self._mean_velocity = float(value)
        else:
            raise TypeError("mean_velocity must be a numeric value")
        return

    @property
    def angle(self):
        """
        The impact angle in degrees.

        Returns
        -------
        float
        """
        return self._angle

    @angle.setter
    def angle(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("angle must be a numeric value")
        if value < 0:
            raise ValueError("angle must be a positive number")
        self._angle = float(value)

    @property
    def direction(self):
        """
        The impact direction in degrees.

        Returns
        -------
        float
        """
        return self._direction

    @direction.setter
    def direction(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("direction must be a numeric value")
        if value < 0:
            raise ValueError("direction must be a positive number")
        self._direction = float(value)

    @property
    def location(self) -> tuple[float, float]:
        """
        The location of the projectile on the target body.

        Returns
        -------
        tuple
            A tuple containing the (lon,lat) coordinates of the projectile in degrees.
        """
        return self._location

    @location.setter
    def location(self, value: tuple[float, float]):
        self._location = validate_and_normalize_location(value)

    @property
    def velocity(self):
        """
        The impact velocity in m/s.

        Returns
        -------
        float
        """
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("velocity must be a numeric value")
            if value < 0:
                raise ValueError("velocity must be a positive number")
            self._velocity = float(value)
        else:
            self._velocity = None
        return

    @parameter
    def density(self):
        """
        The density of the projectile in kg/m^3.

        Returns
        -------
        float
        """
        return self._density

    @density.setter
    def density(self, value):
        """
        Sets the density of the projectile in kg/m^3.

        Parameters
        ----------
        value : float | None
            The density in kg/m^3 or None to use a default value.

        Raises
        ------
        ValueError
            If the provided value is negative or None.
        TypeError
            If the provided value is not a numeric value or None.
        """
        if not isinstance(value, FloatLike):
            raise TypeError("density must be a numeric value")
        if value < 0:
            raise ValueError("density must be a positive number")
        self._density = float(value)

    @property
    def vertical_velocity(self):
        """
        The vertical component of the impact velocity in m/s.

        Returns
        -------
        float
        """
        return self.velocity * math.sin(math.radians(self.angle))

    @property
    def population(self):
        """
        The name of the population of the projectile model.

        Returns
        -------
        int
        """
        return self._component_name

    @property
    def target(self):
        """
        The target object for the projectile model.

        Returns
        -------
        Target
        """
        return self._target

    @target.setter
    def target(self, value):
        from cratermaker.components.target import Target

        self._target = Target.maker(value)

    @property
    def catalogue(self) -> str:
        from cratermaker.utils.general_utils import format_large_units

        if self.__class__._catalogue is None:
            return "This Projectile component does not have a catalogue."
        lines = []
        header1 = "Target"
        lines.append(f"\n{header1:<11}| Mean Velocity")
        lines.append(26 * "-")
        for name, velocity in self.__class__._catalogue.items():
            formatted_velocity = format_large_units(velocity, quantity="velocity")
            lines.append(f"{name:<11}|     {formatted_velocity}")
        return "\n".join(lines)


import_components(__name__, __path__, ignore_private=True)
