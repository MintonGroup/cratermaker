from typing import Any
from warnings import warn

from numpy.random import Generator

from cratermaker.components.projectile import Projectile
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike


@Projectile.register("asteroids")
class AsteroidProjectiles(Projectile):
    _catalogue = {
        "Mercury": 41100.0,
        "Venus": 29100.0,
        "Earth": 24600.0,
        "Moon": 22100.0,
        "Mars": 10700.0,
        "Ceres": 5300.0,
        "Vesta": 5300.0,
        "MBA": 5300.0,
    }

    def __init__(
        self,
        monte_carlo: bool = True,
        projectile_density: FloatLike | None = None,
        projectile_angle: FloatLike | None = None,
        projectile_direction: FloatLike | None = None,
        projectile_location: tuple[FloatLike, FloatLike] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        An operations class for computing the projectile properties of an asteroid source population.

        Parameters
        ----------
        monte_carlo : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, the projectile projectile_velocity will be the mean projectile_velocity for the given target body.
        projectile_density : float, optional
            The projectile_density of the projectile in kg/m³. Default is 2250 kg/m³.
        projectile_angle : float, optional
            The impact projectile_angle in degrees. Default is 90.0 degrees (vertical impact) if `monte_carlo` is False. If `monte_carlo` is True, this value is ignored.
        projectile_direction : float | None
            The impact projectile_direction in degrees. Default is 0.0 degrees (due North) if `monte_carlo` is False. If `monte_carlo` is True, this value is ignored.`
        projectile_location : tuple[float, float] | None
            The impact projectile_location as a tuple of (longitude, latitude) in degrees. Default is (0.0, 0.0) (the equator and prime meridian) if `monte_carlo` is False. If `monte_carlo` is True, this value is ignored.
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
        projectile_density = 2250.0 if projectile_density is None else projectile_density
        kwargs.pop("projectile_mean_velocity", None)
        self._target = Target.maker(target, **kwargs)
        projectile_mean_velocity = self._set_mean_velocity()
        super().__init__(
            monte_carlo=monte_carlo,
            projectile_mean_velocity=projectile_mean_velocity,
            projectile_angle=projectile_angle,
            projectile_direction=projectile_direction,
            projectile_location=projectile_location,
            projectile_density=projectile_density,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += f"Target: {self.target.name}\n"
        return str_repr

    def _set_mean_velocity(self):
        """
        Sets the mean projectile_velocity of the projectile in m/s based on the target body.

        Returns
        -------
        float
        """
        if self.target.name in self.__class__._catalogue:
            pmv = float(self.__class__._catalogue[self.target.name])
        else:
            warn(
                f"Target {self.target.name} not found in known targets. Known targets include {list(self.__class__._catalogue.keys())}. Defaulting to the Moon.",
                stacklevel=2,
            )
            pmv = float(self.__class__._catalogue["Moon"])
        return pmv
