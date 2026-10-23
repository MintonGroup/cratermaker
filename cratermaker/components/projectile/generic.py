from typing import Any

from numpy.random import Generator

from cratermaker.components.projectile import Projectile
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike


@Projectile.register("generic")
class GenericProjectiles(Projectile):
    def __init__(
        self,
        monte_carlo: bool = True,
        projectile_mean_velocity: FloatLike | None = None,
        projectile_velocity: FloatLike | None = None,
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
        An operations class for computing the projectile properties of a generic source population.

        The user is required to either set `projectile_mean_velocity` (if `monte_carlo==True`) or `projectile_velocity` (if `monte_carlo==False`).  If `monte_carlo==True`, the impact velocities, angles, and directions are sampled from distributions. If `monte_carlo==False`, the impact velocities are set to the mean projectile_velocity, and `projectile_angle` and `projectile_direction` can also be set, but will default to set 90 degrees (vertical impact) and 0 degrees (due North`).

        Parameters
        ----------
        monte_carlo : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to True, the `projectile_mean_velocity` argument is required. If set to False, the `projectile_velocity` argument is required.
        projectile_mean_velocity : float, optional
            The mean projectile_velocity of the projectile in m/s. Required if `monte_carlo` is True, ignored if `monte_carlo` is False.
        projectile_velocity : float | None
            The impact projectile_velocity in m/s. If `monte_carlo` is True, this value is ignored. If `monte_carlo` is False, this value is required.
        projectile_density : float, optional
            The projectile_density of the projectile in kg/m³. Default is 1000 kg/m³.
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
        # set the defaults
        projectile_density = 1000.0 if projectile_density is None else projectile_density
        projectile_angle = 90.0 if projectile_angle is None else projectile_angle
        projectile_direction = 0.0 if projectile_direction is None else projectile_direction
        projectile_location = (0.0, 0.0) if projectile_location is None else projectile_location
        if projectile_mean_velocity is None:
            if monte_carlo:
                projectile_mean_velocity = 20.0e3
            else:
                projectile_mean_velocity = projectile_velocity
        super().__init__(
            monte_carlo=monte_carlo,
            projectile_mean_velocity=projectile_mean_velocity,
            projectile_velocity=projectile_velocity,
            projectile_density=projectile_density,
            projectile_angle=projectile_angle,
            projectile_direction=projectile_direction,
            projectile_location=projectile_location,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )
        return
