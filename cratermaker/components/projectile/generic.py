from typing import Any

from numpy.random import Generator

from cratermaker.components.projectile import Projectile
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike


@Projectile.register("generic")
class GenericProjectiles(Projectile):
    def __init__(
        self,
        sample: bool = True,
        mean_velocity: FloatLike | None = None,
        velocity: FloatLike | None = None,
        density: FloatLike | None = None,
        angle: FloatLike | None = None,
        direction: FloatLike | None = None,
        location: tuple[FloatLike, FloatLike] | None = None,
        target: Target | str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        An operations class for computing the projectile properties of a generic source population. The user is required to either set `mean_velocity` (if `sample==True`) or `velocity` (if `sample==False`).
        If `sample==True`, the impact velocities, angles, and directions are sampled from distributions. If `sample==False`, the impact velocities are set to the mean velocity, and `angle` and `direction` can also be set, but will default to set 90 degrees (vertical impact) and 0 degrees (due North`).

        Parameters
        ----------
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to True, the `mean_velocity` argument is required. If set to False, the `velocity` argument is required.
        mean_velocity : float, optional
            The mean velocity of the projectile in m/s. Required if `sample` is True, ignored if `sample` is False.
        velocity : float | None
            The impact velocity in m/s. If `sample` is True, this value is ignored. If `sample` is False, this value is required.
        density : float, optional
            The density of the projectile in kg/m^3. Default is 1000 kg/m^3.
        angle : float, optional
            The impact angle in degrees. Default is 90.0 degrees (vertical impact) if `sample` is False. If `sample` is True, this value is ignored.
        direction : float | None
            The impact direction in degrees. Default is 0.0 degrees (due North) if `sample` is False. If `sample` is True, this value is ignored.`
        location : tuple[float, float] | None
            The impact location as a tuple of (longitude, latitude) in degrees. Default is (0.0, 0.0) (the equator and prime meridian) if `sample` is False. If `sample` is True, this value is ignored.
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

        # set the defaults
        density = 1000.0 if density is None else density
        angle = 90.0 if angle is None else angle
        direction = 0.0 if direction is None else direction
        location = (0.0, 0.0) if location is None else location
        super().__init__(
            sample=sample,
            mean_velocity=mean_velocity,
            velocity=velocity,
            density=density,
            angle=angle,
            direction=direction,
            location=location,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )
        return
