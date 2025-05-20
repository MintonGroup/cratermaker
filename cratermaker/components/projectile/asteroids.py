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
        sample: bool = True,
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
        An operations class for computing the projectile properties of an asteroid source population.

        Parameters
        ----------
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, the projectile velocity will be the mean velocity for the given target body.
        density : float, optional
            The density of the projectile in kg/m^3. Default is 2250 kg/m^3.
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
        density = 2250.0 if density is None else density
        kwargs.pop("mean_velocity", None)
        kwargs.pop("velocity", None)
        self._target = Target.maker(target, **kwargs)
        mean_velocity = self._set_mean_velocity()
        super().__init__(
            sample=sample,
            mean_velocity=mean_velocity,
            velocity=mean_velocity,
            angle=angle,
            direction=direction,
            location=location,
            density=density,
            target=target,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nTarget: {self.target.name}\n"

    def _set_mean_velocity(self):
        """
        Sets the mean velocity of the projectile in m/s based on the target body.

        Returns
        -------
        float
        """

        if self.target.name in self.__class__._catalogue:
            pmv = float(self.__class__._catalogue[self.target.name])
        else:
            warn(
                f"Target {self.target.name} not found in known targets. Known targets include {list(self.__class__._catalogue.keys())}. Defaulting to the Moon."
            )
            pmv = float(self.__class__._catalogue["Moon"])
        return pmv
