from typing import Any
from warnings import warn

from numpy.random import Generator

from cratermaker.components.projectile import Projectile
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import parameter


@Projectile.register("asteroids")
class AsteroidProjectiles(Projectile):
    def __init__(
        self,
        target: Target | str | None = None,
        sample: bool = True,
        density: FloatLike | None = None,
        angle: FloatLike | None = None,
        direction: FloatLike | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        An operations class for computing the projectile properties of an asteroid source population.

        Parameters
        ----------
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, the projectile velocity will be the mean velocity for the given target body.
        density : float, optional
            The density of the projectile in kg/m^3. Default is 2250 kg/m^3.
        angle : float, optional
            The impact angle in degrees. Default is 90.0 degrees (vertical impact) if `sample` is False. If `sample` is True, this value is ignored.
        direction : float | None
            The impact direction in degrees. Default is 0.0 degrees (due North) if `sample` is False. If `sample` is True, this value is ignored.`
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """
        self._target = Target.maker(target, **kwargs)
        if density is None:
            density = 2250.0

        mean_velocity = self._set_mean_velocity()
        super().__init__(
            sample=sample,
            mean_velocity=mean_velocity,
            velocity=mean_velocity,
            angle=angle,
            direction=direction,
            density=density,
            rng=rng,
            rng_seed=rng_seed,
            rng_state=rng_state,
            **kwargs,
        )

    def __repr__(self) -> str:
        base = super().__repr__()
        return f"{base}\nTarget: {self.target.name}\n"

    def _set_mean_velocity(self):
        """
        Sets the mean velocity of the projectile in m/s based on the target body.

        Returns
        -------
        float
        """
        known_targets = [
            "Mercury",
            "Venus",
            "Earth",
            "Moon",
            "Mars",
            "Ceres",
            "Vesta",
            "MBA",
        ]
        target_velocities = [
            41100.0,
            29100.0,
            24600.0,
            22100.0,
            10700.0,
            5300.0,
            5300.0,
            5300.0,
        ]
        catalogue = dict(zip(known_targets, target_velocities))
        if self.target_name in known_targets:
            pmv = float(catalogue[self.target_name])
        else:
            warn(
                f"Target {self.target_name} not found in known targets. Known targets include {known_targets}. Defaulting to the Moon."
            )
            pmv = float(catalogue["Moon"])
        return pmv

    @parameter
    def target_name(self):
        """
        The name of the target body.

        Returns
        -------
        str
        """
        return self._target.name

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
