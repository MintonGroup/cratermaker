from numpy.random import Generator
from typing import Any
from cratermaker.components.target import Target
from cratermaker.components.projectile import Projectile
from cratermaker.utils.custom_types import FloatLike
from warnings import warn

@Projectile.register("asteroids")
class AsteroidProjectiles(Projectile):
    def __init__(self, 
                 target : Target | str | None = None,
                 density : FloatLike | None = None,
                 rng: Generator | None = None,
                 rng_seed : int | None = None,
                 rng_state : dict | None = None, 
                 **kwargs: Any):
        """
        An operations class for computing the projectile properties of an asteroid source population.

        Parameters
        ----------
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        density : float
            The density of the projectile in kg/m^3. Default is 2250.0 kg/m^3.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """

        # This model always samples velocities, angles, and directions, so override any values that may have been passed.
        if density is None:
            density = 2250.0
        super().__init__(target=target, density=density, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        self.mean_velocity = self._set_mean_velocity()


    def _set_mean_velocity(self):
        """
        Sets the mean velocity of the projectile in m/s based on the target body.
        
        Returns
        -------
        float 
        """
        known_targets = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Ceres', 'Vesta', 'MBA']
        target_velocities = [41100.0, 29100.0, 24600.0, 22100.0, 10700.0, 5300.0, 5300.0, 5300.0]
        catalogue = dict(zip(known_targets, target_velocities))
        if self.target_name in known_targets:
            pmv = float(catalogue[self.target_name])
        else:
            warn(f"Target {self.target_name} not found in known targets. Known targets include {known_targets}. Defaulting to the Moon.")
            pmv = float(catalogue["Moon"])
        return pmv
