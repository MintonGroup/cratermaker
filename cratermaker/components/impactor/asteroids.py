from numpy.random import Generator
from typing import Any
from cratermaker.components.impactor import ImpactorModel, register_impactor_model
from warnings import warn

@register_impactor_model("asteroids")
class AsteroidImpactors(ImpactorModel):
    def __init__(self, 
                 target_name : str = "Moon",
                 density : float = 2000.0,
                 rng: Generator | None = None,
                 **kwargs: Any):
        """
        An operations class for computing the impactor properties of an asteroid source population.

        Parameters
        ----------
        target_name : str
            The name of the target body for the impact.
        density : float
            The density of the impactor in kg/m^3. Default is 2000.0 kg/m^3.
        rng : Generator | None
            A random number generator for Monte Carlo simulations. If None, a default generator will be used.

        **kwargs : Any
            Additional keyword arguments to be passed to internal functions.
        """

        # This model always samples velocities, angles, and directions, so override any values that may have been passed.
        kwargs["sample_velocities"] = True
        kwargs["sample_angles"] = True
        kwargs["sample_directions"] = True
        kwargs["mean_velocity"] = self._set_mean_velocity()
        super().__init__(target_name=target_name, density=density, rng=rng, **kwargs)


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
            pmv = float(catalogue[self.target.name])
        else:
            warn(f"Target {self.target_name} not found in known targets. Known targets include {known_targets}. Defaulting to the Moon.")
            pmv = float(catalogue["Moon"])
        return pmv
