from numpy.random import Generator
from typing import Any
from cratermaker.components.projectile import Projectile
from warnings import warn
@Projectile.register("comets")
class CometProjectiles(Projectile):
    def __init__(self, 
                 target_name : str | None = None,
                 density : float | None = None,
                 rng: Generator | None = None,
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs: Any):
        """
        An operations class for computing the projectile properties of a comet source population.

        Parameters
        ----------
        target_name : str
            The name of the target body for the impact.
        density : float
            The density of the projectile in kg/m^3. Default is 500.0 kg/m^3.
        rng : Generator | None
            A random number generator for Monte Carlo simulations. If None, a default generator will be used.
        rng_seed : int | None
            The random rng_seed for the simulation if rng is not provided. If None, a random rng_seed is used.    
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        Notes
        -----
        The mean impact velocities for outer solar system bodies come from Table 1 of Zahnle et al. [1]_. For inner solar system bodies, from Table 2 of Borin et al. [2]_. 

        References
        ----------
        .. [1] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263-289. https://doi.org/10.1016/S0019-1035(03)00048-4
        .. [2] Borin, P., Cremonese, G., Marzari, F., Lucchetti, A., 2017. Asteroidal and cometary dust flux in the inner solar system. A&A 605, A94. https://doi.org/10.1051/0004-6361/201730617


        **kwargs : Any
            Additional keyword arguments to be passed to internal functions.
        """

        # This model always samples velocities, angles, and directions, so override any values that may have been passed.
        if density is None:
            density = 500.0
        super().__init__(target_name=target_name, density=density, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        self.mean_velocity = self._set_mean_velocity()


    def _set_mean_velocity(self):
        """
        Sets the mean velocity of the projectile in m/s based on the target body.
        
        Returns
        -------
        float 
        """

        catalogue = {
            "Mercury" : 18400.0,
            "Venus" : 26300.0,
            "Earth" : 17670.0,
            "Moon" : 17670.0,
            "Mars" : 14210.0,
            "Metis": 59000.0,
            "Amalthea": 50000.0,
            "Thebe": 45000.0,
            "Io": 32000.0,
            "Europa": 26000.0,
            "Ganymede": 20000.0,
            "Callisto": 15000.0,
            "Himalia": 6100.0,
            "Prometheus": 32000.0,
            "Pandora": 31000.0,
            "Epimetheus": 30000.0,
            "Janus": 30000.0,
            "Mimas": 27000.0,
            "Enceladus": 24000.0,
            "Tethys": 21000.0,
            "Telesto": 21000.0,
            "Calypso": 21000.0,
            "Dione": 19000.0,
            "Helene": 19000.0,
            "Rhea": 16000.0,
            "Titan": 10500.0,
            "Hyperion": 9400.0,
            "Iapetus": 6100.0,
            "Phoebe": 3200.0,
            "Cordelia": 20000.0,
            "Ophelia": 19000.0,
            "Bianca": 18000.0,
            "Cressida": 18000.0,
            "Desdemona": 18000.0,
            "Juliet": 18000.0,
            "Portia": 18000.0,
            "Rosalind": 17000.0,
            "Belinda": 16000.0,
            "Puck": 15000.0,
            "Miranda": 12500.0,
            "Ariel": 10300.0,
            "Umbriel": 8700.0,
            "Titania": 6800.0,
            "Oberon": 5900.0,
            "Naiad": 22000.0,
            "Thalassa": 22000.0,
            "Despina": 21000.0,
            "Galatea": 20000.0,
            "Larissa": 18000.0,
            "Proteus": 14000.0,
            "Triton": 8200.0,
            "Nereid": 2800.0,
            "Pluto": 1900.0,
            "Charon": 1800.0,
            "KBO" : 1800.0,
        }
        if self.target_name in catalogue:
            pmv = float(catalogue[self.target_name])
        else:
            warn(f"Target {self.target_name} not found in known targets. Known targets include {list(catalogue.keys())}. Defaulting to KBO.")
            pmv = float(catalogue["KBO"])
        return pmv
