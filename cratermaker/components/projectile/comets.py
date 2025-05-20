from typing import Any
from warnings import warn

from numpy.random import Generator

from cratermaker.components.projectile import Projectile
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike


@Projectile.register("comets")
class CometProjectiles(Projectile):
    _catalogue = {
        "Mercury": 18400.0,
        "Venus": 26300.0,
        "Earth": 17670.0,
        "Moon": 17670.0,
        "Mars": 14210.0,
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
        "KBO": 1800.0,
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
        An operations class for computing the projectile properties of a comet source population.

        Parameters
        ----------
        sample : bool, default True
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, the projectile velocity will be the mean velocity for the given target body.
        density : float
            The density of the projectile in kg/m^3. Default is 500.0 kg/m^3.
        angle : float, optional
            The impact angle in degrees. Default is 90.0 degrees (vertical impact) if `sample` is False. If `sample` is True, this value is ignored.
        direction : float | None
            The impact direction in degrees. Default is 0.0 degrees (due North) if `sample` is False. If `sample` is True, this value is ignored.`
        location : tuple[float, float] | None
            The impact location as a tuple of (longitude, latitude) in degrees. Default is (0.0, 0.0) (the equator and prime meridian) if `sample` is False. If `sample` is True, this value is ignored.
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        rng : Generator | None
            A random number generator for Monte Carlo simulations. If None, a default generator will be used.
        rng_seed : int | None
            The random rng_seed for the simulation if rng is not provided. If None, a random rng_seed is used.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.

        Notes
        -----
        The comets are assumed to be Jupiter-family Comets. The mean impact velocities for outer solar system bodies come from Table 1 of Zahnle et al. [#]_, and for inner solar system bodies, Table 2 of Borin et al. [#]_.

        References
        ----------
        .. [#] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263-289. https://doi.org/10.1016/S0019-1035(03)00048-4
        .. [#] Borin, P., Cremonese, G., Marzari, F., Lucchetti, A., 2017. Asteroidal and cometary dust flux in the inner solar system. A&A 605, A94. https://doi.org/10.1051/0004-6361/201730617


        **kwargs : Any
            Additional keyword arguments to be passed to internal functions.
        """
        density = 500.0 if density is None else density
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
                f"Target {self.target.name} not found in known targets. Known targets include {list(self.__class__._catalogue.keys())}. Defaulting to KBO."
            )
            pmv = float(self.__class__._catalogue["KBO"])
        return pmv
