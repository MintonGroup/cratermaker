"""
Define a custom production function
=====================================

.. rubric:: By David Minton

This example demonstrates how to create a custom production function by defining a class derived from the base Production class. The only requirements are that the following two abstract methods are given concrete definitions:

.. code-block:: python

    @abstractmethod
    def _chronology(
        self,
        time_start: FloatLike | ArrayLike,
        time_end: FloatLike | ArrayLike,
        **kwargs: Any,
    ) -> FloatLike | ArrayLike: ...


    @abstractmethod
    def _csfd(self, diameter: FloatLike | ArrayLike, **kwargs: Any) -> FloatLike | ArrayLike: ...


These are called by the public :py:meth:`~cratermaker.components.production.Production.chronology` and :py:meth:`~cratermaker.components.production.Production.csfd` methods after they have validated the inputs, so your implementations do not need to do any input checking. In particular, ``time_end`` will never be ``None`` when ``_chronology`` is called.

You then register a name for your production function and then can apply it to a Simulation.

For this example, we will generate an outer solar system production function based on Case A of Zahnle et al. (2003) [#]_. This particular productdion function is defined in terms of projectiles and the impact rates are scaled relative to Jupiter. Therefore a conversion must be done using values given in Table 1 of Zahnle et al. (2003) for outer solar system satellites.


References
----------

.. [#] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263-289. `doi: 10.1016/S0019-1035(03)00048-4 <https://doi.org/10.1016/S0019-1035(03)00048-4>`_
"""

import numpy as np

from cratermaker import Production, Simulation


@Production.register("zahnle_case_a")
class ZahnleCaseA(Production):
    # Ecliptic comet impact probability relative to Jupiter (Zahnle et al. 2003, Table 1)
    _P_EC = {
        # Jupiter system
        "Jupiter": 1.0,
        "Metis": 2.8e-7,
        "Amalthea": 7.7e-7,
        "Thebe": 2.9e-8,
        "Io": 1.4e-4,
        "Europa": 6.6e-5,
        "Ganymede": 1.2e-4,
        "Callisto": 6.1e-5,
        "Himalia": 1.4e-8,
        # Saturn system
        "Saturn": 0.42,
        "Prometheus": 1.7e-7,
        "Pandora": 1.0e-7,
        "Epimetheus": 1.8e-7,
        "Janus": 4.5e-7,
        "Mimas": 1.7e-6,
        "Enceladus": 2.2e-6,
        "Tethys": 7.9e-6,
        "Telesto": 3.5e-9,  # Table 1 lists Telesto/Calypso as a single row
        "Calypso": 3.5e-9,
        "Dione": 7.1e-6,
        "Helene": 5.8e-9,
        "Rhea": 9.6e-6,
        "Titan": 5.4e-5,
        "Hyperion": 1.0e-7,
        "Iapetus": 1.4e-6,
        "Phoebe": 8.7e-9,
        # Uranus system
        "Uranus": 0.25,
        "Cordelia": 3.3e-8,
        "Ophelia": 4.2e-8,
        "Bianca": 8.2e-8,
        "Cressida": 1.6e-7,
        "Desdemona": 1.2e-7,
        "Juliet": 3.0e-7,
        "Portia": 4.7e-7,
        "Rosalind": 1.3e-7,
        "Belinda": 1.4e-7,
        "Puck": 6.6e-7,
        "Miranda": 5.2e-6,
        "Ariel": 2.1e-5,
        "Umbriel": 1.6e-5,
        "Titania": 1.8e-5,
        "Oberon": 1.3e-5,
        # Neptune system
        "Neptune": 0.27,  # Halved from the 0.54 given by Levison et al. (see Zahnle et al. text)
        "Naiad": 2.0e-7,
        "Thalassa": 3.8e-7,
        "Despina": 1.2e-6,
        "Galatea": 1.2e-6,
        "Larissa": 1.6e-6,
        "Proteus": 4.8e-6,
        "Triton": 7.5e-5,
        "Nereid": 1.2e-7,
        # Pluto system
        "Pluto": 2.0e-4,
        "Charon": 3.2e-5,
    }

    def _csfd(self, diameter, **kwargs):
        """
        Implements the Zahnle et al. Case A projectile CSFD for Jupiter given by eq. (14).

        Return the cumulative size frequency distribution of  projectiles larger than diameter in m relative to age = 1 My ago per m².

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            units of meter
        **kwargs : Any
            Keyword arguments that are passed to the base class (required in order to pass argument validation flags).

        Returns
        -------
        Numpy array
           Cumulative number density of per square meter greater than the input diameter.
        """
        if self.target.name not in self.__class__._P_EC:
            raise ValueError(f"Target {self.target.name} is not valid for this production function")
        p_ec = self.__class__._P_EC[self.target.name]

        # Zahnle et al. (2003) calibrate Case A to 0.005 impacts per year on Jupiter by comets with d > 1.5 km (p. 274). P_EC is the ratio
        # of the total number of impacts on the target body to the number on Jupiter, so dividing by the target's surface area gives the
        # impact rate per unit area. The factor of 1e6 converts from per year to per My.
        n_target = p_ec * 0.005 * 1e6 / self.target.surface_area

        d = diameter * 1e-3  # Convert from m to km to define "d" as in the equation
        return n_target * np.where(
            d < 1.5,
            (d / 1.5) ** (-1.0),
            np.where(
                d < 5.0,
                (d / 1.5) ** (-1.7),
                np.where(
                    d < 30.0,
                    (1.5 / 5.0) ** (1.7) * (d / 5.0) ** (-2.5),
                    (1.5 / 5.0) ** (1.7) * (5 / 30.0) ** (2.5) * (d / 30) ** (-3.2),
                ),
            ),
        )

    def _chronology(self, time_start, time_end, **kwargs):
        """
        Return the number of impacts between time_start and time_end relative to the number produced in the last 1 My at the present-day rate.

        Zahnle et al. (2003) assume that the impact flux declines as 1/t, where t is the time since the origin of the Solar System 4.56 Gy ago
        (Section 5). The rate at time τ My before present is therefore proportional to T0 / (T0 - τ), with T0 = 4560 My, normalized so that the
        present-day rate is the one given by Eq. (14). Integrating from time_end to time_start gives the expression below.

        Parameters
        ----------
        time_start : FloatLike or ArrayLike
            Age in the past relative to the present day in units of My.
        time_end : FloatLike or ArrayLike
            The ending age in the past relative to the present day in units of My.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        FloatLike or ArrayLike
            The number of impacts relative to the amount produced in 1 My at the present-day rate.
        """
        T0 = 4560.0
        return T0 * np.log((T0 - time_end) / (T0 - time_start))

    @property
    def valid_time(self) -> tuple[float, float]:
        """
        The range of ages over which the production function is valid, in My.

        The 1/t chronology diverges at the origin of the Solar System (4560 My ago), so ages are restricted to 4550 My.
        """
        return (0.0, 4550.0)

    @property
    def generator_type(self) -> str:
        """
        Get the generator type of the production function.

        Returns
        -------
        str
            The generator type of the production function.
        """
        return "projectile"

    @generator_type.setter
    def generator_type(self, value: str) -> None:
        """
        The generator type for this production function is determined by the version.
        """
        raise NotImplementedError("The generator type cannot be changed for this production function.")


if __name__ == "__main__":
    # Sample the production function for Juliet
    # ------------------------------------------
    #
    # Here we make a version of Fig. 8 of Zahnle et al. (2003) for Case A on the uranian satellite Juliet. Unlike Fig. 8, which shows crater
    # diameters, our production function generates projectiles, so the x-axis is crater diameter. We sample the population of projectiles
    # that strike Juliet over 100 My, repeat this 100 times, and plot each sample against the expected value given by
    # :py:meth:`Production.function`. The right-hand axis converts the cumulative number into the equivalent present-day cumulative impact rate
    # per 10⁶ km² per year, which is the quantity plotted in Fig. 8.

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    simdir = "simdata-2.9"

    sim = Simulation(
        simdir=simdir, gridlevel=6, target="Juliet", production="zahnle_case_a", rng_seed=788563223, ask_overwrite=False, reset=True
    )
    target = sim.target
    area = target.surface_area
    age = 100.0
    diameter_range = (1.0, 2 * target.radius)  # Range of projectile diameters to generate in m (up to the size of Juliet itself)
    Nevaluations = 100

    x_min = diameter_range[0] * 1e-3
    x_max = diameter_range[1] * 1e-3
    fig, ax = plt.subplots(figsize=(6, 6), layout="constrained")
    ax.set_title(f"Zahnle et al. (2003) Case A: {target.name}, {age:.0f} My")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Crater Diameter (km)")
    ax.set_ylabel(f"$\\mathregular{{N_{{>D}}}}$ on {target.name}")
    ax.set_ylim(1, 100)
    ax.set_xlim(1, 100)
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Plot each sampled population Monte Carlo style
    for _ in range(Nevaluations):
        craters = sim.populate(age=age, emplace=False)
        Dsampled = [c.diameter for c in craters]
        Dsampled = np.sort(Dsampled)[::-1]
        Nsampled = range(1, len(Dsampled) + 1)
        ax.plot(Dsampled * 1e-3, Nsampled, "-", color="cornflowerblue", linewidth=2.0, zorder=50, alpha=0.2)

    # Secondary axis: equivalent present-day cumulative impact rate in units of per 10⁶ km² per year, as in Fig. 8 of Zahnle et al. (2003)
    rate_factor = 1.0 / (area * 1e-12) / (sim.production.chronology(age=age) * 1e6)
    secax = ax.secondary_yaxis("right", functions=(lambda N: N * rate_factor, lambda R: R / rate_factor))
    secax.set_ylabel("Present-day cratering rate [per 10⁶ km² yr⁻¹]")

    plt.show()
