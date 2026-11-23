from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray

from cratermaker.components.counting import R_to_CSFD
from cratermaker.components.production import Production
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter


@Production.register("neukum")
class NeukumProduction(Production):
    """
    An operations class for computing the the Neukum production function for the Moon and Mars.

    Notes
    -----
    The CSFD is computed using the model of Ivanov, Neukum, and Hartmann (2001) for the Moon and Mars, with
    minor changes. Notably, there is a typo in the table of coefficients in the original paper. The value of a0 should be -3.0768 and not -3.0876, as the
    linear term of the chronology function is 10^(a0) = 8.38e-4. We change this due to the fact that we compute the linear term of the chronology function from the coefficients
    in our implementation of the chronology function.

    The Mars production function is based on Ivanov (2001). The projectile production function is from Ivanov, Neukum, and Wagner (2001).

    References
    ----------
    - Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to
      the Lunar Reference System. *Space Science Reviews*, 96, 55-86. `doi:10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`_
    - Ivanov, B.A., 2001. Mars/Moon Cratering Rate Ratio Estimates. *Space Science Reviews*, 96, 87-104. `doi:10.1023/A:1011941121102 <https://doi.org/10.1023/A:1011941121102>`_
    - Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters
      and Asteroids. In *Collisional Processes in the Solar System*, Springer Netherlands, Dordrecht, pp. 1-34.
      https://doi.org/10.1007/978-94-010-0712-2_1
    - Neukum, G., 1983. Meteorite bombardment and planetary surfaces dating. Habilitation Dissertation for Faculty Membership, University of Munich.

    """

    def __init__(
        self,
        version: str | None = None,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        version : {"Moon", "Mars", "Projectile"}, optional
            The specific model to use for the production function. "Moon" and "Mars" are both crater production functions, and
            "Projectile" is a projectile function. Default is "Moon".
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        object.__setattr__(self, "_Cexp", None)
        object.__setattr__(self, "_Clin", None)
        object.__setattr__(self, "_tau", None)

        self.version = version or "Moon"
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    def __str__(self) -> str:
        str_repr = super().__str__()
        timelo = format_large_units(self.valid_time[0], quantity="time")
        timehi = format_large_units(self.valid_time[1], quantity="time")
        dlo = format_large_units(self.sfd_range[0], quantity="length")
        dhi = format_large_units(self.sfd_range[1], quantity="length")
        str_repr += f"Version: {self.version}\n"
        str_repr += f"Valid Time Range: {timelo} - {timehi}\n"
        str_repr += f"Valid Diameter Range: {dlo} - {dhi}"
        return str_repr

    def chronology(
        self,
        time_start: ArrayLike = (1000.0),
        time_end: ArrayLike | None = None,
        validate_inputs: bool = True,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Returns the relative number of craters produced over a given age range.

        This implements the chronology function given in Eq. 5 of Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86, but takes in the age argument in the Cratermaker unit system of My instead of Gy. The result is scaled to the Cratermaker default unit system to give the cumulative number of craters per m² relative to the number in 1 My.

        Parameters
        ----------
        time_start : ArrayLike, default=1000.0
            Age in the past relative to the present day to compute cumulative SFD in units of My.
        time_end : ArrayLike, default=0.0
            The ending age in the past relative to the present day to compute cumulative SFD in units of My. The default is 0 (present day).
        validate_inputs : bool, default=True
            If True, the function will check that the age is within the valid range of ages for this production function. If False, no check is performed.

        Returns
        -------
        NDArray
            The number of craters relative to the amount produced in the last 1 My between age and ange_end.

        """
        if "age" in kwargs:
            time_start = kwargs.pop("age")
        if validate_inputs:
            time_start, time_end = self._validate_age(time_start, time_end)

        def _N1km(
            age: Sequence[FloatLike] | ArrayLike,
        ) -> FloatLike | ArrayLike:
            """
            Return the cumulative number of 1 km craters per km² as a function of age in Gy.

            This is a direct implementation of Eq. 5 in Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86 (with corrected coefficient for the linear term).

            Parameters
            ----------
            age : FloatLike or numpy array
                Time ago in units of My

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than 1 km in diameter
            """
            N1 = self.Cexp * (np.exp(age * 1e-3 / self.tau) - 1.0) + self.Clin * age * 1e-3
            return N1

        _N1km_norm = _N1km(1.0)

        if time_end is None:
            return _N1km(age=time_start) / _N1km_norm
        else:
            return (_N1km(age=time_start) - _N1km(age=time_end)) / _N1km_norm

    @property
    def valid_versions(self) -> list[str]:
        """
        The valid versions of the production function.

        Returns
        -------
        list
            The valid versions of the production function.
        """
        return ["Moon", "Mars", "Projectile"]

    @property
    def generator_type(self) -> str:
        """
        Get the generator type of the production function.

        Returns
        -------
        str
            The generator type of the production function.
        """
        return self._generator_type

    @generator_type.setter
    def generator_type(self, value: str) -> None:
        """
        The generator type for this production function is determined by the version.
        """
        raise NotImplementedError("The generator type cannot be changed for this production function.")

    @parameter
    def version(self) -> str:
        """
        Get the version of the production function.

        Returns
        -------
        str
            The version of the production function.
        """
        return self._version

    @version.setter
    def version(self, value: str) -> None:
        """
        Set the version of the production function.

        Parameters
        ----------
        value : str
            The version of the production function.
        """
        if value.title() not in self.valid_versions:
            raise ValueError(f"Invalid version '{value}'. Valid options are {self.valid_versions}")
        self._version = value.title()

        if self._version == "Projectile":
            self._generator_type = "projectile"
        else:
            self._generator_type = "crater"

    @property
    def sfd_tables(self) -> dict[str, list[float]]:
        """
        Get the table of coefficients for the production function.

        Returns
        -------
        dict
            The table of coefficients for the production function.
        """
        return {
            "Moon": [
                -3.0768,  # Typo corrected
                -3.557528,
                +0.781027,
                +1.021521,
                -0.156012,
                -0.444058,
                +0.019977,
                +0.086850,
                -0.005874,
                -0.006809,
                +8.25e-4,
                +5.54e-5,
            ],
            "Mars": [
                -3.384,
                -3.197,
                +1.257,
                +0.7915,
                -0.4861,
                -0.3630,
                +0.1016,
                +6.756e-2,
                -1.181e-2,
                -4.753e-3,
                +6.233e-4,
                +5.805e-5,
            ],
            "Projectile": [
                0,
                +1.375458,
                +1.272521e-1,
                -1.282166,
                -3.074558e-1,
                +4.149280e-1,
                +1.910668e-1,
                -4.260980e-2,
                -3.976305e-2,
                -3.180179e-3,
                +2.799369e-3,
                +6.892223e-4,
                +2.614385e-6,
                -1.416178e-5,
                -1.191124e-6,
            ],
        }

    @property
    def sfd_coef(self) -> list[float]:
        """
        Get the table of coefficients for the production function.

        Returns
        -------
        list
            The table of SFD coefficients for the production function.
        """
        return self.sfd_tables[self.version]

    @property
    def Clin(self) -> float:
        """
        Get the linear coefficient for the production function.

        Returns
        -------
        float
            The linear coefficient for the production function.
        """
        if self._Clin is None:
            self._Clin = 10 ** (self.sfd_coef[0])
        return self._Clin

    @property
    def Cexp(self) -> float:
        """
        Get the exponential coefficient for the production function.

        Returns
        -------
        float
            The exponential coefficient for the production function.
        """
        if self._Cexp is None:
            Cexp_moon = 5.44e-14
            if self.version == "Moon":
                self._Cexp = Cexp_moon
            else:
                Clin_moon = 10 ** (self.sfd_tables["Moon"][0])
                self._Cexp = Cexp_moon * self.Clin / Clin_moon
        return self._Cexp

    @property
    def tau(self) -> float:
        """
        Get the time constant for the production function.

        Returns
        -------
        float
            The time constant for the production function.
        """
        if self._tau is None:
            self._tau = 1.0 / 6.93
        return self._tau

    @property
    def sfd_range(self) -> tuple[float, float]:
        """
        The range of diameters over which the SFD is valid. The range is given in m.

        Returns
        -------
        tuple
            The lower and upper bounds of the SFD range in m.
        """
        sfd_range = {
            "Moon": (10.0, 1000.0e3),
            "Mars": (15.0, 362.0e3),
            "Projectile": (
                0.1,
                200.0e3,
            ),  # Estimated based on Fig. 16 of Ivanov et al. (2001)
        }

        return sfd_range[self.version]

    @property
    def valid_time(self) -> tuple[float, float]:
        """
        The range of ages over which the production function is valid. The range is given in My.

        Returns
        -------
        tuple
            The lower and upper bounds of the valid time range in My.
        """
        return (0, 4500)

    def csfd(
        self,
        diameter: FloatLike | ArrayLike,
        validate_inputs: bool = True,
        **kwargs,
    ) -> FloatLike | ArrayLike:
        """
        Return the cumulative size frequency distribution of craters larger than diameter in m relative to age = 1 My ago per m².

        Parameters
        ----------
        diameter : FloatLike or ArrayLike
            units of meter
        validate_inputs: bool, default=True
            If True, the function will check that the validity of age, time_end, and diameter arguments. If False, no check is performed, and the arguments are assumed to be valid. This can be used to speed up the function, particularly if it is called as part of a solver or optimization routine where the inputs are already known to be valid.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        FloatLike or numpy array
           Cumulative number density of per square meter greater than the input diameter.
        """
        # Convert m to km for internal functions
        Dkm_lo = self.sfd_range[0] * 1e-3
        Dkm_hi = self.sfd_range[1] * 1e-3

        def _extrapolate_sfd(side: str = "lo") -> FloatLike | ArrayLike:
            """
            Return the exponent, p, and and proportionality constant, A, for  the extrapolated CSFD in the form N(D) = A * D**-p.

            Parameters
            ----------
            side : str
                The side of the range to extrapolate. Valid values are "lo" and "hi"

            Returns
            -------
            A, p
            """
            if side == "lo":
                Dkm = Dkm_lo
            elif side == "hi":
                Dkm = Dkm_hi
            else:
                raise ValueError("side must be 'lo' or 'hi'")
            p = _dNdD(Dkm)
            A = _CSFD(Dkm)
            return A, p

        def _dNdD(
            Dkm: FloatLike | Sequence[FloatLike] | ArrayLike,
        ) -> FloatLike | ArrayLike:
            """
            Return the derivative of the cumulative size-frequency distribution as a function of diameter.

            For diameter values outside the range of the NPF, the derivative is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The differential number of craters (dN/dD) per square kilometer greater than Dkm in diameter at age = 1 Gy ago.
            """

            def _dNdD_scalar(Dkm):
                dcoef = self.sfd_coef[1:]
                if Dkm < Dkm_lo:
                    A, p = _extrapolate_sfd(side="lo")
                    return A * (p / Dkm) * (Dkm / Dkm_lo) ** p
                elif Dkm > Dkm_hi:
                    A, p = _extrapolate_sfd(side="hi")
                    return A * (p / Dkm) * (Dkm / Dkm_hi) ** p
                else:
                    return sum(co * np.log10(Dkm) ** i for i, co in enumerate(dcoef))

            return _dNdD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_dNdD_scalar)(Dkm)

        def _CSFD(
            Dkm: FloatLike | Sequence[FloatLike] | ArrayLike,
        ) -> FloatLike | ArrayLike:
            """
            Return the cumulative size-frequency distribution at the reference age of 1 Gy ago.

            For diameter values outside the range of the NPF, the CSFD is extrapolated using a power law.

            Parameters
            ----------
            Dkm : FloatLike or numpy array
                diameter in units of km

            Returns
            -------
            FloatLike or numpy array
                The number of craters per square kilometer greater than Dkm in diameter at age=1 Gy ago.
            """

            def _CSFD_scalar(Dkm):
                if Dkm < Dkm_lo:
                    A, p = _extrapolate_sfd(side="lo")
                    return A * (Dkm / Dkm_lo) ** p
                elif Dkm > Dkm_hi:
                    A, p = _extrapolate_sfd(side="hi")
                    p -= 5.0  # Steepen the upper branch of the SFD to prevent anomolously large craters from forming
                    return A * (Dkm / Dkm_hi) ** p
                else:
                    logCSFD = sum(co * np.log10(Dkm) ** i for i, co in enumerate(self.sfd_coef))
                    return 10**logCSFD

            return _CSFD_scalar(Dkm) if np.isscalar(Dkm) else np.vectorize(_CSFD_scalar)(Dkm)

        if validate_inputs:
            diameter, _ = self._validate_csfd(diameter=diameter)

        Dkm = diameter * 1e-3  # Convert m to km for internal functions

        # Cratermaker uses 1 My as the reference age instead of 1 Gy, so this conversion factor ensures that the value of the CSFD is exactly the NPF value at 1 Gy
        age_norm = 1000.0 / self.chronology(age=1000.0)

        if self.version == "Projectile":
            # This is a normalization factor that gets the projectile CSFD to approximately match the lunar crater CSFD
            Rnorm = 2.94e-5
            Ncumulative = R_to_CSFD(R=_CSFD, D=Dkm) * Rnorm
        else:
            Ncumulative = _CSFD(Dkm)

        return Ncumulative * age_norm * 1e-9  # convert from Gy^-1 km^-2 to My^-1 m^-2
