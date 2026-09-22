from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import InitVar, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy import fft

from cratermaker.bindings import basicmoon_bindings, realmoon_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import (
    BasicMoonCrater,
    BasicMoonCraterFixed,
    BasicMoonMorphology,
)
from cratermaker.constants import FloatLike
from cratermaker.core.base import CratermakerBase
from cratermaker.utils.general_utils import format_large_units, parameter

_PSD1D_COEF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_PSD1D_MIN_POINTS = 64


class PSD1D(CratermakerBase):
    def __init__(
        self,
        mean: FloatLike,
        pix: FloatLike,
        npoints: int | None = None,
        control_points: dict[str, np.float64] | None = None,
        power: NDArray[np.float64] | None = None,
        wavelength: NDArray[np.float64] | None = None,
        phase: NDArray[np.float64] | None = None,
        add_noise: bool = True,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        Contains data and methods for building and using a 1D power spectral density model for profile functions.

        Parameters
        ----------
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_mean", None)
        object.__setattr__(self, "_pix", None)
        object.__setattr__(self, "_npoints", None)
        object.__setattr__(self, "_control_points", None)
        object.__setattr__(self, "_power", None)
        object.__setattr__(self, "_wavelength", None)
        object.__setattr__(self, "_phase", None)
        object.__setattr__(self, "_add_noise", None)

        self.mean = mean
        self.pix = pix
        self.npoints = npoints
        self.control_points = control_points
        self.wavelength = wavelength
        self.power = power
        self.phase = phase
        self.add_noise = add_noise

    @property
    def mean(self) -> np.float64:
        """
        The mean value of the profile.
        """
        return self._mean

    @mean.setter
    def mean(self, value):
        if value > 0.0:
            self._mean = np.float64(value)
        else:
            raise ValueError("mean must be a positive number")

    @property
    def pix(self) -> np.float64:
        """
        The resolution of the profile in meters.
        """
        return self._pix

    @pix.setter
    def pix(self, value):
        if value > 0.0:
            self._pix = np.float64(value)
        else:
            raise ValueError("pix must be a positive number")

    @property
    def npoints(self) -> int:
        """
        The number of power spectral density points in the signal.
        """
        if self._npoints is None:
            if self._wavelength is not None:
                return len(self._wavelength)
            elif self.pix is not None:
                return max(int(4 * math.pi * self.mean / self.pix), _PSD1D_MIN_POINTS)
        return self._npoints

    @npoints.setter
    def npoints(self, value):
        if value is not None and (not isinstance(value, int) or value < _PSD1D_MIN_POINTS):
            raise ValueError(f"npoints must be a positive integer larger than {_PSD1D_MIN_POINTS}")
        else:
            self._npoints = value

    @property
    def wavelength(self) -> NDArray[np.float64]:
        if self._wavelength is None:
            psd_arrs = realmoon_bindings.get_1d_psd_from_control_points(
                control_points=self.control_points,
                npoints=self.npoints,
                add_noise=self.add_noise,
                rng_seed=self.rng_seed,
            )
            return psd_arrs[0]
        else:
            return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        if value is not None:
            value = np.asarray(value)
            if self._npoints is None:
                self.npoints = len(value)
            else:
                if len(value) != self.npoints:
                    raise ValueError(f"Size of wavelength must be {self.npoints}")
        self._wavelength = value

    @property
    def power(self) -> NDArray[np.float64]:
        if self._power is None:
            psd_arrs = realmoon_bindings.get_1d_psd_from_control_points(
                control_points=self.control_points,
                npoints=self.npoints,
                add_noise=self.add_noise,
                rng_seed=self.rng_seed,
            )
            return psd_arrs[1]
        else:
            return self._power

    @power.setter
    def power(self, value):
        if value is not None:
            value = np.asarray(value)
            if self._npoints is None:
                self.npoints = len(value)
            else:
                if len(value) != self.npoints:
                    raise ValueError(f"Size of power must be {self.npoints}")
        self._power = value

    @property
    def phase(self) -> NDArray[np.float64]:
        if self._phase is None:
            psd_arrs = realmoon_bindings.get_1d_psd_from_control_points(
                control_points=self.control_points,
                npoints=self.npoints,
                add_noise=self.add_noise,
                rng_seed=self.rng_seed,
            )
            return psd_arrs[2]
        else:
            return self._phase

    @phase.setter
    def phase(self, value):
        if value is not None:
            value = np.asarray(value)
            if self._npoints is None:
                self.npoints = len(value)
            else:
                if len(value) != self.npoints:
                    raise ValueError(f"Size of phase must be {self.npoints}")
            self._phase = value

    @property
    def control_points(self) -> dict[str, np.float64]:
        """
        Control points for a procedurally-generated PSD model.
        """
        return self._control_points

    @control_points.setter
    def control_points(self, value):
        if not isinstance(value, dict):
            raise TypeError("control_points must be a dict")
        required_keys = ["sn", "yn", "y1", "y2", "y3", "y4", "y5"]
        if any(k not in value for k in required_keys):
            raise ValueError(f"control_points must have all required keys {required_keys}")
        self._control_points = value

    def set_from_profile(self, profile):
        """
        Sets the 1D power spectral density (PSD) of a given 1D profile.

        This function will set the following properties:
            wavelength : NDArray[np.float64]
                An array of wavelengths corresponding to the frequencies of the PSD, and the second column contains the corresponding PSD values. The array is sorted in descending order of wavelength (i.e., ascending order of frequency).
            power : NDArray[np.float64]
                An array of power spectral density values corresponding to the frequencies of the input signal. The array is sorted in descending order of wavelength (i.e., ascending order of frequency). The PSD values are normalized by the interval and the number of points in the input signal, such that they represent the power per unit wavelength. The normalization is done by multiplying the squared magnitude of the Fourier coefficients by 2 and dividing by the product of the interval and the number of points in the input signal. The factor of 2 accounts for the fact that we are using a one-sided PSD (i.e., only considering positive frequencies).
            phases : NDArray[np.float64]
                An array of phase values corresponding to the frequencies of the input signal. The phase values are computed from the angle of the Fourier coefficients and are normalized by the frequency to represent the phase shift in terms of spatial units (e.g., meters). The phase values are sorted in descending order of wavelength (i.e., ascending order of frequency).

        Parameters
        ----------
        profile : ArrayLike
            The input signal for which to compute the 1D PSD. This should be a 1D array of values representing the signal in the spatial domain.

        """
        y = np.asarray(profile, dtype=np.float64)
        n = len(y)
        ymean = np.mean(y)
        dfft = fft.rfft(y - ymean) / n

        interval = 2 * math.pi / n
        power = (2 * np.abs(dfft)) ** 2 / (interval * n)

        if n % 2 == 0:
            index_end = n // 2
        else:
            index_end = n // 2 + 1
        freq = fft.fftfreq(n, interval)
        wavelength = 1 / freq[1:index_end]
        phase = np.angle(dfft[1:index_end])

        self.npoints = n
        self.wavelength = wavelength
        self.power = power[1:index_end]
        self.phase = phase
        return

    def to_xarray(self) -> xr.Dataset:
        return xr.Dataset(
            coords={
                "wavelength": self.wavelength,
            },
            data_vars={
                "power": ("wavelength", self.power),
                "phase": ("wavelength", self.phase),
            },
        )

    @parameter
    def add_noise(self) -> bool:
        """
        Whether to add noise to the control points and PSD spectra based on the standard deviations of the PSD fits, both in the control points and in the PSD itself.
        """
        return self._add_noise

    @add_noise.setter
    def add_noise(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"add_noise must be a boolean value. Got {value} of type {type(value)}.")
        self._add_noise = value


@dataclass(frozen=True, slots=True)
class RealMoonCraterFixed(BasicMoonCraterFixed):
    rim_radius_rng_seed: int | None = None
    """The random seed used to generate the rim radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    rim_height_rng_seed: int | None = None
    """The random seed used to generate the rim elevation PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    floor_radius_rng_seed: int | None = None
    """The random seed used to generate the floor radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    wall_texture_rng_seed: int | None = None
    """The random seed used to generate the wall texture PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    ejecta_texture_rng_seed: int | None = None
    """The random seed used to generate the ejecta texture PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    floor_texture_rng_seed: int | None = None
    """The random seed used to generate the floor texture PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""


class RealMoonCraterVariable(MorphologyCraterVariable):
    def __init__(
        self,
        rim_radius_control: np.ndarray | None = None,
        rim_height_control: np.ndarray | None = None,
        floor_radius_control: np.ndarray | None = None,
        wall_texture_control: np.ndarray | None = None,
        ejecta_texture_control: np.ndarray | None = None,
        rim_radius_psd: PSD1D | None = None,
        floor_radius_psd: PSD1D | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_rim_radius_control", rim_radius_control)
        object.__setattr__(self, "_rim_height_control", rim_height_control)
        object.__setattr__(self, "_floor_radius_control", floor_radius_control)
        object.__setattr__(self, "_wall_texture_control", wall_texture_control)
        object.__setattr__(self, "_ejecta_texture_control", ejecta_texture_control)
        object.__setattr__(self, "_rim_radius_psd", rim_radius_psd)
        object.__setattr__(self, "_floor_radius_psd", floor_radius_psd)
        return

    @property
    def wall_texture_control(self) -> np.ndarray | None:
        """
        The control points for the wall texture PSD.
        """
        return self._wall_texture_control

    @property
    def ejecta_texture_control(self) -> np.ndarray | None:
        """
        The control points for the ejecta texture PSD.
        """
        return self._ejecta_texture_control

    @property
    def rim_radius_psd(self) -> PSD1D:
        return self._rim_radius_psd

    @property
    def floor_radius_psd(self) -> PSD1D:
        return self._floor_radius_psd


@Crater.register("realmooncrater")
class RealMoonCrater(BasicMoonCrater):
    def __init__(self, crater: Crater | None = None, fixed_cls=RealMoonCraterFixed, variable_cls=RealMoonCraterVariable, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        rim_radius_psd: PSD1D | None = None,
        floor_radius_psd: PSD1D | None = None,
        rim_radius_control: dict[str, np.float64] | None = None,
        floor_radius_control: dict[str, np.float64] | None = None,
        **kwargs: Any,
    ) -> RealMoonCrater:
        """
        Initialize a RealMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with parameters used to generate realistic craters as defined in Du et al. (2024)a. [#]_ and Du et al. (2024)b [#]_

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        rim_radius_psd : PSD1D, optional
            Rim radius PSD. If None, then it will be computed from control points.
        floor_radius_psd : PSD1D, optional
            Floor radius PSD. If None, then it will be computed from control points.
        rim_height_psd : PSD1D, optional
            Rim height PSD. If None, then it will be computed from control points.
        **kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        References
        ----------
        .. [#] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2024). Spectral analysis of the morphology of fresh lunar craters I: Rim crest, floor, and rim flank outlines. Journal of Geophysical Research: Planets, 129(11), e2024JE008357. `doi: 10.1029/2024JE008357 <https://doi.org/10.1029/2024JE008357>`_
        .. [#] Du, J., Minton, D.A., Blevins, A.M., Fassett, C.I., Huang, Y.-H., 2025. Spectral Analysis of the Morphology of Fresh Lunar Craters II: Two-Dimensional Surface Elevations of the Continuous Ejecta, Wall, and Floor. Journal of Geophysical Research: Planets 130, e2024JE008890. `doi: 10.1029/2024JE008890 <https://doi.org/10.1029/2024JE008890>`_
        """
        from cratermaker.components.morphology import Morphology
        from cratermaker.utils.montecarlo_utils import bounded_norm, sample_logfit_heteroskedastic, sample_pikefit

        input_args = locals()

        # This is a copy operation, to use old values for any un-specified arguments
        if crater is not None and isinstance(crater, RealMoonCrater):
            rim_radius_psd = crater.rim_radius_psd if rim_radius_psd is None else rim_radius_psd
            rim_radius_control = crater.rim_radius_control if rim_radius_control is None else rim_radius_control
            floor_radius_psd = crater.floor_radius_psd if floor_radius_psd is None else floor_radius_psd
            floor_radius_control = crater.floor_radius_control if floor_radius_control is None else floor_radius_control

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        args = {}

        for var in [
            "rim_radius",
            "floor_radius",
            "rim_height",
            "wall_texture",
            "ejecta_texture",
            "floor_texture",
        ]:
            argname = f"{var}_rng_seed"
            args[argname] = morphology.rng.integers(0, 2**32 - 1)
            argname = f"{var}_control"
            args[argname] = input_args.get(argname)

        if crater.nrings > 0:
            for i, ring in enumerate(crater.rings):
                crater.rings[i] = cls(crater=ring, morphology=morphology, isring=True, **args)

        return cls(
            crater=crater,
            morphology=morphology,
            rim_radius_psd=rim_radius_psd,
            floor_radius_psd=floor_radius_psd,
            **args,
        )

    def rim_radius_profile(self, bearings: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the rim radius profile of the crater based on the rim radius PSD.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater for which to compute the rim radius profile.
        bearings : ArrayLike
            The bearings (in degrees) at which to compute the rim radius profile. This is used to compute the azimuthal variation in the rim radius based on the 2D PSD model.

        Returns
        -------
        rim_radius_profile : NDArray[np.float64]
            The computed rim radius profile at each bearing.
        """
        theta = np.radians(bearings)
        return realmoon_bindings.profile_from_psd(
            crater_radius=self.radius,
            ymean=self.radius,
            psd=self.rim_radius_psd,
            theta=theta,
        )

    def floor_radius_profile(self, bearings: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the floor radius profile of the crater based on the floor radius PSD.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater for which to compute the floor radius profile.
        bearings : ArrayLike
            The bearings (in degrees) at which to compute the floor radius profile. This is used to compute the azimuthal variation in the floor radius based on the 2D PSD model.

        Returns
        -------
        floor_radius_profile : NDArray[np.float64]
            The computed floor radius profile at each bearing.
        """
        theta = np.radians(bearings)
        return realmoon_bindings.profile_from_psd(
            crater_radius=self.crater_radius,
            ymean=self.floor_radius,
            psd=self.floor_radius_psd,
            theta=theta,
        )

    def rim_height_profile(self, bearings: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the rim elevation profile of the crater based on the rim elevation PSD.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater for which to compute the rim elevation profile.
        bearings : ArrayLike
            The bearings (in degrees) at which to compute the rim elevation profile. This is used to compute the azimuthal variation in the rim elevation based on the 2D PSD model.

        Returns
        -------
        rim_height_profile : NDArray[np.float64]
            The computed rim elevation profile at each bearing.
        """
        theta = np.radians(bearings)
        return realmoon_bindings.profile_from_psd(
            crater_radius=self.radius,
            ymean=self.rim_height,
            psd=self.rim_radius_psd,
            theta=theta,
        )

    @property
    def rim_radius_psd(self) -> dict[str, NDArray]:
        """
        The power spectral density distribution of the rim radius outline.
        """
        if self._var._rim_radius_psd is None:
            self._var._rim_radius_psd = PSD1D(
                mean=self.radius,
                pix=self.morphology.surface.pix,
                control_points=self.rim_radius_control,
                rng_seed=self.rim_radius_rng_seed,
                add_noise=self.morphology.add_noise,
            )
        return self._var._rim_radius_psd

    @property
    def floor_radius_psd(self) -> dict[str, NDArray]:
        """
        The power spectral density distribution of the floor radius outline.
        """
        if self._var._floor_radius_psd is None:
            self._var._floor_radius_psd = PSD1D(
                mean=self.floor_radius,
                pix=self.morphology.surface.pix,
                control_points=self.floor_radius_control,
                rng_seed=self.floor_radius_rng_seed,
                add_noise=self.morphology.add_noise,
            )
        return self._var._floor_radius_psd

    @property
    def rim_radius_control(self) -> np.ndarray | None:
        if self._var._rim_radius_control is None:
            self._var._rim_radius_control = self.morphology.get_control_points(
                crater=self, psd1d_coef=self.morphology.psd1d_coef["rim"]
            )
        return self._var._rim_radius_control

    @property
    def floor_radius_control(self) -> np.ndarray | None:
        if self._var._floor_radius_control is None:
            self._var._floor_radius_control = self.morphology.get_control_points(
                crater=self,
                psd1d_coef=self.morphology.psd1d_coef["floor"],
            )
        return self._var._floor_radius_control


@Morphology.register("realmoon")
class RealmoonMorphology(BasicMoonMorphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh. This uses the morphology model of Du et al. 2025a,b.

    This uses the morphology model of Du et al. 2025a,b.

    Parameters
    ----------
    crater : Crater, optional
        The crater object to be converted into a RealMoonCrater. If None, then a new crater is created using the provided parameters.
    fixed_cls : type[RealMoonCraterFixed], optional
        The class definition for the fixed parameters of the RealMoonCrater. Default is RealMoonCraterFixed.
    variable_cls : type[MorphologyCraterVariable], optional
        The class definition for the variable parameters of the RealMoonCrater. Default is MorphologyCraterVariable.
    psd1d_coef_file : str or Path, optional
        The file path for the 1D power spectral density coefficients. If None, then it defaults to the default internal file.
    psd2d_coef_file : str or Path, optional
        The file path for the 2D power spectral density coefficients. If None, then it defaults to the default internal file.
    **kwargs : Any
        |kwargs|

    """

    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=RealMoonCraterFixed,
        variable_cls=RealMoonCraterVariable,
        psd1d_coef_file: str | Path = _PSD1D_COEF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEF_FILE,
        **kwargs,
    ):
        object.__setattr__(self, "_add_noise", None)
        object.__setattr__(self, "_psd1d_coef", None)
        object.__setattr__(self, "_psd2d_coef", None)
        self.add_noise = kwargs.pop("add_noise", True)  # Can be disabled for testing

        psd1d_coef_file = Path(psd1d_coef_file)
        if not psd1d_coef_file.exists():
            raise FileNotFoundError(
                f"1D power spectral density coefficient file not found at {psd1d_coef_file}. Please provide a valid file path."
            )
        psd2d_coef_file = Path(psd2d_coef_file)
        if not psd2d_coef_file.exists():
            raise FileNotFoundError(
                f"2D power spectral density coefficient file not found at {psd2d_coef_file}. Please provide a valid file path."
            )

        self._psd1d_coef = xr.open_dataset(psd1d_coef_file, engine="h5netcdf")
        self._psd2d_coef = xr.open_dataset(psd2d_coef_file, engine="h5netcdf")
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def crater_profile(
        self,
        crater: RealMoonCrater,
        radial_distances: ArrayLike,
        bearings: ArrayLike,
        reference_elevations: ArrayLike | None = None,
        crater_cls: type[Crater] = RealMoonCrater,
        profile_func: Callable = realmoon_bindings.realmoon_profile,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater object containing the parameters for the crater profile.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings : ArrayLike
            Bearings (in degrees) corresponding to the radial distances.
        reference_elevations : ArrayLike, optional
            Reference elevation values to be modified by the crater profile.
        crater_cls : type[Crater], optional
            The class of the crater type used. If the crater object doesn't match, then it is cast as this type. Default is BasicMoonCrater.
        profile_func: Callable, optional
            The backend function used to draw the crater profile. Default is bhe basicmoon_profile from the basicmoon_bindings Rust library.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed crater elevation profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        return super().crater_profile(
            crater=crater,
            radial_distances=radial_distances,
            bearings=bearings,
            reference_elevations=reference_elevations,
            crater_cls=crater_cls,
            profile_func=profile_func,
            **kwargs,
        )

    def ejecta_profile(
        self,
        crater: RealMoonCrater,
        radial_distances: ArrayLike,
        bearings: ArrayLike,
        crater_cls: type[Crater] = RealMoonCrater,
        profile_func: Callable = realmoon_bindings.realmoon_profile,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta profile.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings : ArrayLike
            Bearings (in degrees) corresponding to the radial distances.
        crater_cls : type[Crater], optional
            The class of the crater type used. If the crater object doesn't match, then it is cast as this type. Default is BasicMoonCrater.
        profile_func: Callable, optional
            The backend function used to draw the crater profile. Default is bhe basicmoon_profile from the basicmoon_bindings Rust library.

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed ejecta profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        return super().ejecta_profile(
            crater=crater,
            radial_distances=radial_distances,
            bearings=bearings,
            crater_cls=crater_cls,
            profile_func=profile_func,
            **kwargs,
        )

    @parameter
    def add_noise(self) -> bool:
        """
        Whether to add noise to the control points and PSD spectra based on the standard deviations of the PSD fits, both in the control points and in the PSD itself.
        """
        return self._add_noise

    @add_noise.setter
    def add_noise(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"add_noise must be a boolean value. Got {value} of type {type(value)}.")
        self._add_noise = value

    def get_control_points(self, crater: Crater, psd1d_coef: xr.DataArray):
        """
        Get the control points for the PSD model based on the crater diameter and the provided coefficient and sigma values.

        Parameters
        ----------
        crater : Crater
            The crater for which to compute the control points.
        psd1d_coef : xr.DataArray
            A DataArray containing the coefficients and sigma values for the control points. The expected dimensions are "index" and "term", where "index" corresponds to the different coefficients (e.g., mean and sigma for the control points) and "term" corresponds to the different control points (e.g., "s12", "x2", "y2", etc.).
        add_noise : bool
            Whether to add noise to the control points based on the sigma values in the coef_sigma DataArray. If True, then the control points will be sampled from a normal distribution with mean given by the coefficients and standard deviation given by the sigma values. If False, then the control points will be set to the mean values given by the coefficients without any noise. The default value is True.
        """
        diameter_km = crater.diameter * 1e-3

        control_points = {}
        sigma = {}
        if crater.morphology_type == "ring":
            morphology_type = "complex"
        else:
            morphology_type = crater.morphology_type
        coef = psd1d_coef.sel(morphology_type=morphology_type)
        for term in coef.term:
            c = coef.sel(term=term)
            control_points[str(term.data)] = c.sel(param="m") * diameter_km + c.sel(param="b")
            if self.add_noise:
                sigma = (
                    c.sel(param="sigma") / len(coef.term)
                )  # There are too many correlations in the breakpoints for the control points to be random. Dividing by the number of terms suppresses the noise in the control points a bit so that the PSDs don't become extreme
                cmid = control_points[str(term.data)]
                control_points[str(term.data)] = self.rng.normal(cmid, sigma)

        return control_points

    @property
    def psd1d_coef(self) -> xr.Dataset:
        """
        The coefficients for the 1D PSD models used in the RealmoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd1d_coef

    @property
    def psd2d_coef(self) -> xr.Dataset:
        """
        The coefficients for the 2D PSD models used in the RealmoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd2d_coef

    @property
    def _CraterType(self) -> type[RealMoonCrater]:
        """
        The class definition of the associated Crater type.
        """
        return RealMoonCrater
