from __future__ import annotations

import copy
import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray

from cratermaker.bindings import realmoon_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import (
    BasicMoonCrater,
    BasicMoonCraterFixed,
    BasicMoonMorphology,
)
from cratermaker.constants import FloatLike
from cratermaker.core.base import CratermakerBase
from cratermaker.utils.general_utils import parameter

_PSD1D_COEF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_PSD1D_MIN_POINTS = 64
# Profile values below this fraction of the profile mean are treated as missing when computing a PSD from a profile
_PSD1D_LOG_FLOOR = 0.05
# Highest wavenumber constrained by the measured profiles (360 bearings). Must match PSD_DATA_KMAX in realmoon.rs.
_PSD1D_DATA_KMAX = 180
_PSD1D_PARAMETERS = (
    "radius",
    "floor_radius",
    "rim_width",
    "rim_height",
    "wall_curvature",
    "floor_blend",
    "frac_ejrim",
    "ejprofile",
)
# Parameters whose PSD describes the negative of the parameter, so that it is a positive quantity that can be modeled in log
# space like the others (the ejecta profile exponent is negative). Must match the negation in realmoon.rs (at_bearing).
_PSD1D_NEGATED = ("ejprofile",)
# Control points for parameters that have no PSD model in the coefficient file. Their PSDs are given a roughness scale
# (log_scale) of zero, so these values do not matter: the profile is constant at the crater's value.
_PSD1D_CONSTANT_CONTROL_POINTS = {"y1": 0.0, "y2": 0.0, "y3": 0.0, "y4": 0.0, "y5": 0.0, "yn": 0.0, "sn": 2.0}


def _psd_sign(feature: str) -> float:
    """The sign that converts a parameter to the (positive) quantity its PSD describes."""
    return -1.0 if feature in _PSD1D_NEGATED else 1.0


def fill_collapsed_log_profile(profile: ArrayLike) -> tuple[NDArray[np.float64], int]:
    """
    The natural log of a periodic profile, with collapsed values treated as missing.

    Values below _PSD1D_LOG_FLOOR times the mean of the profile are fitter collapses in the measured profiles: isolated
    bearings that are numerically zero or pinned at the fitter's lower bound between neighbors with ordinary values. They are
    replaced by periodic linear interpolation of ln(profile) between the neighboring points, assuming uniform spacing.

    Parameters
    ----------
    profile : ArrayLike
        The profile values, uniformly spaced around a full circle. The mean must be positive.

    Returns
    -------
    logy : NDArray[np.float64]
        ln(profile) with the collapsed values filled in.
    nmissing : int
        The number of values that were filled in.
    """
    y = np.asarray(profile, dtype=np.float64)
    missing = y < _PSD1D_LOG_FLOOR * y.mean()
    logy = np.log(np.where(missing, 1.0, y))
    if missing.any():
        index = np.arange(len(y), dtype=np.float64)
        logy[missing] = np.interp(index[missing], index[~missing], logy[~missing], period=len(y))
    return logy, int(missing.sum())


class PSD1D(CratermakerBase):
    def __init__(
        self,
        normalization_length: FloatLike,
        profile: NDArray[np.float64] | None = None,
        control_points: dict[str, np.float64] | None = None,
        amplitude: NDArray[np.float64] | None = None,
        phase: NDArray[np.float64] | None = None,
        mean: FloatLike | None = None,
        pix: FloatLike | None = None,
        nprofile: int | None = None,
        add_noise: bool = False,
        skewness: FloatLike = 0.0,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        Contains data and methods for building and using a 1D power spectral density model for profile functions.

        Parameters
        ----------
        normalization_length : FloatLike
            The normalization factor used when defining the profile units (typically this is crater radius)
        profile : NDArray[np.float64], optional
            If set, then the PSD model will be generated using the provided linear profile.
        control_points: dict[str, np.float64], optional
            If set, then the PSD model will be generated using a control point model. Dict must contain the following keys: "sn", "yn", "y1", "y2", "y3", "y4", "y5", where "sn" is the log slope of the PSD from yn to y5 and y1...y5, yn are the power values of the first 5 points and the nth (resolution limit) point
        amplitude: NDArray[np.float64], optional
            If set, then the PSD model will ne generated using the provided amplitude values directly. Reuires phase to also be set.
        phase : NDArray[np.float64], optional
            If set, then the PSD model will ne generated using the provided phase values directly. Reuires amplitude to also be set.
        mean : FloaLike, optional
            A mean value of the signal. This is only used for contol points models and will be ignored if either a profile, or amplitude and phase is provided.
        pix : FloatLike, optional
            Defines the limiting resolution of the signal. Ignored when profile input is used.
        nprofile : int, optional
            Defines the number of profile points to define the PSD. Only used for control_points method
        add_noise : bool, optional
            Whether to add Rayleigh noise to the amplitudes of a PSD generated from control points. Default is False.
        skewness : FloatLike, optional
            Target skewness of the log profile synthesized from the PSD. Zero (the default) leaves the phases as they are, which
            for random phases gives a symmetric distribution of log values.
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        # Set a unique rng seed for this object so that the same profile gets generated every time.
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_normalization_length", None)
        object.__setattr__(self, "_mean", None)
        object.__setattr__(self, "_pix", None)
        object.__setattr__(self, "_nprofile", None)
        object.__setattr__(self, "_nfreq", None)
        object.__setattr__(self, "_control_points", None)
        object.__setattr__(self, "_amplitude", None)
        object.__setattr__(self, "_phase", None)
        object.__setattr__(self, "_wavelength", None)
        object.__setattr__(self, "_add_noise", None)
        object.__setattr__(self, "_skewness", 0.0)
        # Defaults that leave the synthesized profile unchanged. RealMoonMorphology sets these from its coefficient file.
        object.__setattr__(self, "log_scale", 1.0)
        object.__setattr__(self, "log_ratio_min", -np.inf)
        object.__setattr__(self, "log_ratio_max", np.inf)
        object.__setattr__(self, "max_radius_fraction", np.inf)

        self.add_noise = add_noise
        self.skewness = skewness
        self.normalization_length = normalization_length
        if profile is not None:
            min_pix = 2 * math.pi * np.abs(profile.mean().item()) / len(profile)
            self.pix = max(np.abs(np.diff(profile)).max() / self.normalization_length, min_pix)
            self.set_from_profile(profile)
        elif (amplitude is not None and phase is None) or (amplitude is None and phase is not None):
            raise ValueError("Both amplitude and phase must be set to use this method of constructing a PSD")
        elif amplitude is not None and phase is not None:
            self.mean = amplitude[0]
            self.amplitude = amplitude
            self.phase = phase
            if pix is None:
                self.pix = 2 * math.pi * self.mean / self.nprofile
            else:
                self.pix = pix / self.normalization_length
        else:
            self.mean = mean / self.normalization_length
            self.pix = pix / self.normalization_length
            self.nprofile = nprofile
            self.control_points = control_points

    def set_from_profile(self, profile):
        """
        Sets the 1D power spectral density (PSD) of a given 1D profile.

        This function will set the following properties:
            wavelength : NDArray[np.float64]
                An array of corresponding to the frequencies of the PSD, and the second column contains the corresponding PSD values. The array is sorted in descending order of wavelength (i.e., ascending order of frequency).
            amplitude : NDArray[np.float64]
                An array of amplitude values corresponding to the frequencies of the input signal. The array is sorted in descending order of wavelength (i.e., ascending order of frequency). The first (DC) element is the arithmetic mean of the profile in units of normalization_length. The remaining (AC) elements are the one-sided Fourier amplitudes of ln(profile), which are dimensionless. Profile values below _PSD1D_LOG_FLOOR times the mean are treated as missing (see `fill_collapsed_log_profile`). Synthesizing a profile from these values exponentiates the log profile and rescales it to the mean, so generated profiles are always positive.
            phases : NDArray[np.float64]
                An array of phase values corresponding to the frequencies of the input signal. The phase values are computed from the angle of the Fourier coefficients and are normalized by the frequency to represent the phase shift in terms of spatial units (e.g., meters). The phase values are sorted in descending order of wavelength (i.e., ascending order of frequency).

        Parameters
        ----------
        profile : ArrayLike
            The input signal for which to compute the 1D PSD. This should be a 1D array of values representing the signal in the spatial domain in its original units. Normalization is done internally using the value of the normalization_length attribute.

        """
        y = np.asarray(profile, dtype=np.float64) / self.normalization_length
        self.nprofile = len(y)
        self.mean = np.mean(y)
        if self.mean <= 0.0:
            raise ValueError("The mean of the profile must be positive to compute its log-space PSD")
        # The AC components describe the spectrum of ln(profile), with collapsed values in the measured profiles filled in
        logy, _ = fill_collapsed_log_profile(y)
        amplitude, phase = realmoon_bindings.psd_from_profile(logy)
        amplitude = np.asarray(amplitude)
        phase = np.asarray(phase)
        # The DC component holds the arithmetic mean of the profile itself, in units of normalization_length
        amplitude[0] = self.mean
        phase[0] = 0.0
        self.amplitude = amplitude[0 : self.nfreq]
        self.phase = phase[0 : self.nfreq]
        return

    def log_sd(self, kmax: int = _PSD1D_DATA_KMAX) -> float:
        """
        Standard deviation of the log profile from the wavenumbers up to kmax, including log_scale.

        By default this includes only the wavenumbers resolved by the measured profiles, so that it is comparable to the
        roughness measured from them.

        Parameters
        ----------
        kmax : int, optional
            Highest wavenumber to include.

        Returns
        -------
        float
        """
        amplitude = np.asarray(self.amplitude, dtype=np.float64)
        nyquist = self.nprofile // 2 if self.nprofile % 2 == 0 else None
        k = np.arange(1, min(kmax, len(amplitude) - 1) + 1)
        variance = np.where(k == nyquist, 1.0, 0.5) * amplitude[k] ** 2
        return float(self.log_scale * np.sqrt(variance.sum()))

    @property
    def skewness(self) -> float:
        """
        Target skewness of the log profile synthesized from the PSD.
        """
        return self._skewness

    @skewness.setter
    def skewness(self, value: FloatLike):
        if not np.isfinite(value):
            raise ValueError("skewness must be finite")
        self._skewness = float(value)

    def to_xarray(self) -> xr.Dataset:
        """
        Exports the amplitude and phase values to an xarray Dataset.
        """
        return xr.Dataset(
            coords={"wavelength": self.wavelength, "bearing": self.bearing},
            data_vars={
                "amplitude": ("wavelength", self.amplitude),
                "phase": ("wavelength", self.phase),
                "power": ("wavelength", self.power),
                "profile": ("bearing", self.profile()),
            },
        )

    @property
    def normalization_length(self) -> np.float64:
        """
        The radius of the crater associated with this profile, which is used for unit normalization.
        """
        return self._normalization_length

    @normalization_length.setter
    def normalization_length(self, value: FloatLike):
        if value > 0.0:
            self._normalization_length = np.float64(value)
        else:
            raise ValueError("normalization_length must be a positive number")

    @property
    def mean(self) -> np.float64:
        """
        The mean value of the profile.
        """
        return self._mean

    @mean.setter
    def mean(self, value):
        self._mean = np.float64(value)

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

    def profile(self, bearing: NDArray[np.float64] | None = None) -> NDArray[np.float64]:
        """
        Generates a 1D profile from the PSD.

        Parameters
        ----------
        bearing : array of floats, optional
            The bearing angles to at which to compute the profile. If not passed, then the default will be `nprofile`` equal polar angle intervals from 0 to 2π.
        """
        if bearing is None:
            bearing = self.bearing
        return realmoon_bindings.profile_from_psd(self, bearing)

    @property
    def bearing(self) -> NDArray[np.float64]:
        """
        Bearing angles of the profile generated by the PSD.
        """
        return np.linspace(0, 2 * np.pi, self.nprofile)

    @property
    def nprofile(self) -> int:
        """
        The number of points in the linear profile that defines the PSD, which is 2x the number of points in the spectral density (nfreq).

        By default this samples the circle of radius normalization_length (the crater radius) at two points per pixel, so
        every parameter of a crater gets a spectrum of the same length regardless of the size of the parameter itself.
        """
        if self._nprofile is None:
            if self._amplitude is not None:
                return 2 * (len(self._amplitude) - 1)
            elif self.pix is not None:
                # pix is in units of normalization_length, so this is 4π * radius / pix
                return max(int(4 * math.pi / self.pix), _PSD1D_MIN_POINTS)
        return self._nprofile

    @nprofile.setter
    def nprofile(self, value):
        if value is not None and (not isinstance(value, int) or value < _PSD1D_MIN_POINTS):
            raise ValueError(f"nprofile must be a positive integer larger than {_PSD1D_MIN_POINTS}")
        else:
            self._nprofile = value

    @property
    def nfreq(self) -> int:
        """
        The number of power spectra points, which is half the number of points in the profile (nprofile).
        """
        if self._nfreq is None:
            return self.nprofile // 2 + 1
        else:
            return self._nfreq

    @nfreq.setter
    def nfreq(self, value):
        if value is not None and (not isinstance(value, int) or value < _PSD1D_MIN_POINTS // 2):
            raise ValueError(f"nfreq must be a positive integer larger than {_PSD1D_MIN_POINTS // 2}")
        else:
            self._nfreq = value

    @property
    def wavelength(self) -> NDArray[np.float64]:
        """
        Wavelength values of the PSD.
        """
        if self._wavelength is None:
            acvals = 2 * np.pi / np.arange(1.0, self.nfreq)
            self._wavelength = np.insert(acvals, 0, np.nan)
        return self._wavelength

    @property
    def amplitude(self) -> NDArray[np.float64]:
        if self._amplitude is None:
            psd_arrs = realmoon_bindings.get_1d_psd_from_control_points(
                control_points=self.control_points,
                mean=self.mean,
                nprofile=self.nprofile,
                add_noise=self.add_noise,
                rng_seed=self.rng_seed,
            )
            return psd_arrs[1]
        else:
            return self._amplitude

    @amplitude.setter
    def amplitude(self, value):
        """
        Amplitude values of the PSD.
        """
        if value is not None:
            value = np.asarray(value)
            if self._nfreq is None:
                self._nfreq = len(value)
            else:
                if len(value) != self.nfreq:
                    raise ValueError(f"Size of amplitude must be {self.nfreq}")
        self._amplitude = value

    @property
    def power(self) -> NDArray[np.float64]:
        """
        Power values of the PSD.
        """
        return self.amplitude**2 / (2 * np.pi)

    @property
    def phase(self) -> NDArray[np.float64]:
        """
        Phase values of the PSD.
        """
        if self._phase is None:
            psd_arrs = realmoon_bindings.get_1d_psd_from_control_points(
                control_points=self.control_points,
                mean=self.mean,
                nprofile=self.nprofile,
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
            if self._nfreq is None:
                self._nfreq = len(value)
            else:
                if len(value) != self.nfreq:
                    raise ValueError(f"Size of phase must be {self.nfreq}")
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

    @parameter
    def add_noise(self) -> bool:
        """
        Whether to add noise to the control points and PSD spectra based on the standard deviations of the PSD fits, both in the control points and in the PSD itself.
        """
        if self._add_noise is None:
            return False
        return self._add_noise

    @add_noise.setter
    def add_noise(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(f"add_noise must be a boolean value. Got {value} of type {type(value)}.")
        self._add_noise = value


def register_psd1d_parameters(*features: str):
    """
    Class decorator to auto-inject _psd, _control, and _profile functions.
    """

    def decorator(cls):
        for name in features:

            def make_control_prop(feat_name=name):
                """
                Defines the [feature]_control property, which defines a 1D profile using control points.
                """
                private_attr = f"_{feat_name}_control"

                @property
                def control_prop(self):
                    if getattr(self._var, private_attr, None) is None:
                        if feat_name in self.morphology.psd1d_coef:
                            val = self.morphology.get_control_points(crater=self, psd1d_coef=self.morphology.psd1d_coef[feat_name])
                        else:
                            val = dict(_PSD1D_CONSTANT_CONTROL_POINTS)
                        setattr(self._var, private_attr, val)
                    return getattr(self._var, private_attr)

                return control_prop

            def make_psd_prop(feat_name=name):
                """
                Defines the [feature]_psd property, which contains the object that is used to store and retrieve PSD information.
                """
                private_attr = f"_{feat_name}_psd"

                @property
                def psd_prop(self):
                    if getattr(self._var, private_attr, None) is None:
                        val = PSD1D(
                            normalization_length=self.radius,
                            mean=_psd_sign(feat_name) * getattr(self, feat_name),
                            pix=self.morphology.psd_pix,
                            control_points=getattr(self, f"{feat_name}_control"),
                            add_noise=self.morphology.add_noise,
                            **{**vars(self.morphology.common_args), **{"rng_seed": self.morphology.rng.integers(0, 2**32 - 1)}},
                        )
                        self.morphology.configure_psd(crater=self, feature=feat_name, psd=val)
                        setattr(self._var, private_attr, val)
                    return getattr(self._var, private_attr)

                return psd_prop

            def make_profile_method(feat_name=name):
                """
                Defines the [feature]_profile method, which returns the profile value at a set of bearing angles.
                """

                def profile_method(self, bearing: ArrayLike) -> NDArray[np.float64]:
                    bearing = np.radians(bearing)
                    return _psd_sign(feat_name) * realmoon_bindings.profile_from_psd(
                        psd=getattr(self, f"{feat_name}_psd"),
                        bearing=bearing,
                    )

                return profile_method

            # Attach them dynamically to the class definition
            setattr(cls, f"{name}_control", make_control_prop())
            setattr(cls, f"{name}_psd", make_psd_prop())
            setattr(cls, f"{name}_profile", make_profile_method())

        return cls

    return decorator


# Keyword arguments that set a crater's PSD data directly
_PSD_DATA_KEYWORDS = tuple(f"{name}_{kind}" for name in _PSD1D_PARAMETERS for kind in ("control", "psd"))


def _copy_psd_data(source: Crater, target: Crater) -> None:
    """
    Copy the PSD control points and PSD objects that `source` has already generated to `target`.

    The data for each feature is generated from the crater radius, the morphology type, and the feature's own value, so it is
    copied only for the features where these are the same in both craters. Data that `source` has not generated yet is left to
    be generated when it is first accessed, so copying a crater does not draw any random numbers.

    Parameters
    ----------
    source : Crater
        The crater being copied.
    target : Crater
        The new crater. Any PSD data it already has is kept.
    """
    if not isinstance(source, RealMoonCrater):
        return
    if source.radius != target.radius or source.morphology_type != target.morphology_type:
        return
    for name in _PSD1D_PARAMETERS:
        if getattr(source, name) != getattr(target, name):
            continue
        for kind in ("control", "psd"):
            attr = f"_{name}_{kind}"
            value = getattr(source._var, attr, None)
            if value is not None and getattr(target._var, attr, None) is None:
                setattr(target._var, attr, copy.deepcopy(value))


@Crater.register("realmooncrater")
@register_psd1d_parameters(*_PSD1D_PARAMETERS)
class RealMoonCrater(BasicMoonCrater):
    def __init__(
        self, crater: Crater | None = None, fixed_cls=BasicMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):
        # PSD data passed explicitly, as <feature>_psd or <feature>_control
        psd_data = {key: kwargs.pop(key) for key in _PSD_DATA_KEYWORDS if kwargs.get(key) is not None}
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        for key, value in psd_data.items():
            setattr(self._var, f"_{key}", value)
        if crater is not None:
            _copy_psd_data(crater, self)

    def __str__(self) -> str:
        str_repr = super().__str__()
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        **kwargs: Any,
    ) -> RealMoonCrater:
        """
        Initialize a RealMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with parameters used to generate realistic craters as defined in Du et al. (2024)a. [#]_ and Du et al. (2024)b [#]_ with updates from Minton et al. (2026) [#]_. Each of the morphometric parameters ``radius``, ``floor_radius``, ``wall_curvature``, ``floor_blend``, ``rim_width``, ``rim_height``, ``frac_ejrim``, ``ejprofile``, has an associated _psd and _control property, and a _profile method that can be accessed or initialized (e.g. ``crater.radius_psd``, ``crater.wall_curvature_profile``, etc.). If the properties are initialized using keyword arguments, then the _psd parameter will be set to a default value using default control points.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        **kwargs : Any
            The keyword arguments provided can be used to set the _control passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        References
        ----------
        .. [#] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2024). Spectral analysis of the morphology of fresh lunar craters I: Rim crest, floor, and rim flank outlines. Journal of Geophysical Research: Planets, 129(11), e2024JE008357. `doi: 10.1029/2024JE008357 <https://doi.org/10.1029/2024JE008357>`_
        .. [#] Du, J., Minton, D.A., Blevins, A.M., Fassett, C.I., Huang, Y.-H., 2025. Spectral Analysis of the Morphology of Fresh Lunar Craters II: Two-Dimensional Surface Elevations of the Continuous Ejecta, Wall, and Floor. Journal of Geophysical Research: Planets 130, e2024JE008890. `doi: 10.1029/2024JE008890 <https://doi.org/10.1029/2024JE008890>`_
        .. [#] Minton, D.A. et al. 2026. Morphometry of Lunar Craters Great and Small. In prep.
        """
        from cratermaker.components.morphology import Morphology

        source = crater
        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        # Make sure rings are the correct type
        for i in range(crater.nrings):
            ring = crater.rings[i]
            if not isinstance(ring, cls):
                ring = cls(crater=ring, morphology=morphology)
                crater.rings[i] = ring

        if not isinstance(crater, cls):
            crater = cls(crater=crater, morphology=morphology, **kwargs)
        # Keep the PSD data the source crater has already generated, where it still applies
        if source is not None:
            _copy_psd_data(source, crater)
        return crater


@Morphology.register("realmoon")
class RealMoonMorphology(BasicMoonMorphology):
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
    psd_pix : float, optional
        Pixel size value used for generating PSDs. If None, then the pix value of the associated morphology.surface object will be used. Default is None.
    **kwargs : Any
        |kwargs|

    """

    # The compiled function that forms a batch of RealMoonCrater objects
    _form_craters_binding = staticmethod(realmoon_bindings.form_craters)

    def _prepare_crater(self, crater: RealMoonCrater) -> None:
        """
        Generate the PSD objects of the crater and its rings.

        Each PSD object draws a seed from the morphology's random number generator on first access. They are accessed in the
        order the compiled `form_craters` reads them (each PSD parameter of the crater, then of each ring), so random numbers
        are drawn in the same order as when the profile is computed directly.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater about to be formed.
        """
        for c in [crater, *(crater.rings or [])]:
            for name in _PSD1D_PARAMETERS:
                getattr(c, f"{name}_psd")

    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=BasicMoonCraterFixed,
        variable_cls=MorphologyCraterVariable,
        psd1d_coef_file: str | Path = _PSD1D_COEF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEF_FILE,
        psd_pix: float | None = None,
        **kwargs,
    ):
        object.__setattr__(self, "_add_noise", None)
        object.__setattr__(self, "_psd1d_coef", None)
        object.__setattr__(self, "_psd2d_coef", None)
        object.__setattr__(self, "_psd_pix", None)
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

        ds = xr.open_dataset(psd1d_coef_file, engine="h5netcdf")
        fit_categories = [str(f) for f in ds.fit_category.data]
        terms = [str(t) for t in ds.term.data]
        params = [str(p) for p in ds.param.data]
        self._psd1d_coef = {}
        for var in ds.data_vars:
            if set(ds[var].dims) != {"term", "fit_category", "param"}:  # Skewness and anchor variables are read below
                continue
            self._psd1d_coef[var] = {}
            for fit_category in fit_categories:
                self._psd1d_coef[var][fit_category] = {}
                for term in terms:
                    self._psd1d_coef[var][fit_category][term] = {}
                    for param in params:
                        self._psd1d_coef[var][fit_category][term][param] = (
                            ds[var].sel(fit_category=fit_category, term=term, param=param).data.item()
                        )
            # Correlations between the control points of a crater
            correlation_var = f"{var}_correlation"
            if correlation_var in ds.data_vars:
                self._psd1d_coef[var]["correlation_terms"] = [str(t) for t in ds.term2.data]
                self._psd1d_coef[var]["correlation"] = {
                    fit_category: ds[correlation_var].sel(fit_category=fit_category).to_numpy().astype(np.float64)
                    for fit_category in fit_categories
                }
            # Model of the skewness of the log profile
            # Guards against extreme values
            if f"{var}_log_ratio_bounds" in ds.data_vars:
                bounds = ds[f"{var}_log_ratio_bounds"]
                self._psd1d_coef[var]["log_ratio_bounds"] = (float(bounds.sel(bound="min")), float(bounds.sel(bound="max")))
            if f"{var}_max_radius_fraction" in ds.data_vars:
                self._psd1d_coef[var]["max_radius_fraction"] = float(ds[f"{var}_max_radius_fraction"])
            skewness_var = f"{var}_skewness"
            if skewness_var in ds.data_vars:
                self._psd1d_coef[var]["skewness"] = {
                    fit_category: {
                        str(param): ds[skewness_var].sel(fit_category=fit_category, skew_param=param).data.item()
                        for param in ds.skew_param.data
                    }
                    for fit_category in fit_categories
                }
            # Measured control points of large craters, used in place of the fit for large non-simple craters
            anchor_var = f"{var}_anchor"
            if anchor_var in ds.data_vars:
                order = np.argsort(ds["anchor_radius"].to_numpy())
                self._psd1d_coef[var]["anchors"] = {
                    "transition_radius": float(ds.attrs["anchor_transition_radius"]),
                    "radius": ds["anchor_radius"].to_numpy().astype(np.float64)[order],
                    "terms": {term: ds[anchor_var].sel(term=term).to_numpy().astype(np.float64)[order] for term in terms},
                }
                anchor_skewness_var = f"{var}_anchor_skewness"
                if anchor_skewness_var in ds.data_vars:
                    self._psd1d_coef[var]["anchors"]["skewness"] = ds[anchor_skewness_var].to_numpy().astype(np.float64)[order]
                anchor_sdln_var = f"{var}_anchor_sdln"
                if anchor_sdln_var in ds.data_vars:
                    self._psd1d_coef[var]["anchors"]["sdln"] = ds[anchor_sdln_var].to_numpy().astype(np.float64)[order]

        self._psd2d_coef = xr.open_dataset(psd2d_coef_file, engine="h5netcdf")
        self.psd_pix = psd_pix
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
            Bearing values (in degrees) corresponding to the radial distances.
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

    def get_control_points(self, crater: Crater, psd1d_coef: dict):
        """
        Get the control points for the PSD model based on the crater diameter and the provided coefficient and sigma values.

        Parameters
        ----------
        crater : Crater
            The crater for which to compute the control points.
        psd1d_coef : dict
            A nested dict containing the coefficients and sigma values for the control points. The expected structure is psd1d_coef[fit_category][term][param], where fit_category is "simple" or "non-simple", term is ["y1","y2","y3","y4","y5","yn","sn"], and param is ["m","b","sigma"] (e.g., mean and sigma for the control points). It may also contain an "anchors" entry; see `control_point_means`.

        Notes
        -----
        If the morphology's add_noise is True, the control points are drawn jointly from a multivariate normal distribution centered on their mean values, with the intrinsic crater-to-crater sigma of each fit and the correlations between the measured control points (scaled by `control_point_noise_scale`). Otherwise the mean values are returned.
        """
        fit_category = "simple" if crater.morphology_type == "simple" else "non-simple"
        radius_m = 0.5 * min(max(crater.diameter, self.valid_diameter_range[0]), self.valid_diameter_range[1])
        control_points = self.control_point_means(psd1d_coef, fit_category, radius_m)
        if self.add_noise:
            noise_scale = self.control_point_noise_scale(psd1d_coef, fit_category, radius_m)
            terms = list(psd1d_coef[fit_category].keys())
            # The intrinsic sigma of the fit reproduces the crater-to-crater scatter of each measured control point
            sigma = np.array([noise_scale * psd1d_coef[fit_category][term]["sigma"] for term in terms])
            correlation = psd1d_coef.get("correlation", {}).get(fit_category)
            if correlation is None:
                z = self.rng.standard_normal(len(terms))
            else:
                # The control points of a crater are correlated (a crater tends to be smooth or rough at all wavelengths)
                order = [psd1d_coef["correlation_terms"].index(term) for term in terms]
                z = self.rng.multivariate_normal(np.zeros(len(terms)), correlation[np.ix_(order, order)], method="cholesky")
            for term, s, zi in zip(terms, sigma, z, strict=True):
                control_points[term] += s * zi

        return control_points

    def configure_psd(self, crater: Crater, feature: str, psd: PSD1D) -> None:
        """
        Apply the morphology's model of the profile shape to a PSD generated from control points.

        This sets the guards against extreme values (log_ratio_min, log_ratio_max, max_radius_fraction), then the roughness
        scale (log_scale), then the target skewness, which depends on the roughness.

        Parameters
        ----------
        crater : Crater
            The crater the profile belongs to.
        feature : str
            The name of the feature (e.g. "rim_width").
        psd : PSD1D
            The PSD of the feature, generated from the crater's control points.
        """
        if feature not in self.psd1d_coef:
            # No model of this feature's variation around the crater: its profile is constant at the crater's value
            psd.log_scale = 0.0
            psd.skewness = 0.0
            return
        coef = self.psd1d_coef[feature]
        psd.log_ratio_min, psd.log_ratio_max = coef.get("log_ratio_bounds", (-np.inf, np.inf))
        psd.max_radius_fraction = coef.get("max_radius_fraction", np.inf)
        psd.log_scale = 1.0
        psd.log_scale = self.get_log_scale(crater=crater, feature=feature, psd=psd)
        psd.skewness = self.get_skewness(crater=crater, feature=feature, psd=psd)

    def get_log_scale(self, crater: Crater, feature: str, psd: PSD1D) -> float:
        """
        The factor applied to the log profile of a crater feature to set its roughness.

        For non-simple craters larger than the anchor transition radius, the roughness (the standard deviation of the log
        profile over the wavenumbers resolved by the measured profiles) is interpolated in ln(radius) from the crater's own
        roughness at the transition radius to the measured roughness of the anchor craters, and held at that of the largest
        one beyond it, consistent with the control points and the skewness. Otherwise the factor is 1.

        Parameters
        ----------
        crater : Crater
            The crater the profile belongs to.
        feature : str
            The name of the feature (e.g. "rim_width").
        psd : PSD1D
            The PSD of the feature, with its amplitudes defined and log_scale equal to 1.

        Returns
        -------
        float
        """
        anchors = self.psd1d_coef[feature].get("anchors")
        if crater.morphology_type == "simple" or anchors is None or "sdln" not in anchors:
            return 1.0
        radius_m = 0.5 * min(max(crater.diameter, self.valid_diameter_range[0]), self.valid_diameter_range[1])
        rt = anchors["transition_radius"]
        own = psd.log_sd()
        if radius_m <= rt or own <= 0.0:
            return 1.0
        lnr = np.concatenate(([np.log(rt)], np.log(anchors["radius"])))
        target = np.exp(np.interp(np.log(radius_m), lnr, np.concatenate(([np.log(own)], np.log(anchors["sdln"])))))
        return float(target / own)

    def get_skewness(self, crater: Crater, feature: str, psd: PSD1D) -> float:
        """
        Draw the target skewness of the log profile of a crater feature.

        The skewness of measured profiles is modeled as a + b * ln(sdln) + c * ln(radius) with normally distributed residuals
        of standard deviation sigma, where sdln is the standard deviation of the log profile (its roughness). The roughness
        is computed from the PSD over the wavenumbers resolved by the measured profiles. As for the control points (see
        `control_point_means` and `control_point_noise_scale`), for non-simple craters larger than the anchor transition
        radius the mean is instead interpolated in ln(radius) from the model at the transition radius to the measured
        skewness of the anchor craters, and the noise is tapered to zero at the largest anchor crater. If add_noise is
        False, the mean value is returned.

        Parameters
        ----------
        crater : Crater
            The crater the profile belongs to.
        feature : str
            The name of the feature (e.g. "rim_height").
        psd : PSD1D
            The PSD of the feature, with its amplitudes defined.

        Returns
        -------
        float
            The target skewness, or 0 if the coefficient file has no skewness model.
        """
        model = self.psd1d_coef[feature].get("skewness")
        if model is None:
            return 0.0
        fit_category = "simple" if crater.morphology_type == "simple" else "non-simple"
        c = model[fit_category]
        sdln = psd.log_sd()
        if sdln <= 0.0:
            return 0.0
        radius_m = 0.5 * min(max(crater.diameter, self.valid_diameter_range[0]), self.valid_diameter_range[1])

        def model_value(r):
            return c["a"] + c["b_ln_sdln"] * np.log(sdln) + c["c_ln_radius"] * np.log(r)

        # As for the control points, large non-simple craters are based on the measured anchor craters instead of the model
        anchors = self.psd1d_coef[feature].get("anchors")
        if (
            fit_category != "simple"
            and anchors is not None
            and "skewness" in anchors
            and radius_m > anchors["transition_radius"]
        ):
            rt = anchors["transition_radius"]
            lnr = np.concatenate(([np.log(rt)], np.log(anchors["radius"])))
            values = np.concatenate(([model_value(rt)], anchors["skewness"]))
            skewness = float(np.interp(np.log(radius_m), lnr, values))
        else:
            skewness = model_value(radius_m)
        if self.add_noise:
            noise_scale = self.control_point_noise_scale(self.psd1d_coef[feature], fit_category, radius_m)
            skewness = self.rng.normal(skewness, noise_scale * c["sigma"])
        return float(skewness)

    @staticmethod
    def control_point_noise_scale(psd1d_coef: dict, fit_category: str, radius: float) -> float:
        """
        The factor applied to the sigma of the control point noise for a crater of a given radius.

        The sigmas are calibrated on the whole population of a fit category, which for non-simple craters is dominated by
        complex craters tens of km across. The few large craters cannot show how much they vary from one another, so across
        the anchor range the noise is tapered linearly in ln(radius), from its full value at the anchor transition radius to
        zero at the largest anchor crater, so that large craters vary about the anchors through the amplitude noise alone.

        Parameters
        ----------
        psd1d_coef : dict
            The coefficients of one feature, as stored in `psd1d_coef[feature]`.
        fit_category : str
            "simple" or "non-simple"
        radius : float
            Crater radius in meters.

        Returns
        -------
        float
            A factor between 0 and 1.
        """
        anchors = psd1d_coef.get("anchors")
        if fit_category == "simple" or anchors is None:
            return 1.0
        lnr0 = np.log(anchors["transition_radius"])
        lnr1 = np.log(anchors["radius"][-1])
        return float(np.clip((lnr1 - np.log(radius)) / (lnr1 - lnr0), 0.0, 1.0))

    @staticmethod
    def control_point_means(psd1d_coef: dict, fit_category: str, radius: float) -> dict[str, float]:
        """
        The mean value of each PSD control point for a crater of a given radius.

        This is the fit function m / ln(radius)^2 + b. There are too few fresh large craters to constrain the fit at the large
        end, so for non-simple craters larger than the anchor transition radius, the values are instead interpolated in
        ln(radius) from the fit at the transition radius to the measured control points of the anchor craters (e.g.
        Schrödinger and Orientale), holding the values of the largest anchor crater beyond it.

        Parameters
        ----------
        psd1d_coef : dict
            The coefficients of one feature, as stored in `psd1d_coef[feature]`.
        fit_category : str
            "simple" or "non-simple"
        radius : float
            Crater radius in meters.

        Returns
        -------
        dict[str, float]
            The mean value of each control point.
        """

        def fit_value(c, r):
            return c["m"] / np.log(r) ** 2 + c["b"]

        coef = psd1d_coef[fit_category]
        anchors = psd1d_coef.get("anchors")
        use_anchors = fit_category != "simple" and anchors is not None and radius > anchors["transition_radius"]
        means = {}
        for term, c in coef.items():
            if use_anchors:
                rt = anchors["transition_radius"]
                lnr = np.concatenate(([np.log(rt)], np.log(anchors["radius"])))
                values = np.concatenate(([fit_value(c, rt)], anchors["terms"][term]))
                means[term] = float(np.interp(np.log(radius), lnr, values))
            else:
                means[term] = fit_value(c, radius)
        return means

    @property
    def psd1d_coef(self) -> dict:
        """
        The coefficients for the 1D PSD models used in the RealMoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd1d_coef

    @property
    def psd2d_coef(self) -> dict:
        """
        The coefficients for the 2D PSD models used in the RealMoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd2d_coef

    @property
    def psd_pix(self) -> float:
        if self._psd_pix is None:
            return self.surface.pix
        else:
            return self._psd_pix

    @psd_pix.setter
    def psd_pix(self, value):
        if not isinstance(value, FloatLike | None):
            raise TypeError(f"psd_pix must be a a float or None. Got {value} of type {type(value)}.")
        self._psd_pix = value

    @property
    def _CraterType(self) -> type[RealMoonCrater]:
        """
        The class definition of the associated Crater type.
        """
        return RealMoonCrater
