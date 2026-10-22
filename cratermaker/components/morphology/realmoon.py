from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy import fft

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
_PSD1D_PARAMETERS = (
    "radius",
    "floor_radius",
    "rim_width",
    "rim_height",
)


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
        object.__setattr__(self, "_control_points", None)
        object.__setattr__(self, "_amplitude", None)
        object.__setattr__(self, "_phase", None)
        object.__setattr__(self, "_wavelength", None)
        object.__setattr__(self, "_add_noise", None)

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
                An array of amplitude values corresponding to the frequencies of the input signal. The array is sorted in descending order of wavelength (i.e., ascending order of frequency). The PSD values are normalized by the interval and the number of points in the input signal, so amplitude is in units such that the mean=1.0. The normalization is done by multiplying the squared magnitude of the Fourier coefficients by 2 and dividing by the product of the interval and the number of points in the input signal. The factor of 2 accounts for the fact that we are using a one-sided PSD (i.e., only considering positive frequencies).
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
        amplitude, phase = realmoon_bindings.psd_from_profile(y)
        self.amplitude = amplitude[0 : self.nfreq]
        self.phase = phase[0 : self.nfreq]
        return

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
        """
        if self._nprofile is None:
            if self._amplitude is not None:
                return 2 * (len(self._amplitude) - 1)
            elif self.pix is not None:
                return max(int(4 * math.pi * self.mean / self.pix), _PSD1D_MIN_POINTS)
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
        return self.nprofile // 2 + 1

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
            if self._nprofile is None:
                self.nprofile = len(value) * 2
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
                        val = self.morphology.get_control_points(crater=self, psd1d_coef=self.morphology.psd1d_coef[feat_name])
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
                            mean=getattr(self, feat_name),
                            pix=self.morphology.surface.pix,
                            control_points=getattr(self, f"{feat_name}_control"),
                            add_noise=self.morphology.add_noise,
                            **{**vars(self.morphology.common_args), **{"rng_seed": self.morphology.rng.integers(0, 2**32 - 1)}},
                        )
                        setattr(self._var, private_attr, val)
                    return getattr(self._var, private_attr)

                return psd_prop

            def make_profile_method(feat_name=name):
                """
                Defines the [feature]_profile method, which returns the profile value at a set of bearing angles.
                """

                def profile_method(self, bearing: ArrayLike) -> NDArray[np.float64]:
                    bearing = np.radians(bearing)
                    return realmoon_bindings.profile_from_psd(
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


@Crater.register("realmooncrater")
@register_psd1d_parameters(*_PSD1D_PARAMETERS)
class RealMoonCrater(BasicMoonCrater):
    def __init__(
        self, crater: Crater | None = None, fixed_cls=BasicMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)

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

        # This is a copy operation, to use old values for any un-specified arguments
        if crater is not None and isinstance(crater, RealMoonCrater):
            for param in _PSD1D_PARAMETERS:
                var = f"{param}_psd"
                if var not in kwargs and hasattr(crater, var):
                    kwargs[var] = getattr(crater, var)

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        # Make sure rings are the correct type
        for i in range(crater.nrings):
            ring = crater.rings[i]
            if not isinstance(ring, cls):
                ring = cls(crater=ring, morphology=morphology)
                crater.rings[i] = ring

        return cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )


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
    **kwargs : Any
        |kwargs|

    """

    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=BasicMoonCraterFixed,
        variable_cls=MorphologyCraterVariable,
        psd1d_coef_file: str | Path = _PSD1D_COEF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEF_FILE,
        **kwargs,
    ):
        object.__setattr__(self, "_add_noise", None)
        object.__setattr__(self, "_psd1d_coef", None)
        object.__setattr__(self, "_psd2d_coef", None)
        self.add_noise = kwargs.pop("add_noise", False)  # Can be disabled for testing

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
        bearing: ArrayLike,
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
        bearing : ArrayLike
            bearing (in degrees) corresponding to the radial distances.
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
            bearing=bearing,
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
        fit_category = "simple" if crater.morphology_type == "simple" else "non-simple"

        coef = psd1d_coef.sel(fit_category=fit_category)
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
        The coefficients for the 1D PSD models used in the RealMoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd1d_coef

    @property
    def psd2d_coef(self) -> xr.Dataset:
        """
        The coefficients for the 2D PSD models used in the RealMoonMorphology. This is loaded from the specified file during initialization and stored as an attribute.
        """
        return self._psd2d_coef

    @property
    def _CraterType(self) -> type[RealMoonCrater]:
        """
        The class definition of the associated Crater type.
        """
        return RealMoonCrater
