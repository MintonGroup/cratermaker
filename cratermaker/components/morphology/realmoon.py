from __future__ import annotations

import math
from dataclasses import InitVar, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray
from scipy import fft

from cratermaker.bindings import morphology_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import BasicMoonCrater, BasicMoonCraterFixed, BasicMoonMorphology
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.utils.general_utils import format_large_units, parameter

_PSD1D_COEFF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEFF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_PSD1D_NUM_POINTS = 5000  # Number of points used in the construction of the 1D PSD


@dataclass(frozen=True, slots=True)
class RealmoonCraterFixed(BasicMoonCraterFixed):
    rim_radius_control: dict[str, float] | None = None
    """Control points for the rim crest distance PSD model"""
    rim_flank_radius_control: dict[str, float] | None = None
    """Control points for the rim flank radius PSD model"""
    rim_elevation_control: dict[str, float] | None = None
    """Control points for the rim elevation PSD model"""
    floor_radius_control: dict[str, float] | None = None
    """Control points for the floor radius PSD model"""
    wall_texture_control: dict[str, float] | None = None
    """Control points for the wall texture PSD model"""
    ejecta_texture_control: dict[str, float] | None = None
    """Control points for the ejecta texture PSD model"""
    floor_texture_control: dict[str, float] | None = None
    """Control points for the floor texture PSD model"""


@Crater.register("realmooncrater")
class RealmoonCrater(BasicMoonCrater):
    def __init__(
        self, crater: Crater | None = None, fixed_cls=RealmoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):

        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += f"Rim flank radius: {format_large_units(self.rim_flank_radius, quantity='length')}\n"
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        rim_radius_control: np.ndarray | None = None,
        rim_elevation_control: np.ndarray | None = None,
        rim_flank_radius_control: np.ndarray | None = None,
        floor_radius_control: np.ndarray | None = None,
        wall_texture_control: np.ndarray | None = None,
        ejecta_texture_control: np.ndarray | None = None,
        psd1d_coef_file: str | Path = _PSD1D_COEFF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEFF_FILE,
        **kwargs: Any,
    ) -> RealmoonCrater:
        """
        Initialize a RealmoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with parameters used to generate realistic craters as defined in Du et al. (2024)a. [#]_ and Du et al. (2024)b [#]_

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        rim_radius_control : np.ndarray, optional
            Control points for the rim crest radius PSD. If None, then it will be computed
        rim_elevation_control : np.ndarray, optional
            Conntrol points for the rim elevation PSD. If None then it will be computed.
        rim_flank_radius_control : np.ndarray, optional
            Control points for the rim flank radius PDS. If None then it will be computed.
        floor_radius_control : np.ndarray, optional
            Control points for the floor radius profile. If None then it will be computed.
        wall_texture_control : np.ndarray, optional
            Control points for the wall texture. If None then it will be computed.
        ejecta_texture_control : np.ndarray, optional
            Control points for the ejecta texture. If None then it will be computed.
        psd1d_coef_file : str or Path, optional
            The file path for the 1D power spectral density coefficients. If None, then it defaults to the default internal file.
        psd2d_coef_file : str or Path, optional
            The file path for the 2D power spectral density coefficients. If None, then it defaults to the default internal file.
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

        psd1d_coeff = xr.open_dataset(psd1d_coef_file)
        psd2d_coeff = xr.open_dataset(psd2d_coef_file)

        # This is a copy operation, to use old values for any un-specified arguments
        if crater is not None and isinstance(crater, RealmoonCrater):
            rim_radius_control = crater.rim_radius_control if rim_radius_control is None else rim_radius_control
            rim_flank_radius_control = (
                crater.rim_flank_radius_control if rim_flank_radius_control is None else rim_flank_radius_control
            )
            floor_radius_control = crater.floor_radius_control if floor_radius_control is None else floor_radius_control
            rim_elevation_control = crater.rim_elevation_control if rim_elevation_control is None else rim_elevation_control
            wall_texture_control = crater.wall_texture_control if wall_texture_control is None else wall_texture_control
            ejecta_texture_control = crater.ejecta_texture_control if ejecta_texture_control is None else ejecta_texture_control

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        args = {}

        for var in psd1d_coeff.data_vars:
            argname = f"{var}_control"
            if argname in input_args and input_args[argname] is None:
                args[argname] = morphology.get_control_points(
                    crater=crater, coef_sigma=psd1d_coeff[var], add_noise=morphology.add_noise
                )

        kwargs = {**args, **kwargs}

        return cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )

    @property
    def rim_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim radius outline.
        """
        return self.morphology.calculate_target_1D_PSD_from_breakpoint_slope(self.rim_radius_control)

    @property
    def rim_flank_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim flank radius outline.
        """
        return self.morphology.calculate_target_1D_PSD_from_breakpoint_slope(self.rim_flank_radius_control)

    @property
    def floor_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the floor radius outline.
        """
        return self.morphology.calculate_target_1D_PSD_from_breakpoint_slope(self.floor_radius_control)

    @property
    def rim_elevation_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim elevation profile.
        """
        return self.morphology.calculate_target_1D_PSD_from_breakpoint_slope(self.rim_elevation_control)


@Morphology.register("realmoon")
class RealmoonMorphology(BasicMoonMorphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh. This uses the morphology model of Du et al. 2025a,b.

    This uses the morphology model of Du et al. 2025a,b.

    Parameters
    ----------
    crater : Crater, optional
        The crater object to be converted into a RealmoonCrater. If None, then a new crater is created using the provided parameters.
    fixed_cls : type[RealmoonCraterFixed], optional
        The class definition for the fixed parameters of the RealmoonCrater. Default is RealmoonCraterFixed.
    variable_cls : type[MorphologyCraterVariable], optional
        The class definition for the variable parameters of the RealmoonCrater. Default is MorphologyCraterVariable.
    add_noise : bool, optional
        Whether to add noise to the control points and PSD spectra based on the standard deviations of the PSD fits (both in the control points
        and in the PSD itself). Default is True.
    **kwargs : Any
        |kwargs|

    """

    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=RealmoonCraterFixed,
        variable_cls=MorphologyCraterVariable,
        add_noise: bool = True,
        **kwargs,
    ):
        object.__setattr__(self, "_add_noise", None)
        self.add_noise = add_noise
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def crater_profile(
        self,
        crater: BasicMoonCrater,
        r: ArrayLike,
        bearing: ArrayLike,
        r_ref: ArrayLike | None = None,
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the crater profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        bearing : ArrayLike
            Bearings (in degrees) corresponding to the radial distances. This is used to compute the azimuthal variation in the crater profile based on the 2D PSD models.
        r_ref : ArrayLike, optional
            Reference elevation values to be modified by the crater profile.
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

        @dataclass(frozen=True, slots=True)
        class DummyCrater:
            """
            A simple dataclass to hold crater morphology parameters for use in the profile function. This is necessary because the morphology_bindings.basicmoon_profile function expects a crater object with attributes corresponding to the crater morphology parameters, and we want to be able to construct this crater object from the optimized parameter vector returned by curve_fit without having to create a full cratermaker.Crater object with all of its associated methods and properties.
            """

            crater: RealmoonCrater
            diameter: float
            radius: float
            rim_elevation: float
            floor_radius: float
            id: np.uint32 = field(default=None, init=False)
            semimajor_axis: float | None = field(default=None, init=False)
            semiminor_axis: float | None = field(default=None, init=False)
            orientation: float | None = field(default=None, init=False)
            transient_diameter: float | None = field(default=None, init=False)
            projectile_diameter: float | None = field(default=None, init=False)
            projectile_velocity: float | None = field(default=None, init=False)
            projectile_angle: float | None = field(default=None, init=False)
            projectile_density: float | None = field(default=None, init=False)
            location: tuple[float, float] | None = field(default=None, init=False)
            morphology_type: str | None = field(default=None, init=False)
            measured_semimajor_axis: float | None = field(default=None, init=False)
            measured_semiminor_axis: float | None = field(default=None, init=False)
            measured_orientation: float | None = field(default=None, init=False)
            measured_diameter: float | None = field(default=None, init=False)
            measured_radius: float | None = field(default=None, init=False)
            measured_location: tuple[float, float] | None = field(default=None, init=False)
            time: float | None = field(default=None, init=False)
            floor_elevation: float | None = field(default=None, init=False)
            wall_curvature: float | None = field(default=None, init=False)
            rim_width: float | None = field(default=None, init=False)
            rimdrop: float | None = field(default=None, init=False)
            ejrim: float | None = field(default=None, init=False)
            ejprofile: float | None = field(default=None, init=False)
            peak_height: float | None = field(default=None, init=False)
            peak_width: float | None = field(default=None, init=False)
            peak_offset: float | None = field(default=None, init=False)

            def __post_init__(self):
                for f in self.__dataclass_fields__.values():
                    if getattr(self, f.name) is None:
                        object.__setattr__(self, f.name, getattr(self.crater, f.name))

        if not isinstance(crater, RealmoonCrater):
            crater = RealmoonCrater.maker(crater, morphology=self)
        if r_ref is None:
            r_ref = np.zeros_like(r)

        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)

        # flatten r to 1D array
        rflat = np.ravel(r)
        r_ref_flat = np.ravel(r_ref)
        bflat = np.ravel(bearing)

        rim_radius_profile = self.profile_from_psd(
            crater_radius=crater.radius, ymean=crater.radius, psd=crater.rim_radius_psd, bearings=bflat
        )
        rim_elevation_profile = self.profile_from_psd(
            crater_radius=crater.radius, ymean=crater.rim_elevation, psd=crater.rim_elevation_psd, bearings=bflat
        )
        floor_radius_profile = self.profile_from_psd(
            crater_radius=crater.floor_radius, ymean=crater.floor_radius, psd=crater.floor_radius_psd, bearings=bflat
        )

        elevation = np.empty_like(rflat, dtype=np.float64)
        # I need to implement a Rust binding function that can pass the profiles in as arrays
        for i in range(len(rflat)):
            tmp_crater = DummyCrater(
                crater=crater,
                diameter=2 * rim_radius_profile[i],
                radius=rim_radius_profile[i],
                floor_radius=floor_radius_profile[i],
                rim_elevation=rim_elevation_profile[i],
            )
            elevation[i] = morphology_bindings.basicmoon_profile(
                radial_distances=rflat[i : i + 1],
                reference_elevations=r_ref_flat[i : i + 1],
                crater=tmp_crater,
                include_crater=True,
                include_ejecta=False,
            )[0]
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)

        return elevation

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

    def get_control_points(self, crater: Crater, coef_sigma: xr.DataArray, add_noise: bool = True):
        """
        Get the control points for the PSD model based on the crater diameter and the provided coefficient and sigma values.

        Parameters
        ----------
        crater : Crater
            The crater for which to compute the control points.
        coef_sigma : xr.DataArray
            A DataArray containing the coefficients and sigma values for the control points. The expected dimensions are "index" and "term", where "index" corresponds to the different coefficients (e.g., mean and sigma for the control points) and "term" corresponds to the different control points (e.g., "Slope_12", "Breakpoint_2_x", "Breakpoint_2_y", "Breakpoint_3_y", "Breakpoint_4_y").
        add_noise : bool
            Whether to add noise to the control points based on the sigma values in the coef_sigma DataArray. If True, then the control points will be sampled from a normal distribution with mean given by the coefficients and standard deviation given by the sigma values. If False, then the control points will be set to the mean values given by the coefficients without any noise. The default value is True.
        """
        # Do pre-processing
        if crater.morphology_type == "simple":
            index = 0
        else:
            index = 3
        diameter_km = crater.diameter * 1e-3

        control_points = {}
        sigma = {}
        for term in coef_sigma.term:
            control_points[str(term.data)] = coef_sigma.sel(index=index, term=term) * diameter_km + coef_sigma.sel(
                index=index + 1, term=term
            )

        # ------------------------------------------------------------------------------------------------------------------
        if add_noise:
            for term in coef_sigma.term:
                sigma = coef_sigma.sel(index=index + 2, term=term)
                cmid = control_points[str(term.data)]
                control_points[str(term.data)] = self.rng.normal(cmid, sigma)

        return control_points

    def calculate_target_1D_PSD_from_breakpoint_slope(
        self, control_points: dict[str, float], npoints: int = _PSD1D_NUM_POINTS
    ) -> np.ndarray:
        """
        Given a set of control points, this will compute the power spectral density distribution of a linear feature.

        Parameters
        ----------
        control_points : dict[str, float]
            A dictionary containing the control points for the PSD model. The expected keys are "Slope_12", "Breakpoint_2_x", "Breakpoint_2_y", "Breakpoint_3_y", "Breakpoint_4_y". The breakpoint x values are in log10(wavelength) space and the breakpoint y values are in log10(PSD) space. The slope is in log-log space.
        npoints : int
            The number of points to use in the construction of the PSD. This determines the frequency resolution of the PSD and should be sufficiently high to capture the features of the PSD. The default value is 5000, which provides a good balance between resolution and computational efficiency for typical crater sizes.
        """
        slope_12 = control_points["Slope_12"]
        bp2_x = control_points["Breakpoint_2_x"]
        bp2_y = control_points["Breakpoint_2_y"]
        bp3_y = control_points["Breakpoint_3_y"]
        bp4_y = control_points["Breakpoint_4_y"]
        bp4_x = math.log10(2 * math.pi)
        bp3_x = math.log10(10**bp4_x / 2)
        interval = 10**bp4_x / npoints
        dfft = fft.rfft(np.ones(npoints))
        iend = dfft.size - 1
        freq = fft.fftfreq(npoints, interval)
        wavelength = 1 / freq[1:iend]
        psd = [[0 for i in range(2)] for i in range(iend - 1)]
        psd = np.array(psd)
        psd = np.float64(psd)
        psd[:, 0] = wavelength
        psd[:, 1] = np.abs(dfft[1:iend])
        # ----------------------------------------------------------------------------------------------------------------------
        for i in range(len(psd)):
            if psd[i, 0] < 10**bp2_x:
                bp2_x_index = i - 1
                break
        k_23 = (bp3_y - bp2_y) / (bp3_x - bp2_x)
        b_23 = bp3_y - k_23 * bp3_x
        k_12 = slope_12
        b_12 = bp2_y - k_12 * bp2_x
        # -----------------------------------------------------------------------------------------------------------------------
        psd[0, 1] = 10**bp4_y
        psd[1, 1] = 10**bp3_y
        psd[2 : bp2_x_index + 1, 1] = 10 ** (k_23 * np.log10(psd[2 : bp2_x_index + 1, 0]) + b_23)
        psd[bp2_x_index + 1 :, 1] = 10 ** (k_12 * np.log10(psd[bp2_x_index + 1 :, 0]) + b_12)
        # ---------------------------------------------------------------------------------------------------------------------------------
        psd = np.flipud(psd)
        if self.add_noise:
            noise_stddev = 0.55
            power_target_log = np.log10(psd[:, 1])
            power_target_log_noise = np.random.normal(0, noise_stddev, power_target_log.shape) + power_target_log
            psd[:, 1] = 10**power_target_log_noise
        return psd

    def profile_from_psd(self, crater_radius: float, ymean: float, psd: np.ndarray, bearings: np.ndarray, phases: None = None):
        """
        Generate a profile based on a given PSD. This is done by summing sinusoidal components with frequencies and amplitudes determined by the PSD.

        Parameters
        ----------
        crater_radius : float
            The radius of the crater for which to generate the profile in meters. This is used to determine the scale of the profile.
        ymean : float
            The mean value of the profile. The generated profile will be centered around this mean value.
        psd : np.ndarray
            A 2D array containing the frequencies and corresponding power values of the PSD. The first
            column should contain the frequencies (or wavelengths) and the second column should contain the power values of the PSD at those frequencies.
        bearings : np.ndarray
            A 1D array containing the bearings (in degrees) at which to compute the profile values. The profile will be computed at these bearings and returned as an array of the same length.
        phases : np.ndarray or None
            An optional 1D array containing the phase values (in degrees) for each frequency component in the PSD. If None, then the phases will be randomly generated from a uniform distribution between 0 and 360 degrees. The length of this array should match the number of frequency components in the PSD (i.e., the number of rows in the psd array).
        """
        period_total = psd[-1, 0]
        npoints = psd.size
        nfreq = psd.shape[0]
        thetavals = np.radians(bearings)
        if phases is None:
            phases = self.rng.uniform(size=nfreq) * period_total
        delta_y = np.zeros(npoints)
        amplitude = np.sqrt(psd[:, 1] * period_total / (npoints**2))
        y_ind = amplitude[:, np.newaxis] * np.sin(2 * np.pi * (1 / psd[:, 0][:, np.newaxis]) * (thetavals + phases[:, np.newaxis]))
        delta_y = np.sum(y_ind, axis=0)

        return delta_y * crater_radius + ymean

    @property
    def _CraterType(self) -> type[RealmoonCrater]:
        """
        The class definition of the associated Crater type.
        """
        return RealmoonCrater
