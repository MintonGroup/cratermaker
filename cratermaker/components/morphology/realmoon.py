from __future__ import annotations

import math
from dataclasses import InitVar, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray
from scipy import fft

from cratermaker.bindings import basicmoon_bindings, realmoon_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import BasicMoonCrater, BasicMoonCraterFixed, BasicMoonMorphology
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter

_PSD1D_COEF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_PSD1D_MIN_POINTS = 8


@dataclass(frozen=True, slots=True)
class RealMoonCraterFixed(BasicMoonCraterFixed):
    rim_radius_rng_seed: int | None = None
    """The random seed used to generate the rim radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    rim_flank_radius_rng_seed: int | None = None
    """The random seed used to generate the rim flank radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    rim_elevation_rng_seed: int | None = None
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
        rim_elevation_control: np.ndarray | None = None,
        rim_flank_radius_control: np.ndarray | None = None,
        floor_radius_control: np.ndarray | None = None,
        wall_texture_control: np.ndarray | None = None,
        ejecta_texture_control: np.ndarray | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_rim_radius_control", rim_radius_control)
        object.__setattr__(self, "_rim_elevation_control", rim_elevation_control)
        object.__setattr__(self, "_rim_flank_radius_control", rim_flank_radius_control)
        object.__setattr__(self, "_floor_radius_control", floor_radius_control)
        object.__setattr__(self, "_wall_texture_control", wall_texture_control)
        object.__setattr__(self, "_ejecta_texture_control", ejecta_texture_control)
        return

    @property
    def rim_radius_control(self) -> np.ndarray | None:
        """
        The control points for the rim radius PSD.
        """
        return self._rim_radius_control

    @property
    def rim_elevation_control(self) -> np.ndarray | None:
        """
        The control points for the rim elevation PSD.
        """
        return self._rim_elevation_control

    @property
    def rim_flank_radius_control(self) -> np.ndarray | None:
        """
        The control points for the rim flank radius PSD.
        """
        return self._rim_flank_radius_control

    @property
    def floor_radius_control(self) -> np.ndarray | None:
        """
        The control points for the floor radius PSD.
        """
        return self._floor_radius_control

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


@Crater.register("realmooncrater")
class RealMoonCrater(BasicMoonCrater):
    def __init__(self, crater: Crater | None = None, fixed_cls=RealMoonCraterFixed, variable_cls=RealMoonCraterVariable, **kwargs):
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

        for var in [
            "rim_radius",
            "rim_flank_radius",
            "floor_radius",
            "rim_elevation",
            "wall_texture",
            "ejecta_texture",
            "floor_texture",
        ]:
            argname = f"{var}_rng_seed"
            args[argname] = morphology.rng.integers(0, 2**32 - 1)
            argname = f"{var}_control"
            args[argname] = input_args.get(argname)

        kwargs = {**args, **kwargs}

        return cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
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
            phases=None,
            rng_seed=self.rim_radius_rng_seed,
        )

    def rim_flank_radius_profile(self, bearings: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the rim flank radius profile of the crater based on the rim flank radius PSD.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater for which to compute the rim flank radius profile.
        bearings : ArrayLike
            The bearings (in degrees) at which to compute the rim flank radius profile. This is used to compute the azimuthal variation in the rim flank radius based on the 2D PSD model.

        Returns
        -------
        rim_flank_radius_profile : NDArray[np.float64]
            The computed rim flank radius profile at each bearing.
        """
        theta = np.radians(bearings)
        return realmoon_bindings.profile_from_psd(
            crater_radius=self.radius,
            ymean=self.rim_flank_radius,
            psd=self.rim_flank_radius_psd,
            theta=theta,
            phases=None,
            rng_seed=self.rim_flank_radius_rng_seed,
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
            crater_radius=self.floor_radius,
            ymean=self.floor_radius,
            psd=self.floor_radius_psd,
            theta=theta,
            phases=None,
            rng_seed=self.floor_radius_rng_seed,
        )

    def rim_elevation_profile(self, bearings: ArrayLike) -> NDArray[np.float64]:
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
        rim_elevation_profile : NDArray[np.float64]
            The computed rim elevation profile at each bearing.
        """
        theta = np.radians(bearings)
        return realmoon_bindings.profile_from_psd(
            crater_radius=self.radius,
            ymean=self.rim_elevation,
            psd=self.rim_elevation_psd,
            theta=theta,
            phases=None,
            rng_seed=self.rim_elevation_rng_seed,
        )

    @property
    def rim_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim radius outline.
        """
        npoints = max(int(2 * math.pi * self.radius / self.morphology.surface.pix), _PSD1D_MIN_POINTS)
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_radius_control,
            npoints=npoints,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_radius_rng_seed,
        )

    @property
    def rim_flank_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim flank radius outline.
        """
        npoints = int(2 * math.pi * self.rim_flank_radius / self.morphology.surface.pix)
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_flank_radius_control,
            npoints=npoints,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_flank_radius_rng_seed,
        )

    @property
    def floor_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the floor radius outline.
        """
        npoints = max(int(2 * math.pi * self.floor_radius / self.morphology.surface.pix), _PSD1D_MIN_POINTS)
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.floor_radius_control,
            npoints=npoints,
            add_noise=self.morphology.add_noise,
            rng_seed=self.floor_radius_rng_seed,
        )

    @property
    def rim_elevation_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim elevation profile.
        """
        npoints = max(int(2 * math.pi * self.radius / self.morphology.surface.pix), _PSD1D_MIN_POINTS)
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_elevation_control,
            npoints=npoints,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_elevation_rng_seed,
        )

    @property
    def rim_radius_control(self) -> np.ndarray | None:
        if self._var._rim_radius_control is None:
            self._var._rim_radius_control = self.morphology.get_control_points(
                crater=self, coef_sigma=self.morphology.psd1d_coef["rim_radius"]
            )
        return self._var._rim_radius_control

    @property
    def rim_flank_radius_control(self) -> np.ndarray | None:
        if self._var._rim_flank_radius_control is None:
            self._var._rim_flank_radius_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["rim_flank_radius"],
            )
        return self._var._rim_flank_radius_control

    @property
    def floor_radius_control(self) -> np.ndarray | None:
        if self._var._floor_radius_control is None:
            self._var._floor_radius_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["floor_radius"],
            )
        return self._var._floor_radius_control

    @property
    def rim_elevation_control(self) -> np.ndarray | None:
        if self._var._rim_elevation_control is None:
            self._var._rim_elevation_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["rim_elevation"],
            )
        return self._var._rim_elevation_control


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
        variable_cls=MorphologyCraterVariable,
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

        self._psd1d_coef = xr.open_dataset(psd1d_coef_file)
        self._psd2d_coef = xr.open_dataset(psd2d_coef_file)
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def crater_profile(
        self,
        crater: RealMoonCrater,
        r: ArrayLike,
        bearing: ArrayLike,
        r_ref: ArrayLike | None = None,
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater object containing the parameters for the crater profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        bearing : ArrayLike
            Bearings (in degrees) corresponding to the radial distances.
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
        if not isinstance(crater, RealMoonCrater):
            crater = RealMoonCrater.maker(crater, morphology=self)
        if r_ref is None:
            r_ref = np.zeros_like(r)

        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)

        # flatten r to 1D array
        rflat = np.ravel(r)
        r_ref_flat = np.ravel(r_ref)
        bflat = np.ravel(np.radians(bearing))

        elevation = realmoon_bindings.realmoon_profile(
            radial_distances=rflat,
            bearings=bflat,
            reference_elevations=r_ref_flat,
            crater=crater,
            include_crater=True,
            include_ejecta=False,
        )

        elevation = np.reshape(elevation, r.shape)

        return elevation

    def ejecta_profile(self, crater: RealMoonCrater, r: ArrayLike, bearing: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        bearing : ArrayLike
            Bearings (in degrees) corresponding to the radial distances.

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed ejecta profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)
        # flatten r to 1D array
        rflat = np.ravel(r)
        bflat = np.ravel(np.radians(bearing))
        elevation = realmoon_bindings.realmoon_profile(
            radial_distances=rflat,
            bearings=bflat,
            reference_elevations=np.zeros_like(rflat),
            crater=crater,
            include_crater=False,
            include_ejecta=True,
        )
        elevation = np.array(elevation, dtype=np.float64)
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)
        return elevation

    # These will use the BasicMoon profiles for inversions
    def _profile_invert_ejecta(self, r, crater, minimum_thickness):
        ans = super().ejecta_profile(crater, r) - minimum_thickness
        return ans[0]

    def _profile_invert_crater(self, r, crater, minimum_thickness):
        ans = super().crater_profile(crater, r, np.zeros(1)) - minimum_thickness
        return ans[0]

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

    def get_control_points(self, crater: Crater, coef_sigma: xr.DataArray):
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
        if self.add_noise:
            for term in coef_sigma.term:
                sigma = coef_sigma.sel(index=index + 2, term=term)
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
