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
from cratermaker.utils.general_utils import format_large_units, parameter

_PSD1D_COEF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_PSD1D_NUM_POINTS = 5000  # Number of points used in the construction of the 1D PSD


@dataclass(frozen=True, slots=True)
class RealMoonCraterFixed(BasicMoonCraterFixed):
    rim_radius_psd_seed: int | None = None
    """The random seed used to generate the rim radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    rim_flank_radius_psd_seed: int | None = None
    """The random seed used to generate the rim flank radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    rim_elevation_psd_seed: int | None = None
    """The random seed used to generate the rim elevation PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    floor_radius_psd_seed: int | None = None
    """The random seed used to generate the floor radius PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    wall_texture_psd_seed: int | None = None
    """The random seed used to generate the wall texture PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    ejecta_texture_psd_seed: int | None = None
    """The random seed used to generate the ejecta texture PSD so that they can be computed on the fly from the control points without having to store the full PSD in memory."""
    floor_texture_psd_seed: int | None = None
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
        object.__setattr__(self, "_rim_radius_control", None)
        object.__setattr__(self, "_rim_elevation_control", None)
        object.__setattr__(self, "_rim_flank_radius_control", None)
        object.__setattr__(self, "_floor_radius_control", None)
        object.__setattr__(self, "_wall_texture_control", None)
        object.__setattr__(self, "_ejecta_texture_control", None)
        return

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the crater variable properties.
        """
        dict_repr = super().as_dict()
        keys = (
            "rim_radius_control",
            "rim_elevation_control",
            "rim_flank_radius_control",
            "floor_radius_control",
            "wall_texture_control",
            "ejecta_texture_control",
        )
        for key in keys:
            dict_repr[key] = getattr(self, key)

        return dict_repr

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

        for var in morphology.psd1d_coef.data_vars:
            argname = f"{var}_psd_seed"
            args[argname] = morphology.rng.integers(0, 2**32 - 1)

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
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_radius_control,
            npoints=_PSD1D_NUM_POINTS,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_radius_psd_seed,
        )

    @property
    def rim_flank_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim flank radius outline.
        """
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_flank_radius_control,
            npoints=_PSD1D_NUM_POINTS,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_flank_radius_psd_seed,
        )

    @property
    def floor_radius_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the floor radius outline.
        """
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.floor_radius_control,
            npoints=_PSD1D_NUM_POINTS,
            add_noise=self.morphology.add_noise,
            rng_seed=self.floor_radius_psd_seed,
        )

    @property
    def rim_elevation_psd(self) -> np.ndarray:
        """
        The power spectral density distribution of the rim elevation profile.
        """
        return realmoon_bindings.get_1d_psd_from_control_points(
            control_points=self.rim_elevation_control,
            npoints=_PSD1D_NUM_POINTS,
            add_noise=self.morphology.add_noise,
            rng_seed=self.rim_elevation_psd_seed,
        )

    @property
    def rim_radius_control(self) -> np.ndarray | None:
        if self._var._rim_radius_control is None:
            self._var._rim_radius_control = self.morphology.get_control_points(
                crater=self, coef_sigma=self.morphology.psd1d_coef["rim_radius"], add_noise=self.morphology.add_noise
            )
        return self._var._rim_radius_control

    @property
    def rim_flank_radius_control(self) -> np.ndarray | None:
        if self._var._rim_flank_radius_control is None:
            self._var._rim_flank_radius_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["rim_flank_radius"],
                add_noise=self.morphology.add_noise,
            )
        return self._var._rim_flank_radius_control

    @property
    def floor_radius_control(self) -> np.ndarray | None:
        if self._var._floor_radius_control is None:
            self._var._floor_radius_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["floor_radius"],
                add_noise=self.morphology.add_noise,
            )
        return self._var._floor_radius_control

    @property
    def rim_elevation_control(self) -> np.ndarray | None:
        if self._var._rim_elevation_control is None:
            self._var._rim_elevation_control = self.morphology.get_control_points(
                crater=self,
                coef_sigma=self.morphology.psd1d_coef["rim_elevation"],
                add_noise=self.morphology.add_noise,
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
    add_noise : bool, optional
        Whether to add noise to the control points and PSD spectra based on the standard deviations of the PSD fits (both in the control points
        and in the PSD itself). Default is True.
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
        add_noise: bool = True,
        psd1d_coef_file: str | Path = _PSD1D_COEF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEF_FILE,
        **kwargs,
    ):
        object.__setattr__(self, "_add_noise", None)
        object.__setattr__(self, "_psd1d_coef", None)
        object.__setattr__(self, "_psd2d_coef", None)
        self.add_noise = add_noise

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
        # for i in range(len(rflat)):
        #     tmp_crater = DummyCrater(
        #         crater=crater,
        #         diameter=2 * rim_radius_profile[i],
        #         radius=rim_radius_profile[i],
        #         floor_radius=floor_radius_profile[i],
        #         rim_elevation=rim_elevation_profile[i],
        #     )
        #     elevation[i] = basicmoon_bindings.basicmoon_profile(
        #         radial_distances=rflat[i : i + 1],
        #         reference_elevations=r_ref_flat[i : i + 1],
        #         crater=tmp_crater,
        #         include_crater=True,
        #         include_ejecta=False,
        #     )[0]
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
