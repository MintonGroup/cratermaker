from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import BasicMoonCrater, BasicMoonCraterFixed, BasicMoonMorphology
from cratermaker.utils.general_utils import format_large_units, parameter

_PSD1D_COEFF_FILE = Path(__file__).resolve().parent / "psd1d_coeffs.nc"
_PSD2D_COEFF_FILE = Path(__file__).resolve().parent / "psd2d_coeffs.nc"
_BREAKPOINT_DIAMETER_IN_FIT = 20e3  # Diameter where the PSD breakpoints occur


@dataclass(frozen=True, slots=True)
class RealisticMoonCraterFixed(BasicMoonCraterFixed):
    rim_radius_control: dict[str, float] | None = None
    """Control points for the rim crest distance PSD model"""
    rim_flank_radius_control: dict[str, float] | None = None
    """Control points for the rim flank radius PSD model"""
    rim_elevation_control: dict[str, float] | None = None
    """Control points for the rim elevation PSD model"""
    floor_radius_control: dict[str, float] | None = None
    """Control points for the floor radius PSD model"""
    floor_elevation_control: dict[str, float] | None = None
    """Control points for the floor radius PSD model"""
    wall_texture_control: dict[str, float] | None = None
    """Control points for the wall texture PSD model"""
    ejecta_texture_control: dict[str, float] | None = None
    """Control points for the ejecta texture PSD model"""
    floor_texture_control: dict[str, float] | None = None
    """Control points for the floor texture PSD model"""


@Crater.register("realisticmooncrater")
class RealisticMoonCrater(BasicMoonCrater):
    def __init__(
        self, crater: Crater | None = None, fixed_cls=RealisticMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
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
        floor_elevation_control: np.ndarray | None = None,
        wall_texture_control: np.ndarray | None = None,
        ejecta_texture_control: np.ndarray | None = None,
        psd1d_coef_file: str | Path = _PSD1D_COEFF_FILE,
        psd2d_coef_file: str | Path = _PSD2D_COEFF_FILE,
        flag_bp_sigma: bool = True,
        **kwargs: Any,
    ) -> RealisticMoonCrater:
        """
        Initialize a RealisticMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with parameters used to generate realistic craters as defined in Du et al. (2024)a. [#]_ and Du et al. (2024)b [#]_

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        rim_radius_control : np.ndarray, optional
            Control points for the rim crest distance profile. If None, then it is generated using the morphology model.
        rim_flank_radius_control : np.ndarray, optional
            Control points for the rim flank radius profile. If None, then it is generated using the morphology model.
        floor_radius_control : np.ndarray, optional
            Control points for the floor radius profile. If None, then it is generated using the morphology model.
        floor_elevation_control : np.ndarray, optional
            Control points for the floor elevation profile. If None, then it is generated using the morphology model.
        wall_texture_control : np.ndarray, optional
            Control points for the wall texture. If None, then it is generated using the morphology model.
        ejecta_texture_control : np.ndarray, optional
            Control points for the ejecta texture. If None, then it is generated using the morphology model.
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
        if crater is not None and isinstance(crater, RealisticMoonCrater):
            rim_radius_control = crater.rim_radius_control if rim_radius_control is None else rim_radius_control
            rim_flank_radius_control = (
                crater.rim_flank_radius_control if rim_flank_radius_control is None else rim_flank_radius_control
            )
            floor_radius_control = crater.floor_radius_control if floor_radius_control is None else floor_radius_control
            rim_elevation_control = crater.rim_elevation_control if rim_elevation_control is None else rim_elevation_control
            floor_elevation_control = crater.floor_elevation_control if floor_elevation_control is None else floor_elevation_control
            wall_texture_control = crater.wall_texture_control if wall_texture_control is None else wall_texture_control
            ejecta_texture_control = crater.ejecta_texture_control if ejecta_texture_control is None else ejecta_texture_control

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        args = {}

        for var in psd1d_coeff.data_vars:
            argname = f"{var}_control"
            if argname in input_args and input_args[argname] is None:
                args[argname] = morphology.get_control_points(
                    diameter=crater.diameter, coef_sigma=psd1d_coeff[var], flag_bp_sigma=flag_bp_sigma
                )

        kwargs = {**args, **kwargs}

        return cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )


@Morphology.register("realisticmoon")
class RealisticMoonMorphology(BasicMoonMorphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh. This uses the morphology model of Du et al. 2025a,b.

    This uses the morphology model of Du et al. 2025a,b.

    Parameters
    ----------
    **kwargs : Any
        |kwargs|

    """

    def __init__(
        self, crater: Crater | None = None, fixed_cls=RealisticMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def _get_1D_power_spectral_density(self, feature, psd_coef, num_psd_component_effec=100) -> NDArray:
        """
        Construct a 1D power spectral density.

        Coeffcients are from Du et al. (2024) [#]_ and  Du et al. (2025) [#]_.

        Parameters
        ----------
        feature : string
            For a 1D feature, choose from ejecta, rim, and floor
        psd_coef : dict (.json)
            Coeffcients used to constract a 1D power spectral density
        num_psd_component_effec : int
            Only reconstrut the sine waves with wavelengths smaller than 2pi/num_psd_component_effec to improve computational efficiency
        **kwargs : Any
            |kwargs|

        References
        ----------
        .. [#] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2024). Spectral analysis of the morphology of fresh lunar craters I: Rim crest, floor, and rim flank outlines. Journal of Geophysical Research: Planets, 129(11), e2024JE008357. `doi: 10.1029/2024JE008357 <https://doi.org/10.1029/2024JE008357>`_
        .. [#] Du, J., Minton, D.A., Blevins, A.M., Fassett, C.I., Huang, Y.-H., 2025. Spectral Analysis of the Morphology of Fresh Lunar Craters II: Two-Dimensional Surface Elevations of the Continuous Ejecta, Wall, and Floor. Journal of Geophysical Research: Planets 130, e2024JE008890. `doi: 10.1029/2024JE008890 <https://doi.org/10.1029/2024JE008890>`_

        """
        # ------------------------------------------------------------------------------------------------------------------
        num_psd_component = 5000
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.diameter < psd_coef["1D"][feature]["slope_12"]["D_tie"]:
            slope_12 = psd_coef["1D"][feature]["slope_12"]["k1"] * self.crater.diameter + psd_coef["1D"][feature]["slope_12"]["b1"]
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_1"]
        else:
            slope_12 = psd_coef["1D"][feature]["slope_12"]["k2"] * self.crater.diameter + psd_coef["1D"][feature]["slope_12"]["b2"]
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_2"]
        if self.crater.diameter < psd_coef["1D"][feature]["bp2_x"]["D_tie"]:
            bp2_x = psd_coef["1D"][feature]["bp2_x"]["k1"] * self.crater.diameter + psd_coef["1D"][feature]["bp2_x"]["b1"]
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_1"]
        else:
            bp2_x = psd_coef["1D"][feature]["bp2_x"]["k2"] * self.crater.diameter + psd_coef["1D"][feature]["bp2_x"]["b2"]
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_2"]
        if self.crater.diameter < psd_coef["1D"][feature]["bp2_y"]["D_tie"]:
            bp2_y = psd_coef["1D"][feature]["bp2_y"]["k1"] * self.crater.diameter + psd_coef["1D"][feature]["bp2_y"]["b1"]
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_1"]
        else:
            bp2_y = psd_coef["1D"][feature]["bp2_y"]["k2"] * self.crater.diameter + psd_coef["1D"][feature]["bp2_y"]["b2"]
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_2"]
        if self.crater.diameter < psd_coef["1D"][feature]["bp3_y"]["D_tie"]:
            bp3_y = psd_coef["1D"][feature]["bp3_y"]["k1"] * self.crater.diameter + psd_coef["1D"][feature]["bp3_y"]["b1"]
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_1"]
        else:
            bp3_y = psd_coef["1D"][feature]["bp3_y"]["k2"] * self.crater.diameter + psd_coef["1D"][feature]["bp3_y"]["b2"]
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_2"]
        if self.crater.diameter < psd_coef["1D"][feature]["bp4_y"]["D_tie"]:
            bp4_y = psd_coef["1D"][feature]["bp4_y"]["k1"] * self.crater.diameter + psd_coef["1D"][feature]["bp4_y"]["b1"]
            bp4_y_sigma = psd_coef["1D"][feature]["bp4_y"]["sigma_1"]
        else:
            bp4_y = psd_coef["1D"][feature]["bp4_y"]["k2"] * self.crater.diameter + psd_coef["1D"][feature]["bp4_y"]["b2"]
            bp4_y_sigma = psd_coef["1D"][feature]["bp4_y"]["sigma_2"]
        slope_12 = self.rng.normal(slope_12, slope_12_sigma)
        bp2_x = self.rng.normal(bp2_x, bp2_x_sigma)
        bp2_y = self.rng.normal(bp2_y, bp2_y_sigma)
        bp3_y = self.rng.normal(bp3_y, bp3_y_sigma)
        bp4_y = self.rng.normal(bp4_y, bp4_y_sigma)
        psd_sigma = psd_coef["1D"][feature]["psd_sigma"]
        # ------------------------------------------------------------------------------------------------------------------
        bp4_x = math.log10(2 * math.pi)
        bp3_x = math.log10(10**bp4_x / 2)
        interval = 2 * math.pi / num_psd_component
        freq = fft.fftfreq(num_psd_component, interval)
        wavelength = 1 / freq[1 : num_psd_component // 2]
        psd = np.zeros((num_psd_component // 2 - 1, 2))
        psd[:, 0] = wavelength
        # ----------------------------------------------------------------------------------------------------------------------
        bp2_x_index = int(10**bp4_x / 10**bp2_x - 1)
        k_23 = (bp3_y - bp2_y) / (bp3_x - bp2_x)
        b_23 = bp3_y - k_23 * bp3_x
        k_12 = slope_12
        b_12 = bp2_y - k_12 * bp2_x
        # -----------------------------------------------------------------------------------------------------------------------
        psd[0, 1] = 10**bp4_x
        psd[1, 1] = 10**bp3_y
        psd[2 : bp2_x_index + 1, 1] = 10 ** (k_23 * np.log10(psd[2 : bp2_x_index + 1, 0]) + b_23)
        psd[bp2_x_index + 1 :, 1] = 10 ** (k_12 * np.log10(psd[bp2_x_index + 1 :, 0]) + b_12)
        # ----------------------------------------------------------------------------------------------------------------------
        psd = psd[:num_psd_component_effec]
        psd_log = np.log10(psd[:, 1])
        psd_log += self.rng.normal(0, psd_sigma, psd_log.shape)
        psd[:, 1] = 10**psd_log
        psd = np.flipud(psd)
        return psd

    def _get_2D_power_spectral_density(
        self, feature, psd_coef, ejecta_radius_norm, floor_radius_norm, max_effec_freq
    ) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
        """
        Constructs a 2D power spectral density.

        Coeffcients are from [#]_.

        Parameters
        ----------
        feature : string
            For a 2D feature, choose from ejecta, wall, and floor
        psd_coef : dict (.json)
            Coeffcients used to constract a 2D power spectral density
        ejecta_radius_norm : float
            The radius of the continuous ejecta normalized by the crater radius
        floor_radius_norm : float
            The radius of the floor normalized by the crater radius
        max_effec_freq : int
            Only reconstrut the sine waves with frequencies smaller than a maxmium effective frequency to improve computational efficiency

        References
        ----------
        .. [#] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2025). Spectral Analysis of the Morphology of Fresh Lunar Craters II: Two-Dimensional Surface Elevations of the Continuous Ejecta, Wall, and Floor. Journal of Geophysical Research: Planets

        """
        # read the psd_coef outside this function. this is temporary until we find a better way to read the psd_coef
        psd_file = Path(__file__).parent / "psd_coef.json"
        with open(psd_file) as f:
            psd_coef = json.load(f)  # read in from _init_?
        # ------------------------------------------------------------------------------------------------------------------
        if feature == "ejecta":
            max_effec_freq = 40
        if feature == "wall":
            max_effec_freq = 40
        if feature == "floor":
            max_effec_freq = 100
        # ------------------------------------------------------------------------------------------------------------------
        if feature == "ejecta":
            theta_min = 0
            theta_max = 2 * math.pi
            r_min = 1
            r_max = ejecta_radius_norm
            area = (r_max - r_min) * (theta_max - theta_min)
            theta_number = 2000
            r_number = 400
        if feature == "wall":
            theta_min = 0
            theta_max = 2 * math.pi
            r_min = floor_radius_norm
            r_max = 1
            area = (r_max - r_min) * (theta_max - theta_min)
            theta_number = 1000
            r_number = 200
        if feature == "floor":
            r_min = -floor_radius_norm / math.sqrt(2)
            r_max = floor_radius_norm / math.sqrt(2)
            theta_min = r_min
            theta_max = r_max
            area = (r_max - r_min) ** 2
            theta_number = 512
            r_number = 512
        # ------------------------------------------------------------------------------------------------------------------
        side_r = np.linspace(0, r_max - r_min, num=r_number, endpoint=False)
        side_theta = np.linspace(0, theta_max - theta_min, num=theta_number, endpoint=False)
        grid_theta, grid_r = np.meshgrid(side_theta, side_r)
        side_freq_theta = np.fft.fftfreq(theta_number, side_theta[1] - side_theta[0])
        side_freq_r = np.fft.fftfreq(r_number, side_r[1] - side_r[0])
        freq_theta, freq_r = np.meshgrid(np.fft.fftshift(side_freq_theta), np.fft.fftshift(side_freq_r))
        freq_theta_quadrant = freq_theta[r_number // 2 :, theta_number // 2 :]
        freq_r_quadrant = freq_r[r_number // 2 :, theta_number // 2 :]
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.diameter < psd_coef["2D"][feature]["p_max"]["D_tie"]:
            p_max = psd_coef["2D"][feature]["p_max"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["p_max"]["b1"]
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_1"]
        else:
            p_max = psd_coef["2D"][feature]["p_max"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["p_max"]["b2"]
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_2"]
        if self.crater.diameter < psd_coef["2D"][feature]["p_diff"]["D_tie"]:
            p_diff = psd_coef["2D"][feature]["p_diff"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["p_diff"]["b1"]
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_1"]
        else:
            p_diff = psd_coef["2D"][feature]["p_diff"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["p_diff"]["b2"]
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_2"]
        if self.crater.diameter < psd_coef["2D"][feature]["nu_fall"]["D_tie"]:
            nu_fall = psd_coef["2D"][feature]["nu_fall"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["nu_fall"]["b1"]
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_1"]
        else:
            nu_fall = psd_coef["2D"][feature]["nu_fall"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["nu_fall"]["b2"]
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_2"]
        if self.crater.diameter < psd_coef["2D"][feature]["psd_sigma"]["D_tie"]:
            psd_sigma = (
                psd_coef["2D"][feature]["psd_sigma"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["psd_sigma"]["b1"]
            )
            psd_sigma_sigma = psd_coef["2D"][feature]["psd_sigma"]["sigma_1"]
        else:
            psd_sigma = (
                psd_coef["2D"][feature]["psd_sigma"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["psd_sigma"]["b2"]
            )
            psd_sigma_sigma = psd_coef["2D"][feature]["psd_sigma"]["sigma_2"]
        p_max = self.rng.normal(p_max, p_max_sigma)
        p_diff = self.rng.normal(p_diff, p_diff_sigma)
        nu_fall = self.rng.normal(nu_fall, nu_fall_sigma)
        psd_sigma = self.rng.normal(psd_sigma, psd_sigma_sigma)
        # ------------------------------------------------------------------------------------------------------------------
        freq_r_matrix = np.sqrt(freq_theta_quadrant**2 + freq_r_quadrant**2)
        freq_theta_matrix = np.arctan2(freq_r_quadrant, freq_theta_quadrant)
        if feature == "ejecta":
            if self.crater.diameter < psd_coef["2D"][feature]["E_rad"]["D_tie"]:
                E_rad = psd_coef["2D"][feature]["E_rad"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["E_rad"]["b1"]
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_1"]
            else:
                E_rad = psd_coef["2D"][feature]["E_rad"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["E_rad"]["b2"]
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_2"]
            E_rad = self.rng.normal(E_rad, E_rad_sigma)
            psd_log_quadrant = (
                p_diff / np.sqrt(1 - (E_rad * np.sin(freq_theta_matrix)) ** 2) * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1)
                + p_max
            )
        if feature == "wall":
            if self.crater.diameter < psd_coef["2D"][feature]["E_circ"]["D_tie"]:
                E_circ = psd_coef["2D"][feature]["E_circ"]["k1"] * self.crater.diameter + psd_coef["2D"][feature]["E_circ"]["b1"]
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_1"]
            else:
                E_circ = psd_coef["2D"][feature]["E_circ"]["k2"] * self.crater.diameter + psd_coef["2D"][feature]["E_circ"]["b2"]
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_2"]
            E_circ = self.rng.normal(E_circ, E_circ_sigma)
            psd_log_quadrant = (
                p_diff / np.sqrt(1 - (E_circ * np.cos(freq_theta_matrix)) ** 2) * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1)
                + p_max
            )
        if feature == "floor":
            psd_log_quadrant = p_diff * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1) + p_max
        psd_log_quadrant_1st = self.rng.normal(psd_log_quadrant, psd_sigma)
        psd_log_quadrant_2nd = self.rng.normal(psd_log_quadrant, psd_sigma)
        # ------------------------------------------------------------------------------------------------------------------
        psd_log = np.zeros((r_number, theta_number))
        psd_log[r_number // 2 :, theta_number // 2 :] = psd_log_quadrant_1st
        psd_log[r_number // 2 + 1 :, 1 : theta_number // 2] = np.fliplr(psd_log_quadrant_2nd[1:, 1:])
        psd_log[1 : r_number // 2, theta_number // 2 + 1 :] = np.flipud(psd_log_quadrant_2nd[1:, 1:])
        psd_log[1 : r_number // 2, 1 : theta_number // 2] = np.fliplr(np.flipud(psd_log_quadrant_1st[1:, 1:]))
        # ---------------------------------------------------------------------------------------------------------------
        psd_log[1 : r_number // 2, theta_number // 2] = np.flip(psd_log[r_number // 2 + 1 :, theta_number // 2])
        psd_log[r_number // 2, 1 : theta_number // 2] = np.flip(psd_log[r_number // 2, theta_number // 2 + 1 :])
        # ---------------------------------------------------------------------------------------------------------------
        psd_log[0, theta_number // 2 + 1 :] = psd_log[1, theta_number // 2 + 1 :]
        psd_log[r_number // 2 + 1 :, 0] = psd_log[r_number // 2 + 1 :, 1]
        # ---------------------------------------------------------------------------------------------------------------
        psd_log[0, 1 : theta_number // 2] = np.flip(psd_log[0, theta_number // 2 + 1 :])
        psd_log[1 : r_number // 2 :, 0] = np.flip(psd_log[r_number // 2 + 1 :, 0])
        # --------------------------------------------------------------------------------------------------------------
        psd_log[0, 0] = psd_log[1, 1]
        psd_log[0, theta_number // 2] = psd_log[1, theta_number // 2]
        psd_log[r_number // 2, 0] = psd_log[r_number // 2, 1]
        # ------------------------------------------------------------------------------------------------------------------
        psd = 10**psd_log
        phase = self.rng.uniform(-math.pi, math.pi, size=(r_number, theta_number))
        if feature == "ejecta" or feature == "wall":
            psd[:, theta_number // 2] = np.zeros(r_number)
        # ------------------------------------------------------------------------------------------------------------------
        freq_theta_delta = (freq_theta.max() - freq_theta.min()) / theta_number
        freq_r_delta = (freq_r.max() - freq_r.min()) / r_number
        freq_theta_num = int(max_effec_freq / freq_theta_delta)
        freq_r_num = int(max_effec_freq / freq_r_delta)
        freq_theta = freq_theta[
            r_number // 2 - freq_r_num : r_number // 2 + freq_r_num,
            theta_number // 2 - freq_theta_num : theta_number // 2 + freq_theta_num,
        ]
        freq_r = freq_r[
            r_number // 2 - freq_r_num : r_number // 2 + freq_r_num,
            theta_number // 2 - freq_theta_num : theta_number // 2 + freq_theta_num,
        ]
        psd = psd[
            r_number // 2 - freq_r_num : r_number // 2 + freq_r_num,
            theta_number // 2 - freq_theta_num : theta_number // 2 + freq_theta_num,
        ]
        phase = phase[
            r_number // 2 - freq_r_num : r_number // 2 + freq_r_num,
            theta_number // 2 - freq_theta_num : theta_number // 2 + freq_theta_num,
        ]
        # ------------------------------------------------------------------------------------------------------------------
        psd = np.sqrt(psd * area)
        return psd, phase, freq_theta, freq_r, grid_theta, grid_r

    def get_control_points(self, diameter: float, coef_sigma: xr.DataArray, flag_bp_sigma: bool = True):
        # Do pre-processing
        if diameter < _BREAKPOINT_DIAMETER_IN_FIT:
            index = 0
        else:
            index = 3
        diameter_km = diameter * 1e-3

        control_points = {}
        sigma = {}
        for term in coef_sigma.term:
            control_points[str(term.data)] = coef_sigma.sel(index=index, term=term) * diameter_km + coef_sigma.sel(
                index=index + 1, term=term
            )

        # ------------------------------------------------------------------------------------------------------------------
        if flag_bp_sigma:
            for term in coef_sigma.term:
                sigma = coef_sigma.sel(index=index + 2, term=term)
                cmid = control_points[str(term.data)]
                control_points[str(term.data)] = self.rng.normal(cmid, sigma)

        return control_points

    def calculate_target_PSD_from_breakpoint_slope(diameter_temp, control_points, num_vertices):
        slope_12 = control_points[0]
        bp2_x = control_points[1]
        bp2_y = control_points[2]
        bp3_y = control_points[3]
        bp4_y = control_points[4]
        bp4_x = math.log10(2 * math.pi)
        bp3_x = math.log10(10**bp4_x / 2)
        interval = 10**bp4_x / num_vertices
        dfft = fft.rfft(np.ones(num_vertices))
        iend = dfft.size - 1
        freq = fft.fftfreq(num_vertices, interval)
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
        for idx in range(len(psd)):
            if idx == 0:
                psd[idx, 1] = 10**bp4_y
            if idx == 1:
                psd[idx, 1] = 10**bp3_y
            if idx > 1 and idx <= bp2_x_index:
                psd[idx, 1] = 10 ** (k_23 * math.log10(psd[idx, 0]) + b_23)
            if idx > bp2_x_index:
                psd[idx, 1] = 10 ** (k_12 * math.log10(psd[idx, 0]) + b_12)
        # ---------------------------------------------------------------------------------------------------------------------------------
        psd = np.flipud(psd)
        return psd

    def add_noise_to_array(psd_target, noise):
        psd_target_noise = psd_target
        power_target_log = np.log10(psd_target[:, 1])
        power_target_log_noise = np.random.normal(0, noise, power_target_log.shape) + power_target_log
        power_target_noise = 10**power_target_log_noise
        psd_target_noise[:, 1] = power_target_noise
        return psd_target_noise

    @property
    def _CraterType(self) -> type[RealisticMoonCrater]:
        return RealisticMoonCrater
