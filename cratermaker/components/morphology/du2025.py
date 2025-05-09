import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.morphology import Morphology
from cratermaker.components.morphology.simplemoon import SimpleMoon


@Morphology.register("du2025")
class Du2025(SimpleMoon):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh. This uses the morphology model of Du et al. 2025a,b

    Parameters
    ----------
    **kwargs : Any
        Additional keyword arguments to be passed to internal functions.
    """

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)

    def _get_1D_power_spectral_density(
        self, feature, psd_coef, num_psd_component_effec=100
    ) -> NDArray:
        """
        This method construct a 1D power spectral density.
        Coeffcients are from [#]_.

        Parameters
        ----------
        feature : string
            For a 1D feature, choose from ejecta, rim, and floor
        psd_coef : dict (.json)
            Coeffcients used to constract a 1D power spectral density
        num_psd_component_effec : int
            Only reconstrut the sine waves with wavelengths smaller than 2pi/num_psd_component_effec to improve computational efficiency

        References
        ----------
        .. [#] Du, J., Minton, D. A., Blevins, A. M., Fassett, C. I., & Huang, Y. H. (2024). Spectral analysis of the morphology of fresh lunar craters I: Rim crest, floor, and rim flank outlines. Journal of Geophysical Research: Planets, 129(11), e2024JE008357. https://doi.org/10.1029/2024JE008357

        """
        # read the psd_coef outside this function. this is temporary until we find a better way to read the psd_coef
        psd_file = Path(__file__).parent / "psd_coef.json"
        with open(psd_file) as f:
            psd_coef = json.load(f)
        # ------------------------------------------------------------------------------------------------------------------
        num_psd_component = 5000
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.final_diameter < psd_coef["1D"][feature]["slope_12"]["D_tie"]:
            slope_12 = (
                psd_coef["1D"][feature]["slope_12"]["k1"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["slope_12"]["b1"]
            )
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_1"]
        else:
            slope_12 = (
                psd_coef["1D"][feature]["slope_12"]["k2"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["slope_12"]["b2"]
            )
            slope_12_sigma = psd_coef["1D"][feature]["slope_12"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["1D"][feature]["bp2_x"]["D_tie"]:
            bp2_x = (
                psd_coef["1D"][feature]["bp2_x"]["k1"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp2_x"]["b1"]
            )
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_1"]
        else:
            bp2_x = (
                psd_coef["1D"][feature]["bp2_x"]["k2"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp2_x"]["b2"]
            )
            bp2_x_sigma = psd_coef["1D"][feature]["bp2_x"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["1D"][feature]["bp2_y"]["D_tie"]:
            bp2_y = (
                psd_coef["1D"][feature]["bp2_y"]["k1"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp2_y"]["b1"]
            )
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_1"]
        else:
            bp2_y = (
                psd_coef["1D"][feature]["bp2_y"]["k2"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp2_y"]["b2"]
            )
            bp2_y_sigma = psd_coef["1D"][feature]["bp2_y"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["1D"][feature]["bp3_y"]["D_tie"]:
            bp3_y = (
                psd_coef["1D"][feature]["bp3_y"]["k1"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp3_y"]["b1"]
            )
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_1"]
        else:
            bp3_y = (
                psd_coef["1D"][feature]["bp3_y"]["k2"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp3_y"]["b2"]
            )
            bp3_y_sigma = psd_coef["1D"][feature]["bp3_y"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["1D"][feature]["bp4_y"]["D_tie"]:
            bp4_y = (
                psd_coef["1D"][feature]["bp4_y"]["k1"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp4_y"]["b1"]
            )
            bp4_y_sigma = psd_coef["1D"][feature]["bp4_y"]["sigma_1"]
        else:
            bp4_y = (
                psd_coef["1D"][feature]["bp4_y"]["k2"] * self.crater.final_diameter
                + psd_coef["1D"][feature]["bp4_y"]["b2"]
            )
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
        psd[2 : bp2_x_index + 1, 1] = 10 ** (
            k_23 * np.log10(psd[2 : bp2_x_index + 1, 0]) + b_23
        )
        psd[bp2_x_index + 1 :, 1] = 10 ** (
            k_12 * np.log10(psd[bp2_x_index + 1 :, 0]) + b_12
        )
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
        This method constructs a 2D power spectral density.
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
        side_theta = np.linspace(
            0, theta_max - theta_min, num=theta_number, endpoint=False
        )
        grid_theta, grid_r = np.meshgrid(side_theta, side_r)
        side_freq_theta = np.fft.fftfreq(theta_number, side_theta[1] - side_theta[0])
        side_freq_r = np.fft.fftfreq(r_number, side_r[1] - side_r[0])
        freq_theta, freq_r = np.meshgrid(
            np.fft.fftshift(side_freq_theta), np.fft.fftshift(side_freq_r)
        )
        freq_theta_quadrant = freq_theta[r_number // 2 :, theta_number // 2 :]
        freq_r_quadrant = freq_r[r_number // 2 :, theta_number // 2 :]
        # ------------------------------------------------------------------------------------------------------------------
        if self.crater.final_diameter < psd_coef["2D"][feature]["p_max"]["D_tie"]:
            p_max = (
                psd_coef["2D"][feature]["p_max"]["k1"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["p_max"]["b1"]
            )
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_1"]
        else:
            p_max = (
                psd_coef["2D"][feature]["p_max"]["k2"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["p_max"]["b2"]
            )
            p_max_sigma = psd_coef["2D"][feature]["p_max"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["2D"][feature]["p_diff"]["D_tie"]:
            p_diff = (
                psd_coef["2D"][feature]["p_diff"]["k1"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["p_diff"]["b1"]
            )
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_1"]
        else:
            p_diff = (
                psd_coef["2D"][feature]["p_diff"]["k2"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["p_diff"]["b2"]
            )
            p_diff_sigma = psd_coef["2D"][feature]["p_diff"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["2D"][feature]["nu_fall"]["D_tie"]:
            nu_fall = (
                psd_coef["2D"][feature]["nu_fall"]["k1"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["nu_fall"]["b1"]
            )
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_1"]
        else:
            nu_fall = (
                psd_coef["2D"][feature]["nu_fall"]["k2"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["nu_fall"]["b2"]
            )
            nu_fall_sigma = psd_coef["2D"][feature]["nu_fall"]["sigma_2"]
        if self.crater.final_diameter < psd_coef["2D"][feature]["psd_sigma"]["D_tie"]:
            psd_sigma = (
                psd_coef["2D"][feature]["psd_sigma"]["k1"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["psd_sigma"]["b1"]
            )
            psd_sigma_sigma = psd_coef["2D"][feature]["psd_sigma"]["sigma_1"]
        else:
            psd_sigma = (
                psd_coef["2D"][feature]["psd_sigma"]["k2"] * self.crater.final_diameter
                + psd_coef["2D"][feature]["psd_sigma"]["b2"]
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
            if self.crater.final_diameter < psd_coef["2D"][feature]["E_rad"]["D_tie"]:
                E_rad = (
                    psd_coef["2D"][feature]["E_rad"]["k1"] * self.crater.final_diameter
                    + psd_coef["2D"][feature]["E_rad"]["b1"]
                )
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_1"]
            else:
                E_rad = (
                    psd_coef["2D"][feature]["E_rad"]["k2"] * self.crater.final_diameter
                    + psd_coef["2D"][feature]["E_rad"]["b2"]
                )
                E_rad_sigma = psd_coef["2D"][feature]["E_rad"]["sigma_2"]
            E_rad = self.rng.normal(E_rad, E_rad_sigma)
            psd_log_quadrant = (
                p_diff
                / np.sqrt(1 - (E_rad * np.sin(freq_theta_matrix)) ** 2)
                * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1)
                + p_max
            )
        if feature == "wall":
            if self.crater.final_diameter < psd_coef["2D"][feature]["E_circ"]["D_tie"]:
                E_circ = (
                    psd_coef["2D"][feature]["E_circ"]["k1"] * self.crater.final_diameter
                    + psd_coef["2D"][feature]["E_circ"]["b1"]
                )
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_1"]
            else:
                E_circ = (
                    psd_coef["2D"][feature]["E_circ"]["k2"] * self.crater.final_diameter
                    + psd_coef["2D"][feature]["E_circ"]["b2"]
                )
                E_circ_sigma = psd_coef["2D"][feature]["E_circ"]["sigma_2"]
            E_circ = self.rng.normal(E_circ, E_circ_sigma)
            psd_log_quadrant = (
                p_diff
                / np.sqrt(1 - (E_circ * np.cos(freq_theta_matrix)) ** 2)
                * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1)
                + p_max
            )
        if feature == "floor":
            psd_log_quadrant = (
                p_diff * (np.exp(-np.sqrt(freq_r_matrix / nu_fall)) - 1) + p_max
            )
        psd_log_quadrant_1st = self.rng.normal(psd_log_quadrant, psd_sigma)
        psd_log_quadrant_2nd = self.rng.normal(psd_log_quadrant, psd_sigma)
        # ------------------------------------------------------------------------------------------------------------------
        psd_log = np.zeros((r_number, theta_number))
        psd_log[r_number // 2 :, theta_number // 2 :] = psd_log_quadrant_1st
        psd_log[r_number // 2 + 1 :, 1 : theta_number // 2] = np.fliplr(
            psd_log_quadrant_2nd[1:, 1:]
        )
        psd_log[1 : r_number // 2, theta_number // 2 + 1 :] = np.flipud(
            psd_log_quadrant_2nd[1:, 1:]
        )
        psd_log[1 : r_number // 2, 1 : theta_number // 2] = np.fliplr(
            np.flipud(psd_log_quadrant_1st[1:, 1:])
        )
        # ---------------------------------------------------------------------------------------------------------------
        psd_log[1 : r_number // 2, theta_number // 2] = np.flip(
            psd_log[r_number // 2 + 1 :, theta_number // 2]
        )
        psd_log[r_number // 2, 1 : theta_number // 2] = np.flip(
            psd_log[r_number // 2, theta_number // 2 + 1 :]
        )
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
            psd[:, theta_number // 2] = np.zeros((r_number))
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
