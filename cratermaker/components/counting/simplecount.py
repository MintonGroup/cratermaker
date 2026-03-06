import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from cratermaker._cratermaker import counting_bindings
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.counting import _MIN_FACE_FOR_COUNTING, _N_LAYER, Counting
from cratermaker.components.crater import Crater
from cratermaker.components.surface import LocalSurface, Surface


@Counting.register("simplecount")
class SimpleCount(Counting):
    """
    A basic crater counting model that uses depth-to-diameter values to estimate degradation states of craters for counting.

    Parameters
    ----------
    surface : Surface | LocalSurface
        The surface or local surface view to be counted.
    **kwargs : Any
        |kwargs|
    """

    def __init__(self, surface, **kwargs: Any):
        super().__init__(surface=surface, **kwargs)
        self._component_name = "simplecount"

    def measure_degradation_state(self, crater: Crater, **kwargs: Any) -> float:
        """
        Measure the degradation state of a crater by using a variation of the depth-to-diameter relationship from Minton et al. (2019) [#]_ and Riedel et al. (2020) [#]_.

        Parameters
        ----------
        crater : Crater
            The crater to measure.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        float
           The estimated degradation state of the crater in m². The crater object will also have its degradation state value updated in place

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021
        .. [#] Riedel, C., Minton, D.A., Michael, G., Orgel, C., Bogert, C.H. van der, Hiesinger, H., 2020. Degradation of Small Simple and Large Complex Lunar Craters: Not a Simple Scale Dependence. Journal of Geophysical Research: Planets 125, e2019JE006273. https://doi.org/10.1029/2019JE006273


        """
        from cratermaker.constants import _VSMALL

        # Recalibrated parameters based on Cratermaker's depth/diam calculations
        a = 0.13984926122036828
        diam_correction = 20e3  # depth/diameter correction transition diameter from Riedel et al. (2020)
        correction_factor = 2.0e-7
        if crater.measured_depth_to_diameter is None:
            return 0.0
        depth_over_depth_orig = crater.measured_depth_to_diameter / crater.depth_to_diameter
        if depth_over_depth_orig < _VSMALL:
            depth_over_depth_orig = _VSMALL
        # if crater.measured_diameter > diam_correction:
        #     depth_over_depth_orig += (crater.measured_diameter - diam_correction) * correction_factor
        K = a * (1.0 / np.sqrt(depth_over_depth_orig) - 1.0) * crater.measured_radius**2
        crater.degradation_state = K

        return K

    def visibility_function(self, crater: Crater, Kv1: float = 0.30, gamma: float = 2.0, **kwargs: Any) -> float:
        """
        Calculate the visibility function for a crater using eq. 7 from Minton et al. (2019) [#]_.

        Parameters
        ----------
        crater : Crater
            The crater to calculate the visibility function for.
        Kv1: float
            The visibility function parameter Kv1 from Minton et al. (2019). Default value is 0.30 based on recalibration with Cratermaker.
        gamma: float
            The visibility function parameter gamma from Minton et al. (2019). Default value is 2.0.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        float
            The visibility function value for the crater in units of m².

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        return Kv1 * crater.measured_radius**gamma
