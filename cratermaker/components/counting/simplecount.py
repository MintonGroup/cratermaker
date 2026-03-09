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
        Measure the degradation state of a crater by using a variation of the depth-to-diameter proxy for degradation state from Minton et al. (2019) [#]_ and Riedel et al. (2020) [#]_.

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


        Notes
        -----
        The technique used here is a major improvement over the one used previously in CTEM. Rather than using the depth-to-diameter directly, it instead uses the ratio of the measured depth-to-diameter to the initial computed depth-to-diameter of the crater. Fitting functions for both simple and complex craters were developed using this metric.
        """
        from cratermaker.constants import _VSMALL

        # A simplified transition diameter for the two morphology regimes. A more comprehensive study including other planetary bodies is needed to improve this
        _DTR = 10.0e3

        # The following are empirically-derived constants determined through numerical experiments
        c0 = 0.2771649498066961
        c1 = 0.1032775007149389
        d_scale = 178.40182304735092
        s0 = 0.11511000581387373
        s1 = 0.017078961587273955
        p0 = 2.0700104291102317
        p1 = 0.0920770504148557

        def a_vs_diameter(diameter):
            if diameter < _DTR:
                return s0 + s1 * np.log(diameter)
            else:
                dkm = diameter * 1e-3
                dtrkm = _DTR * 1e-3
                return c0 + c1 * (1.0 - np.exp((dtrkm - dkm) / d_scale))

        def p_vs_diameter(diameter):
            if diameter > _DTR:  # In the fit, we fixed the exponent p to be constant for complex craters
                diameter = _DTR
            return p0 + p1 * np.log(diameter)

        def K_vs_depth_over_depth_orig(depth_over_depth_orig, a, p):
            return a * (depth_over_depth_orig ** (-1 / p) - 1.0) * crater.measured_radius**2

        if crater.measured_depth_to_diameter is None:
            return 0.0
        depth_over_depth_orig = crater.measured_depth_to_diameter / crater.depth_to_diameter
        if depth_over_depth_orig < _VSMALL:
            depth_over_depth_orig = _VSMALL
        a = a_vs_diameter(crater.measured_diameter)
        p = p_vs_diameter(crater.measured_diameter)
        crater.degradation_state = K_vs_depth_over_depth_orig(depth_over_depth_orig, a, p)

        return crater.degradation_state

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
