import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from cratermaker._cratermaker import counting_bindings
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.counting import _MIN_FACE_FOR_COUNTING, _N_LAYER, _TALLY_VARIABLE_NAME, Counting
from cratermaker.components.crater import Crater
from cratermaker.components.surface import LocalSurface, Surface


@Counting.register("minton2019")
class Minton2019Counting(Counting):
    """
    Minton 2019 crater counting model.

    Parameters
    ----------
    surface : Surface | LocalSurface
        The surface or local surface view to be counted.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(self, surface, **kwargs: Any):
        super().__init__(surface=surface, **kwargs)
        self._component_name = "minton2019"

    def measure_degradation_state(self, crater: Crater, **kwargs: Any) -> Crater:
        """
        Measure the degradation state of a crater by measuring its depth-to-diameter ratio and using eq. 9 from Minton et al. (2019) [#]_ with a correction factor for complex craters from Riedel et al. (2020) [#]_.

        Parameters
        ----------
        crater : Crater
            The crater to measure.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        crater : Crater
            The updated Crater object with the measured degradation state and updated depth properties.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021
        .. [#] Riedel, C., Minton, D.A., Michael, G., Orgel, C., Bogert, C.H. van der, Hiesinger, H., 2020. Degradation of Small Simple and Large Complex Lunar Craters: Not a Simple Scale Dependence. Journal of Geophysical Research: Planets 125, e2019JE006273. https://doi.org/10.1029/2019JE006273


        """
        from cratermaker.constants import _VSMALL

        a = 0.07
        b = 0.15
        diam_correction = 20e3  # depth/diameter correction transition diameter from Riedel et al. (2020)
        correction_factor = 2.0e-7
        crater = self.measure_crater_depth(crater)
        depth = crater.measured_rim_height - crater.measured_floor_depth
        depth_diam = depth / crater.measured_diameter
        if depth_diam < _VSMALL:
            depth_diam = _VSMALL
        if crater.measured_diameter > diam_correction:
            depth_diam += (crater.measured_diameter - diam_correction) * correction_factor
        K = (a / np.sqrt(depth_diam) - b) * crater.measured_radius**2
        crater = Crater.maker(crater, measured_degradation_state=K)

        return crater

    def visibility_function(self, crater: Crater, Kv1: float = 0.17, gamma: float = 2.0, **kwargs: Any) -> float:
        """
        Calculate the visibility function for a crater using eq. 7 from Minton et al. (2019) [#]_.

        Parameters
        ----------
        crater : Crater
            The crater to calculate the visibility function for.
        Kv1: float
            The visibility function parameter Kv1 from Minton et al. (2019).
        gamma: float
            The visibility function parameter gamma from Minton et al. (2019).
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        float
            The visibility function value for the crater.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        return Kv1 * crater.measured_radius**gamma
