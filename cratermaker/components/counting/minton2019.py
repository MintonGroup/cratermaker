import json
import math
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from cratermaker._cratermaker import counting_bindings
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.counting import _MIN_FACE_FOR_COUNTING, _N_LAYER, _TALLY_ID, Counting
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

    def tally(self, region: LocalSurface | None = None, Kv1: float = 0.17, gamma: float = 2.0) -> dict[int:Crater]:
        """
        Tally the craters on the surface using the counting model of Minton et al. (2019) [#]_.

        This method works well for simple craters, but does not work well for complex craters. Use with caution.

        Parameters
        ----------
        region : LocalSurface, optional
            A LocalSurface region to count. If not supplied, then the associated surface property is used.
        Kv1: float
            The visibility function parameter Kv1 from eq. 7 in Minton et al. (2019).
        gamma: float
            The visibility function parameter gamma from eq. 7 in Minton et al. (2019).

        Returns
        -------
        dict[int:Crater]
            A dictionary of observed craters indexed by their ID.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        if region is None:
            region = self.surface
            id_array = self.surface.crater_id
        elif isinstance(region, LocalSurface):
            id_array = region.crater_id
        else:
            raise TypeError(f"Expected a LocalSurface, but got {type(region).__name__}.")

        unique_ids = np.unique(id_array[id_array > 0])
        remove_ids = []
        for id in unique_ids:
            # Check if we have orphaned crater ids for some reason and remove them
            if id not in self.observed:
                self.remove(id)
                continue

            # Update the crater size measurement before computing the degradation and visibility functions
            crater = self.observed[id]
            crater = self.fit_rim(crater=crater, fit_center=False, fit_ellipse=False)
            Kd = self.measure_degradation_state(crater)
            Kv = self.visibility_function(crater, Kv1=Kv1, gamma=gamma)
            if Kd >= Kv:
                remove_ids.append(id)
            else:
                self.observed[id] = crater  # Save the updated measurements to the observed tally

        if len(remove_ids) > 0:
            print(f"Removing {len(remove_ids)} craters from the tally.")
            for id in remove_ids:
                self.remove(id)

        return

    def measure_degradation_state(self, crater: Crater) -> float:
        """
        Measure the degradation state of a crater using eq. 9 from Minton et al. (2019) [#]_ with a correction factor for complex craters from Riedel et al. (2020) [#]_.

        Parameters
        ----------
        crater : Crater
            The crater to measure.

        Returns
        -------
        float
            The degradation state of the crater.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021
        .. [#] Riedel, C., Minton, D.A., Michael, G., Orgel, C., Bogert, C.H. van der, Hiesinger, H., 2020. Degradation of Small Simple and Large Complex Lunar Craters: Not a Simple Scale Dependence. Journal of Geophysical Research: Planets 125, e2019JE006273. https://doi.org/10.1029/2019JE006273


        """
        from cratermaker.constants import _VSMALL

        a = 0.07
        b = 0.15
        diam_correction = 20e3
        correction_factor = 2.0e-7
        depth = self.measure_crater_depth(crater)
        depth_diam = depth / crater.measured_diameter
        if depth_diam < _VSMALL:
            depth_diam = _VSMALL
        if crater.measured_diameter > diam_correction:
            depth_diam += (crater.measured_diameter - diam_correction) * correction_factor
        K = (a / np.sqrt(depth_diam) - b) * crater.measured_radius**2

        return K

    def visibility_function(self, crater: Crater, Kv1: float = 0.17, gamma: float = 2.0) -> float:
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

        Returns
        -------
        float
            The visibility function value for the crater.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        return Kv1 * crater.radius**gamma
