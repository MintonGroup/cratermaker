import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.counting import _MIN_FACE_FOR_COUNTING, _N_LAYER, _TALLY_NAME, Counting
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.core.crater import Crater


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

    def tally(self, region: LocalSurface | None = None) -> dict[int:Crater]:
        """
        Tally the craters on the surface using the counting model of Minton et al. (2019) [#]_.

        This method works well for simple craters, but does not work well for complex craters. Use with caution.

        Parameters
        ----------
        region : LocalSurface, optional
            A LocalSurface region to count. If not supplied, then the associated surface property is used.

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
        elif not isinstance(region, LocalSurface):
            raise TypeError(f"Expected a LocalSurface, but got {type(region).__name__}.")

        tally_data = self.surface.uxds[_TALLY_NAME].data[region.face_indices, :]
        unique_ids, counts = np.unique(tally_data, return_counts=True, equal_nan=False)
        counts = counts[unique_ids > 0]
        unique_ids = unique_ids[unique_ids > 0]

        return
