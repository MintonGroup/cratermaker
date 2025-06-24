import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy import fft

from cratermaker.components.counting import Counting


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
