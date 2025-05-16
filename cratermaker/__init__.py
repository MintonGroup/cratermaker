"""
Copyright 2025 - David Minton
This file is part of Cratermaker.
Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Cratermaker.
If not, see: https://www.gnu.org/licenses.
"""

import warnings

from ._version import version as __version__
from .components.morphology import Morphology
from .components.production import Production
from .components.projectile import Projectile
from .components.scaling import Scaling
from .components.surface import Surface
from .components.target import Target
from .core.crater import Crater
from .core.simulation import Simulation

__all__ = [
    "Crater",
    "Morphology",
    "Production",
    "Projectile",
    "Scaling",
    "Surface",
    "Target",
    "Simulation",
    "__version__",
]


# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
warnings.filterwarnings("ignore", category=FutureWarning, module="xarray")
warnings.filterwarnings("ignore", category=FutureWarning, module="uxarray")
warnings.filterwarnings(
    "ignore", category=DeprecationWarning, module="geopandas._compat"
)
warnings.filterwarnings(
    "ignore", category=RuntimeWarning, message="numpy.ndarray size changed"
)
warnings.filterwarnings(
    "ignore", category=UserWarning, module="uxarray.grid.coordinates"
)
warnings.filterwarnings("ignore", category=DeprecationWarning, module="tqdm.std")
