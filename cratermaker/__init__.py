"""
Copyright 2026 - David Minton.

This file is part of Cratermaker.

Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  You should have received a copy of the GNU General Public License along with Cratermaker.  If not, see: https://www.gnu.org/licenses.
"""

import warnings

from cratermaker._version import version as __version__
from cratermaker.components.counting import Counting
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology
from cratermaker.components.production import Production
from cratermaker.components.projectile import Projectile
from cratermaker.components.scaling import Scaling
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.components.target import Target
from cratermaker.core.simulation import Simulation
from cratermaker.utils.general_utils import cleanup

_COMPONENT_NAMES = ["crater", "counting", "morphology", "production", "projectile", "scaling", "surface", "target"]

__all__ = ["Simulation", "__version__"] + [n.capitalize() for n in _COMPONENT_NAMES]


# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
warnings.filterwarnings("ignore", category=FutureWarning, module="xarray")
warnings.filterwarnings("ignore", category=FutureWarning, module="uxarray")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="geopandas._compat")
warnings.filterwarnings("ignore", category=RuntimeWarning, message="numpy.ndarray size changed")
warnings.filterwarnings("ignore", category=UserWarning, module="uxarray.grid.coordinates")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="tqdm.std")
