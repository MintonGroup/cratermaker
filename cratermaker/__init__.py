"""
 Copyright 2023 - David Minton
 This file is part of Cratermaker.
 Cratermaker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Cratermaker is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Cratermaker. 
 If not, see: https://www.gnu.org/licenses. 
"""

# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
import warnings
from ._version import version as __version__
warnings.filterwarnings("ignore",category=FutureWarning,module="xarray")
warnings.filterwarnings("ignore",category=FutureWarning,module="uxarray")

from .core.simulation import Simulation
from .core.target import Target
from .core.crater import Crater
from .core.surface import Surface
from .components.scaling import available_scaling_models, get_scaling_model
from .components.production import available_production_models, get_production_model
from .components.morphology import available_morphology_models, get_morphology_model
from .components.grid import available_grid_types, get_grid_type
from .utils.general_utils import validate_and_convert_location, normalize_coords, R_to_CSFD
from .utils.montecarlo import get_random_location, get_random_location_on_face, get_random_impact_angle, get_random_velocity, get_random_size, bounded_norm
from . import realistic, crater, ejecta
