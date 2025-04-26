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

# This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
# The API change does not affect the functionality of the code, so we can safely ignore the warning
import warnings
from ._version import version as __version__
warnings.filterwarnings("ignore",category=FutureWarning,module="xarray")
warnings.filterwarnings("ignore",category=FutureWarning,module="uxarray")

from .core.simulation import Simulation
from .core.target import Target, make_target
from .core.crater import Crater, make_crater
from .core.surface import Surface
from .components.scaling import ScalingModel, available_scaling_models, get_scaling_model, make_scaling
from .components.production import ProductionModel, available_production_models, get_production_model, make_production
from .components.morphology import MorphologyModel, available_morphology_models, get_morphology_model, make_morphology
from .components.impactor import ImpactorModel, available_impactor_models, get_impactor_model, make_impactor
from .components.grid import GridMaker, available_grid_types, get_grid_type, make_grid
