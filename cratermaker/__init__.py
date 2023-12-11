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

from .core.simulation import Simulation
from .core.target import Target, Material
from .core.crater import Crater, Projectile
from .core.scale import Scale
from .core.surface import Surface, initialize_surface, generate_grid, generate_data, elevation_to_cartesian
from .core.morphology import Morphology
from .core.production import Production, NeukumProductionFunction, R_to_CSFD
from .utils.general_utils import to_config, set_properties, check_properties, create_catalogue, validate_and_convert_location, float_like, normalize_coords
from .utils.montecarlo import get_random_location, get_random_impact_angle, get_random_velocity, get_random_size, bounded_norm
from .perlin import apply_noise
