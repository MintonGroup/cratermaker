import numpy as np
from cratermaker.plugins.target_catalogue import register_target_catalogue, TargetCataloguePlugin
from cratermaker.utils.general_utils import create_catalogue

@register_target_catalogue("default")
class DefaultTargetCatalogue(TargetCataloguePlugin):
    """Default target catalogue providing properties for celestial bodies."""

    def get_targets(self):
        
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = np.float64(9.80665) # 1 g in SI units
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "transition_scale_type"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", "silicate"),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", "silicate"),
            ("Earth",   6371.01e3, 1.000 * gEarth, "Wet Soil" , "silicate"),
            ("Moon",    1737.53e3, 0.1657* gEarth, "Soft Rock", "silicate"),
            ("Mars",    3389.92e3, 0.379 * gEarth, "Soft Rock", "silicate"),
            ("Ceres",   469.7e3,   0.029 * gEarth, "Ice"      , "ice"),
            ("Vesta",   262.7e3,   0.025 * gEarth, "Soft Rock", "silicate"),
        ]
        
        return create_catalogue(body_properties, body_values)

    
