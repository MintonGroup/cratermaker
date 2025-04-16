import numpy as np
from cratermaker.plugins.material_catalogue import register_material_catalogue, MaterialCataloguePlugin
from cratermaker.utils.general_utils import create_catalogue

@register_material_catalogue("default")
class DefaultMaterialCatalogue(MaterialCataloguePlugin):
    """Default material catalogue providing properties for celestial bodies."""

    def get_materials(self):
        
        # Define some built-in catalogue values for known solar system materials of interest
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "density" 
        ]
        material_values = [
            ("Water",     2.30,     0.55,   0.0,        1000.0),
            ("Sand",      0.24,     0.41,   0.0,        1750.0),
            ("Dry Soil",  0.24,     0.41,   0.18e6,     1500.0),
            ("Wet Soil",  0.20,     0.55,   1.14e6,     2000.0),
            ("Soft Rock", 0.20,     0.55,   7.60e6,     2250.0),
            ("Hard Rock", 0.20,     0.55,   18.0e6,     2500.0),
            ("Ice",       15.625,   0.48,   0.0,        900.0), 
        ]        
        return create_catalogue(material_properties, material_values)

    
