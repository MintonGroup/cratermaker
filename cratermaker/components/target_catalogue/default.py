from cratermaker.components.target_catalogue import register_target_catalogue, TargetCataloguePlugin
from cratermaker.utils.general_utils import _create_catalogue

@register_target_catalogue("default")
class DefaultTargetCatalogue(TargetCataloguePlugin):
    """Default target catalogue providing properties for celestial bodies."""

    def get_targets(self):
        
        # Define some built-in catalogue values for known solar system targets of interest
        target_properties = [
            "name",    "radius",   "mass",      "material_name", "transition_scale_type"
        ]
        # The catalogue was created with Swiftest
        target_values = [
            ("Mercury", 2439.40e3, 3.301001e+23, "Soft Rock", "silicate"),
            ("Venus", 6051.84e3, 4.867306e+24, "Hard Rock", "silicate"),
            ("Earth", 6371.01e3, 5.972168e+24, "Wet Soil", "silicate"),
            ("Moon", 1737.53e3, 7.345789e+22, "Soft Rock", "silicate"),
            ("Mars", 3389.92e3, 6.416909e+23, "Soft Rock", "silicate"),
            ("Phobos", 11.17e3, 1.080000e+16, "Soft Rock", "silicate"),
            ("Deimos", 6.30e3, 1.800000e+15, "Soft Rock", "silicate"),
            ("Ceres", 469.70e3, 9.383516e+20, "Ice", "ice"),
            ("Vesta", 262.70e3, 2.590270e+20, "Soft Rock", "silicate"),
            ("Io", 1821.49e3, 8.929649e+22, "Hard Rock", "silicate"),
            ("Europa", 1560.80e3, 4.798574e+22, "Ice", "ice"),
            ("Ganymede", 2631.20e3, 1.481479e+23, "Ice", "ice"),
            ("Callisto", 2410.30e3, 1.075661e+23, "Ice", "ice"),
            ("Titan", 2575.50e3, 1.345181e+23, "Ice", "ice"),
            ("Rhea", 764.50e3, 2.306459e+21, "Ice", "ice"),
            ("Dione", 562.50e3, 1.095486e+21, "Ice", "ice"),
            ("Tethys", 536.30e3, 6.174430e+20, "Ice", "ice"),
            ("Enceladus", 252.30e3, 1.080318e+20, "Ice", "ice"),
            ("Mimas", 198.80e3, 3.750939e+19, "Ice", "ice"),
            ("Ariel", 578.90e3, 1.250019e+21, "Ice", "ice"),
            ("Umbriel", 584.70e3, 1.279535e+21, "Ice", "ice"),
            ("Titania", 788.90e3, 3.338178e+21, "Ice", "ice"),
            ("Oberon", 761.40e3, 3.076577e+21, "Ice", "ice"),
            ("Miranda", 235.70e3, 6.442623e+19, "Ice", "ice"),
            ("Triton", 1352.60e3, 2.140292e+22, "Ice", "ice"),
            ("Charon", 606.00e3, 1.589680e+21, "Ice", "ice"),
            ("Pluto", 1188.30e3, 1.302498e+22, "Ice", "ice"),
            ("Arrokoth", 9.13e3, 7.485000e+14, "Ice", "ice"),
        ]
        
        return _create_catalogue(target_properties, target_values)

    
