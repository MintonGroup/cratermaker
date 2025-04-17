from cratermaker.plugins.target_catalogue import register_target_catalogue, TargetCataloguePlugin
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
           ("Mercury", 2439.40, 3.301001e+23, "Soft Rock", "silicate"),
           ("Venus", 6051.84, 4.867306e+24, "Hard Rock", "silicate"),
           ("Earth", 6371.01, 5.972168e+24, "Wet Soil", "silicate"),
           ("Moon", 1737.53, 7.345789e+22, "Soft Rock", "silicate"),
           ("Mars", 3389.92, 6.416909e+23, "Soft Rock", "silicate"),
           ("Phobos", 11.17, 1.080000e+16, "Soft Rock", "silicate"),
           ("Deimos", 6.30, 1.800000e+15, "Soft Rock", "silicate"),
           ("Ceres", 469.70, 9.383516e+20, "Ice", "ice"),
           ("Vesta", 262.70, 2.590270e+20, "Soft Rock", "silicate"),
           ("Io", 1821.49, 8.929649e+22, "Hard Rock", "silicate"),
           ("Europa", 1560.80, 4.798574e+22, "Ice", "ice"),
           ("Ganymede", 2631.20, 1.481479e+23, "Ice", "ice"),
           ("Callisto", 2410.30, 1.075661e+23, "Ice", "ice"),
           ("Titan", 2575.50, 1.345181e+23, "Ice", "ice"),
           ("Rhea", 764.50, 2.306459e+21, "Ice", "ice"),
           ("Dione", 562.50, 1.095486e+21, "Ice", "ice"),
           ("Tethys", 536.30, 6.174430e+20, "Ice", "ice"),
           ("Enceladus", 252.30, 1.080318e+20, "Ice", "ice"),
           ("Mimas", 198.80, 3.750939e+19, "Ice", "ice"),
           ("Ariel", 578.90, 1.250019e+21, "Ice", "ice"),
           ("Umbriel", 584.70, 1.279535e+21, "Ice", "ice"),
           ("Titania", 788.90, 3.338178e+21, "Ice", "ice"),
           ("Oberon", 761.40, 3.076577e+21, "Ice", "ice"),
           ("Miranda", 235.70, 6.442623e+19, "Ice", "ice"),
           ("Triton", 1352.60, 2.140292e+22, "Ice", "ice"),
           ("Charon", 606.00, 1.589680e+21, "Ice", "ice"),
           ("Pluto", 1188.30, 1.302498e+22, "Ice", "ice"),
        ]
        
        return _create_catalogue(target_properties, target_values)

    
