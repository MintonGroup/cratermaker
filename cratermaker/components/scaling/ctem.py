from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.target import Target
from cratermaker.components.scaling import Scaling
from cratermaker.components.scaling.default import DefaultScaling
from cratermaker.components.projectile import Projectile
@Scaling.register("ctem")
class CTEMScaling(DefaultScaling):
    """
    This is an operations class for computing the scaling relationships between projectiles and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties, 
    as well as determining crater morphology based on size and target propertiesImplements the scaling laws described in Richardson (2009) [1]_ that were implemented in CTEM.
        
    Parameters
    ----------
    target : Target | str, default="Moon"
        The target body for the impact. Can be a Target object or a string representing the target name.
    projectile : Projectile | str, default="asteroids"
        The projectile model for the impact. Can be an Projectile object or a string representing the projectile name.
    material_name : str or None
        Name of the target material composition of the target body to look up from the built-in catalogue. Options include "water", "sand", "dry soil", "wet soil", "soft rock", "hard rock", and "ice".
    K1 : FloatLike, optional
        Variable used in crater scaling (see _[1])
    mu : FloatLike, optional
        Variable used in crater scaling (see _[1])
    Ybar : FloatLike, optional
        The strength of the target material, (Pa)
    density : FloatLike, optional
        Volumentric density of target material, (kg/m^3)
    **kwargs : Any
        Additional keyword arguments.
    Notes
    -----
    - The `target` parameter is required and must be an instance of the `Target` class.
    - The `material_name` parameter is optional. If not provided, it will be retrieved from `target`. Setting it explicitly will override the value in `target`.
    - The `K1`, `mu`, `Ybar`, and `density` parameters are optional. If not provided, they will be retrieved from the material catalogue based on the `material_name`. Setting them explicitly will override the values in the catalogue.
    - The built-in material property values are from Holsapple (1993) [2]_ and Kraus et al. (2011) [3]_. 
    
    References
    ----------
    .. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
    .. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
    .. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016        
    """  

    def __init__(self, 
                 target: Target | str | None = None,
                 projectile: Projectile | str | None = None,
                 material_name: str | None = None,
                 K1: FloatLike | None = None,
                 mu: FloatLike | None = None,
                 Ybar: FloatLike | None = None,
                 density: FloatLike | None = None,
                 **kwargs):
        super().__init__(target=target, projectile=projectile, material_name=material_name, K1=K1, mu=mu, Ybar=Ybar, density=density, **kwargs)
        object.__setattr__(self, "_montecarlo_scaling", False)
        self._compute_simple_to_complex_transition_factors() # Recreate this with MC turned off 
