from cratermaker.components.projectile import Projectile
from cratermaker.components.scaling import Scaling
from cratermaker.components.scaling.montecarlo import MonteCarloScaling
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike


@Scaling.register("ctem")
class CTEMScaling(MonteCarloScaling):
    """
    This is an operations class for computing the scaling relationships between projectiles and craters.

    This class encapsulates the logic for converting between projectile properties and crater properties,
    as well as determining crater morphology based on size and target propertiesImplements the scaling laws described in Richardson (2009) that were implemented in CTEM.

    Parameters
    ----------
    target : Target | str, default="Moon"
        The target body for the impact. Can be a Target object or a string representing the target name.
    projectile : Projectile | str, default="asteroids"
        The projectile model for the impact. Can be an Projectile object or a string representing the projectile name.
    material : str or None
        Name of the target material composition of the target body to look up from the built-in catalogue. Options include "water", "sand", "dry soil", "wet soil", "soft rock", "hard rock", and "ice".
    K1 : FloatLike, optional
        Variable used in crater scaling (see Richardson 2009)
    mu : FloatLike, optional
        Variable used in crater scaling (see Richardson 2009)
    Ybar : FloatLike, optional
        The strength of the target material, (Pa)
    density : FloatLike, optional
        Volumentric density of target material, (kg/m^3)
    **kwargs : Any
        Additional keyword arguments.

    Notes
    -----
    - The `target` parameter is required and must be an instance of the `Target` class.
    - The `material` parameter is optional. If not provided, it will be retrieved from `target`. Setting it explicitly will override the value in `target`.
    - The `K1`, `mu`, `Ybar`, and `density` parameters are optional. If not provided, they will be retrieved from the material catalogue based on the `material`. Setting them explicitly will override the values in the catalogue.
    - The built-in material property values are from Holsapple (1993) and Kraus et al. (2011).
    - Complex craater scaling parameters are a synthesis of Pike (1980), Croft (1985), and Schenk et al. (2004).

    References
    ----------
    - Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
    - Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
    - Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016
    - Pike, R.J., 1980. Control of crater morphology by gravity and target type - Mars, earth, moon. In: Lunar and Planetary Science Conference 11, 2159-2189.
    - Croft, S.K., 1985. The scaling of complex craters. Proceedings of the Fifteenth Lunar and Planetary Science Conference, Part 2 Journal of Geophysical Research 90, Supplement, C828-C842.
    - Schenk, P.M., Chapman, C.R., Zahnle, K., Moore, J.M., 2004. Ages and interiors: the cratering record of the Galilean satellites, Cambridge University Press. Cambridge University Press, Cambridge, UK.
    """

    def __init__(
        self,
        target: Target | str | None = None,
        projectile: Projectile | str | None = None,
        material: str | None = None,
        K1: FloatLike | None = None,
        mu: FloatLike | None = None,
        Ybar: FloatLike | None = None,
        density: FloatLike | None = None,
        **kwargs,
    ):
        super().__init__(
            target=target,
            projectile=projectile,
            material=material,
            K1=K1,
            mu=mu,
            Ybar=Ybar,
            density=density,
            **kwargs,
        )
        object.__setattr__(self, "_montecarlo_scaling", False)
        self.recompute()  # Recreate this with MC turned off

    def recompute(self):
        """
        Computes and sets the internal attributes for transition factors between simple and complex craters.
        """
        # Constants from CTEM
        CXEXPS = (
            1 / 0.885 - 1.0
        )  # Complex crater scaling explonent for silicate rock (Croft 1985)
        SIMCOMKS = 16533.8  # Simple-to-complex transition scaling coefficient for silicate rock
        SIMCOMPS = (
            -1.0303
        )  # Simple-to-complex transition scaling exponent for silicate rock
        CXEXPI = 0.155  # Complex crater scaling exponent for ice
        SIMCOMKI = 3081.39  # Simple-to-complex transition scaling coefficient for ice
        SIMCOMPI = -1.22486  # Simple-to-complex transition scaling exponent for ice

        # These terms are used to compute the ratio of the transient crater to simple crater size
        simple_enlargement = 0.84  # See Melosh (1989) pg. 129 just following eq. 8.2.1

        # These terms are used in the exponent in the final rim radius/ simple crater radius vs  final radius / transition radius relationship
        # See Holsapple (1993) eq. 28
        complex_enlargement_factor = 1.02

        # These terms are used to compute the transition diameter as a function of gravity
        # The transition values come from CTEM and are a synthesis of Pike (1980), Croft (1985), Schenk et al. (2004).
        if self.target.transition_scale_type == "silicate":
            simple_complex_exp = SIMCOMPS
            simple_complex_fac = SIMCOMKS
            final_exp = CXEXPS
        elif self.target.transition_scale_type == "ice":
            simple_complex_exp = SIMCOMPI
            simple_complex_fac = SIMCOMKI
            final_exp = CXEXPI

        # The nominal value will be used for determining the range of the "transitional" morphology type
        self._transition_nominal = float(
            simple_complex_fac * self.target.gravity**simple_complex_exp
        )

        simple_enlargement_factor = 1.0 / simple_enlargement
        self._simple_enlargement_factor = float(simple_enlargement_factor)
        self._complex_enlargement_factor = float(complex_enlargement_factor)
        self._final_exp = float(final_exp)
        self._transition_diameter = float(self._transition_nominal)
        return
