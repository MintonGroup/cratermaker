from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad
from scipy.optimize import root_scalar
from tqdm import tqdm

from cratermaker.bindings import morphology_bindings
from cratermaker.components.crater import Crater, CraterFixed
from cratermaker.components.morphology import Morphology, MorphologyCrater, MorphologyCraterVariable
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface

_RIMDROP = -6.0
_EJPROFILE = -3.0


@dataclass(frozen=True, slots=True)
class BasicMoonCraterFixed(CraterFixed):
    rim_elevation: float | None = None
    """Original rim height of the crater in meters relative to the reference surface."""
    floor_elevation: float | None = None
    """Original floor depth of the crater in meters relative to the reference surface."""
    floor_radius: float | None = None
    """Original floor diameter of the crater in meters."""
    wall_curvature: float | None = None
    """The curvature of the crater walls."""
    rim_width: float | None = None
    """The width of the crater rim in meters."""
    rim_flank_radius: float | None = None
    """Mean radius of the rim flank of the crater in meters."""
    rimdrop: float | None = None
    """The power law exponent for the structural uplift underneath the ejecta"""
    ejrim: float | None = None
    """Original ejecta rim thickness of the crater in meters."""
    ejprofile: float | None = None
    """Power law exponent for the ejecta thickness profile of the crater."""
    peak_height: float | None = None
    """Original central peak height of the crater in meters relative to the reference surface. 0 for simple craters."""
    peak_width: float | None = None
    """Original central peak width of the crater in meters. 0 for simple craters."""
    peak_offset: float | None = None
    """Original central peak offset of the crater in meters. 0 for simple craters."""

    @property
    def depth_to_diameter(self) -> float | None:
        """
        The depth to diameter ratio of the crater.

        This is computed from `rim_elevation`-`floor_elevation`.
        """
        floor_elevation = self.floor_elevation
        rim_elevation = self.rim_elevation
        if floor_elevation is not None and rim_elevation is not None:
            return (rim_elevation - floor_elevation) / self.diameter
        else:
            return None


@Crater.register("basicmooncrater")
class BasicMoonCrater(MorphologyCrater):
    def __init__(
        self, crater: Crater | None = None, fixed_cls=BasicMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += (
            f"Rim height: {format_large_units(self.rim_elevation, quantity='length')}\n"
            f"Floor depth: {format_large_units(self.floor_elevation, quantity='length')}\n"
            f"Floor diameter: {format_large_units(self.floor_radius, quantity='length')}\n"
            f"Central peak height: {format_large_units(self.peak_height, quantity='length') if self.peak_height else 'None'}\n"
            f"Ejecta rim thickness: {format_large_units(self.ejrim, quantity='length')}\n"
        )
        if self.morphology_subtype is not None:
            str_repr += f"Morphology subtype: {self.morphology_subtype}\n"
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        floor_elevation: float | None = None,
        floor_radius: float | None = None,
        wall_curvature: float | None = None,
        rim_width: float | None = None,
        rim_elevation: float | None = None,
        rim_flank_radius: float | None = None,
        rimdrop: float | None = None,
        ejrim: float | None = None,
        ejprofile: float | None = None,
        peak_height: float | None = None,
        peak_width: float | None = None,
        peak_offset: float | None = None,
        **kwargs: Any,
    ) -> BasicMoonCrater:
        """
        Initialize a BasicMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters. The morphometric parameters are mostly taken from Pike (1977) [#]_ for D>5 km craters with a higher value of d/D and floor_radius from Fassett and Thomson (2014) [#]_, Yang et al. (2021) [#]_ for D<50 m craters, and a random weighted mixture of the two models using the d/D vs D trend seen in Hoover at al. (2024) [#]_.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        rim_elevation : float, optional
            Original rim height of the crater in meters relative to the reference surface. If None, it will be computed.
        rim_flank_radius : float, optional
            The mean radius of the rim flank in meters. If None, it will be computed.
        floor_elevation : float, optional
            Original floor depth of the crater in meters relative to the reference surface. If None, it will be computed.
        floor_radius : float, optional
            Original floor radius of the crater in meters. If None, it will be computed.
        wall_curvature : float, optional
            The curvature of the crater walls. If None, it will be computed based on the morphology type and diameter.
        ejrim : float, optional
            Original ejecta rim thickness of the crater in meters. If None, it will be computed.
        peak_height : float, optional
            Original central peak height of the crater in meters relative to the reference surface. If None, it will be computed for complex craters and set to 0 for simple craters.
        peak_width : float, optional
            Original central peak width of the crater in meters. If None, it will be computed for complex craters and set to 0 for simple craters.
        peak_offset : float, optional
            Original central peak offset of the crater in meters. If None, it will be computed for complex craters and set to 0for simple craters.
        **kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        References
        ----------
        .. [#] Pike, R.J., 1977. Size-dependence in the shape of fresh impact craters on the moon. Presented at the In: Impact and explosion cratering: Planetary and terrestrial implications; Proceedings of the Symposium on Planetary Cratering Mechanics, pp. 489-509.
        .. [#] Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. `doi:10.1002/2014JE004698 <https://doi.org/10.1002/2014JE004698>`_
        .. [#] Yang, X., Fa, W., Du, J., Xie, M., Liu, T., 2021. Effect of Topographic Degradation on Small Lunar Craters: Implications for Regolith Thickness Estimation. Geophysical Research Letters 48, e2021GL095537. `doi:10.1029/2021GL095537 <https://doi.org/10.1029/2021GL095537>`_
        .. [#] Hoover, R.H., Robbins, S.J., Hynek, B.M., Hayne, P.O., 2024. Depth-to-diameter Ratios of Fresh Craters on the Moon and Implications for Surface Age Estimates. Planet. Sci. J. 5, 26. `doi:10.3847/PSJ/ad18d4 <https://doi.org/10.3847/PSJ/ad18d4>`_
        """
        from cratermaker.components.morphology import Morphology
        from cratermaker.utils.montecarlo_utils import bounded_norm, sample_logfit_heteroskedastic, sample_pikefit

        # This is a copy operation, to use old values for any un-specified arguments
        if crater is not None and isinstance(crater, BasicMoonCrater):
            rim_elevation = crater.rim_elevation if rim_elevation is None else rim_elevation
            rim_flank_radius = crater.rim_flank_radius if rim_flank_radius is None else rim_flank_radius
            floor_elevation = crater.floor_elevation if floor_elevation is None else floor_elevation
            floor_radius = crater.floor_radius if floor_radius is None else floor_radius
            rimdrop = crater.rimdrop if rimdrop is None else rimdrop
            ejrim = crater.ejrim if ejrim is None else ejrim
            ejprofile = crater.ejprofile if ejprofile is None else ejprofile
            peak_height = crater.peak_height if peak_height is None else peak_height
            peak_width = crater.peak_width if peak_width is None else peak_width
            peak_offset = crater.peak_offset if peak_offset is None else peak_offset

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)

        depth_params = {
            "simple_sub500m": {
                "coefficients": [-1.910001619942053, 0.7831875571105978, 0.042836046225245894],
                "c": -2.975325517042965,
                "alpha": 1.4136290732609256,
            },
            "simple": {
                "coefficients": [-4.521021365770051, 1.7347527856394551, -0.04489022602475602],
                "c": -5.26570819908019,
                "alpha": 1.8641293440183668,
            },
            "transitional": {
                "coefficients": [22.18782422399827, -3.203597438379582, 0.17631543814150283],
                "c": 63.24942043364388,
                "alpha": -6.715217581955981,
            },
            "complex": {
                "coefficients": [-1.6315393638746256, 1.5032727501015142, -0.055610926613342146],
                "c": -27.654075266623785,
                "alpha": 4.694621309150626,
            },
        }
        rim_elevation_params = {
            "simple": {
                "coefficients": [-9.087441060124089, 2.3723804170387552, -0.07722891774004866],
                "c": -1.3414979152290676,
                "alpha": 1.3146593627541205,
            },
            "transitional": {
                "coefficients": [-10.564003835229924, 2.9763623205780294, -0.12704846222410626],
                "c": -53.730085877961216,
                "alpha": 9.537258205109097,
            },
            "complex": {
                "coefficients": [-12.791362969209594, 3.1852601403185417, -0.12566215471641845],
                "c": -4.245272638928989,
                "alpha": 1.9263059567962026,
            },
        }
        floor_radius_params = {
            "simple": {
                "coefficients": [-19.235109460529202, 4.840583566677377, -0.21450124590854683],
                "c": -8.422750876479917,
                "alpha": 2.569365299227605,
            },
            "transitional": {
                "coefficients": [-47.25759607077153, 10.01636202766342, -0.4444138708218017],
                "c": 8.073586072170954,
                "alpha": 0.41963094367654863,
            },
            "complex": {
                "coefficients": [4.571914820366517, -0.21608677058312598, 0.061233414548021586],
                "c": 0.7771019915975753,
                "alpha": 1.295875187989571,
            },
        }
        rim_width_params = {
            "simple": {
                "coefficients": [-24.82061509600679, 6.470961147998662, -0.32200704888815185],
                "c": 0.9641271074363086,
                "alpha": 1.619163211456247,
            },
            "transitional": {
                "coefficients": [-129.98317617103737, 27.374040577016952, -1.3529946148839573],
                "c": 13.424721979662507,
                "alpha": 0.01699281614062642,
            },
            "complex": {
                "coefficients": [-6.156942924353989, 1.9832268288802843, -0.055861398599379546],
                "c": -4.278121587075394,
                "alpha": 2.348577422207449,
            },
        }
        args = {}
        diameter_m = crater.diameter
        diameter_km = diameter_m * 1e-3
        if crater.morphology_type in ["basin", "multiring", "peakring"]:
            morphology_type = "complex"
        else:
            morphology_type = crater.morphology_type

        if rim_elevation is None:
            rim_elevation = max(sample_logfit_heteroskedastic(diameter_m, **rim_elevation_params[morphology_type])[0], 0.0)
        args["rim_elevation"] = rim_elevation

        if floor_elevation is None:
            if crater.diameter < 500.0:
                floor_elevation = -sample_logfit_heteroskedastic(diameter_m, **depth_params["simple_sub500m"])[0] + rim_elevation
            else:
                floor_elevation = -sample_logfit_heteroskedastic(diameter_m, **depth_params[morphology_type])[0] + rim_elevation
            floor_elevation = min(floor_elevation, 0.0)
        args["floor_elevation"] = floor_elevation

        if floor_radius is None:
            floor_radius = max(sample_logfit_heteroskedastic(diameter_m, **floor_radius_params[morphology_type])[0], 0.0)
        args["floor_radius"] = floor_radius

        if peak_height is None:
            if morphology_type != "complex":
                peak_height = 0.0
            else:
                peak_height = sample_pikefit(diameter_km, a=0.900, b=0.032, errhi=0.0011, errlo=-0.008, n=22)[0] * 1e3
        args["peak_height"] = peak_height

        if wall_curvature is None:
            wall_curvature = morphology.rng.uniform(low=0, high=6, size=1)[0]
        args["wall_curvature"] = wall_curvature

        if rim_width is None:
            rim_width = max(sample_logfit_heteroskedastic(diameter_m, **rim_width_params[morphology_type])[0], 0.0)
        args["rim_width"] = rim_width

        if rim_flank_radius is None:
            # Values from Jun Du
            a = 2.8740358071397156
            b = -0.12400776345642424
            sigma = 0.18025709610312335
            mean = a * diameter_km**b
            rim_flank_radius = morphology.rng.normal(mean, sigma, 1)[0] * diameter_km / 2 * 1e3
        args["rim_flank_radius"] = rim_flank_radius

        args["peak_width"] = args["peak_height"] * 2 if peak_width is None else peak_width
        args["peak_offset"] = 0.0 if peak_offset is None else peak_offset

        # Try to approximately conserve volume when setting the ejecta thickness at the rim value
        args["ejprofile"] = _EJPROFILE if ejprofile is None else ejprofile
        args["rimdrop"] = _RIMDROP if rimdrop is None else rimdrop
        if ejrim is None:
            ejrim = 0.14 * (diameter_m / 2) ** 0.74
            if diameter_km > 300:
                hf = -args["floor_elevation"]
                hr = args["rim_elevation"]
                fr = args["floor_radius"] / (0.5 * diameter_m)
                prd = args["rimdrop"]
                pej = args["ejprofile"]
                ejrim = max(ejrim, (hf * fr**2 - 2 * hr / (2 - prd)) / (2 * (1.0 / (2 - pej) - 1.0 / (2 - prd))))
                ejrim = min(ejrim, args["rim_elevation"])
        args["ejrim"] = ejrim

        kwargs = {**args, **kwargs}

        return cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )


@Morphology.register("basicmoon")
class BasicMoonMorphology(Morphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh.
    """

    def __init__(
        self,
        surface: Surface | str | None = None,
        ejecta_truncation: FloatLike | None = None,
        dorays: bool = False,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
        ejecta_truncation : float, optional
            The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a
            truncation distance based on where the ejecta thickness reaches a small value.
        dorays : bool, optional
            A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket, default is False.
        rng : numpy.random.Generator | None
            |rng|
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            |rng_seed|
        rng_state : dict, optional
            |rng_state|
        **kwargs : Any
            |kwargs|
        """
        object.__setattr__(self, "_ejecta_truncation", None)
        object.__setattr__(self, "_node", None)
        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays
        super().__init__(surface=surface, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        if self.ejecta_truncation is not None:
            str_repr += f"Ejecta Trunction: {self.ejecta_truncation:.2f} * crater.radius\n"
        else:
            str_repr += "Ejecta Truncation: Off\n"
        str_repr += f"Ejecta Rays: {self.dorays}\n"
        return str_repr

    def emplace(self, craters: Crater | list[Crater] | None = None, **kwargs: Any) -> list[BasicMoonCrater]:
        """
        Convenience method to immediately emplace a crater onto the surface.

        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        crater : Crater | list[Crater] | None
            The crater or list of craters to be emplaced. If None, then a crater will be created using the provided parameters in kwargs and emplaced.
        kwargs : Any
            |kwargs|

        Returns
        -------
        list[BasicMoonCrater]
            The list of BasicMoonCrater objects that were emplaced.
        """
        if craters is None:
            craters = [BasicMoonCrater.maker(morphology=self, **kwargs)]
        elif isinstance(craters, (list | tuple)) and len(craters) > 0:
            processed_craters = []
            for c in tqdm(
                craters,
                total=len(craters),
                desc="Preparing craters for emplacement",
                unit="crater",
                position=0,
                leave=False,
            ):
                if isinstance(c, BasicMoonCrater):
                    processed_craters.append(c)
                else:
                    processed_craters.append(BasicMoonCrater.maker(c, morphology=self, **kwargs))
            craters = processed_craters
        elif isinstance(craters, BasicMoonCrater):
            craters = [craters]
        elif isinstance(craters, Crater):
            craters = [BasicMoonCrater.maker(craters, morphology=self, **kwargs)]
        else:
            raise ValueError(
                "Invalid input for crater emplacement. Must be a Crater object, a list of Crater objects, or None with additional arguments for Crater.maker()."
            )

        return super().emplace(craters, **kwargs)

    def form_crater(self, crater: Crater | BasicMoonCrater, **kwargs: Any) -> None:
        """
        Form a crater or list of craters on the surface.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object to be formed.
        kwargs : Any
            |kwargs|
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)

        if crater.emplaceable:
            super().form_crater(crater, **kwargs)
        return

    def form_ejecta(self, crater: Crater | BasicMoonCrater, **kwargs: Any) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Form the ejecta blanket of the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        crater : Crater
            The crater object whose ejecta is to be formed.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        tuple[NDArray[np.float64], NDArray[np.float64]]
            The computed ejecta thickness and intensity at the face and node elevations.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        ejecta_thickness, ejecta_intensity = super().form_ejecta(crater, **kwargs)
        k_ej = self.ejecta_burial_degradation(ejecta_thickness[: crater.ejecta_region.n_face], ejecta_soften_factor=1.50)
        crater.ejecta_region.apply_diffusion(k_ej)

        k_deg = self.degradation_function(crater.diameter, fe=100) * ejecta_intensity[: crater.ejecta_region.n_face]
        crater.ejecta_region.apply_diffusion(k_deg)
        return ejecta_thickness, ejecta_intensity

    def crater_shape(self, crater: BasicMoonCrater, region: LocalSurface, **kwargs: Any) -> NDArray[np.float64]:
        """
        Compute the crater shape based on the region view and surface.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the crater shape.
        region:  LocalSurface
            The region view of the surface mesh centered at the crater center.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        NDArray[np.float64]
            The computed crater shape at the face and node elevations.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        reference_elevation = region.get_reference_surface(reference_radius=crater.radius)

        # Combine distances and references for nodes and faces
        r = np.concatenate([region.face_distance, region.node_distance])
        bearing = np.concatenate([region.face_bearing, region.node_bearing])

        original_elevation = np.concatenate([region.face_elevation, region.node_elevation])

        new_elevation = self.crater_profile(crater=crater, r=r, bearing=bearing, r_ref=reference_elevation)
        elevation_change = new_elevation - original_elevation

        return elevation_change

    def crater_profile(
        self, crater: BasicMoonCrater, r: ArrayLike, r_ref: ArrayLike | None = None, **kwargs: Any
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the crater profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        r_ref : ArrayLike, optional
            Reference elevation values to be modified by the crater profile.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed crater elevation profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        if r_ref is None:
            r_ref = np.zeros_like(r)

        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)

        # flatten r to 1D array
        rflat = np.ravel(r)
        r_ref_flat = np.ravel(r_ref)
        elevation = morphology_bindings.basicmoon_profile(
            radial_distances=rflat, reference_elevations=r_ref_flat, crater=crater, include_crater=True, include_ejecta=False
        )
        # reshape elevation to match the shape of r
        elevation = np.array(elevation, dtype=np.float64)
        elevation = np.reshape(elevation, r.shape)

        return elevation

    def ejecta_shape(
        self, crater: BasicMoonCrater, region: LocalSurface, **kwargs: Any
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the ejecta shape based on the region view and surface.

        Parameters
        ----------
        region : LocalSurface
            The region view of the surface mesh centered at the crat er center.
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta shape.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        NDArray[np.float64]
            The computed ejecta shape at the face and node elevations
        NDArray[np.float64]
            The computed ejecta intensity at the face and node elevations which is used to compute the degradation function. When `dorays` is True, this is the ray intensity function. When `dorays` is False, this is just a constant 1 over the ejecta region.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)

        distance = np.concatenate([region.face_distance, region.node_distance])
        if self.dorays:
            bearing = np.concatenate([region.face_bearing, region.node_bearing])
            thickness, intensity = self.ejecta_distribution(crater, distance, bearing)
        else:
            thickness = self.ejecta_profile(crater, distance)
            intensity = np.ones_like(thickness)

        return thickness, intensity

    def ejecta_profile(self, crater: BasicMoonCrater, r: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed ejecta profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)
        # flatten r to 1D array
        rflat = np.ravel(r)
        elevation = morphology_bindings.basicmoon_profile(
            radial_distances=rflat,
            reference_elevations=np.zeros_like(rflat),
            crater=crater,
            include_crater=False,
            include_ejecta=True,
        )
        elevation = np.array(elevation, dtype=np.float64)
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)
        return elevation

    def ejecta_distribution(self, crater, r: ArrayLike, theta: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the ejecta thickness distribution modulated by ray patterns.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta distribution.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        theta : ArrayLike
            Angular bearings from the crater center (in degrees).

        Returns
        -------
        thickness : NDArray[np.float64]
            The computed ejecta thickness for each (r, theta) pair.
        ray_intensity : NDArray[np.float64]
            The computed ray intensity for each (r, theta) pair.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        # flatten r and theta to 1D arrays
        thickness = self.ejecta_profile(crater, r)
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        intensity = morphology_bindings.ray_intensity(
            radial_distances=rflat,
            initial_bearing=np.radians(theta_flat),
            crater_diameter=crater.diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        thickness = np.array(thickness, dtype=np.float64)
        intensity = np.array(intensity, dtype=np.float64)
        thickness *= intensity
        # reshape thickness to match the shape of r and theta
        return np.reshape(thickness, r.shape), np.reshape(intensity, r.shape)

    def ray_intensity(self, crater: BasicMoonCrater, r: ArrayLike, theta: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ray pattern intensity modulation at each (r, theta) pair.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ray intensity.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        theta : ArrayLike
            Angular bearings from the crater center (in degrees).

        Returns
        -------
        intensity : NDArray[np.float64]
            The computed ray intensity values.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        # flatten r and theta to 1D arrays
        rflat = np.ravel(r)
        theta_flat = np.radians(np.ravel(theta))
        intensity = morphology_bindings.ray_intensity(
            rflat,
            theta_flat,
            crater.diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        intensity = np.array(intensity, dtype=np.float64)
        # reshape intensity to match the shape of r and theta
        intensity = np.reshape(intensity, r.shape)
        return intensity

    def ejecta_burial_degradation(self, ejecta_thickness, ejecta_soften_factor=1.50) -> NDArray[np.float64]:
        """
        Computes the change in degradation state due to ejecta burial.

        This function implements a combination of the model by Minton et al. (2019) [#]_.

        Parameters
        ----------
        region : LocalSurface
            The region view of the surface mesh centered at the crater center.
        ejecta_thickness : NDArray[np.float64]
            The computed ejecta thickness at the face and node elevations.

        Returns
        -------
        NDArray[np.float64]
            The computed change in degradation state for all faces in the

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. `doi: 10.1016/j.icarus.2019.02.021 <https://doi.org/10.1016/j.icarus.2019.02.021>`_
        """
        return ejecta_soften_factor * ejecta_thickness**2

    def compute_subpixel_degradation(
        self,
        time_start: float,
        time_end: float,
        **kwargs,
    ) -> None:
        """
        Performs the subpixel degradation.

        This models the combined degradation of the part of the production population that is below the resolution of the mesh on each face. It is called between batches of craters by the `emplace` method.

        Parameters
        ----------
        time_start : float
            The time of the surface at the start of the degradation.
        time_end : float
            The time of the surface at the end of the degradation.
            |kwargs|
        """
        dc_min = 1e-8  # Minimum crater size for subpixel degradation calculation.

        if time_end >= time_start:
            raise ValueError("time_end must be less than time_start.")
        if self.production is None:
            raise RuntimeError("Production model must be set in the Morphology object to compute subpixel degradation.")

        if not hasattr(self, "_Kdiff"):
            self._Kdiff = np.zeros_like(self.surface.face_elevation)

        def _subpixel_degradation(diameter):
            fe = 100.0
            k = self.degradation_function(diameter, fe)
            n = self.production.function(
                diameter=diameter,
                time_start=time_start,
                time_end=time_end,
                validate_inputs=False,
            ).item()
            degradation_region_area = np.pi * (diameter / 2) * fe
            return k * n * degradation_region_area

        for face_indices, dc_max in zip(self.surface.face_bin_indices, self.surface.face_bin_max_sizes, strict=False):
            delta_kdiff, _ = quad(_subpixel_degradation, dc_min, dc_max)
            self._Kdiff[face_indices] += delta_kdiff

        # If any Kdiff values reaches a threshold where a meaningful amount of diffusion will occur on the surface, then go ahead and apply it.
        # Otherwise, degradation will continue to accumulate until the next batch of craters is processed.
        if np.any(self._Kdiff / self.surface.face_area > 1):
            self.apply_subpixel_degradation()
            if self.do_counting:
                measure_rim = kwargs.pop("measure_rim", False)
                self.counting.tally(measure_rim=measure_rim, **kwargs)

        return

    def apply_subpixel_degradation(self, **kwargs) -> None:
        """
        Apply subpixel degradation to the surface using the current Kdiff values.

        This method is called after all craters have been processed and is used to
        apply the accumulated degradation effects.
        """
        if not hasattr(self, "_Kdiff"):
            return

        self.surface.apply_diffusion(self._Kdiff)
        self._Kdiff = np.zeros_like(self.surface.face_elevation)
        return

    def rmax(
        self,
        crater: Crater,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
    ) -> float:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor, whichever is smaller.

        Parameters
        ----------
        crater : Crater
            The crater object to be used.
        minimum_thickness : FloatLike
            The minimum thickness of the feature blanket in meters.
        feature : str, optional, default = "ejecta"
            The feature to compute the maximum extent. Either "crater" or "ejecta". If "crater" is chosen, the rmax is based
            on where the raised rim is smaller than minimum thickness.

        Returns
        -------
        float
            The maximum extent of the crater or ejecta blanket in meters.
        """
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)

        def _profile_invert_ejecta(r):
            ans = self.ejecta_profile(crater, r) - minimum_thickness
            return ans[0]

        def _profile_invert_crater(r):
            ans = self.crater_profile(crater, r, np.zeros(1)) - minimum_thickness
            return ans[0]

        if feature == "ejecta":
            _profile_invert = _profile_invert_ejecta
        elif feature == "crater":
            _profile_invert = _profile_invert_crater
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")

        # Get the maximum extent
        lower_limit = crater.radius * 1.0001
        upper_limit = self.ejecta_truncation * crater.radius if self.ejecta_truncation else np.pi * self.surface.target.radius

        if _profile_invert(lower_limit) < 0:
            ans = lower_limit
        elif _profile_invert(upper_limit) > 0:
            ans = upper_limit
        else:
            sol = root_scalar(
                _profile_invert,
                bracket=[lower_limit, upper_limit],
                method="brentq",
            )
            ans = sol.root if sol.converged else crater.radius

        return float(ans)

    def degradation_function(
        self,
        diameter: FloatLike,
        fe: FloatLike = 100.0,
    ) -> float:
        """
        Computes the degradation function, which defines the topographic degradation that each crater contributes to the surface.

        This function implements a combination of the model by Minton et al. (2019) [#]_ for small craters and Riedel et al. (2020) [#]_ for large craters. It is currently not well-constrained, so may change in the future.

        Parameters
        ----------
        diameter : FloatLike
            The final diameter of the crater in meters.
        fe : FloatLike, optional
            The degradation function size factor, which is a scaling factor for the degradation function. Default is 100.0.

        Returns
        -------
        float
            The computed degradation function


        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. `doi:10.1016/j.icarus.2019.02.021 <https://doi.org/10.1016/j.icarus.2019.02.021>`_
        .. [#] Riedel, C., Minton, D.A., Michael, G., Orgel, C., Bogert, C.H. van der, Hiesinger, H., 2020. Degradation of Small Simple and Large Complex Lunar Craters: Not a Simple Scale Dependence. Journal of Geophysical Research: Planets 125, e2019JE006273. `doi:10.1029/2019JE006273 <https://doi.org/10.1029/2019JE006273>`_
        """

        def _kdmare(r, fe, psi):
            """
            The mare-scale degradation function from Minton et al. (2019). See eq. (32).
            """
            kv1 = 0.30
            neq1 = 0.0084
            eta = 3.2
            gamma = 2.0
            beta = 2.0
            kd1 = kv1 * (math.pi * fe**2 * neq1 * (gamma * beta / ((eta - 2.0) * (beta + gamma - eta)))) ** (-gamma / (eta - beta))
            psi = gamma * ((eta - 2.0) / (eta - beta))
            return kd1 * r**psi

        def _smooth_broken(x, a, x_break, alpha_1, alpha_2, delta):
            return a * (x / x_break) ** (alpha_1) * (0.5 * (1.0 + (x / x_break) ** (1.0 / delta))) ** ((alpha_2 - alpha_1) / delta)

        def _kd(r, fe):
            psi_1 = 2.0  # Mare scale power law exponent
            psi_2 = 1.2  # Highlands scale power law exponent
            rb = 0.5e3  # breakpoint radius
            delta = 1.0e0  # Smoothing function

            kd1 = _kdmare(rb, fe, psi_1) / (1 + (psi_1 - psi_2) / psi_1) ** 2
            return _smooth_broken(r, kd1, rb, psi_1, psi_2, delta)

        return float(_kd(diameter / 2, fe))

    @staticmethod
    def overlap_function(crater: BasicMoonCrater) -> tuple[set[int], set[int]]:
        """
        Get the affected node and face indices for the crater.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object to be used.

        Returns
        -------
        tuple[set[int], set[int]]
            The affected node and face indices for the crater.
        """
        if not crater.emplaceable or crater.affected_node_indices is None or crater.affected_face_indices is None:
            return set(), set()
        else:
            return crater.affected_node_indices, crater.affected_face_indices

    @parameter
    def ejecta_truncation(self) -> float:
        """
        The radius at which the crater is truncated relative to the crater radius.

        Returns
        -------
        float or None
        """
        return self._ejecta_truncation

    @ejecta_truncation.setter
    def ejecta_truncation(self, value: FloatLike | None):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("truction_radius must be of type FloatLike")
            self._ejecta_truncation = float(value)
        else:
            self._ejecta_truncation = None

    @parameter
    def dorays(self) -> bool:
        """
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket.

        Returns
        -------
        bool
        """
        return self._dorays

    @dorays.setter
    def dorays(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("dorays must be of type bool")
        self._dorays = value

    @property
    def _CraterType(self) -> type[BasicMoonCrater]:
        """
        The class definition of the associated Crater type.
        """
        return BasicMoonCrater
