from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad

from cratermaker.bindings import basicmoon_bindings
from cratermaker.components.crater import Crater, CraterFixed
from cratermaker.components.morphology import Morphology, MorphologyCrater, MorphologyCraterVariable
from cratermaker.components.surface import Surface
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter

_EJPROFILE = -3.0

# The crater properties that define the axisymmetric profile, as read by the compiled profile and volume functions
_PROFILE_PARAMETERS = (
    "diameter",
    "radius",
    "floor_elevation",
    "floor_radius",
    "wall_curvature",
    "floor_blend",
    "rim_width",
    "rim_height",
    "frac_ejrim",
    "ejprofile",
    "peak_height",
    "peak_width",
    "peak_ring_radius",
    "peak_center_distance",
    "peak_center_bearing",
)

# Nodes and weights of the 8-point Gauss-Legendre rule on [-1, 1], used for the fixed-order composite quadratures below.
_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(8)


def _gauss_legendre_points(edges: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Return the abscissas and weights of a composite 8-point Gauss-Legendre rule.

    Parameters
    ----------
    edges : NDArray[np.float64]
        Sorted, strictly increasing boundaries of the subintervals. Each subinterval gets its own 8-point rule.

    Returns
    -------
    x : NDArray[np.float64]
        The abscissas of all subintervals, flattened.
    w : NDArray[np.float64]
        The corresponding weights, such that ``np.sum(w * f(x))`` approximates the integral of f from edges[0] to edges[-1].
    """
    a = edges[:-1, np.newaxis]
    b = edges[1:, np.newaxis]
    half_width = 0.5 * (b - a)
    x = (half_width * _GL_NODES + 0.5 * (a + b)).ravel()
    w = (half_width * _GL_WEIGHTS).ravel()
    return x, w


@dataclass(frozen=True, slots=True)
class BasicMoonCraterFixed(CraterFixed):
    rim_height: float | None = None
    """Original rim height of the crater in meters relative to the reference surface."""
    floor_elevation: float | None = None
    """Original floor depth of the crater in meters relative to the reference surface."""
    floor_radius: float | None = None
    """Original floor diameter of the crater in meters."""
    wall_curvature: float | None = None
    """A factor that controls the curvature of the crater wall (0 for straight walls, and 1.0 for very curvy walls)."""
    floor_blend: float | None = None
    """A factor that controls how sharply the wall and floor blend (0 for sharp floor-to-wall transition, 1.0 for gentle floor-to-wall transition)."""
    rim_width: float | None = None
    """The width of the crater rim in meters."""
    ejprofile: float | None = None
    """Power law exponent for the ejecta thickness profile of the crater."""
    peak_height: float | None = None
    """Central peak height of the crater in meters relative to the reference surface. 0 for simple craters."""
    peak_width: float | None = None
    """Central peak width of the crater in meters. 0 for simple craters."""
    peak_ring_radius: float | None = None
    """Radius of peak ring in meters. 0 for central peaks."""
    peak_center_distance: float | None = None
    """Distance of central peak/peak ring from crater center in meters."""
    peak_center_bearing: float | None = None
    """Bearing angle of central peak/peak ring in degrees."""
    isring: bool | None = field(default=False, init=True)
    """Flag that indicates that this is a ring rather than a crater."""
    parent: np.uint64 | None = field(default=np.uint64(0), init=True)

    @property
    def depth_to_diameter(self) -> float | None:
        """
        The depth to diameter ratio of the crater.

        This is computed from `rim_height`-`floor_elevation`.
        """
        floor_elevation = self.floor_elevation
        rim_height = self.rim_height
        if floor_elevation is not None and rim_height is not None:
            return (rim_height - floor_elevation) / self.diameter
        else:
            return None


@Crater.register("basicmooncrater")
class BasicMoonCrater(MorphologyCrater):
    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=BasicMoonCraterFixed,
        variable_cls=MorphologyCraterVariable,
        **kwargs,
    ):
        if hasattr(crater, "isring"):
            isring = crater.isring
        else:
            isring = False
        isring = kwargs.pop("isring", isring)
        super().__init__(
            crater=crater,
            fixed_cls=fixed_cls,
            variable_cls=variable_cls,
            isring=isring,
            **kwargs,
        )

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += (
            f"Rim height: {format_large_units(self.rim_height, quantity='length')}\n"
            f"Ejecta fraction of rim: {self.frac_ejrim}\n"
            f"Ejecta profile: {self.ejprofile}\n"
            f"Rim width: {format_large_units(self.rim_width, quantity='length')}\n"
            f"Floor elevation: {format_large_units(self.floor_elevation, quantity='length')}\n"
            f"Floor radius: {format_large_units(self.floor_radius, quantity='length')}\n"
            f"Floor blend factor: {self.floor_blend}\n"
            f"Wall curvature factor: {self.wall_curvature}\n"
            f"Central peak height: {format_large_units(self.peak_height, quantity='length') if self.peak_height else 'None'}\n"
            f"Central peak width: {format_large_units(self.peak_width, quantity='length') if self.peak_width else 'None'}\n"
            f"Central peak offset: {format_large_units(self.peak_ring_radius, quantity='length') if self.peak_ring_radius else 'None'}\n"
        )
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        floor_elevation: float | None = None,
        floor_radius: float | None = None,
        wall_curvature: float | None = None,
        floor_blend: float | None = None,
        rim_width: float | None = None,
        rim_height: float | None = None,
        frac_ejrim: float | None = None,
        ejprofile: float | None = None,
        peak_height: float | None = None,
        peak_width: float | None = None,
        peak_ring_radius: float | None = None,
        peak_center_distance: float | None = None,
        peak_center_bearing: float | None = None,
        conserve_volume: bool = True,
        monte_carlo: bool | None = None,
        make_rings: bool | None = None,
        **kwargs: Any,
    ) -> BasicMoonCrater:
        """
        Initialize a BasicMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        rim_height : float, optional
            Original rim height of the crater in meters relative to the reference surface. If None, it will be computed.
        floor_elevation : float, optional
            Original floor depth of the crater in meters relative to the reference surface. If None, it will be computed.
        floor_radius : float, optional
            Original floor radius of the crater in meters. If None, it will be computed.
        wall_curvature : float, optional
            A factor that controls the curvature of the crater wall (0 for straight walls, and 1.0 for very curvy walls). If None, it will be computed based on the morphology type and diameter.
        floor_blend: float, optional
            A factor that controls how sharply the wall and floor blend (0 for sharp floor-to-wall transition, 1.0 for gentle floor-to-wall transition). If None, it will be set to 0.5.
        frac_ejrim : float, optional
            Original ejecta rim thickness of the crater in meters. If None, it will be computed.
        peak_height : float, optional
            Original central peak height of the crater in meters relative to the reference surface. If None, it will be computed for complex craters and set to 0 for simple craters.
        peak_width : float, optional
            Original central peak width of the crater in meters. If None, it will be computed for complex craters and set to 0 for simple craters.
        peak_ring_radius : float, optional
            Original central peak offset of the crater in meters. If None, it will be computed for complex craters and set to 0for simple craters.
        peak_center_distance: float, optional
            Distance of central peak/peak ring from crater center in meters.
        peak_center_bearing: float, optional
            Bearing angle of central peak/peak ring in degrees.
        conserve_volume: bool, optional
            If True, then the value of frac_ejrim will be adjusted in order to attempt to conserve volume between excavation and deposition. Note that even when True, further adjustments may be needed when emplacing the crater due to the local topography. Default is True.
        monte_carlo : bool or None, optional
            If True, then crater parameters will be drawn from random distributions. If False, then nominal values will be used for all parameters. If None, then the value of the ``monte_carlo`` of the associated Scaling object will be used. Default is None.
        make_rings : bool, optional
            If True, then rings will be generated for complex and multiring craters. Default is True for new craters and False for copied craters unless a parameter was changed.
        **kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        References
        ----------
        TODO: Update with Minton et al. (2026) info once done.
        """
        input_args = locals()
        source = crater
        from cratermaker.components.morphology import Morphology
        from cratermaker.utils.montecarlo_utils import sample_logfit_heteroskedastic

        if crater is not None and isinstance(crater, BasicMoonCrater):
            for k in ["__class__", "kwargs", "morphology", "crater", "conserve_volume", "cls", "make_rings"]:
                input_args.pop(k, None)
            conserve_volume = conserve_volume and any(list(input_args.values()))
            if make_rings is None:
                make_rings = any(list(input_args.values()))
            # This is a copy operation, to use old values for any un-specified arguments
            floor_elevation = crater.floor_elevation if floor_elevation is None else floor_elevation
            floor_radius = crater.floor_radius if floor_radius is None else floor_radius
            wall_curvature = crater.wall_curvature if wall_curvature is None else wall_curvature
            floor_blend = crater.floor_blend if floor_blend is None else floor_blend
            rim_width = crater.rim_width if rim_width is None else rim_width
            rim_height = crater.rim_height if rim_height is None else rim_height
            frac_ejrim = crater.frac_ejrim if frac_ejrim is None else frac_ejrim
            ejprofile = crater.ejprofile if ejprofile is None else ejprofile
            peak_height = crater.peak_height if peak_height is None else peak_height
            peak_width = crater.peak_width if peak_width is None else peak_width
            peak_ring_radius = crater.peak_ring_radius if peak_ring_radius is None else peak_ring_radius
            peak_center_distance = crater.peak_center_distance if peak_center_distance is None else peak_center_distance
            peak_center_bearing = crater.peak_center_bearing if peak_center_bearing is None else peak_center_bearing
        elif make_rings is None:
            make_rings = True

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, monte_carlo=monte_carlo, **kwargs)
        monte_carlo = monte_carlo if monte_carlo is not None else morphology.scaling.monte_carlo
        compute_nominal = not monte_carlo
        rng = morphology.rng

        # The following are derived based on fitting profiles to fresh craters (Minton et al. 2026 in prep)

        floor_elevation_params = {
            "simple": {
                "coefficients": [-4.167687018033455, 1.6001202454523342, -0.03747552110388699],
                "c": -2.274464155672326,
                "alpha": 1.455516197513936,
                "fit_scale": "loglog",
            },
            "non-simple": {
                "coefficients": [5.4915689253420545, 0.23298682432901896, -0.0022430661802085037],
                "c": -3.7362261125113427,
                "alpha": 1.866447443526078,
                "fit_scale": "loglog",
            },
        }
        rim_height_params = {
            "simple": {
                "coefficients": [-4.461857399058903, 1.2645885691929473, -0.011531387005605066],
                "c": -4.218604205452387,
                "alpha": 1.72797341370465,
                "fit_scale": "loglog",
            },
            "non-simple": {
                "coefficients": [-5.10327300754197, 1.8991299658269363, -0.0711637556510648],
                "c": -6.080795135614891,
                "alpha": 2.0966521781546383,
                "fit_scale": "loglog",
            },
        }
        rim_width_params = {
            "simple": {
                "coefficients": [-3.0877789815769905, 1.4473834897759599, -0.036680803391976054],
                "c": -5.504615471104058,
                "alpha": 2.1143223655406027,
                "fit_scale": "loglog",
            },
            "non-simple": {
                "coefficients": [-0.13959296561846923, 0.6920030607817031, 0.006073402057851526],
                "c": -6.197739434603876,
                "alpha": 2.3569121813546987,
                "fit_scale": "loglog",
            },
        }
        floor_radius_params = {
            "simple": {
                "coefficients": [-0.900242255950584, 0.455179909117242, 0.04575920054549337],
                "c": -0.3616562884185268,
                "alpha": 1.4684589969693498,
                "fit_scale": "loglog",
            },
            "non-simple": {
                "coefficients": [-12.34730111609687, 2.8374446646267857, -0.07447326061389413],
                "c": -2.377257077688099,
                "alpha": 1.752512234751084,
                "fit_scale": "loglog",
            },
        }
        wall_curvature_params = {
            "simple": {
                "coefficients": [-0.25824388040970825, 0.17625334392925013, -0.014915775447250382],
                "c": -0.008556260200970858,
                "alpha": 0.10594060062350721,
                "fit_scale": "logx",
            },
            "non-simple": {
                "coefficients": [0.18987288872317, -0.08418984102324649, 0.007551647630654395],
                "c": 0.001958457521099688,
                "alpha": 0.04951276138707961,
                "fit_scale": "logx",
            },
        }
        floor_blend_params = {
            "simple": {
                "coefficients": [1.059340070020194, -0.13182446192525693, 0.00902946332988967],
                "c": 0.1404880414615158,
                "alpha": -0.17880038251603678,
                "fit_scale": "logx",
            },
            "non-simple": {
                "coefficients": [0.5808979959114502, 0.03860715120274549, -0.0037083300328609988],
                "c": -0.05625078483323789,
                "alpha": 0.13327426573917356,
                "fit_scale": "logx",
            },
        }
        peak_height_params = {
            "non-simple": {
                "coefficients": [3446.114633388009, -888.9879125470408, 58.43540553673689],
                "c": 49177.29323139237,
                "alpha": 108.72351648253827,
                "fit_scale": "logx",
            },
        }
        peak_width_params = {
            "non-simple": {
                "coefficients": [370074.37674616196, -73224.92406382009, 3630.2460325568773],
                "c": 3843578.4915611,
                "alpha": 252.2461102371035,
                "fit_scale": "logx",
            },
        }
        ring_radius_params = {
            "non-simple": {
                "coefficients": [-6.613689951394371, 1.976602153637156, -0.042115405413714696],
                "c": -3.610867824558726,
                "alpha": 1.7741307984706642,
                "fit_scale": "loglog",
            },
        }
        ring_rim_height_params = {
            "non-simple": {
                "coefficients": [-31.93167277655391, 5.496001921151941, -0.19131828795098776],
                "c": -0.3766269768654428,
                "alpha": 1.7427748575666637,
                "fit_scale": "loglog",
            },
        }
        ring_floor_radius_params = {
            "non-simple": {
                "coefficients": [-5.551393952731496, 1.6979145956073238, -0.029252203610892827],
                "c": -6.467424320735284,
                "alpha": 2.139551688327142,
                "fit_scale": "loglog",
            },
        }
        ring_rim_width_params = {
            "non-simple": {
                "coefficients": [-15.962694436765481, 3.2135115326621557, -0.08918699188536301],
                "c": -1.0673993551823213,
                "alpha": 1.7610723251290876,
                "fit_scale": "loglog",
            },
        }

        args = {}
        diameter_m = min(max(crater.diameter, morphology.valid_diameter_range[0]), morphology.valid_diameter_range[1])
        fcorrection = crater.diameter / diameter_m

        if crater.morphology_type in ["complex", "transitional", "basin", "multiring", "peakring", "ring"]:
            fit_category = "non-simple"
        else:
            fit_category = "simple"

        if rim_height is None:
            rim_height = (
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_height_params[fit_category]
                )[0]
                * fcorrection
            )
        args["rim_height"] = rim_height

        if rim_width is None:
            rim_width = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_width_params[fit_category]
                )[0]
                * fcorrection,
                0.0,
            )
        args["rim_width"] = rim_width

        # Try to approximately conserve volume when setting the ejecta thickness at the rim value
        args["ejprofile"] = _EJPROFILE if ejprofile is None else ejprofile

        # This is an initial guess of the ejecta rim. We will adjust it later by integrating the volume of the crater and ejecta profiles
        if frac_ejrim is None:
            if rim_height > 0.0:
                if monte_carlo:
                    frac_ejrim = rng.uniform(low=0.3, high=0.7, size=1)[0]  # Temporary until a morphometric analysis ic complete
                else:
                    frac_ejrim = 0.5
            else:
                frac_ejrim = 0.0
        args["frac_ejrim"] = min(max(frac_ejrim, 0.0), 1.0)

        if floor_elevation is None:
            floor_elevation = (
                -sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_elevation_params[fit_category]
                )[0]
                * fcorrection
            )
            floor_elevation = min(floor_elevation, 0.0)
        args["floor_elevation"] = floor_elevation

        if floor_radius is None:
            floor_radius = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_radius_params[fit_category]
                )[0]
                * fcorrection,
                0.0,
            )
        args["floor_radius"] = min(floor_radius, 0.8 * crater.radius)

        if peak_height is None:
            if crater.morphology_type == "complex":
                peak_height = max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **peak_height_params[fit_category]
                    )[0],
                    0.0,
                )
            else:
                peak_height = 0.0
        args["peak_height"] = peak_height

        if peak_width is None:
            if crater.morphology_type == "complex":
                peak_width = max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **peak_width_params[fit_category]
                    )[0],
                    0.0,
                )
            else:
                peak_width = 0.0
        args["peak_width"] = peak_width
        args["peak_ring_radius"] = 0.0 if peak_ring_radius is None else peak_ring_radius
        args["peak_center_distance"] = 0.0 if peak_center_distance is None else peak_center_distance
        args["peak_center_bearing"] = 0.0 if peak_center_bearing is None else peak_center_bearing

        if wall_curvature is None:
            wall_curvature = min(
                max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **wall_curvature_params[fit_category]
                    )[0],
                    0.0,
                ),
                1.0,
            )

        args["wall_curvature"] = wall_curvature

        if floor_blend is None:
            floor_blend = min(
                max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_blend_params[fit_category]
                    )[0],
                    0.0,
                ),
                1.0,
            )

        args["floor_blend"] = floor_blend

        kwargs = {**args, **kwargs}

        crater = cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )

        # A copy that does not make new rings keeps copies of the rings of the crater it was made from
        if not make_rings and crater.nrings == 0 and isinstance(source, MorphologyCrater) and source.nrings > 0:
            object.__setattr__(crater._var, "_rings", [type(ring)(crater=ring, morphology=morphology) for ring in source.rings])

        num_new_rings = 0
        if make_rings and not crater.isring and crater.nrings == 0:
            if crater.morphology_type == "complex":
                num_new_rings = 1
            elif crater.morphology_type == "multiring":
                num_new_rings = 3

        if num_new_rings > 0:
            for i in range(num_new_rings):
                if i == 0:
                    radius_outer = crater.radius
                    floor_radius_outer = crater.floor_radius
                else:
                    radius_outer = crater.rings[i - 1].radius
                    floor_radius_outer = crater.rings[i - 1].floor_radius
                radius = min(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **ring_radius_params[fit_category]
                    )[0]
                    * radius_outer
                    / crater.radius,
                    radius_outer,
                )
                floor_radius = min(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **ring_floor_radius_params[fit_category]
                    )[0]
                    * floor_radius_outer
                    / crater.floor_radius,
                    floor_radius_outer,
                )
                rim_height = max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **ring_rim_height_params[fit_category]
                    )[0],
                    0.0,
                )
                if i == 0:
                    frac_ejrim = crater.frac_ejrim * rim_height / crater.rim_height
                else:
                    frac_ejrim = 0.0
                crater.add_ring(
                    radius=radius,
                    floor_radius=floor_radius,
                    rim_height=rim_height,
                    frac_ejrim=frac_ejrim,
                )
        elif crater.nrings > 0:
            # Make sure rings are the correct type
            for i in range(crater.nrings):
                ring = crater.rings[i]
                if not isinstance(ring, cls):
                    ring = cls(crater=ring, morphology=morphology)
                    crater.rings[i] = ring

        # Adjust frac_ejrim value(s) in order to get closer to a volume-conserving solution for the ejecta
        if conserve_volume and not crater.isring:
            # The search only changes the rim height and the ejecta profile exponent, so the volumes are evaluated on a copy of
            # the profile parameters and the crater is rebuilt once with the final values.
            params = SimpleNamespace(**{name: getattr(crater, name) for name in _PROFILE_PARAMETERS})

            def ejvol(p):
                return 2 * np.pi * p.frac_ejrim * p.rim_height * p.radius**2 / (-p.ejprofile - 2)

            for _ in range(10):
                excavated_volume = morphology.estimate_volume(params, include_crater=True, include_ejecta=False)
                while excavated_volume > 0.0:
                    # This occurs when rim_height is too high
                    kwargs["rim_height"] *= 0.9
                    params.rim_height = kwargs["rim_height"]
                    excavated_volume = morphology.estimate_volume(params, include_crater=True, include_ejecta=False)

                conservation_factor = -ejvol(params) / excavated_volume

                ejprofile = -2.0 + conservation_factor * (params.ejprofile + 2)
                if ejprofile < -4:
                    kwargs["rim_height"] *= 0.9
                else:
                    kwargs["ejprofile"] = ejprofile
                params.rim_height = kwargs["rim_height"]
                params.ejprofile = kwargs["ejprofile"]
                if ejprofile > -4:
                    break
            crater = cls(crater=crater, morphology=morphology, **kwargs)

        return crater

    def as_dict(self, ignore_keys: list[str] | tuple[str] = (), skip_complex_data: bool = False, **kwargs) -> dict:
        """
        Return a dictionary representation of the crater properties, including morphology-specific variable properties.

        Parameters
        ----------
        ignore_keys : list[str] or tuple[str], optional
            A list or tuple of property names to ignore when creating the dictionary representation. Default is an empty tuple unless `skip_complex_data` is True, in which case it will be extended to include ("morphology", "affected_face_indices", "affected_node_indices", "ejecta_region", "ejecta_region", "crater_region").
        skip_complex_data : bool, optional
            If True, skip complex data types when creating the dictionary representation. This is useful when serializing the object for saving to a file, as it removes complex data types that may not be serializable. Default is False.
        """
        if skip_complex_data:
            ignore_keys += ("rings",)
        dict_repr = super().as_dict(ignore_keys=ignore_keys, skip_complex_data=skip_complex_data, **kwargs)
        return dict_repr

    def add_ring(
        self,
        ring: BasicMoonCrater | None = None,
        radius: float | None = None,
        floor_radius: float | None = None,
        wall_curvature: float | None = None,
        frac_ejrim: float | None = None,
        rim_width: float | None = None,
        rim_height: float | None = None,
        **kwargs: Any,
    ):
        """
        Add a ring to the crater.

        Parameters
        ----------
        ring: BasicMoonCrater, optional
            An existing Crater object that will be converted into a ring object.
        radius : float, optional
            The radius of the ring in meters.
        floor_radius : float, optional
            The floor radius of the ring in meters.
        wall_curvature : float, optional
            The wall curvature of the ring wall (between 0 for straight walls and 1 for very curvy walls).
        frac_ejrim: float, optional
            The fraction of the ring rim_height that is ejecta (usually 0 except for maybe the outerost ring of a multiring basin).
        rim_width : float, optional
            The rim width of the ring in meters.
        rim_height : float, optional
            The rim elevation of the ring in meters.
        **kwargs : Any
            Additional keyword arguments that are passed to the .maker() method. Any arguments that are valid for a Crater ar valid for a ring. Otherwise the ring properties are copied from its associated crater.
        """
        if ring is not None:
            crater = ring
            radius = ring.radius if radius is None else radius
            floor_radius = ring.floor_radius if floor_radius is None else floor_radius
            wall_curvature = ring.wall_curvature if wall_curvature is None else wall_curvature
            frac_ejrim = ring.frac_ejrim if frac_ejrim is None else frac_ejrim
            rim_width = ring.rim_width if rim_width is None else rim_width
            rim_height = ring.rim_height if rim_height is None else rim_height
        else:
            crater = self
        frac_ejrim = 0.0 if frac_ejrim is None else frac_ejrim

        # Link together parameters that are controlled by the main crater feature by passing them in as argtuments
        newring = self.__class__.maker(
            crater=crater,
            morphology_type="ring",
            parent=self.id,
            morphology=self.morphology,
            radius=radius,
            floor_radius=floor_radius,
            wall_curvature=wall_curvature,
            rim_width=rim_width,
            rim_height=rim_height,
            frac_ejrim=frac_ejrim,
            isring=True,
            conserve_volume=False,
            **kwargs,
        )
        if newring.radius > self.radius:
            raise ValueError("Ring radius cannot be larger than the crater radius!")

        # If there are rings, we need to insert this ring into the correct order. Rings are added in descending order by radius
        if self.nrings == 0:
            self._rings.append(newring)
        else:
            oldtot = self.nrings
            for i in range(oldtot):
                if newring.radius > self.rings[i].radius:
                    self._rings.insert(i, newring)
                    break
                if i == oldtot - 1:
                    self._rings.append(newring)

        return

    @property
    def rim_height(self) -> float | None:
        """Original rim height of the crater in meters relative to the reference surface."""
        return self._fixed.rim_height

    @property
    def floor_elevation(self) -> float | None:
        """Original floor depth of the crater in meters relative to the reference surface."""
        return self._fixed.floor_elevation

    @property
    def floor_radius(self) -> float | None:
        """Original floor diameter of the crater in meters."""
        return self._fixed.floor_radius

    @property
    def wall_curvature(self) -> float | None:
        """A factor that controls the curvature of the crater wall (0 for straight walls, and 1.0 for very curvy walls). If None, it will be computed based on the morphology type and diameter."""
        return self._fixed.wall_curvature

    @property
    def floor_blend(self) -> float | None:
        """A factor that controls how sharply the wall and floor blend (0 for sharp floor-to-wall transition, 1.0 for gentle floor-to-wall transition)."""
        return self._fixed.floor_blend

    @property
    def rim_width(self) -> float | None:
        """The width of the crater rim in meters."""
        return self._fixed.rim_width

    @property
    def ejprofile(self) -> float | None:
        """Power law exponent for the ejecta thickness profile of the crater."""
        return self._fixed.ejprofile

    @property
    def peak_height(self) -> float | None:
        """Central peak height of the crater in meters relative to the reference surface. 0 for simple craters."""
        return self._fixed.peak_height

    @property
    def peak_width(self) -> float | None:
        """Central peak width of the crater in meters. 0 for simple craters."""
        return self._fixed.peak_width

    @property
    def peak_ring_radius(self) -> float | None:
        """Radius of peak ring in meters. 0 for central peaks."""
        return self._fixed.peak_ring_radius

    @property
    def peak_center_distance(self) -> float | None:
        """Distance of central peak/peak ring from crater center in meters."""
        return self._fixed.peak_center_distance

    @property
    def peak_center_bearing(self) -> float | None:
        """Bearing angle of central peak/peak ring in degrees."""
        return self._fixed.peak_center_bearing

    @property
    def isring(self) -> bool | None:
        """Flag that indicates that this is a ring rather than a crater."""
        return self._fixed.isring

    @property
    def frac_ejrim(self) -> float | None:
        """Ejecta rim thickness of the crater in meters."""
        return self._var._frac_ejrim

    @property
    def rings(self) -> list[BasicMoonCrater]:
        return self._var.rings

    @property
    def nrings(self) -> int:
        """
        Returns the number of rings associated with this Crater.
        """
        return self._var.nrings


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

    # The compiled function that forms a batch of craters of this model's crater type. Models with a different crater type
    # override this.
    _form_craters_binding = staticmethod(basicmoon_bindings.form_craters)

    # The angle of repose used for slope collapse after a crater forms, in degrees
    _CRITICAL_SLOPE_ANGLE = 35.0

    def form_crater(self, crater: Crater, **kwargs: Any) -> None:
        """
        Form a crater and its ejecta blanket on the surface.

        This alters the face and node elevations and the ejecta thickness of the surface around the crater: the crater is drawn
        relative to a reference plane fitted to the pre-impact surface, the ejecta is scaled to conserve volume (if
        `conserve_volume` is set), the pre-existing surface is degraded by ejecta burial and by the crater's degradation function,
        and steep slopes are collapsed (if `do_slope_collapse` is set).

        A crater whose shape removes no net volume from the surface mesh is not formed, and is not recorded by the counting
        component. This happens to some craters that are only a few faces across, where the sampled rim uplift outweighs the
        sampled bowl, and keeps the craters and their ejecta from adding volume to the surface.

        Parameters
        ----------
        crater : Crater
            The crater object to be formed.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        This is a wrapper for a compiled Rust function that runs the whole formation sequence in a single call.
        """
        self.form_craters([crater], **kwargs)
        return

    def form_craters(self, craters: list[Crater], **kwargs: Any) -> None:
        """
        Form a batch of craters whose regions do not overlap.

        Each crater is formed as described in `form_crater`. The bookkeeping that must happen in order (drawing each crater's
        random numbers, and then recording the craters that were formed with the counting component) is done one crater at a
        time in Python, and all of the craters are formed in parallel by a single compiled call in between. Because their
        regions do not overlap, the result is the same as forming them one after another.

        Parameters
        ----------
        craters : list of Crater
            The craters to be formed. Their regions (see `overlap_function`) must not share any faces or nodes.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        This is a wrapper for a compiled Rust function that forms the craters in parallel.
        """
        formed = []
        regions = []
        ray_seeds = []
        for crater in craters:
            if not isinstance(crater, Crater):
                raise TypeError("crater must be an instance of Crater")
            if not isinstance(crater, self._CraterType):
                crater = self.Crater.maker(crater=crater, **kwargs)
            if not crater.emplaceable:
                continue

            region = crater.crater_region
            crater_area = math.pi * region.radius**2

            # Skip craters that are smaller than the face they are centered on
            if crater_area <= self.surface.face_area[crater.face_index]:
                crater.reset()
                continue

            # The crater and ejecta are computed on the same set of faces and nodes
            ejecta_region = crater.ejecta_region
            if ejecta_region is not None and ejecta_region is not region:
                if ejecta_region.n_face != region.n_face or ejecta_region.n_node != region.n_node:
                    raise ValueError("The crater and ejecta regions of a crater must contain the same faces and nodes")

            # Any crater data that is generated lazily from the random number generator must exist before the ray seed is drawn,
            # so that random numbers are drawn in the same order as when the crater shape was computed in Python
            self._prepare_crater(crater)

            formed.append(crater)
            regions.append(region)
            ray_seeds.append(int(self.rng.integers(0, 2**32 - 1)) if self.dorays else None)

        if not formed:
            return

        surface = regions[0].surface
        if "ejecta_thickness" not in surface.uxds.data_vars:
            surface._add_new_data("ejecta_thickness", long_name="ejecta thickness", units="m", data=0.0, isfacedata=True)
        variables = surface.uxds.variables

        def _global_indices(indices, n_total):
            if isinstance(indices, slice):
                indices = np.arange(n_total)[indices]
            return np.ascontiguousarray(indices, dtype=np.int64)

        results = self._form_craters_binding(
            craters=formed,
            rings=[crater.rings for crater in formed],
            regions=regions,
            face_indices=[_global_indices(region.face_indices, surface.n_face) for region in regions],
            node_indices=[_global_indices(region.node_indices, surface.n_node) for region in regions],
            ray_seeds=ray_seeds,
            face_elevation=variables["face_elevation"].data,
            node_elevation=variables["node_elevation"].data,
            ejecta_thickness=variables["ejecta_thickness"].data,
            conserve_volume=kwargs.get("conserve_volume", self.conserve_volume),
            critical_slope=float(np.tan(np.deg2rad(self._CRITICAL_SLOPE_ANGLE))) if self.do_slope_collapse else None,
        )

        for crater, (_, conservation_factor) in zip(formed, results, strict=True):
            # A crater that removes no net volume on the mesh (one the mesh barely resolves) is not formed, and is not recorded
            if conservation_factor is not None:
                self._scale_ejecta_rim(crater, conservation_factor)
                # Record the crater to the counting layer, in emplacement order. This does not depend on the surface elevation
                # or ejecta thickness, so it can be done after the batch is formed.
                if self.counting is not None:
                    self.counting.add(crater, add_to_observed=self.save_observed, add_to_emplaced=self.save_emplaced, **kwargs)
            # Remove the saved regions now that the crater has been processed
            crater.reset()
        return

    def _prepare_crater(self, crater: BasicMoonCrater) -> None:
        """
        Generate any crater data that is created lazily from the morphology's random number generator.

        This is called by `form_craters` before it draws any random numbers of its own. Models whose craters draw random numbers
        on first access to some of their data override this to access that data, in the order the compiled `form_craters` reads
        it. The basic model has no such data.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater about to be formed.
        """
        return

    def crater_profile(
        self,
        crater: BasicMoonCrater,
        radial_distances: ArrayLike,
        bearings: ArrayLike | None = None,
        reference_elevations: ArrayLike | None = None,
        crater_cls: type[Crater] = BasicMoonCrater,
        profile_func: Callable = basicmoon_bindings.basicmoon_profile,
        include_crater: bool = True,
        include_ejecta: bool = False,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the crater profile.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings : ArrayLike, optional.
            Bearings (in degrees) corresponding to the radial distances for profiles that are non-axisymmetric.
        reference_elevations : ArrayLike, optional
            Reference elevation values to be modified by the crater profile.
        crater_cls : type[Crater], optional
            The class of the crater type used. If the crater object doesn't match, then it is cast as this type. Default is BasicMoonCrater.
        profile_func: Callable, optional
            The backend function used to draw the crater profile. Default is bhe basicmoon_profile from the basicmoon_bindings Rust library.
        include_crater : bool, optional
            Draws the crater profile without the ejecta. Default is True.
        include_ejecta : bool, optional
            Draws the ejecta profile without the crater. Default is False.
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
        if not isinstance(crater, crater_cls):
            crater = crater_cls.maker(crater, morphology=self)
        if reference_elevations is None:
            reference_elevations = np.zeros_like(radial_distances)

        if np.isscalar(radial_distances):
            radial_distances = np.array([radial_distances], dtype=np.float64)
        elif isinstance(radial_distances, (list | tuple)):
            radial_distances = np.array(radial_distances, dtype=np.float64)

        # flatten r to 1D array
        orig_shape = radial_distances.shape
        radial_distances = np.ravel(radial_distances)
        reference_elevations = np.ravel(reference_elevations)
        if bearings is None:
            bearings = np.zeros_like(radial_distances)
        else:
            bearings = np.ravel(np.radians(bearings))

        elevation = profile_func(
            radial_distances=radial_distances,
            bearings=bearings,
            reference_elevations=reference_elevations,
            crater=crater,
            rings=crater.rings,
            include_crater=include_crater,
            include_ejecta=include_ejecta,
        )
        elevation = np.array(elevation, dtype=np.float64)
        elevation = np.reshape(elevation, orig_shape)

        return elevation

    def ejecta_profile(
        self,
        crater: BasicMoonCrater,
        radial_distances: ArrayLike,
        bearings: ArrayLike | None = None,
        crater_cls: type[Crater] = BasicMoonCrater,
        profile_func: Callable = basicmoon_bindings.basicmoon_profile,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta profile.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings : ArrayLike, optional.
            Bearings (in degrees) corresponding to the radial distances for profiles that are non-axisymmetric.
        crater_cls : type[Crater], optional
            The class of the crater type used. If the crater object doesn't match, then it is cast as this type. Default is BasicMoonCrater.
        profile_func: Callable, optional
            The backend function used to draw the crater profile. Default is bhe basicmoon_profile from the basicmoon_bindings Rust library.
        kwargs : Any
            |kwargs|

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed ejecta profile at each radial point.
        **kwargs : Any
            |kwargs|

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        return self.crater_profile(
            crater=crater,
            radial_distances=radial_distances,
            bearings=bearings,
            crater_cls=crater_cls,
            profile_func=profile_func,
            include_crater=False,
            include_ejecta=True,
            **kwargs,
        )

    def ejecta_distribution(
        self, crater, radial_distances: ArrayLike, bearings: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the ejecta thickness distribution modulated by ray patterns.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ejecta distribution.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings : ArrayLike
            Angular bearingss from the crater center (in degrees).

        Returns
        -------
        thickness : NDArray[np.float64]
            The computed ejecta thickness for each (radial_distance, bearing) pair.
        ray_intensity : NDArray[np.float64]
            The computed ray intensity for each (radial_distance, bearing) pair.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, self._CraterType):
            crater = self.Crater.maker(crater=crater)
        thickness = self.ejecta_profile(crater, radial_distances, bearings=bearings)
        radial_distances = np.ravel(radial_distances)
        bearings = np.ravel(np.radians(bearings))
        intensity = basicmoon_bindings.ray_intensity(
            radial_distances=radial_distances,
            bearings=bearings,
            crater_diameter=crater.diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        thickness = np.array(thickness, dtype=np.float64)
        intensity = np.array(intensity, dtype=np.float64)
        thickness *= intensity
        return np.reshape(thickness, radial_distances.shape), np.reshape(intensity, radial_distances.shape)

    def ray_intensity(self, crater: BasicMoonCrater, radial_distances: ArrayLike, bearings: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ray pattern intensity modulation at each (r, theta) pair.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object containing the parameters for the ray intensity.
        radial_distances : ArrayLike
            Radial distances from the crater center (in meters).
        bearings: ArrayLike
            Angular bearings from the crater center (in degrees).

        Returns
        -------
        intensity : NDArray[np.float64]
            The computed ray intensity values.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, self._CraterType):
            crater = self.Crater.maker(crater=crater)
        radial_distances = np.ravel(radial_distances)
        bearings = np.ravel(np.radians(bearings))
        intensity = basicmoon_bindings.ray_intensity(
            radial_distances,
            bearings,
            crater.diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        intensity = np.array(intensity, dtype=np.float64)
        intensity = np.reshape(intensity, radial_distances.shape)
        return intensity

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
            n = np.asarray(
                self.production.function(
                    diameter=diameter,
                    time_start=time_start,
                    time_end=time_end,
                    validate_inputs=False,
                )
            ).reshape(diameter.shape)
            degradation_region_area = np.pi * (diameter / 2) ** 2 * fe
            return k * n * degradation_region_area

        for face_indices, dc_max in zip(self.surface.face_bin_indices, self.surface.face_bin_max_sizes, strict=False):
            # The integrand spans many orders of magnitude in diameter, so integrate in log(D) with a composite Gauss-Legendre rule
            # (two 8-point subintervals per decade). All points are evaluated in a single vectorized call.
            n_intervals = max(1, math.ceil(2 * math.log10(dc_max / dc_min)))
            log_d, weights = _gauss_legendre_points(np.linspace(math.log(dc_min), math.log(dc_max), n_intervals + 1))
            diameters = np.exp(log_d)
            delta_kdiff = np.sum(weights * diameters * _subpixel_degradation(diameters))
            self._Kdiff[face_indices] += delta_kdiff

        # If any Kdiff values reaches a threshold where a meaningful amount of diffusion will occur on the surface, then go ahead and apply it.
        # Otherwise, degradation will continue to accumulate until the next batch of craters is processed.
        if np.any(self._Kdiff / self.surface.face_area > 1):
            self.apply_subpixel_degradation()
            if self.save_observed:
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
        if not isinstance(crater, self._CraterType):
            crater = self.Crater.maker(crater=crater)

        if feature not in ("crater", "ejecta"):
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")

        # Get the maximum extent
        lower_limit = crater.radius * 1.0001
        upper_limit = self.ejecta_truncation * crater.radius if self.ejecta_truncation else np.pi * self.surface.target.radius

        # The profile crossing is found with a compiled port of scipy.optimize.brentq
        return float(
            basicmoon_bindings.rmax(
                crater=crater,
                minimum_thickness=float(minimum_thickness),
                feature=feature,
                lower_limit=float(lower_limit),
                upper_limit=float(upper_limit),
            )
        )

    def estimate_volume(
        self,
        crater: BasicMoonCrater,
        include_crater: bool,
        include_ejecta: bool,
        crater_profile_function: Callable = basicmoon_bindings.crater_profile_function,
        ejecta_profile_function: Callable = basicmoon_bindings.ejecta_profile_function,
        **kwargs,
    ) -> np.float64:
        """
        Estimates the volume change of a crater by integrating its axisymmetric profile.

        With the default profile functions, the integral is computed with a composite Gauss-Legendre rule whose subintervals are
        split at the radii where the profile changes character (floor blend, rim, central peak), and whose integrand is evaluated
        in a single vectorized call. If custom profile functions or quadrature keyword arguments are given, scipy.integrate.quad
        is used instead.
        """
        if not (include_crater or include_ejecta):
            return np.float64(0.0)

        if self.ejecta_truncation is None:
            rmax = 20 * crater.radius
        else:
            rmax = self.ejecta_truncation * crater.radius

        if (
            crater_profile_function is basicmoon_bindings.crater_profile_function
            and ejecta_profile_function is basicmoon_bindings.ejecta_profile_function
            and not kwargs
        ):
            edges = self._volume_quadrature_edges(crater, include_crater, include_ejecta, rmax)
            if edges is not None:
                r, w = _gauss_legendre_points(edges)
                integrand = basicmoon_bindings.volume_integrand(r, crater, include_crater, include_ejecta)
                return np.float64(2 * np.pi * np.sum(w * integrand))

        def _crater_func(r):
            hc = crater_profile_function(crater, r)
            he = ejecta_profile_function(crater, r)

            return r * (hc - he)

        def _ejecta_func(r):
            h = ejecta_profile_function(crater, r)
            return r * h

        def _combo_func(r):
            hc = crater_profile_function(crater, r)
            return r * hc

        if include_crater and not include_ejecta:
            func = _crater_func
        elif include_ejecta and not include_crater:
            func = _ejecta_func
        else:
            func = _combo_func
        v = quad(func, 0.0, rmax, **{"limit": 100, "epsrel": 1e-3, **kwargs}, full_output=1)[0]
        return 2 * np.pi * v

    @staticmethod
    def _volume_quadrature_edges(
        crater: BasicMoonCrater, include_crater: bool, include_ejecta: bool, rmax: float
    ) -> NDArray[np.float64] | None:
        """
        Subinterval boundaries for the composite Gauss-Legendre volume integral in `estimate_volume`.

        The boundaries are placed where the profile (see `crater_profile_function` in the Rust crate) switches between pieces or
        blends between them, so that the integrand is smooth within each subinterval:

        - the floor-to-wall blend, centered on the floor radius with half-width ``rfw``
        - the wall-to-rim blend and the start of the ejecta, centered on the rim radius with half-width ``rim_width``
        - the central peak, if present

        Beyond ``radius + rim_width``, the crater-only integrand is purely the Gaussian rim, so it is truncated at 12 rim widths
        (where it is below 1e-31 of its peak). When the ejecta is included, the power-law tail is covered with logarithmically
        spaced subintervals out to `rmax`. Each subinterval is then split in two.

        Returns None if the crater parameters do not produce a valid set of boundaries.
        """
        radius = crater.radius
        floor_radius = crater.floor_radius
        rim_width = crater.rim_width
        floor_blend_width = (floor_radius if floor_radius < 0.5 * radius else radius - floor_radius) * crater.floor_blend
        rim_outer = radius + rim_width

        points = [0.0, floor_radius - floor_blend_width, floor_radius, floor_radius + floor_blend_width, radius - rim_width, radius]
        if crater.peak_width > 0:
            points += [crater.peak_ring_radius + k * crater.peak_width for k in (-3.0, 0.0, 3.0)]

        if include_crater and not include_ejecta:
            top = min(rmax, radius + 12 * rim_width)
            if rim_outer < top:
                points += list(np.linspace(rim_outer, top, 4))
        else:
            top = rmax
            if 0 < rim_outer < top:
                points += list(np.geomspace(rim_outer, top, 12))
        points.append(top)

        points = np.asarray(points, dtype=np.float64)
        if not (np.all(np.isfinite(points)) and top > 0):
            return None
        edges = np.unique(np.clip(points, 0.0, top))
        return np.sort(np.concatenate((edges, 0.5 * (edges[1:] + edges[:-1]))))

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

        if np.ndim(diameter) == 0:
            return basicmoon_bindings.degradation_function(float(diameter), float(fe))
        diameter = np.asarray(diameter, dtype=np.float64)
        degradation = basicmoon_bindings.degradation_function_array(np.ascontiguousarray(diameter.ravel()), float(fe))
        return degradation.reshape(diameter.shape)

    @staticmethod
    def overlap_function(crater: BasicMoonCrater) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        """
        Get the affected node and face indices for the crater.

        Parameters
        ----------
        crater : BasicMoonCrater
            The crater object to be used.

        Returns
        -------
        tuple[NDArray[np.int64], NDArray[np.int64]]
            The affected global node and face indices for the crater.
        """
        if not crater.emplaceable or crater.affected_node_indices is None or crater.affected_face_indices is None:
            empty = np.empty(0, dtype=np.int64)
            return empty, empty
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

    @property
    def valid_diameter_range(self) -> tuple[float, float]:
        """
        Valid diameter range for the morphometric property models.

        Craters with diameters outside the range will use the morphometric property models at the nearest boundary diameter.
        """
        return (100.0, 1000.0e3)
