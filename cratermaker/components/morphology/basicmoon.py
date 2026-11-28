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
_FRAC_EJRIM = 0.3
# The profiles of the largest craters. Above the range of the fresh craters the non-simple size relations are fitted to (x_max,
# 138 km), the profile parameters (in proportion to the diameter or radius) move from those of the relations at x_max to those
# of the peak-ring basin Schrodinger (340 km) and the multiring basin Orientale (989 km), interpolated in ln D, and keep those of
# Orientale beyond it. The values are the medians of fits of the basic profile (with a Gaussian peak ring, and for Orientale its
# ring scarp, see _SCARP) to narrow wedges of their H-field extracts out to 1.5 rim radii, avoiding younger craters and ray
# trenches, each wedge with its own background level (Lunar_Morphometry_Revisited, hfield_wedge_fit.py). The relief (rim height
# plus floor depth) of a wedge is well determined, but not how it divides between rim and floor: the relief is the median of the
# wedges and the rim height is 0.31 of it, the median of the wedges of Orientale (that of Schrodinger is poorly constrained).
# Both basins are degraded (freshness class D) and their floors hold later mare and volcanic deposits, so their
# floor depths and peak-ring heights are lower bounds.
_BASIN_ANCHORS = {
    340.2e3: {  # Schrodinger: relief 0.0134 D
        "rim_height_over_D": 0.00412,
        "floor_depth_over_D": 0.00927,
        "floor_radius_over_R": 0.735,
        "rim_width_over_R": 0.144,
        "wall_curvature": 0.609,
        "floor_blend": 0.721,
    },
    988.6e3: {  # Orientale: relief 0.00588 D
        "rim_height_over_D": 0.00181,
        "floor_depth_over_D": 0.00407,
        "floor_radius_over_R": 0.666,
        "rim_width_over_R": 0.109,
        "wall_curvature": 0.0,
        "floor_blend": 0.939,
    },
}

# The ejecta of multiring basins: that of the profile fit of Orientale out to 2 rim radii. A basin's final surface keeps little of
# its excavation (Orientale has 0.36e15 m^3 below the surrounding level, against 2.5e15 m^3 above it out to two radii: the floor
# rebounds and is compensated), so its ejecta cannot be balanced against it, and craters with a ring scarp (multiring basins) are
# not held to volume conservation.
_MULTIRING_EJECTA = {"frac_ejrim": 0.4047612450358902, "ejprofile": -2.8211551267661767}

# The ring scarp of multiring basins: the nested ring fault (Orientale's Outer Rook ring), an asymmetric ridge with an
# inward-facing face (see basicmoon.rs, scarp_function). Its radius, height and face width are the medians of the wedge fits of
# Orientale (see _BASIN_ANCHORS); as Orientale is degraded, the height is a lower bound. Its radius (0.682 of the rim radius) has
# the lognormal scatter of the middle rings of the multiring basins of Neumann et al. (2015) (0.68-0.75 of the main ring in their
# frame, where Orientale's Outer Rook ring is at 0.682). Multiring basins appear on the Moon from 500 to 650 km (Hertzsprung,
# 571 km, is the smallest; Neumann et al. 2015, Johnson et al. 2018): the height of the scarp rises from 0 to its full value over
# that range (a smooth step in ln D).
_SCARP = {
    "radius_over_R": 0.682,
    "radius_sigma_ln": 0.04,
    "height_over_D": 0.00210,
    "width_over_R": 0.074,
    "onset": (500e3, 650e3),
}

# Peak rings. Central peaks open into peak rings as craters grow from 200 to 300 km (a smooth step in ln D). The ring radius, as
# a fraction of the rim radius, is 0.460 (D / 340 km)^0.107 with lognormal scatter 0.086: the slope and scatter are those of the
# inner rings of the peak-ring and multiring basins of Neumann et al. (2015), whose ratio is nearly constant (the innermost ring of
# a multiring basin is its peak ring), and the normalization is the geometric mean of those implied by the peak rings of
# Schrodinger (0.476 at 340 km) and Orientale (its Inner Rook ring, 0.498 at 989 km) in the wedge fits of _BASIN_ANCHORS (their
# ratio, 1.046, is below the 1.121 of the slope; the slope of the population is kept).
_PEAK_RING_ONSET = (200e3, 300e3)
_PEAK_RING_RATIO = 0.460
_PEAK_RING_D_REF = 340.2e3
_PEAK_RING_EXPONENT = 0.107
_PEAK_RING_SIGMA_LN = 0.086
# The width (standard deviation) of the Gaussian peak ring of the basic profile, as a fraction of the rim radius: that of the
# peak ring of Schrodinger (Lunar_Morphometry_Revisited, peakring_measure.py)
_PEAK_RING_WIDTH = 0.084


# The largest crater whose rings are measured (Orientale). Beyond it, the ring ratios are held at their values for Orientale
# rather than extrapolated, and the scatter of the ring radii is held at its size for Orientale in m rather than in proportion
# to the radius, so that megabasins (South Pole-Aitken) do not get ring positions that vary by tens of km.
_RING_DATA_MAX_DIAMETER = 988.6e3


def _ring_scatter(diameter: float, sigma_ln: float, rng: Generator | None) -> float:
    """The lognormal scatter factor of a ring radius, whose spread in m stops growing beyond _RING_DATA_MAX_DIAMETER."""
    if rng is None:
        return 1.0
    return math.exp(sigma_ln * min(1.0, _RING_DATA_MAX_DIAMETER / diameter) * rng.standard_normal())


def _smoothstep(x: float) -> float:
    """The smooth step 3x^2 - 2x^3 between 0 (x <= 0) and 1 (x >= 1)."""
    x = min(max(x, 0.0), 1.0)
    return x * x * (3.0 - 2.0 * x)


def _log_step(diameter: float, bounds: tuple[float, float]) -> float:
    """A smooth step in ln D from 0 at bounds[0] to 1 at bounds[1]."""
    return _smoothstep(math.log(diameter / bounds[0]) / math.log(bounds[1] / bounds[0]))


def peak_ring_layout(diameter: float, rng: Generator | None = None) -> tuple[float, float]:
    """
    The layout of the central peak or peak ring of a non-simple crater (see _PEAK_RING_RATIO).

    Beyond the diameter of Orientale, the ring ratio is held at its value for Orientale and its scatter in m is held at its size
    for Orientale.

    Parameters
    ----------
    diameter : float
        The diameter of the crater (m).
    rng : numpy.random.Generator, optional
        The random number generator for the scatter of the ring radius. If None, the nominal radius is returned.

    Returns
    -------
    ring_fraction : float
        How far the central peak has opened into a ring (0 for a central peak, 1 for a peak ring).
    ring_radius : float
        The radius of the peak ring (m); 0 for a central peak.
    """
    ring_fraction = _log_step(diameter, _PEAK_RING_ONSET)
    if ring_fraction <= 0.0:
        return 0.0, 0.0
    ratio = _PEAK_RING_RATIO * (min(diameter, _RING_DATA_MAX_DIAMETER) / _PEAK_RING_D_REF) ** _PEAK_RING_EXPONENT
    ratio *= _ring_scatter(diameter, _PEAK_RING_SIGMA_LN, rng)
    return ring_fraction, ring_fraction * ratio * diameter / 2


def scarp_layout(diameter: float, morphology_type: str, rng: Generator | None = None) -> tuple[float, float, float]:
    """
    The ring scarp of a crater (see _SCARP): only multiring basins have one.

    Beyond the diameter of Orientale, the scatter of its radius in m is held at its size for Orientale.

    Parameters
    ----------
    diameter : float
        The diameter of the crater (m).
    morphology_type : str
        The morphology type of the crater.
    rng : numpy.random.Generator, optional
        The random number generator for the scatter of the scarp radius. If None, the nominal radius is returned.

    Returns
    -------
    scarp_radius, scarp_height, scarp_width : float
        The radius of the crest of the scarp, its height and the width of its face (m); all 0 for a crater without a scarp.
    """
    if morphology_type != "multiring":
        return 0.0, 0.0, 0.0
    ratio = _SCARP["radius_over_R"] * _ring_scatter(diameter, _SCARP["radius_sigma_ln"], rng)
    height = _SCARP["height_over_D"] * diameter * _log_step(diameter, _SCARP["onset"])
    return ratio * diameter / 2, height, _SCARP["width_over_R"] * diameter / 2


def _basin_ratio(diameter: float, key: str, top_diameter: float, top_value: float) -> float:
    """
    A profile ratio of a crater larger than the range of the non-simple size relations (see _BASIN_ANCHORS).

    Parameters
    ----------
    diameter : float
        The diameter of the crater (m).
    key : str
        The ratio (a key of the _BASIN_ANCHORS entries).
    top_diameter, top_value : float
        The upper end of the range of the size relations, and the nominal ratio of the relation there.

    Returns
    -------
    float
        The ratio, interpolated in ln D between the anchors.
    """
    xs = [math.log(top_diameter), *(math.log(d) for d in _BASIN_ANCHORS)]
    ys = [top_value, *(anchor[key] for anchor in _BASIN_ANCHORS.values())]
    return float(np.interp(math.log(diameter), xs, ys))


# Range of the ejecta profile exponent allowed by volume conservation (the volume budget of fresh craters gives -2.7 to -3.9)
_EJPROFILE_RANGE = (-4.0, -2.5)
# The steepest exponent allowed when no balance is possible within _EJPROFILE_RANGE
_EJPROFILE_LIMIT = -6.0

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
    "scarp_radius",
    "scarp_height",
    "scarp_width",
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
    """The fixed (immutable) parameters of a basicmoon crater: its morphometric parameters, set when the crater is made."""

    rim_height: float | None = None
    """Original rim height of the crater in meters relative to the reference surface."""
    floor_elevation: float | None = None
    """Original floor elevation of the crater in meters relative to the reference surface."""
    floor_radius: float | None = None
    """Original floor radius of the crater in meters."""
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
    scarp_radius: float | None = None
    """Radius of the crest of the ring scarp (the nested ring fault of a multiring basin) in meters. 0 for no scarp."""
    scarp_height: float | None = None
    """Height of the crest of the ring scarp above the rest of the profile in meters."""
    scarp_width: float | None = None
    """Width of the inward-facing face of the ring scarp in meters."""

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
    """The basic lunar crater, with an axisymmetric profile and its morphometric parameters (rim height, floor, peak, etc.)."""

    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=BasicMoonCraterFixed,
        variable_cls=MorphologyCraterVariable,
        **kwargs,
    ):
        super().__init__(
            crater=crater,
            fixed_cls=fixed_cls,
            variable_cls=variable_cls,
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
            f"Ring scarp radius: {format_large_units(self.scarp_radius, quantity='length') if self.scarp_radius else 'None'}\n"
            f"Ring scarp height: {format_large_units(self.scarp_height, quantity='length') if self.scarp_radius else 'None'}\n"
            f"Ring scarp width: {format_large_units(self.scarp_width, quantity='length') if self.scarp_radius else 'None'}\n"
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
        scarp_radius: float | None = None,
        scarp_height: float | None = None,
        scarp_width: float | None = None,
        conserve_volume: bool = True,
        monte_carlo: bool | None = None,
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
            Original floor elevation of the crater in meters relative to the reference surface. If None, it will be computed.
        floor_radius : float, optional
            Original floor radius of the crater in meters. If None, it will be computed.
        wall_curvature : float, optional
            A factor that controls the curvature of the crater wall (0 for straight walls, and 1.0 for very curvy walls). If None, it will be computed based on the morphology type and diameter.
        floor_blend: float, optional
            A factor that controls how sharply the wall and floor blend (0 for sharp floor-to-wall transition, 1.0 for gentle floor-to-wall transition). If None, it will be computed based on the morphology type and diameter.
        rim_width : float, optional
            Original width of the crater rim in meters. If None, it will be computed.
        frac_ejrim : float, optional
            The fraction of the rim height that is ejecta (the ejecta thickness at the rim divided by the rim height). If None, it is 0.3 (or that of Orientale for multiring basins), subject to volume conservation.
        ejprofile : float, optional
            The power law exponent of the ejecta thickness profile. If None, it is -3 (or that of Orientale for multiring basins), subject to volume conservation.
        peak_height : float, optional
            Original central peak height of the crater in meters relative to the reference surface. If None, it will be computed for complex craters and multiring basins and set to 0 for all other types.
        peak_width : float, optional
            Original central peak width (standard deviation) of the crater in meters. If None, it will be computed for complex craters and multiring basins (widening to that of the peak ring of Schrodinger as the peak opens into a ring) and set to 0 for all other types.
        peak_ring_radius : float, optional
            Original radius of the peak ring of the crater in meters. If None, it follows `peak_ring_layout`: 0 (a central peak) below 200 km, opening into a peak ring between 200 and 300 km, at about 0.4 of the rim radius.
        peak_center_distance : float, optional
            Distance of the center of the central peak or peak ring from the crater center (m). The basicmoon profile is
            axisymmetric and ignores it (default 0); realmoon builds its peaks around this offset center, and if it is None draws
            a random offset (see realmoon._PEAK_OFFSET_SIGMA).
        peak_center_bearing : float, optional
            Bearing of the center of the central peak or peak ring from the crater center (degrees, clockwise from north). The
            basicmoon profile ignores it (default 0); realmoon draws a random bearing if it is None.
        scarp_radius : float, optional
            Radius of the crest of the ring scarp (the nested ring fault of a multiring basin) in meters; 0 for no scarp. If None, it follows `scarp_layout`: multiring basins have the Outer Rook ring of Orientale, at 0.67 of the rim radius, and other craters none.
        scarp_height : float, optional
            Height of the crest of the ring scarp above the rest of the profile in meters. If None, it follows `scarp_layout` (0.0025 of the diameter, rising from 0 as multiring basins grow from 500 to 650 km).
        scarp_width : float, optional
            Width of the inward-facing face of the ring scarp in meters; its back slope is three times as wide. If None, it follows `scarp_layout` (0.074 of the rim radius).
        conserve_volume: bool, optional
            If True, then the ejecta parameters (ejprofile and frac_ejrim) will be adjusted in order to attempt to conserve volume between excavation and deposition; the rim height is never changed. Craters with a ring scarp (multiring basins) are exempt. Note that even when True, further adjustments may be needed when emplacing the crater due to the local topography. Default is True.
        monte_carlo : bool or None, optional
            If True, then crater parameters will be drawn from random distributions. If False, then nominal values will be used for all parameters. If None, then the value of the ``monte_carlo`` of the associated Scaling object will be used. Default is None.
        **kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`MorphologyCrater.maker() <cratermaker.components.morphology.MorphologyCrater.maker>`.  Refer to its documentation for a detailed description of valid keyword arguments.

        Returns
        -------
        BasicMoonCrater
            The new crater.

        References
        ----------
        TODO: Update with Minton et al. (2026) info once done.
        """
        input_args = locals()
        source = crater
        from cratermaker.components.morphology import Morphology
        from cratermaker.utils.montecarlo_utils import sample_logfit_heteroskedastic

        if crater is not None and isinstance(crater, BasicMoonCrater):
            for k in ["__class__", "kwargs", "morphology", "crater", "conserve_volume", "cls"]:
                input_args.pop(k, None)
            conserve_volume = conserve_volume and any(list(input_args.values()))
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
            scarp_radius = crater.scarp_radius if scarp_radius is None else scarp_radius
            scarp_height = crater.scarp_height if scarp_height is None else scarp_height
            scarp_width = crater.scarp_width if scarp_width is None else scarp_width

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, monte_carlo=monte_carlo, **kwargs)
        monte_carlo = monte_carlo if monte_carlo is not None else morphology.scaling.monte_carlo
        compute_nominal = not monte_carlo
        rng = morphology.rng

        # The following are derived based on fitting profiles to fresh craters (Minton et al. 2026 in prep). The floor elevation,
        # rim height, rim width, floor radius, wall curvature and floor blend relations are the fits of the ringless baseline
        # profiles of the fresh craters (Lunar_Morphometry_Revisited: scripts/hfield_baseline.py and heteroskedastic_fit.py,
        # data/morphometry_params.json), with the scatter between craters separated from their measurement uncertainty and
        # modeled as normal or lognormal, whichever fits better. Outside the range of diameters of the data (x_min, x_max) they are
        # extrapolated by geometric similarity. Wall curvature and floor blend only smooth the profile. Multiring basins are not
        # part of these fits.

        floor_elevation_params = {
            "simple": {
                "coefficients": [-3.979752385199293, 1.54381759221194, -0.03334186751664677],
                "c": -3.9931109902489275,
                "alpha": -0.26071155480588115,
                "fit_scale": "loglog",
                "scatter_model": "lognormal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
                "x_ref": 8.389456184881865,
            },
            "non-simple": {
                "coefficients": [6.228704419081566, 0.043398512731655174, 0.009068374449593998],
                "c": -3.896394812514182,
                "alpha": -0.5546987550121527,
                "fit_scale": "loglog",
                "scatter_model": "lognormal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
                "x_ref": 10.30162133305592,
            },
        }
        rim_height_params = {
            "simple": {
                "coefficients": [-4.796922969834622, 1.451816414155027, -0.028410426419867563],
                "c": -3.050412634099904,
                "alpha": 0.04133920854978147,
                "fit_scale": "loglog",
                "scatter_model": "lognormal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
                "x_ref": 8.389456184881865,
            },
            "non-simple": {
                "coefficients": [-16.448212792700943, 3.9948614194744763, -0.16833861904099015],
                "c": -0.5932181155690884,
                "alpha": 1.6246168522237205,
                "fit_scale": "loglog",
                "scatter_model": "normal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
            },
        }
        rim_width_params = {
            "simple": {
                "coefficients": [-3.1580381994811417, 1.4800170697740265, -0.03936837662995788],
                "c": -3.1229375645622865,
                "alpha": 1.9553060602140493,
                "fit_scale": "loglog",
                "scatter_model": "normal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
            },
            "non-simple": {
                "coefficients": [19.09819682376756, -3.0444968967088917, 0.18658559498130467],
                "c": -4.019524610605341,
                "alpha": 2.1406629022085246,
                "fit_scale": "loglog",
                "scatter_model": "normal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
            },
        }
        # Smallest value / radius measured in the fresh crater data set, used as the floor of the sampled values. The
        # sampled values have additive noise, so without a floor they can occasionally reach zero or become negative.
        min_radius_fraction = {
            "rim_width": {"simple": 0.1489, "non-simple": 0.0520},  # Cauchy, Orientale
            "rim_height": {"simple": 0.0291, "non-simple": 0.0041},  # D123.9m in RUMKERDOME6, Orientale
            "floor_radius": {"simple": 0.0119, "non-simple": 0.2850},  # D51.2m in NEWCRATER6, Mandel'shtam F
        }
        floor_radius_params = {
            "simple": {
                "coefficients": [-1.977448632884101, 0.7739798351891516, 0.023992241532671363],
                "c": -2.917365314699089,
                "alpha": -0.46252594459213675,
                "fit_scale": "loglog",
                "scatter_model": "lognormal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
                "x_ref": 8.389456184881865,
            },
            "non-simple": {
                "coefficients": [-15.93569252908765, 3.5759345794649127, -0.11275862232298753],
                "c": -4.112698583745451,
                "alpha": -1.1089376740320276,
                "fit_scale": "loglog",
                "scatter_model": "lognormal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
                "x_ref": 10.30162133305592,
            },
        }
        wall_curvature_params = {
            "simple": {
                "coefficients": [0.04276154426076039, 0.08274120819561653, -0.008220292162235484],
                "c": -4.9704008369860535,
                "alpha": 0.0,
                "fit_scale": "logx",
                "scatter_model": "normal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
            },
            "non-simple": {
                "coefficients": [8.487150118544074, -1.6751371267568027, 0.08324044638098378],
                "c": -5.749831239704206,
                "alpha": 0.0,
                "fit_scale": "logx",
                "scatter_model": "normal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
            },
        }
        floor_blend_params = {
            "simple": {
                "coefficients": [1.0970431931049336, -0.12376285402087304, 0.007579980711761793],
                "c": -3.3198599210894835,
                "alpha": 0.0,
                "fit_scale": "logx",
                "scatter_model": "normal",
                "x_min": 114.29288745466218,
                "x_max": 15869.094334407944,
            },
            "non-simple": {
                "coefficients": [-8.446081548018446, 1.791107986657214, -0.08645564637214528],
                "c": -3.7768924080345405,
                "alpha": 0.0,
                "fit_scale": "logx",
                "scatter_model": "normal",
                "x_min": 12060.776126394443,
                "x_max": 136569.08247240804,
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
        args = {}
        diameter_m = min(max(crater.diameter, morphology.valid_diameter_range[0]), morphology.valid_diameter_range[1])
        fcorrection = crater.diameter / diameter_m

        if crater.morphology_type in ["complex", "transitional", "basin", "multiring", "peakring", "ring"]:
            fit_category = "non-simple"
        else:
            fit_category = "simple"

        # Craters larger than the range of the non-simple relations move towards the profiles of Schrodinger and Orientale (see
        # _BASIN_ANCHORS): each drawn length is scaled by the ratio of the anchor profile to the nominal relation at the top of the
        # range (so that the scatter of the relations is kept), and the shape factors are shifted by their difference.
        top_diameter = rim_height_params["non-simple"]["x_max"]
        is_basin = fit_category == "non-simple" and crater.diameter > top_diameter

        def nominal_at_top(params: dict) -> float:
            return float(sample_logfit_heteroskedastic(top_diameter, compute_nominal=True, **params["non-simple"])[0])

        def basin_scale(key: str, params: dict, scale_length: float) -> float:
            top_ratio = nominal_at_top(params) / scale_length
            return _basin_ratio(crater.diameter, key, top_diameter, top_ratio) / top_ratio

        def basin_shift(key: str, params: dict) -> float:
            top_value = nominal_at_top(params)
            return _basin_ratio(crater.diameter, key, top_diameter, top_value) - top_value

        if crater.morphology_type == "multiring":
            frac_ejrim = _MULTIRING_EJECTA["frac_ejrim"] if frac_ejrim is None else frac_ejrim
            ejprofile = _MULTIRING_EJECTA["ejprofile"] if ejprofile is None else ejprofile

        if rim_height is None:
            rim_height = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_height_params[fit_category]
                )[0]
                * fcorrection,
                min_radius_fraction["rim_height"][fit_category] * 0.5 * crater.diameter,
            )
            if is_basin:
                rim_height *= basin_scale("rim_height_over_D", rim_height_params, top_diameter)
        args["rim_height"] = rim_height

        if rim_width is None:
            rim_width = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_width_params[fit_category]
                )[0]
                * fcorrection,
                min_radius_fraction["rim_width"][fit_category] * 0.5 * crater.diameter,
            )
            if is_basin:
                rim_width *= basin_scale("rim_width_over_R", rim_width_params, top_diameter / 2)
        args["rim_width"] = rim_width

        # Try to approximately conserve volume when setting the ejecta thickness at the rim value
        args["ejprofile"] = _EJPROFILE if ejprofile is None else ejprofile

        # This is an initial guess of the ejecta rim. We will adjust it later by integrating the volume of the crater and ejecta profiles
        # The ejecta parameters are not well determined by crater profiles (the topography near the rim mixes structural uplift,
        # ejecta and inherited terrain); the volume budget of fresh craters supports an ejecta thickness at the rim of ~0.3 of the
        # rim height with a -3 power law fall-off (_EJPROFILE), which volume conservation then adjusts
        if frac_ejrim is None:
            frac_ejrim = _FRAC_EJRIM if rim_height > 0.0 else 0.0
        args["frac_ejrim"] = min(max(frac_ejrim, 0.0), 1.0)

        if floor_elevation is None:
            floor_elevation = (
                -sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_elevation_params[fit_category]
                )[0]
                * fcorrection
            )
            floor_elevation = min(floor_elevation, 0.0)
            if is_basin:
                floor_elevation *= basin_scale("floor_depth_over_D", floor_elevation_params, top_diameter)
        args["floor_elevation"] = floor_elevation

        if floor_radius is None:
            floor_radius = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_radius_params[fit_category]
                )[0]
                * fcorrection,
                min_radius_fraction["floor_radius"][fit_category] * 0.5 * crater.diameter,
            )
            if is_basin:
                floor_radius *= basin_scale("floor_radius_over_R", floor_radius_params, top_diameter / 2)
        args["floor_radius"] = min(floor_radius, 0.8 * crater.radius)

        # The peak opens into a peak ring as craters grow (see peak_ring_layout). Complex craters and multiring basins take the
        # height of the fit to central peaks (the peak ring of Orientale is degraded, so its height is only a lower bound), and the
        # width of the peak moves from that of the fit to central peaks to that of the peak ring of Schrodinger as it opens.
        peaked = crater.morphology_type in ("complex", "multiring")
        if peak_height is None:
            if peaked:
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
            if peaked:
                central_width = max(
                    sample_logfit_heteroskedastic(
                        diameter_m, compute_nominal=compute_nominal, rng=rng, **peak_width_params[fit_category]
                    )[0],
                    0.0,
                )
                ring_fraction, _ = peak_ring_layout(crater.diameter)
                peak_width = (1.0 - ring_fraction) * central_width + ring_fraction * _PEAK_RING_WIDTH * crater.radius
            else:
                peak_width = 0.0
        args["peak_width"] = peak_width
        if peak_ring_radius is None:
            peak_ring_radius = 0.0
            if args["peak_height"] > 0.0:
                _, peak_ring_radius = peak_ring_layout(crater.diameter, rng=None if compute_nominal else rng)
        args["peak_ring_radius"] = peak_ring_radius
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
            if is_basin:
                wall_curvature = min(max(wall_curvature + basin_shift("wall_curvature", wall_curvature_params), 0.0), 1.0)

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
            if is_basin:
                floor_blend = min(max(floor_blend + basin_shift("floor_blend", floor_blend_params), 0.0), 1.0)

        args["floor_blend"] = floor_blend

        # The ring scarp of multiring basins (see scarp_layout)
        if scarp_radius is None or scarp_height is None or scarp_width is None:
            layout = scarp_layout(crater.diameter, crater.morphology_type, rng=None if compute_nominal else rng)
            scarp_radius = layout[0] if scarp_radius is None else scarp_radius
            scarp_height = layout[1] if scarp_height is None else scarp_height
            scarp_width = layout[2] if scarp_width is None else scarp_width
        args["scarp_radius"] = scarp_radius
        args["scarp_height"] = scarp_height
        args["scarp_width"] = scarp_width

        kwargs = {**args, **kwargs}

        crater = cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )

        # Adjust the ejecta so that its volume matches the volume excavated by the crater. Only the ejecta parameters are changed:
        # first the exponent of the thickness profile, within the range supported by the volume budget of fresh craters
        # (_EJPROFILE_RANGE); if the exponent needed is outside that range, it is held at the nearest end and the rim thickness
        # (frac_ejrim) is solved for instead. The rim height and the rest of the crater are never changed. The structural uplift
        # of the rim, (1 - frac_ejrim) * rim_height, is part of the crater, so its volume depends on frac_ejrim: the solution is
        # iterated.
        # Craters with a ring scarp (multiring basins) are not held to volume conservation (see _MULTIRING_EJECTA).
        if conserve_volume and not crater.has_scarp:
            params = SimpleNamespace(**{name: getattr(crater, name) for name in _PROFILE_PARAMETERS})
            lo, hi = _EJPROFILE_RANGE

            rim_area = 2 * np.pi * params.rim_height * params.radius**2
            frac0 = params.frac_ejrim

            def candidate(frac_ejrim):
                """
                The exponent that balances the volume for a rim thickness, the remaining imbalance, and whether it excavates.

                The exponent is limited to [_EJPROFILE_LIMIT, hi]; the imbalance is the deposit plus the net volume of the
                crater part, relative to the crater part.
                """
                params.frac_ejrim = frac_ejrim
                excavated = morphology.estimate_volume(params, include_crater=True, include_ejecta=False)
                if excavated >= 0.0:
                    return np.nan, np.inf, False
                ejprofile = -2.0 - frac_ejrim * rim_area / -excavated if frac_ejrim > 0 else lo
                ejprofile = min(max(ejprofile, _EJPROFILE_LIMIT), hi)
                deposit = frac_ejrim * rim_area / (-ejprofile - 2)
                return ejprofile, abs(deposit + excavated) / -excavated, True

            if params.rim_height > 0.0:
                # The rim thickness and the exponent are chosen together. For each rim thickness on a grid, the exponent that
                # balances the volume is found; the choice is, in order of preference: (1) an exact balance with the exponent
                # in _EJPROFILE_RANGE, with the rim thickness closest to its initial value; (2) otherwise, the smallest
                # imbalance among those whose crater part removes net volume (the exponent may then go to _EJPROFILE_LIMIT).
                # Raising frac_ejrim moves volume from the structural uplift (part of the crater) to the ejecta, so a crater
                # whose rim holds more volume than its cavity removes can still be balanced.
                grid = np.unique(np.concatenate([[frac0], np.linspace(0.0, 1.0, 21)]))
                results = [(f, *candidate(f)) for f in grid]
                exact = [r for r in results if r[3] and r[2] < 1e-9 and lo <= r[1] <= hi]
                if exact:
                    best = min(exact, key=lambda r: abs(r[0] - frac0))
                else:
                    valid = [r for r in results if r[3]]
                    best = min(valid, key=lambda r: (r[2], abs(r[0] - frac0))) if valid else (frac0, params.ejprofile, 0, False)
                params.frac_ejrim, params.ejprofile = best[0], best[1]
            kwargs["ejprofile"] = params.ejprofile
            kwargs["frac_ejrim"] = params.frac_ejrim
            crater = cls(crater=crater, morphology=morphology, **kwargs)

        return crater

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
    def scarp_radius(self) -> float | None:
        """Radius of the crest of the ring scarp (the nested ring fault of a multiring basin) in meters. 0 for no scarp."""
        return self._fixed.scarp_radius

    @property
    def scarp_height(self) -> float | None:
        """Height of the crest of the ring scarp above the rest of the profile in meters."""
        return self._fixed.scarp_height

    @property
    def scarp_width(self) -> float | None:
        """Width of the inward-facing face of the ring scarp in meters."""
        return self._fixed.scarp_width

    @property
    def has_scarp(self) -> bool:
        """Whether the crater has a ring scarp (a multiring basin)."""
        return bool(self.scarp_radius) and self.scarp_radius > 0.0

    @property
    def frac_ejrim(self) -> float | None:
        """Ejecta rim thickness of the crater in meters."""
        return self._var._frac_ejrim


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
            str_repr += f"Ejecta Truncation: {self.ejecta_truncation:.2f} * crater.radius\n"
        else:
            str_repr += "Ejecta Truncation: Off\n"
        str_repr += f"Ejecta Rays: {self.dorays}\n"
        return str_repr

    # The compiled function that forms a batch of craters of this model's crater type. Models with a different crater type
    # override this.
    _form_craters_binding = staticmethod(basicmoon_bindings.form_craters)

    # The critical angle used for slope collapse after a crater forms, in degrees. The interiors of fresh lunar craters measured
    # with the same mesh gradient (Lunar_Morphometry_Revisited, scripts/interior_slopes.py) sit at the angle of repose, about 35
    # degrees on their walls, but rarely exceed 40 degrees (0.1% of the interior of simple craters at the SLDEM spacing), so
    # collapse removes only slopes steeper than real craters have.
    _CRITICAL_SLOPE_ANGLE = 40.0

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
                self.counting.tally(**kwargs)

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
