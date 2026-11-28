from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray

from cratermaker.bindings import basicmoon_bindings, realmoon_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import (
    BasicMoonCrater,
    BasicMoonCraterFixed,
    BasicMoonMorphology,
    peak_ring_layout,
)

# The H-field model of fresh lunar craters (Lunar_Morphometry_Revisited, data/hfield/model.nc), distributed with the package
# in compact form (the modes packed as 16-bit integers, accurate to ~2e-5). Another model file can be given with the
# hfield_model_file argument.
_HFIELD_MODEL_FILE = Path(__file__).resolve().parent / "hfield_model.nc"

# The polar grid of the model: 720 bearings and warped radii rho = r / R_rim from 0 to 2 in steps of 0.005
_N_THETA = 720
_RHO_MAX = 2.0
_DRHO = 0.005
_KMAX = _N_THETA // 2 - 1
_KMAX_RIM = 8  # the rim outline has azimuthal wavenumbers k <= 8
_N_SCARP = 72  # bearings of the variation of the ring scarp of multiring basins (see RealMoonMorphology.scarp_outline_data)
# The rim outline of craters larger than the largest complex crater of the model (Langrenus, its node at 137 km) does not keep
# growing in proportion to the radius: the shape of a rim outline is set by the collapse of the rim and wall through local
# faulting and fracturing, whose scale is that of the rock rather than of the crater. Its rms wobble (in m) is held at that of
# the model at its largest complex node (about 2.4 km). Schrodinger's rim, traced within 10% of its catalog radius, has lobes of
# at most 3.5-5.6 km rms (an upper bound: the trace still strays onto the terrain and craters around it;
# Lunar_Morphometry_Revisited, rim_spectrum.py).
# The ring scarps of multiring basins are drawn as rims, with the same wobble.
# _FIELD_SCALE_NOTE: beyond the largest complex crater of the model, the field (its terraces, slumps and texture) keeps that
# crater's physical scale around the rim instead of being stretched with the crater: at wavenumber k, a crater of diameter D
# takes the statistics of the model at k' = k d0 / D (d0 the diameter of the largest complex crater; at least those of k' = 1),
# with the variance at each k multiplied by d0 / D, so that the variance per unit length of the rim is that of the model; and
# radially, inside the rim, it is compressed about the rim, rho_model = 1 + (rho - 1) D / d0, so that the features of the wall
# keep their size and shape (at Schrodinger's size the field reaches in to 0.72 R; the floor inside has no field). Outside the
# rim it is not compressed (rho_model = rho): the outer flank of the rim and the ejecta blanket scale with the crater, and the
# rim's texture grades into them (compressed, the field ended abruptly just outside the rim, leaving the blanket smooth).
# Stretched in bearing too, the field of Langrenus made the rims of basins lobate at 30-50 km scales, while the rim of
# Schrodinger is nearly circular with fine scalloping.
# Morphology types without a node of their own use those of complex craters. Multiring basins do too: the model's only basin
# (Orientale) is heavily degraded (freshness class D), and its field (about three times rougher in units of its relief than the
# largest complex crater's on the floor, and 30% rougher on the wall and rim) is mostly the cratering and modification it has
# accumulated, which leaves its rim indistinct. Basins therefore take the field of the largest complex crater of the model
# (Langrenus), scaled to their relief, and keep the ringed profile of Orientale (see basicmoon._MULTIRING).
_TYPE_ALIASES = {"basin": "complex", "peakring": "complex", "ring": "complex", "multiring": "complex"}
# The strength of the lobes of rim outlines varies from crater to crater more than the Gaussian draws of the model give: most
# fresh complex craters are nearly round, a few are strongly lobate. The drawn outline (k <= _KMAX_RIM) is scaled by a lognormal
# factor exp(sigma N - sigma^2), whose mean square is one, with sigma the scatter of the rms rim anomaly of fresh complex craters
# (classes A and B) about its trend with diameter beyond that of the Gaussian draws (Lunar_Morphometry_Revisited,
# rim_lobe_scatter.py).
_RIM_LOBE_SIGMA_LN = 0.32
# The wall of complex craters and basins: stacked rim-like levels (the crest and the terrace scarps below it; realmoon.rs,
# WallLevels; Lunar_Morphometry_Revisited, notes section 19 and scallop_levels.py). The walls of fresh complex craters are made
# of many overlapping scarp arcs at different radii; a single rim outline with lobes, or the H-field's broad structure, moves
# the whole wall smoothly into large lobes, while the real rims are nearly round with small irregular lobes. Each level is a
# continuous, nearly round closed curve with lobes of its own (ln radius a random function of bearing with rms _WALL_WOBBLE and a
# spectrum falling as k^-_WALL_WOBBLE_SLOPE up to _WALL_WOBBLE_KMAX), a scarp whose strength varies along it (a random function
# of rms _WALL_STRENGTH_SD about _WALL_STRENGTH_MEAN up to _WALL_STRENGTH_KMAX, clipped to 0-1) and an elevation offset that
# varies along it as turbulence (octaves of |band-limited noise| from _WALL_OFFSET_K0, rms _WALL_OFFSET_SD of the mean drop
# between levels). The levels go from the crest (1) down to the floor edge (for a multiring basin, to the outer edge of its
# ring scarp), spaced by a lognormal amount (median _WALL_SPACING of the rim radius, ln scatter _WALL_SPACING_SIGMA). They wobble
# into each other, and where they collide their scarps merge, which makes them read as overlapping arcs. A turbulence field of
# rms _WALL_TURBULENCE of the mean drop between levels, with largest wavelength _WALL_TURBULENCE_WAVELENGTH, is added over the
# wall. The per-crater lobe strength scatter (_RIM_LOBE_SIGMA_LN) scales the lobes of all levels. Beyond the largest complex
# crater of the model, the spacing, the lobes and the wavelengths keep that crater's physical scale (see _FIELD_SCALE_NOTE). The
# H-field gives these craters only its texture, above wavenumber _WALL_FIELD_KMIN (at the physical scale of the largest complex
# crater beyond it). The parameters were chosen by comparing renders with Langrenus and Copernicus, and the lobes calibrated so
# that the rim outlines of synthetic craters, traced as the rims of the model are (with the fit extended to k <= 24), have
# close to the rms of those of the fresh complex craters in each band of k (an amplitude of 0.040 matches it within 4-10%; it
# was reduced to 0.035 by eye, about 3-15% below; Lunar_Morphometry_Revisited, rim_shape_spectrum.py).
_WALL_TYPES = ("complex", "multiring")
_WALL_SPACING = 0.08
_WALL_SPACING_SIGMA = 0.3
_WALL_WOBBLE = 0.035
_WALL_WOBBLE_SLOPE = 1.2
_WALL_WOBBLE_KMAX = 60
_WALL_STRENGTH_MEAN = 0.65
_WALL_STRENGTH_SD = 0.45
_WALL_STRENGTH_KMAX = 18
_WALL_OFFSET_SD = 0.6
_WALL_OFFSET_K0 = 2
_WALL_TURBULENCE_OCTAVES = 6
_WALL_TURBULENCE_PERSISTENCE = 0.5
_WALL_TURBULENCE = 0.15
_WALL_TURBULENCE_WAVELENGTH = 0.25
_WALL_FIELD_KMIN = 30
_NO_WALL = (np.zeros(0), np.zeros((3, 0, 0)), np.zeros(0))

# Central peaks (built on the surface by the central peak model of the Rust crate, see central_peak.rs). The height of the
# relaxed peak above the floor, as a fraction of the crater's relief (rim height + floor depth), follows the peaks of the fresh
# craters of 36-139 km measured at the native SLDEM2015 spacing: 0.259 (D / 60 km)^0.89 with lognormal scatter 0.44. Above
# the fitted range it is held at its value at 139 km (it would otherwise exceed the relief by 300 km); below it, it is
# extrapolated, giving the smaller peaks of transitional craters. The envelope of a central peak is a Gaussian of standard
# deviation _PEAK_WIDTH of the rim radius. Central peaks open into peak rings as craters grow from 200 to 300 km (see
# basicmoon.peak_ring_layout), whose envelope is a Gaussian ring of standard deviation _PEAK_RING_WIDTH of the rim radius. The
# widths are calibrated with the rest of the peak model (central_peak.rs, PeakModel) against the central peaks of fresh craters
# and the peak rings of Schrodinger and Orientale (Lunar_Morphometry_Revisited, peak_model_calibrate.py). Multiring basins have a
# peak ring (their innermost ring) with the height of the largest complex craters; the peak ring of Orientale is degraded, so its
# height is a lower bound.
_PEAKED = ("transitional", "complex", "multiring")
_PEAK_HEIGHT_RELIEF = 0.259
_PEAK_HEIGHT_D_REF = 60e3
_PEAK_HEIGHT_EXPONENT = 0.89
_PEAK_HEIGHT_SIGMA_LN = 0.44
_PEAK_HEIGHT_D_MAX = 139e3
_PEAK_WIDTH = 0.108
_PEAK_RING_WIDTH = 0.134


def _interp_rows(values: NDArray[np.float64], x: NDArray[np.float64], x_new: NDArray[np.float64]) -> NDArray[np.float64]:
    """Each row of `values` (given at the increasing points `x`) interpolated linearly at `x_new`, zero outside `x`."""
    j = np.clip(np.searchsorted(x, x_new) - 1, 0, len(x) - 2)
    t = (x_new - x[j]) / (x[j + 1] - x[j])
    out = (1.0 - t)[None, :] * values[:, j] + t[None, :] * values[:, j + 1]
    out[:, (x_new < x[0]) | (x_new > x[-1])] = 0.0
    return out


def _peak_height(diameter: float, relief: float, rng: np.random.Generator | None) -> float:
    """
    Draw the height of a central peak.

    Parameters
    ----------
    diameter : float
        The diameter of the crater (m).
    relief : float
        The relief of the crater (rim height + floor depth, m).
    rng : numpy.random.Generator or None
        The random number generator. If None, the nominal height is returned.

    Returns
    -------
    float
        The height of the relaxed peak above the floor (m).
    """
    d = min(diameter, _PEAK_HEIGHT_D_MAX)
    height = _PEAK_HEIGHT_RELIEF * (d / _PEAK_HEIGHT_D_REF) ** _PEAK_HEIGHT_EXPONENT * relief
    if rng is not None:
        height *= float(np.exp(_PEAK_HEIGHT_SIGMA_LN * rng.standard_normal()))
    return max(height, 0.0)


def _copy_hfield(source: Crater | None, target: Crater) -> None:
    """
    Give the target the field of the source crater, if it has one and the target is the same crater.

    A copy with a different size or type is a different crater, which draws its own field.
    """
    if source is None or source.radius != target.radius or source.morphology_type != target.morphology_type:
        return
    for name in ("_hfield", "_wall", "_scarp_outline"):
        data = getattr(getattr(source, "_var", None), name, None)
        if data is not None:
            object.__setattr__(target._var, name, data)


@Crater.register("realmooncrater")
class RealMoonCrater(BasicMoonCrater):
    """
    A crater of the realistic lunar morphology.

    It is a basic lunar crater plus a rim outline and a field of variation around it (the H-field). The rim outline and field are drawn from the morphology's random number generator on first access (see
    `RealMoonMorphology.hfield_data`) and kept with the crater, including when it is copied without changing its size or type.
    """

    def __init__(
        self, crater: Crater | None = None, fixed_cls=BasicMoonCraterFixed, variable_cls=MorphologyCraterVariable, **kwargs
    ):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        _copy_hfield(crater, self)

    @classmethod
    def maker(cls, crater: Crater | None = None, morphology: Morphology | None = None, **kwargs: Any) -> RealMoonCrater:
        """
        Create a crater (see `BasicMoonCrater.maker`). A copy of a crater that keeps its size and type keeps its field.

        Unlike the analytic peak of basicmoon (complex craters and multiring basins only), a new realmoon crater also gets a
        central peak when it is transitional. This is intended.

        Parameters
        ----------
        crater : Crater, optional
            The crater to copy, if any.
        morphology : Morphology, optional
            The morphology of the crater.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        RealMoonCrater
        """
        new = super().maker(crater=crater, morphology=morphology, **kwargs)
        # A new crater (not a copy of a realmoon crater) gets a central peak from the central peak model, unless one is given
        if kwargs.get("peak_height") is None and not isinstance(crater, RealMoonCrater):
            peaked = new.morphology_type in _PEAKED
            relief = new.rim_height - new.floor_elevation
            rng = getattr(morphology if morphology is not None else new.morphology, "rng", None)
            # The radius of the peak ring is that of the basic crater (see basicmoon.peak_ring_layout); the envelope widens from
            # that of a central peak to that of a peak ring as the peak opens
            ring_fraction, _ = peak_ring_layout(new.diameter)
            width = ((1.0 - ring_fraction) * _PEAK_WIDTH + ring_fraction * _PEAK_RING_WIDTH) * new.radius
            new = super().maker(
                crater=new,
                morphology=morphology,
                peak_height=_peak_height(new.diameter, relief, rng) if peaked else 0.0,
                peak_width=width if peaked else 0.0,
                peak_ring_radius=new.peak_ring_radius if peaked else 0.0,
            )
        _copy_hfield(crater, new)
        return new

    @property
    def hfield(self) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
        """The rim outline R_rim(theta) / R, the field of variation (m, n_theta x n_rho) and its largest warped radius."""
        data = getattr(self._var, "_hfield", None)
        if data is None:
            data = self.morphology.hfield_data(self)
            object.__setattr__(self._var, "_hfield", data)
        return data

    @property
    def wall(self) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """The stacked levels of the wall of a complex crater or basin, drawn with the field (see
        `RealMoonMorphology.hfield_data`): their nominal radii, their radius, scarp strength and elevation offset (3 x n_levels x
        n_theta), and the wall turbulence (standard deviation in m, wavelength in rim radii, seed); all empty for other craters."""
        _ = self.hfield
        return getattr(self._var, "_wall", _NO_WALL)

    @property
    def scarp_outline(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """The variation of the ring scarp around a multiring basin: its radius and height at _N_SCARP bearings, as factors of
        `scarp_radius` and `scarp_height` (empty for a crater without a scarp; see `RealMoonMorphology.scarp_outline_data`)."""
        data = getattr(self._var, "_scarp_outline", None)
        if data is None:
            data = self.morphology.scarp_outline_data(self)
            object.__setattr__(self._var, "_scarp_outline", data)
        return data


@Morphology.register("realmoon")
class RealMoonMorphology(BasicMoonMorphology):
    """
    The realistic lunar crater morphology.

    A crater is the basic lunar crater (its profile and its size relations, with the ejecta and the degradation of the terrain
    under it), warped by an irregular rim outline, plus a random field of variation (terraces, slumps, the texture of the
    floor and walls, and the texture of the ejecta). The outline and field are drawn from the H-field model of fresh lunar
    craters: for each morphology type and size, the covariance between warped radii of the azimuthal Fourier components of
    the variation, decomposed into radial modes (empirical orthogonal functions). Central peaks and peak rings are not part of
    the field: transitional and complex craters get a central peak, and multiring basins a peak ring, built on the surface as a patchy population of columns
    (where noise exceeds a threshold, raised to a Gaussian envelope of standard deviation `peak_width`, offset to radius
    `peak_ring_radius` for a peak ring), with flanks at the critical slope, to the height `peak_height`. The model is described in the
    morphometry of fresh lunar craters by Minton et al. (2026, in prep.).

    Parameters
    ----------
    hfield_model_file : str or Path, optional
        The H-field model file. Default: the model packaged with cratermaker (hfield_model.nc).
    hfield_pix : float, optional
        The resolution (m) the field is drawn to. Default: the size of the smallest face the field of each crater is applied to
        (in its ejecta region), or of the face at its center if it has no region.
    **kwargs : Any
        Arguments of `BasicMoonMorphology`.
    """

    def __init__(self, hfield_model_file: str | Path | None = None, hfield_pix: float | None = None, **kwargs: Any):
        if hfield_model_file is None:
            hfield_model_file = _HFIELD_MODEL_FILE
        hfield_model_file = Path(hfield_model_file)
        if not hfield_model_file.exists():
            raise FileNotFoundError(
                f"H-field model file {hfield_model_file} not found"
            )
        with xr.open_dataset(hfield_model_file) as ds:
            object.__setattr__(
                self,
                "_hfield_model",
                {
                    "type": ds["morphology_type"].values.astype(str),
                    "diameter": ds["diameter_km"].values.astype(np.float64) * 1e3,
                    "k_node": ds["k_node"].values.astype(np.int64),
                    "eigenvalue": ds["eigenvalue"].values.astype(np.float64),
                    "eigenvector": ds["eigenvector"].values.astype(np.float32),
                    "variance": ds["variance"].values.astype(np.float64),
                },
            )
        object.__setattr__(self, "_hfield_model_file", hfield_model_file)
        model = self._hfield_model
        largest = int(np.argmax(np.where(model["type"] == "complex", model["diameter"], -np.inf)))
        object.__setattr__(self, "_rim_wobble_diameter", float(model["diameter"][largest]))
        object.__setattr__(self, "_rim_wobble_ln_rms", self._rim_ln_rms(largest))
        object.__setattr__(self, "_hfield_pix", hfield_pix)
        super().__init__(**kwargs)

    @property
    def _CraterType(self) -> type[RealMoonCrater]:
        """The class definition of the associated Crater type."""
        return RealMoonCrater

    @property
    def hfield_model_file(self) -> Path:
        """The H-field model file."""
        return self._hfield_model_file

    def _nodes(self, morphology_type: str, diameter: float) -> list[tuple[int, float]]:
        """The model nodes of the type that bracket the diameter, with their weights (linear in ln D)."""
        m = self._hfield_model
        mtype = _TYPE_ALIASES.get(morphology_type, morphology_type)
        nodes = np.where(m["type"] == mtype)[0]
        if len(nodes) == 0:
            nodes = np.where(m["type"] == "complex")[0]
        D = m["diameter"][nodes]
        o = np.argsort(D)
        nodes, D = nodes[o], D[o]
        x = np.log(diameter)
        if len(nodes) == 1 or x <= np.log(D[0]):
            return [(int(nodes[0]), 1.0)]
        if x >= np.log(D[-1]):
            return [(int(nodes[-1]), 1.0)]
        j = int(np.searchsorted(np.log(D), x)) - 1
        a = (x - np.log(D[j])) / (np.log(D[j + 1]) - np.log(D[j]))
        return [(int(nodes[j]), 1.0 - a), (int(nodes[j + 1]), a)]

    def hfield_data(self, crater: RealMoonCrater) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
        """
        Draw the rim outline and field of variation of a crater.

        The Fourier coefficients of [ln R_rim, D(rho)] at each azimuthal wavenumber k are drawn from the radial modes of the
        model nodes of the crater's type that bracket its diameter (the modes of the nearest wavenumber node, rescaled to the
        variance of the node at k; the draws of the two nodes are combined with weights sqrt(1 - a) and sqrt(a), so that the
        covariance is interpolated linearly in ln D). Only the wavenumbers and radii the surface resolves are drawn: k up to
        pi rho_max R / pix and steps in rho of at least pix / (2 R), with pix `hfield_pix` or the size of the smallest face the
        field is applied to (that of the crater's ejecta region, or of the face at its center if it has no region).
        The field is in units of the crater's relief (rim height + floor depth), converted to meters. The radial grid always
        ends at rho_max = 2, the extent of the model.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater.

        Returns
        -------
        rim : NDArray[np.float64]
            R_rim(theta) / R at the bearings of the field.
        field : NDArray[np.float64]
            The field (m), n_theta x n_rho, on bearings (i + 1/2) 360 / n_theta degrees and warped radii from 0 to rho_max.
        rho_max : float
            The largest warped radius of the field.
        """
        m = self._hfield_model
        R = crater.radius
        pix = self._hfield_pix
        if pix is None:
            # The finest faces the field is applied to (those of the crater's ejecta region, which extends to the edge of the
            # field), so that a crater centered on coarse faces is still drawn to the resolution of the fine faces it overlaps
            try:
                region = crater.ejecta_region
                if region is not None and region.n_face > 0:
                    pix = float(np.sqrt(np.min(region.face_area)))
                else:
                    pix = float(np.sqrt(self.surface.face_area[crater.face_index]))
            except Exception:  # noqa: BLE001
                pix = 0.0
        kmax = _KMAX if pix <= 0 else int(min(_KMAX, max(_KMAX_RIM, np.ceil(np.pi * _RHO_MAX * R / pix))))
        n_theta = int(min(_N_THETA, max(16, 2 * (kmax + 1))))
        kmax = min(kmax, n_theta // 2 - 1)
        # The radial grid is thinned by a divisor of its number of intervals, so that it always ends at _RHO_MAX
        n_intervals = int(round(_RHO_MAX / _DRHO))
        wanted = 1 if pix <= 0 else max(1, int(np.floor(pix / (2 * R) / _DRHO)))
        step = max(d for d in range(1, n_intervals + 1) if n_intervals % d == 0 and d <= wanted)
        rows = np.arange(0, n_intervals + 1, step)
        rho_max = _RHO_MAX
        elements = np.concatenate(([0], 1 + rows))

        f = min(1.0, self._rim_wobble_diameter / crater.diameter)
        if f < 1.0:
            # Inside the rim, the field keeps the physical scale of the largest complex crater of the model radially too (see
            # _FIELD_SCALE_NOTE): it is drawn at all radii of the model and compressed about the rim
            rim, X = self._draw_hfield(crater.morphology_type, crater.diameter, kmax, n_theta, np.arange(1 + n_intervals + 1))
            x = rows * _DRHO - 1.0
            rho_model = 1.0 + np.where(x <= 0.0, x / f, x)
            X = np.column_stack([X[:, :1], _interp_rows(X[:, 1:], np.arange(n_intervals + 1) * _DRHO, rho_model)])
        else:
            rim, X = self._draw_hfield(crater.morphology_type, crater.diameter, kmax, n_theta, elements)
        relief = max(crater.rim_height - crater.floor_elevation, 0.0)
        wall = _NO_WALL
        if crater.morphology_type in _WALL_TYPES:
            # Only the texture of the field; the wall levels replace its broad structure and the rim outline
            F = np.fft.rfft(X[:, 1:], axis=0)
            F[: int(np.ceil(_WALL_FIELD_KMIN / f)) + 1] = 0.0
            X = np.column_stack([X[:, :1], np.fft.irfft(F, n=n_theta, axis=0)])
            drawn = self._draw_wall(crater, n_theta, f)
            if drawn is not None:
                wall = drawn
                rim = np.ascontiguousarray(wall[1][0].max(axis=0), dtype=np.float64)
        object.__setattr__(crater._var, "_wall", wall)
        field = np.ascontiguousarray(X[:, 1:] * relief, dtype=np.float64)
        return rim, field, float(rho_max)

    def _periodic_noise(self, n: int, rms: float, kmax: float, slope: float = 0.0) -> NDArray[np.float64]:
        """A random periodic function at n bearings with power ~ k^-slope for 1 <= k <= kmax, with the given rms."""
        k = np.arange(n // 2 + 1)
        amp = np.where((k >= 1) & (k <= min(kmax, n // 2 - 1)), np.maximum(k, 1) ** (-slope / 2), 0.0)
        F = amp * (self.rng.standard_normal(len(k)) + 1j * self.rng.standard_normal(len(k)))
        values = np.fft.irfft(F, n=n)
        sd = values.std()
        return values * (rms / sd) if sd > 0 else values

    def _periodic_turbulence(self, n: int, rms: float, k0: float) -> NDArray[np.float64]:
        """Periodic turbulence at n bearings: octaves o of |band-limited noise| over k in [k0 2^o, k0 2^(o+1)), weighted
        _WALL_TURBULENCE_PERSISTENCE^o, with its mean removed and the given rms."""
        k = np.arange(n // 2 + 1)
        values = np.zeros(n)
        for o in range(_WALL_TURBULENCE_OCTAVES):
            band = (k >= k0 * 2**o) & (k < k0 * 2 ** (o + 1)) & (k <= n // 2 - 1)
            noise = self.rng.standard_normal(len(k)) + 1j * self.rng.standard_normal(len(k))
            if not band.any():
                continue
            octave = np.fft.irfft(band * noise, n=n)
            values += _WALL_TURBULENCE_PERSISTENCE**o * np.abs(octave / octave.std())
        values -= values.mean()
        sd = values.std()
        return values * (rms / sd) if sd > 0 else values

    def _draw_wall(
        self, crater: RealMoonCrater, n_theta: int, f: float
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]] | None:
        """
        Draw the stacked levels of the wall of a complex crater or basin (see _WALL_TYPES).

        Parameters
        ----------
        crater : RealMoonCrater
            The crater.
        n_theta : int
            The number of bearings of its field.
        f : float
            The ratio of the diameter of the largest complex crater of the model to the crater's (at most 1).

        Returns
        -------
        tuple or None
            The nominal radii of the levels, their radius, scarp strength and elevation offset (3 x n_levels x n_theta), and the
            wall turbulence (see `RealMoonCrater.wall`); None if the wall holds fewer than two levels.
        """
        R = crater.radius
        base = crater.floor_radius / R
        if crater.has_scarp:
            base = max(base, (crater.scarp_radius + crater.scarp_width) / R)
        levels = [1.0]
        while True:
            nxt = levels[-1] - _WALL_SPACING * f * float(np.exp(_WALL_SPACING_SIGMA * self.rng.standard_normal()))
            if nxt < base:
                break
            levels.append(nxt)
        if len(levels) < 2:
            return None
        nominal = np.array(levels)
        knots = np.asarray(
            basicmoon_bindings.basicmoon_profile(
                radial_distances=nominal * R, bearings=np.zeros(len(nominal)), reference_elevations=np.zeros(len(nominal)),
                crater=crater, include_crater=True, include_ejecta=False,
            ),
            dtype=np.float64,
        )
        drop = max((knots[0] - knots[-1]) / (len(levels) - 1), 0.0)
        lobes = _WALL_WOBBLE * f * float(np.exp(_RIM_LOBE_SIGMA_LN * self.rng.standard_normal() - _RIM_LOBE_SIGMA_LN**2))
        data = np.empty((3, len(levels), n_theta))
        for j, rho in enumerate(levels):
            data[0, j] = rho * np.exp(self._periodic_noise(n_theta, lobes, _WALL_WOBBLE_KMAX / f, _WALL_WOBBLE_SLOPE))
            data[1, j] = np.clip(_WALL_STRENGTH_MEAN + self._periodic_noise(n_theta, _WALL_STRENGTH_SD, _WALL_STRENGTH_KMAX / f), 0, 1)
            data[2, j] = self._periodic_turbulence(n_theta, _WALL_OFFSET_SD * drop, _WALL_OFFSET_K0 / f)
        params = np.array([_WALL_TURBULENCE * drop, _WALL_TURBULENCE_WAVELENGTH * f, float(self.rng.integers(2**31))])
        return nominal, np.ascontiguousarray(data), params

    def _rim_ln_rms(self, node: int) -> float:
        """The rms of the rim anomaly ln(R_rim / R) of the model node (wavenumbers 1 to _KMAX_RIM)."""
        m = self._hfield_model
        ks = np.arange(1, _KMAX_RIM + 1)
        nearest = np.argmin(np.abs(np.log(m["k_node"])[None, :] - np.log(ks)[:, None]), axis=1)
        total = 0.0
        for k, j in zip(ks, nearest):
            lam = m["eigenvalue"][node][j]
            keep = lam > 0
            if keep.any() and m["variance"][node][k] > 0:
                share = lam[keep] / lam[keep].sum()
                total += 2.0 * m["variance"][node][k] * float(np.sum(share * m["eigenvector"][node][j][0, keep].astype(float) ** 2))
        return float(np.sqrt(total))

    def _rim_wobble_scale(self, diameter: float) -> float:
        """The factor applied to the drawn rim anomaly ln(R_rim / R) of a crater of the given diameter: 1 up to the largest
        complex crater of the model, beyond which the rms wobble (in m) is held at that crater's."""
        return min(1.0, self._rim_wobble_diameter / diameter)

    def _draw_hfield(
        self, morphology_type: str, diameter: float, kmax: int, n_theta: int, elements: NDArray[np.int64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Draw a rim outline and field from the model (see `hfield_data`).

        Parameters
        ----------
        morphology_type : str
            The morphology type whose model nodes are used.
        diameter : float
            The diameter the nodes are interpolated to (m).
        kmax : int
            The largest azimuthal wavenumber drawn.
        n_theta : int
            The number of bearings.
        elements : NDArray[np.int64]
            The elements of the model drawn: 0 is ln R_rim, 1 + i is the field at rho = i _DRHO.

        Returns
        -------
        rim : NDArray[np.float64]
            R_rim(theta) / R (k <= _KMAX_RIM, with the wobble held in m beyond the largest complex crater of the model).
        X : NDArray[np.float64]
            The drawn elements at each bearing (n_theta x len(elements)); the field in units of the relief.
        """
        F = np.zeros((n_theta // 2 + 1, len(elements)), dtype=np.complex128)
        ks = np.arange(1, kmax + 1)
        # Beyond the largest complex crater of the model, the field keeps that crater's physical scale around the rim (see
        # _FIELD_SCALE_NOTE); the rim outline (element 0) is drawn at its own wavenumbers
        f = min(1.0, self._rim_wobble_diameter / diameter)
        if f < 1.0:
            self._accumulate_modes(F, morphology_type, diameter, ks, ks.astype(np.float64), 1.0, elements[:1], [0])
            self._accumulate_modes(F, morphology_type, diameter, ks, f * ks, f, elements[1:], list(range(1, len(elements))))
        else:
            self._accumulate_modes(F, morphology_type, diameter, ks, ks.astype(np.float64), 1.0, elements, list(range(len(elements))))
        X = np.fft.irfft(F * n_theta, n=n_theta, axis=0)
        Fr = np.fft.rfft(X[:, 0])
        Fr[_KMAX_RIM + 1 :] = 0.0
        Fr *= self._rim_wobble_scale(diameter) * np.exp(_RIM_LOBE_SIGMA_LN * self.rng.standard_normal() - _RIM_LOBE_SIGMA_LN**2)
        rim = np.exp(np.fft.irfft(Fr, n=n_theta))
        return np.ascontiguousarray(rim, dtype=np.float64), X

    def _accumulate_modes(
        self,
        F: NDArray[np.complex128],
        morphology_type: str,
        diameter: float,
        ks: NDArray[np.int64],
        k_model: NDArray[np.float64],
        variance_factor: float,
        elements: NDArray[np.int64],
        columns: list[int],
    ) -> None:
        """
        Add random draws of the model's modes to the Fourier coefficients of a field (see `_draw_hfield`).

        Parameters
        ----------
        F : NDArray[np.complex128]
            The Fourier coefficients (k x columns), added to in place.
        morphology_type : str
            The morphology type whose model nodes are used.
        diameter : float
            The diameter the nodes are interpolated to (m).
        ks : NDArray[np.int64]
            The wavenumbers drawn.
        k_model : NDArray[np.float64]
            The wavenumber of the model whose modes and variance (interpolated in k, and at least that of k = 1) are used at each
            of `ks`.
        variance_factor : float
            The factor applied to the variance.
        elements : NDArray[np.int64]
            The elements of the model drawn.
        columns : list[int]
            The columns of F they are added to.
        """
        m = self._hfield_model
        k_model = np.maximum(k_model, 1.0)
        nearest = np.argmin(np.abs(np.log(m["k_node"])[None, :] - np.log(k_model)[:, None]), axis=1)
        for node, weight in self._nodes(morphology_type, diameter):
            lam_all, vec_all, var_k = m["eigenvalue"][node], m["eigenvector"][node], m["variance"][node]
            variance = variance_factor * np.interp(k_model, np.arange(len(var_k)), var_k)
            for j in np.unique(nearest):
                lam = lam_all[j]
                keep = np.where(lam > 0)[0]
                if len(keep) == 0:
                    continue
                sel = (nearest == j) & (variance > 0)
                if not sel.any():
                    continue
                kj = ks[sel]
                # Modes of the node, restricted to the drawn radii, scaled to the variance at each k
                V = vec_all[j][elements][:, keep].astype(np.float64)
                scale = np.sqrt(lam[keep] / lam[keep].sum())[None, :] * np.sqrt(variance[sel])[:, None]
                z = (
                    self.rng.standard_normal((len(kj), len(keep))) + 1j * self.rng.standard_normal((len(kj), len(keep)))
                ) / np.sqrt(2.0)
                F[np.ix_(kj, columns)] += np.sqrt(weight) * ((scale * z) @ V.T)

    def scarp_outline_data(self, crater: RealMoonCrater) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Draw the variation of the ring scarp of a multiring basin around the basin.

        The scarp is the rim of the peak-ring basin inside the multiring basin, so it varies as the rim of a complex crater of
        its diameter (2 `scarp_radius`): its radius as that crater's rim outline (with the wobble held in m beyond the largest
        complex crater of the model, as for rims), and its height as that crater's rim height, the field of the model at the
        rim crest (in units of the relief) over the basin's rim fraction of the relief, (rim height) / (rim height + floor
        depth). The height is never negative.

        Parameters
        ----------
        crater : RealMoonCrater
            The basin.

        Returns
        -------
        radius_factor, height_factor : NDArray[np.float64]
            The scarp radius and height at _N_SCARP bearings (i + 1/2) 360 / _N_SCARP degrees, as factors of `scarp_radius`
            and `scarp_height`; both empty for a crater without a scarp.
        """
        if not crater.has_scarp:
            return np.zeros(0), np.zeros(0)
        i_crest = int(round(1.0 / _DRHO))
        rim, X = self._draw_hfield("complex", 2.0 * crater.scarp_radius, _KMAX_RIM, _N_SCARP, np.array([0, 1 + i_crest]))
        relief = crater.rim_height - crater.floor_elevation
        rim_fraction = crater.rim_height / relief if relief > 0 else 1.0
        height = np.clip(1.0 + X[:, 1] / max(rim_fraction, 1e-6), 0.0, None)
        return np.ascontiguousarray(rim, dtype=np.float64), np.ascontiguousarray(height, dtype=np.float64)

    def _prepare_crater(self, crater: RealMoonCrater) -> None:
        """Draw the rim outline and field of the crater before the ray seed, so random numbers are drawn in a fixed order."""
        super()._prepare_crater(crater)
        _ = crater.hfield
        _ = crater.scarp_outline

    def _form_craters_binding(self, craters: list[RealMoonCrater], **kwargs: Any):
        data = [crater.hfield for crater in craters]
        scarps = [crater.scarp_outline for crater in craters]
        walls = [crater.wall for crater in craters]
        return realmoon_bindings.form_craters(
            craters=craters, rims=[d[0] for d in data], fields=[d[1] for d in data], scarp_radii=[s[0] for s in scarps],
            scarp_heights=[s[1] for s in scarps], wall_nominals=[w[0] for w in walls], wall_levels=[w[1] for w in walls],
            wall_params=[w[2] for w in walls], rho_max=_RHO_MAX, **kwargs
        )

    def crater_profile(
        self,
        crater: RealMoonCrater,
        radial_distances: ArrayLike,
        bearings: ArrayLike | None = None,
        reference_elevations: ArrayLike | None = None,
        crater_cls: type[Crater] = RealMoonCrater,
        include_crater: bool = True,
        include_ejecta: bool = False,
        **kwargs: Any,
    ) -> NDArray[np.float64]:
        """
        Compute the elevation of a realistic lunar crater (see `BasicMoonMorphology.crater_profile`).

        The profile does not include the central peak, which is built on the surface when the crater is formed.

        Parameters
        ----------
        crater : RealMoonCrater
            The crater.
        radial_distances : ArrayLike
            Radial distances from the crater center (m).
        bearings : ArrayLike, optional
            Bearings (degrees, clockwise from north) of the points. Default is 0.
        reference_elevations : ArrayLike, optional
            Reference elevations modified by the crater. Default is 0.
        crater_cls : type[Crater], optional
            Ignored: the crater is always cast to a `RealMoonCrater`.
        include_crater, include_ejecta : bool, optional
            Which parts of the profile to include. Where the basic profile stands above the reference surface, the field is shared between the crater and the ejecta in proportion to their parts of its height, so the ejecta thickness is never negative where the topography is above the reference surface.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        NDArray[np.float64]
            The elevation at each point.
        """
        # Always a realistic lunar crater (callers such as ejecta_profile pass the basic crater class)
        if not isinstance(crater, RealMoonCrater):
            crater = RealMoonCrater.maker(crater, morphology=self)
        r = np.asarray(radial_distances, dtype=np.float64)
        shape = r.shape
        r = np.ravel(r)
        b = np.zeros_like(r) if bearings is None else np.ravel(np.radians(np.asarray(bearings, dtype=np.float64)))
        ref = np.zeros_like(r) if reference_elevations is None else np.ravel(np.asarray(reference_elevations, dtype=np.float64))
        rim, field, rho_max = crater.hfield
        h = realmoon_bindings.realmoon_profile(
            radial_distances=np.ascontiguousarray(r),
            bearings=np.ascontiguousarray(np.broadcast_to(b, r.shape)),
            reference_elevations=np.ascontiguousarray(np.broadcast_to(ref, r.shape)),
            crater=crater,
            rim=rim,
            field=field,
            scarp_radius=crater.scarp_outline[0],
            scarp_height=crater.scarp_outline[1],
            wall_nominal=crater.wall[0],
            wall_levels=crater.wall[1],
            wall_params=crater.wall[2],
            rho_max=rho_max,
            include_crater=include_crater,
            include_ejecta=include_ejecta,
        )
        return np.reshape(np.asarray(h, dtype=np.float64), shape)
