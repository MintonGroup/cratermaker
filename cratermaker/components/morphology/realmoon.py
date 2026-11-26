from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray

from cratermaker.bindings import realmoon_bindings
from cratermaker.components.crater import Crater
from cratermaker.components.morphology import Morphology, MorphologyCraterVariable
from cratermaker.components.morphology.basicmoon import BasicMoonCrater, BasicMoonCraterFixed, BasicMoonMorphology

# The H-field model of fresh lunar craters (Lunar_Morphometry_Revisited, data/hfield/model.nc), distributed with the package
# in compact form (the modes packed as 16-bit integers, accurate to ~2e-5). Another model file can be given with the
# hfield_model_file argument or the CRATERMAKER_HFIELD_MODEL environment variable.
_HFIELD_MODEL_FILE = Path(__file__).resolve().parent / "hfield_model.nc"

# The polar grid of the model: 720 bearings and warped radii rho = r / R_rim from 0 to 2 in steps of 0.005
_N_THETA = 720
_RHO_MAX = 2.0
_DRHO = 0.005
_KMAX = _N_THETA // 2 - 1
_KMAX_RIM = 8  # the rim outline has azimuthal wavenumbers k <= 8
# Morphology types without a node of their own use those of complex craters
_TYPE_ALIASES = {"basin": "complex", "peakring": "complex", "ring": "complex"}

# Central peaks (built on the surface by the central peak model of the Rust crate, see central_peak.rs). The height of the
# relaxed peak above the floor, as a fraction of the crater's relief (rim height + floor depth), follows the peaks of the fresh
# craters of 36-139 km measured at the native SLDEM2015 spacing: 0.259 (D / 60 km)^0.89 with lognormal scatter 0.44. Above
# the fitted range it is held at its value at 139 km (it would otherwise exceed the relief by 300 km); below it, it is
# extrapolated, giving the smaller peaks of transitional craters. The envelope of the peak is a Gaussian of standard deviation
# 0.08 R. Multiring basins have no central peak.
_PEAKED = ("transitional", "complex")
_PEAK_HEIGHT_RELIEF = 0.259
_PEAK_HEIGHT_D_REF = 60e3
_PEAK_HEIGHT_EXPONENT = 0.89
_PEAK_HEIGHT_SIGMA_LN = 0.44
_PEAK_HEIGHT_D_MAX = 139e3
_PEAK_WIDTH = 0.08


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
    field = getattr(getattr(source, "_var", None), "_hfield", None) if source is not None else None
    if field is not None and source.radius == target.radius and source.morphology_type == target.morphology_type:
        object.__setattr__(target._var, "_hfield", field)


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
            peaked = new.morphology_type in _PEAKED and not new.isring
            relief = new.rim_height - new.floor_elevation
            rng = getattr(morphology if morphology is not None else new.morphology, "rng", None)
            new = super().maker(
                crater=new,
                morphology=morphology,
                peak_height=_peak_height(new.diameter, relief, rng) if peaked else 0.0,
                peak_width=_PEAK_WIDTH * new.radius if peaked else 0.0,
                peak_ring_radius=kwargs.get("peak_ring_radius") or 0.0,
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


@Morphology.register("realmoon")
class RealMoonMorphology(BasicMoonMorphology):
    """
    The realistic lunar crater morphology.

    A crater is the basic lunar crater (its profile and its size relations, with the ejecta and the degradation of the terrain
    under it), warped by an irregular rim outline, plus a random field of variation (terraces, slumps, the texture of the
    floor and walls, and the texture of the ejecta). The outline and field are drawn from the H-field model of fresh lunar
    craters: for each morphology type and size, the covariance between warped radii of the azimuthal Fourier components of
    the variation, decomposed into radial modes (empirical orthogonal functions). Central peaks and peak rings are not part of
    the field: transitional and complex craters get a central peak built on the surface as a patchy population of columns
    (where noise exceeds a threshold, raised to a Gaussian envelope of standard deviation `peak_width`, offset to radius
    `peak_ring_radius` for a peak ring) relaxed by slope collapse to the height `peak_height`. The model is described in the
    morphometry of fresh lunar craters by Minton et al. (2026, in prep.).

    Parameters
    ----------
    hfield_pix : float, optional
        The resolution (m) the field is drawn to. Default: the size of the face at the center of each crater.
    **kwargs : Any
        Arguments of `BasicMoonMorphology`.
    """

    def __init__(self, hfield_model_file: str | Path | None = None, hfield_pix: float | None = None, **kwargs: Any):
        if hfield_model_file is None:
            hfield_model_file = os.environ.get("CRATERMAKER_HFIELD_MODEL", _HFIELD_MODEL_FILE)
        hfield_model_file = Path(hfield_model_file)
        if not hfield_model_file.exists():
            raise FileNotFoundError(
                f"H-field model file {hfield_model_file} not found: pass hfield_model_file or set CRATERMAKER_HFIELD_MODEL"
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
        pi rho_max R / pix and steps in rho of at least pix / (2 R), with pix `hfield_pix` or the size of the face at the
        crater's center.
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
            try:
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

        F = np.zeros((n_theta // 2 + 1, len(elements)), dtype=np.complex128)
        ks = np.arange(1, kmax + 1)
        # The wavenumber node whose modes are used at each k
        nearest = np.argmin(np.abs(np.log(m["k_node"])[None, :] - np.log(ks)[:, None]), axis=1)
        for node, weight in self._nodes(crater.morphology_type, crater.diameter):
            lam_all, vec_all, var_k = m["eigenvalue"][node], m["eigenvector"][node], m["variance"][node]
            for j in np.unique(nearest):
                kj = ks[nearest == j]
                lam = lam_all[j]
                keep = np.where(lam > 0)[0]
                if len(keep) == 0:
                    continue
                kj = kj[var_k[kj] > 0]
                if len(kj) == 0:
                    continue
                # Modes of the node, restricted to the drawn radii, scaled to the variance of the node at each k
                V = vec_all[j][elements][:, keep].astype(np.float64)
                scale = np.sqrt(lam[keep] / lam[keep].sum())[None, :] * np.sqrt(var_k[kj])[:, None]
                z = (
                    self.rng.standard_normal((len(kj), len(keep))) + 1j * self.rng.standard_normal((len(kj), len(keep)))
                ) / np.sqrt(2.0)
                F[kj] += np.sqrt(weight) * ((scale * z) @ V.T)
        X = np.fft.irfft(F * n_theta, n=n_theta, axis=0)
        Fr = np.fft.rfft(X[:, 0])
        Fr[_KMAX_RIM + 1 :] = 0.0
        rim = np.exp(np.fft.irfft(Fr, n=n_theta))
        relief = max(crater.rim_height - crater.floor_elevation, 0.0)
        field = np.ascontiguousarray(X[:, 1:] * relief, dtype=np.float64)
        return np.ascontiguousarray(rim, dtype=np.float64), field, float(rho_max)

    def _prepare_crater(self, crater: RealMoonCrater) -> None:
        """Draw the rim outline and field of the crater before the ray seed, so random numbers are drawn in a fixed order."""
        super()._prepare_crater(crater)
        _ = crater.hfield

    def _form_craters_binding(self, craters: list[RealMoonCrater], **kwargs: Any):
        data = [crater.hfield for crater in craters]
        return realmoon_bindings.form_craters(
            craters=craters, rims=[d[0] for d in data], fields=[d[1] for d in data], rho_max=_RHO_MAX, **kwargs
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
            Which parts of the profile to include. The field belongs to the crater inside the rim and to the ejecta outside.
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
            rings=crater.rings,
            rim=rim,
            field=field,
            rho_max=rho_max,
            include_crater=include_crater,
            include_ejecta=include_ejecta,
        )
        return np.reshape(np.asarray(h, dtype=np.float64), shape)
