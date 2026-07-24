from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad
from scipy.optimize import root_scalar
from tqdm import tqdm

from cratermaker.bindings import basicmoon_bindings
from cratermaker.components.crater import Crater, CraterFixed
from cratermaker.components.morphology import Morphology, MorphologyCrater, MorphologyCraterVariable
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.constants import _VSMALL, FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter
from cratermaker.utils.montecarlo_utils import bounded_norm

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface

_EJPROFILE = -3.0


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
    parent: np.uint32 | None = field(default=False, init=True)

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


class BasicMoonCraterVariable(MorphologyCraterVariable):
    def __init__(self, frac_ejrim: float | None = None, rings: list[BasicMoonCrater] = None, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_rings", [])
        object.__setattr__(self, "_frac_ejrim", None)
        if frac_ejrim is not None:
            self.frac_ejrim = frac_ejrim
        if rings is not None:
            if isinstance(rings, list):
                self._rings = rings
            elif isinstance(rings, BasicMoonCrater):
                self._ring = [rings]
            else:
                raise ValueError("ring must be a scalar or list of BasicMoonCrater objects")
        return

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the crater variable properties.
        """
        dict_repr = super().as_dict()
        dict_repr["frac_ejrim"] = self.frac_ejrim
        dict_repr["rings"] = self.rings

        return dict_repr

    @property
    def frac_ejrim(self) -> float | None:
        """Ejecta rim thickness of the crater in meters."""
        return self._frac_ejrim

    @frac_ejrim.setter
    def frac_ejrim(self, value: float | None):
        if value is None:
            self._frac_ejrim = None
        else:
            if value < 0.0:
                raise ValueError("frac_ejrim must be positive.")
            self._frac_ejrim = value
        return

    @property
    def rings(self) -> list[BasicMoonCrater]:
        if self.nrings == 0:
            return None
        else:
            return self._rings

    @property
    def nrings(self) -> int:
        """
        Returns the number of rings associated with this Crater.
        """
        if self._rings is None:
            return 0
        else:
            return len(self._rings)


@Crater.register("basicmooncrater")
class BasicMoonCrater(MorphologyCrater):
    def __init__(
        self,
        crater: Crater | None = None,
        fixed_cls=BasicMoonCraterFixed,
        variable_cls=BasicMoonCraterVariable,
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
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += (
            f"Rim height: {format_large_units(self.rim_height, quantity='length')}\n"
            f"Ejecta fraction of rim: {self.frac_ejrim}\n"
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
        **kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        References
        ----------
        TODO: Update with Minton et al. (2026) info once done.
        """
        input_args = locals()
        from cratermaker.components.morphology import Morphology
        from cratermaker.utils.montecarlo_utils import sample_logfit_heteroskedastic, sample_pikefit

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

        morphology = Morphology.maker(morphology, **kwargs)
        crater = super().maker(crater=crater, morphology=morphology, **kwargs)
        monte_carlo_scaling = morphology.scaling.monte_carlo_scaling
        compute_nominal = not monte_carlo_scaling
        rng = morphology.rng
        min_valid_diameter = 500.0  # Cutoff where the model is constrained.
        depth_params = {
            "simple": {
                "coefficients": [-4.12929453954893, 1.6186582785253847, -0.03654228089366738],
                "c": -0.7978977562214146,
                "alpha": 1.230288898576757,
            },
            "transitional": {
                "coefficients": [13.544394609003302, -1.4349772241877992, 0.0864974576815312],
                "c": 29.832704911091195,
                "alpha": -2.4195839537761463,
            },
            "complex": {
                "coefficients": [-8.088308862863402, 2.6975260607735225, -0.11054039080484096],
                "c": 7.993093743747507,
                "alpha": 0.3208726177481036,
            },
        }
        rim_height_params = {
            "simple": {
                "coefficients": [-8.271313414224338, 2.1780660034399393, -0.06590590019301233],
                "c": -2.0063879230179578,
                "alpha": 1.2479301339151092,
            },
            "transitional": {
                "coefficients": [-8.059278782540039, 2.554156289155327, -0.1096155984856503],
                "c": -18.211075653420696,
                "alpha": 4.085584712178312,
            },
            "complex": {
                "coefficients": [-8.891646866181162, 2.5083058440327792, -0.09669451537649532],
                "c": -14.873628790748345,
                "alpha": 3.4661785636393705,
            },
        }
        rim_width_params = {
            "simple": {
                "coefficients": [-9.488293680414406, 2.810408424718434, -0.10907145408865704],
                "c": -11.491659223438303,
                "alpha": 2.914906260269271,
            },
            "transitional": {
                "coefficients": [17.281830217600284, -2.926699529590276, 0.19917534369662218],
                "c": -12.320305204746344,
                "alpha": 3.005434864001613,
            },
            "complex": {
                "coefficients": [-0.5588749331872054, 0.741932527414169, 0.010515195948127332],
                "c": -7.146938235118926,
                "alpha": 2.314670062281234,
            },
        }
        floor_radius_params = {
            "simple": {
                "coefficients": [-22.686946304312624, 5.868957235073512, -0.2859711767854806],
                "c": -6.216199203883163,
                "alpha": 2.255415282938512,
            },
            "transitional": {
                "coefficients": [-38.82514199442406, 8.388394016848952, -0.36490680476970033],
                "c": 4.139305903932077,
                "alpha": 0.9166479383008681,
            },
            "complex": {
                "coefficients": [8.863315887146682, -0.8889562385083289, 0.08746216954176912],
                "c": 18.70159621595362,
                "alpha": -0.4694908752922602,
            },
        }
        args = {}
        diameter_m = crater.diameter if crater.diameter > min_valid_diameter else min_valid_diameter
        diameter_km = diameter_m * 1e-3
        fcorrection = crater.diameter / diameter_m

        if crater.morphology_type in ["basin", "multiring", "peakring", "ring"]:
            morphology_type = "complex"
        else:
            morphology_type = crater.morphology_type

        # Ejecta thickness at the rim nominal value McGetchin, Settle, and Head (1973)
        ejrim = 0.14 * (diameter_m / 2) ** 0.74 * fcorrection

        if rim_height is None:
            rim_height = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_height_params[morphology_type]
                )[0]
                * fcorrection,
                ejrim,
            )
        args["rim_height"] = rim_height

        if rim_width is None:
            rim_width = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **rim_width_params[morphology_type]
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
                frac_ejrim = ejrim / rim_height
                frac_ejrim = min(frac_ejrim, 1.0)
            else:
                frac_ejrim = 0.0
        args["frac_ejrim"] = frac_ejrim

        if floor_elevation is None:
            floor_elevation = (
                -sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **depth_params[morphology_type]
                )[0]
                * fcorrection
                + rim_height
            )
            floor_elevation = min(floor_elevation, 0.0)
        args["floor_elevation"] = floor_elevation

        if floor_radius is None:
            floor_radius = max(
                sample_logfit_heteroskedastic(
                    diameter_m, compute_nominal=compute_nominal, rng=rng, **floor_radius_params[morphology_type]
                )[0]
                * fcorrection,
                0.0,
            )
        args["floor_radius"] = min(floor_radius, 0.8 * crater.radius)

        if peak_height is None:
            if crater.morphology_type == "complex":
                peak_height = (
                    sample_pikefit(
                        diameter_km, compute_nominal=compute_nominal, rng=rng, a=0.900, b=0.032, errhi=0.0011, errlo=-0.008, n=22
                    )[0]
                    * 1e3
                )
            else:
                peak_height = 0.0
        args["peak_height"] = peak_height
        args["peak_width"] = args["peak_height"] * 2 if peak_width is None else peak_width
        args["peak_ring_radius"] = 0.0 if peak_ring_radius is None else peak_ring_radius
        args["peak_center_distance"] = 0.0 if peak_center_distance is None else peak_center_distance
        args["peak_center_bearing"] = 0.0 if peak_center_bearing is None else peak_center_bearing

        if wall_curvature is None:
            if monte_carlo_scaling:
                wall_curvature = rng.uniform(low=0.0, high=0.5, size=1)[0]  # Temporary until a morphometric analysis ic complete
            else:
                wall_curvature = 1.0

        args["wall_curvature"] = wall_curvature

        if floor_blend is None:
            if monte_carlo_scaling:
                floor_blend = rng.uniform(low=0.0, high=0.5, size=1)[0]  # Temporary until a morphometric analysis ic complete
            else:
                floor_blend = 0.25

        args["floor_blend"] = floor_blend

        kwargs = {**args, **kwargs}

        crater = cls(
            crater=crater,
            morphology=morphology,
            **kwargs,
        )

        if crater.morphology_type == "multiring" and crater.nrings == 0 and not crater.isring:
            num_rings = 3  # kwargs.pop("num_rings", rng.integers(low=2, high=4))
            for i in range(num_rings):
                if monte_carlo_scaling:
                    rnd_factor = bounded_norm(loc=1.0, scale=0.1, size=4, lower_bound=0.0, upper_bound=1.0, rng=rng)
                else:
                    rnd_factor = np.ones(4)
                radius = crater.radius * rnd_factor[0] / np.sqrt(2.0) ** (i + 1)
                floor_radius = crater.floor_radius * rnd_factor[1] / np.sqrt(2.0) ** (i + 1)
                rim_height = rnd_factor[3] * crater.rim_height * (0.4) ** (i + 1)
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

        # Make sure rings are the correct type
        for i in range(crater.nrings):
            ring = crater.rings[i]
            if not isinstance(ring, cls):
                ring = cls(crater=ring, morphology=morphology)
                crater.rings[i] = ring

        # Adjust frac_ejrim value(s) in order to get closer to a volume-conserving solution for the ejecta
        if conserve_volume and not crater.isring:

            def _func(frac_ejrim, crater):
                ejrim = frac_ejrim * crater.rim_height
                frac_ejrim_orig = crater.frac_ejrim
                ejrim_orig = crater.frac_ejrim * crater.rim_height
                if ejrim_orig > 0.0:
                    erat = max(ejrim / ejrim_orig, 0.0)
                else:
                    erat = 0.0
                crater.frac_ejrim = max(ejrim / crater.rim_height, 0.0)
                ring_frac_ejrim_orig = []
                if crater.nrings > 0:
                    for ring in crater.rings:
                        if ring.frac_ejrim is not None:
                            ring_frac_ejrim_orig.append(ring.frac_ejrim)
                            ring_ejrim = ring.frac_ejrim * ring.rim_height
                            ring_ejrim *= erat
                            ring.frac_ejrim = ring_ejrim / ring.rim_height
                excavated_volume = morphology.estimate_volume(crater, include_crater=True, include_ejecta=False)
                ejecta_volume = morphology.estimate_volume(crater, include_crater=False, include_ejecta=True)
                result = ejecta_volume + excavated_volume
                crater.frac_ejrim = frac_ejrim_orig
                if len(ring_frac_ejrim_orig) > 0:
                    for ring, ring_frac_ejrim in zip(crater.rings, ring_frac_ejrim_orig, strict=True):
                        ring.frac_ejrim = ring_frac_ejrim
                return result

            for _ in range(10):
                lower_bracket = 0.0
                upper_bracket = 1.0
                lower_bound = _func(lower_bracket, crater)
                upper_bound = _func(upper_bracket, crater)
                if lower_bound < 0.0 and upper_bound > 0.0:
                    break
                while lower_bound > 0.0:
                    # This occurs when rim_height is too high
                    kwargs["rim_height"] *= 0.9
                    crater = cls(
                        crater=crater,
                        morphology=morphology,
                        **kwargs,
                    )
                    lower_bound = _func(lower_bracket, crater)

                upper_bound = _func(upper_bracket, crater)
                while upper_bound < 0.0:
                    # This occurs when the rim_height is too low
                    kwargs["rim_height"] *= 1.1
                    crater = cls(
                        crater=crater,
                        morphology=morphology,
                        **kwargs,
                    )
                    upper_bound = _func(upper_bracket, crater)

            sol = root_scalar(lambda x, crater=crater: _func(x, crater), bracket=[lower_bracket, upper_bracket], method="brentq")
            frac_ejrim = sol.root if sol.converged else crater.frac_ejrim
            if crater.frac_ejrim > _VSMALL:
                conservation_factor = frac_ejrim / crater.frac_ejrim
            else:
                conservation_factor = 1.0
            crater._frac_ejrim *= conservation_factor
            if crater.nrings > 0:
                for ring in crater.rings:
                    if ring.frac_ejrim is not None and ring.frac_ejrim > 0.0:
                        ring.frac_ejrim *= conservation_factor

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
        if ejecta_thickness is None or ejecta_intensity is None:
            return None, None
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
        reference_elevation = region.get_reference_surface(reference_radius=crater.diameter)

        # Combine distances and references for nodes and faces
        radial_distances = np.concatenate([region.face_distance, region.node_distance])
        bearings = np.concatenate([region.face_bearing, region.node_bearing])

        original_elevation = np.concatenate([region.face_elevation, region.node_elevation])
        reference_elevation[radial_distances > crater.radius] = original_elevation[radial_distances > crater.radius]

        new_elevation = self.crater_profile(
            crater=crater, radial_distances=radial_distances, bearings=bearings, reference_elevations=reference_elevation
        )

        elevation_change = new_elevation - original_elevation
        return elevation_change

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

        radial_distances = np.concatenate([region.face_distance, region.node_distance])
        bearings = np.concatenate([region.face_bearing, region.node_bearing])
        if self.dorays:
            thickness, intensity = self.ejecta_distribution(crater, radial_distances, bearings)
        else:
            thickness = self.ejecta_profile(crater, radial_distances, bearings)
            intensity = np.ones_like(radial_distances)

        return thickness, intensity

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
        # reshape thickness to match the shape of r and theta
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
        if not isinstance(crater, BasicMoonCrater):
            crater = BasicMoonCrater.maker(crater, morphology=self)
        # flatten r and theta to 1D arrays
        radial_distances = np.ravel(radial_distances)
        bearings = np.ravel(np.radians(bearings))
        intensity = basicmoon_bindings.ray_intensity(
            radial_distances,
            bearings,
            crater.diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        intensity = np.array(intensity, dtype=np.float64)
        # reshape intensity to match the shape of r and theta
        intensity = np.reshape(intensity, radial_distances.shape)
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

    def _profile_invert_ejecta(self, r, crater, minimum_thickness):
        return basicmoon_bindings.ejecta_profile_function(crater, r) - minimum_thickness

    def _profile_invert_crater(self, r, crater, minimum_thickness):
        return basicmoon_bindings.crater_profile_function(crater, r) - minimum_thickness

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

        def _invert_ejecta(r):
            return self._profile_invert_ejecta(r, crater, minimum_thickness)

        def _invert_crater(r):
            return self._profile_invert_crater(r, crater, minimum_thickness)

        if feature == "ejecta":
            _profile_invert = _invert_ejecta
        elif feature == "crater":
            _profile_invert = _invert_crater
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

    def estimate_volume(
        self,
        crater: BasicMoonCrater,
        include_crater: bool,
        include_ejecta: bool,
        crater_cls: type[Crater] = BasicMoonCrater,
        profile_func: Callable = basicmoon_bindings.basicmoon_profile,
        **kwargs,
    ) -> np.float64:
        """
        Estimates the volume change of a crater by integrating the profile function using the scipy.integrate.quad function.
        """
        if not (include_crater or include_ejecta):
            return np.float64(0.0)

        def _crater_func(r):
            h = basicmoon_bindings.crater_profile_function(crater, r)
            return r * h

        def _ejecta_func(r):
            h = basicmoon_bindings.ejecta_profile_function(crater, r)
            return r * h

        def _combo_func(r):
            hc = basicmoon_bindings.crater_profile_function(crater, r)
            he = basicmoon_bindings.ejecta_profile_function(crater, r)
            return r * (hc + he)

        if include_crater and not include_ejecta:
            func = _crater_func
        elif include_ejecta and not include_crater:
            func = _ejecta_func
        else:
            func = _combo_func
        v = quad(func, 0.0, 20 * crater.radius, **{"limit": 100, **kwargs}, full_output=1)[0]
        return 2 * np.pi * v

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
