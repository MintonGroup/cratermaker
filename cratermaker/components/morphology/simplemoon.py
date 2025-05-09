from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import fsolve

from cratermaker._simplemoon import crater_functions, ejecta_functions
from cratermaker.components.morphology import Morphology
from cratermaker.components.surface import Surface, SurfaceView
from cratermaker.constants import FloatLike
from cratermaker.core.crater import Crater
from cratermaker.utils.general_utils import format_large_units, parameter


@dataclass(frozen=True, slots=True)
class SimpleMoonCrater(Crater):
    rim_height: float | None = None
    rim_width: float | None = None
    floor_depth: float | None = None
    floor_diameter: float | None = None
    peak_height: float | None = None
    ejrim: float | None = None

    def __repr__(self) -> str:
        base = super().__repr__()
        return (
            f"{base}\n"
            f"Rim height: {format_large_units(self.rim_height, quantity='length')}\n"
            f"Rim width: {format_large_units(self.rim_width, quantity='length')}\n"
            f"Floor depth: {format_large_units(self.floor_depth, quantity='length')}\n"
            f"Floor diameter: {format_large_units(self.floor_diameter, quantity='length')}\n"
            f"Central peak height: {format_large_units(self.peak_height, quantity='length') if self.peak_height else 'None'}\n"
            f"Ejecta rim thickness: {format_large_units(self.ejrim, quantity='length')}\n"
        )

    @classmethod
    def maker(cls, crater: Crater | None = None, **kwargs) -> SimpleMoonCrater:
        if crater is None:
            crater = super(cls, cls).maker(**kwargs)
        base_fields = asdict(crater)
        # Remove any SimpleMoonCrater-specific fields to avoid conflicts
        morphology_fields = set(cls.__dataclass_fields__) - set(
            super(cls, cls).__dataclass_fields__
        )
        for key in morphology_fields:
            base_fields.pop(key, None)

        diameter_m = crater.final_diameter
        diameter_km = diameter_m * 1e-3

        if crater.morphology_type in ["simple", "transitional"]:
            rim_height = 0.043 * diameter_km**1.014 * 1e3
            rim_width = 0.257 * diameter_km**1.011 * 1e3
            floor_depth = 0.224 * diameter_km**1.010 * 1e3
            floor_diameter = 0.200 * diameter_km**1.143 * 1e3
            peak_height = None
        elif crater.morphology_type in ["complex", "peakring", "multiring"]:
            rim_height = 0.236 * diameter_km**0.399 * 1e3
            rim_width = 0.467 * diameter_km**0.836 * 1e3
            floor_depth = 1.044 * diameter_km**0.301 * 1e3
            floor_diameter = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m)
            peak_height = 0.032 * diameter_km**0.900 * 1e3
        else:
            raise ValueError(f"Unknown morphology type: {crater.morphology_type}")

        ejrim = 0.14 * (diameter_m * 0.5) ** 0.74

        return cls(
            **base_fields,
            rim_height=rim_height,
            rim_width=rim_width,
            floor_depth=floor_depth,
            floor_diameter=floor_diameter,
            peak_height=peak_height,
            ejrim=ejrim,
        )


@Morphology.register("simplemoon")
class SimpleMoon(Morphology):
    """
    An operations class for computing the morphology of a crater and applying it to a surface mesh.

    Parameters
    ----------
    crater : Crater, optional
        The crater currently attached to the morphology model.
    ejecta_truncation : float, optional
        The relative distance from the rim of the crater to truncate the ejecta blanket, default is None, which will compute a
        truncation distance based on where the ejecta thickness reaches a small value.
    rng : Generator, optional
        A random number generator instance. If not provided, the default numpy RNG will be used.
    dorays : bool, optional
        A flag to determine if the ray pattern should be used instead of the homogeneous ejecta blanket, default is True.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        ejecta_truncation: FloatLike | None = None,
        dorays: bool = True,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_ejecta_truncation", None)
        object.__setattr__(self, "_node", None)
        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    def __repr__(self) -> str:
        base = super().__repr__()
        if self.ejecta_truncation is not None:
            base += f"\nEjecta Trunction: {self.ejecta_truncation:.2f} * crater.final_radius"
        else:
            base += "\nEjecta Truncation: Off"
        return f"{base}\nEjecta Rays: {self.dorays}"

    def emplace(self, crater: Crater, surface: Surface, **kwargs: Any) -> None:
        """
        Convenience method to immediately emplace a crater onto the surface.
        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        crater : Crater
            The crater to be emplaced.
        surface : Surface
            The surface object to modify.
        kwargs : Any
            Additional keyword arguments to pass to the emplace method.
        """
        crater = SimpleMoonCrater.maker(crater)
        super().emplace(crater, surface)
        return

    def crater_shape(
        self, crater: SimpleMoonCrater, region_view: SurfaceView, surface: Surface
    ) -> Surface:
        """
        Compute the crater shape based on the region view and surface.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object containing the parameters for the crater shape.
        region_view : RegionView
            The region view of the surface mesh.
        surface : Surface
            The surface mesh.

        Returns
        -------
        surface : Surface
            The modified surface mesh with the crater shape applied.
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, crater.location
        )
        reference_face_elevation, reference_node_elevation = (
            surface.get_reference_surface(
                region_view,
                face_crater_distance,
                node_crater_distance,
                crater.location,
                crater.final_radius,
            )
        )

        # Combine distances and references for nodes and faces
        combined_distances = np.concatenate(
            [node_crater_distance, face_crater_distance]
        )
        combined_reference = np.concatenate(
            [reference_node_elevation, reference_face_elevation]
        )
        combined_elevation = self.crater_profile(
            crater, combined_distances, combined_reference
        )

        node_elevation = combined_elevation[: len(node_crater_distance)]
        face_elevation = combined_elevation[len(node_crater_distance) :]

        surface.node_elevation[region_view.node_indices] = node_elevation
        surface.face_elevation[region_view.face_indices] = face_elevation
        return surface

    def crater_profile(
        self, crater: SimpleMoonCrater, r: ArrayLike, r_ref: ArrayLike | None = None
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object containing the parameters for the crater profile.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        r_ref : ArrayLike, optional
            Reference elevation values to be modified by the crater profile.

        Returns
        -------
        elevation : NDArray[np.float64]
            The computed crater elevation profile at each radial point.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        if r_ref is None:
            r_ref = np.zeros_like(r)

        # flatten r to 1D array
        rflat = np.ravel(r)
        r_ref_flat = np.ravel(r_ref)
        elevation = crater_functions.profile(
            rflat,
            r_ref_flat,
            crater.final_diameter,
            crater.floor_depth,
            crater.floor_diameter,
            crater.rim_height,
            crater.ejrim,
        )
        # reshape elevation to match the shape of r
        elevation = np.array(elevation, dtype=np.float64)
        elevation = np.reshape(elevation, r.shape)

        return elevation

    def ejecta_shape(
        self, crater: SimpleMoonCrater, region_view: SurfaceView, surface: Surface
    ) -> Surface:
        """
        Compute the ejecta shape based on the region view and surface.

        Parameters
        ----------
        region_view : RegionView
            The region view of the surface mesh.
        surface : Surface
            The surface mesh.

        Returns
        -------
        surface : Surface
            The modified surface mesh with the ejecta shape applied.
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, crater.location
        )
        if self.dorays:
            node_crater_bearings, face_crater_bearings = surface.get_initial_bearing(
                region_view, crater.location
            )
            combined_distances = np.concatenate(
                [node_crater_distance, face_crater_distance]
            )
            combined_bearings = np.concatenate(
                [node_crater_bearings, face_crater_bearings]
            )
            combined_thickness, combined_ray_intensity = self.ejecta_distribution(
                crater, combined_distances, combined_bearings
            )
            surface.ray_intensity[region_view.face_indices] = combined_ray_intensity[
                len(node_crater_distance) :
            ]
        else:
            combined_distances = np.concatenate(
                [node_crater_distance, face_crater_distance]
            )
            combined_thickness = self.ejecta_profile(crater, combined_distances)

        # Slice back the combined thickness without copying
        node_thickness = combined_thickness[: len(node_crater_distance)]
        face_thickness = combined_thickness[len(node_crater_distance) :]

        surface.ejecta_thickness[region_view.face_indices] = face_thickness
        surface.face_elevation[region_view.face_indices] += face_thickness
        surface.node_elevation[region_view.node_indices] += node_thickness

        return surface

    def ejecta_profile(
        self, crater: SimpleMoonCrater, r: ArrayLike
    ) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
        crater : SimpleMoonCrater
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        # flatten r to 1D array
        rflat = np.ravel(r)
        elevation = ejecta_functions.profile(rflat, crater.final_diameter, crater.ejrim)
        elevation = np.array(elevation, dtype=np.float64)
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)
        return elevation

    def ejecta_distribution(
        self, crater, r: ArrayLike, theta: ArrayLike
    ) -> NDArray[np.float64]:
        """
        Compute the ejecta thickness distribution modulated by ray patterns.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object containing the parameters for the ejecta distribution.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        theta : ArrayLike
            Angular bearings from the crater center (in radians).

        Returns
        -------
        thickness : NDArray[np.float64]
            The computed ejecta thickness for each (r, theta) pair.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        # flatten r and theta to 1D arrays
        thickness = self.ejecta_profile(crater, r)
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        ray_intensity = ejecta_functions.ray_intensity(
            rflat,
            theta_flat,
            crater.final_diameter,
        )
        thickness = np.array(thickness, dtype=np.float64)
        ray_intensity = np.array(ray_intensity, dtype=np.float64)
        thickness *= ray_intensity
        # reshape thickness to match the shape of r and theta
        thickness = np.reshape(thickness, r.shape)
        ray_intensity = np.reshape(ray_intensity, r.shape)
        return thickness, ray_intensity

    def ray_intensity(
        self, crater: SimpleMoonCrater, r: ArrayLike, theta: ArrayLike
    ) -> NDArray[np.float64]:
        """
        Compute the ray pattern intensity modulation at each (r, theta) pair.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object containing the parameters for the ray intensity.
        r : ArrayLike
            Radial distances from the crater center (in meters).
        theta : ArrayLike
            Angular bearings from the crater center (in radians).

        Returns
        -------
        intensity : NDArray[np.float64]
            The computed ray intensity values.

        Notes
        -----
        This is a wrapper for a compiled Rust function.
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)
        # flatten r and theta to 1D arrays
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        intensity = ejecta_functions.ray_intensity(
            rflat,
            theta_flat,
            crater.final_diameter,
        )
        intensity = np.array(intensity, dtype=np.float64)
        # reshape intensity to match the shape of r and theta
        intensity = np.reshape(intensity, r.shape)
        return intensity

    def rmax(
        self,
        crater: SimpleMoonCrater,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
    ) -> float:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor,
        whichever is smaller.

        Parameters
        ----------
        crater : Crater
            The crater object to be used. If None, the current crater object is used. If passed, the current crater object is replaced.
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

        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater)

        def _profile_invert_ejecta(r):
            return self.ejecta_profile(crater, r) - minimum_thickness

        def _profile_invert_crater(r):
            return self.crater_profile(crater, r, np.zeros(1)) - minimum_thickness

        if feature == "ejecta":
            _profile_invert = _profile_invert_ejecta
        elif feature == "crater":
            _profile_invert = _profile_invert_crater
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")

        # Get the maximum extent
        rmax = fsolve(_profile_invert, x0=crater.final_radius * 1.01)[0]

        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * crater.final_radius)

        if feature == "crater":
            rmax = max(rmax, crater.final_radius)

        return float(rmax)

    @property
    def node_index(self):
        """
        The index of the node closest to the crater location.

        Returns
        -------
        int
        """
        return self._node_index

    @node_index.setter
    def node_index(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError("node_index must be of type int")
        self._node_index = value

    @property
    def face_index(self):
        """
        The index of the face closest to the crater location.

        Returns
        -------
        int
        """
        return self._face_index

    @face_index.setter
    def face_index(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError("face_index must be of type int")
        self._face_index = value

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
