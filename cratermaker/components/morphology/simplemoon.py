from typing import Any

import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import fsolve

from cratermaker._simplemoon import crater_functions, ejecta_functions
from cratermaker.components.morphology import Morphology
from cratermaker.constants import FloatLike
from cratermaker.core.crater import Crater
from cratermaker.utils.general_utils import format_large_units, parameter


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
        crater: Crater | None = None,
        ejecta_truncation: FloatLike | None = None,
        dorays: bool = True,
        rng: Generator | None = None,
        rng_seed: int | None = None,
        rng_state: dict | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_ejecta_truncation", None)
        object.__setattr__(self, "_rim_height", None)
        object.__setattr__(self, "_rim_width", None)
        object.__setattr__(self, "_peak_height", None)
        object.__setattr__(self, "_floor_diameter", None)
        object.__setattr__(self, "_floor_depth", None)
        object.__setattr__(self, "_ejrim", None)
        object.__setattr__(self, "_node", None)
        self.ejecta_truncation = ejecta_truncation
        self.dorays = dorays
        super().__init__(
            crater=crater, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs
        )

    def __repr__(self) -> str:
        base = super().__repr__()
        if self.ejecta_truncation is not None:
            base += f"\nEjecta Trunction: {self.ejecta_truncation:.2f} * crater.final_radius"
        else:
            base += "\nEjecta Truncation: Off"
        if self.crater is not None:
            base += f"\nRim height: {format_large_units(self.rim_height, quantity='length')}"
            base += (
                f"\nRim width: {format_large_units(self.rim_width, quantity='length')}"
            )
            base += f"\nFloor depth: {format_large_units(self.floor_depth, quantity='length')}"
            base += f"\nFloor diameter: {format_large_units(self.floor_diameter, quantity='length')}"
            if self.peak_height is not None:
                base += f"\nPeak height: {format_large_units(self.peak_height, quantity='length')}"
            base += f"\nEjecta thickness at rim: {format_large_units(self.ejrim, quantity='length')}"
        return f"{base}\nEjecta Rays: {self.dorays}"

    def _set_morphology_config(self) -> None:
        """
        This method adds a crater to the morphology model and sets the parameters for the morphology based on the crater type.
        """
        # Set the morphology based on crater type
        diameter_m = self.crater.final_diameter
        diameter_km = diameter_m * 1e-3  # Convert to km for some models

        if self.crater.morphology_type in ["simple", "transitional"]:
            # A hybrid model between Pike (1977) and Fassett & Thomson (2014)
            self._rim_height = (
                0.043 * diameter_km**1.014 * 1e3
            )  # Closer to Fassett & Thomson
            self._rim_width = 0.257 * diameter_km**1.011 * 1e3  # Pike model
            self._floor_depth = (
                0.224 * diameter_km**1.010 * 1e3
            )  # Closer to Fassett & Thomson
            self._floor_diameter = (
                0.200 * diameter_km**1.143 * 1e3
            )  # Fassett & Thomson for D~1km, Pike for D~20km
            self._peak_height = None
        elif self.crater.morphology_type in ["complex", "peakring", "multiring"]:
            # Following Pike (1977)
            self._rim_height = 0.236 * diameter_km**0.399 * 1e3  # Pike model
            self._rim_width = 0.467 * diameter_km**0.836 * 1e3  # Pike model
            self._floor_depth = 1.044 * diameter_km**0.301 * 1e3  # Pike model
            # Fassett & Thomson for D~1km, Pike for D~20km, but limited to 90% of diameter
            self._floor_diameter = min(
                0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m
            )
            self._peak_height = 0.032 * diameter_km**0.900 * 1e3  # Pike model
        else:
            raise ValueError(f"Unknown morphology type: {self.crater.morphology_type}")

        self._ejrim = 0.14 * (diameter_m * 0.5) ** (
            0.74
        )  # McGetchin et al. (1973) Thickness of ejecta at rim
        return

    def crater_shape(self, region_view, surface):
        """
        Compute the crater shape based on the region view and surface.

        Parameters
        ----------
        region_view : RegionView
            The region view of the surface mesh.
        surface : Surface
            The surface mesh.

        Returns
        -------
        None
        """
        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, self.crater.location
        )
        reference_face_elevation, reference_node_elevation = (
            surface.get_reference_surface(
                region_view,
                face_crater_distance,
                node_crater_distance,
                self.crater.location,
                self.crater.final_radius,
            )
        )

        # Combine distances and references for nodes and faces
        combined_distances = np.concatenate(
            [node_crater_distance, face_crater_distance]
        )
        combined_reference = np.concatenate(
            [reference_node_elevation, reference_face_elevation]
        )
        combined_elevation = self.crater_profile(combined_distances, combined_reference)

        node_elevation = combined_elevation[: len(node_crater_distance)]
        face_elevation = combined_elevation[len(node_crater_distance) :]

        surface.node_elevation[region_view.node_indices] = node_elevation
        surface.face_elevation[region_view.face_indices] = face_elevation
        return surface

    def crater_profile(
        self, r: ArrayLike, r_ref: ArrayLike | None = None
    ) -> NDArray[np.float64]:
        """
        Compute the crater profile elevation at a given radial distance.

        Parameters
        ----------
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
        if r_ref is None:
            r_ref = np.zeros_like(r)

        # flatten r to 1D array
        rflat = np.ravel(r)
        elevation = crater_functions.profile(
            rflat,
            r_ref,
            self.crater.final_diameter,
            self.floor_depth,
            self.floor_diameter,
            self.rim_height,
            self.ejrim,
        )
        # reshape elevation to match the shape of r
        elevation = np.array(elevation, dtype=np.float64)
        elevation = np.reshape(elevation, r.shape)

        return elevation

    def ejecta_shape(self, region_view, surface):
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
        None
        """
        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, self.crater.location
        )
        if self.dorays:
            node_crater_bearings, face_crater_bearings = surface.get_initial_bearing(
                region_view, self.crater.location
            )
            combined_distances = np.concatenate(
                [node_crater_distance, face_crater_distance]
            )
            combined_bearings = np.concatenate(
                [node_crater_bearings, face_crater_bearings]
            )
            combined_thickness, combined_ray_intensity = self.ejecta_distribution(
                combined_distances, combined_bearings
            )
            surface.ray_intensity[region_view.face_indices] = combined_ray_intensity[
                len(node_crater_distance) :
            ]
        else:
            combined_distances = np.concatenate(
                [node_crater_distance, face_crater_distance]
            )
            combined_thickness = self.ejecta_profile(combined_distances)

        # Slice back the combined thickness without copying
        node_thickness = combined_thickness[: len(node_crater_distance)]
        face_thickness = combined_thickness[len(node_crater_distance) :]

        surface.ejecta_thickness[region_view.face_indices] = face_thickness
        surface.face_elevation[region_view.face_indices] += face_thickness
        surface.node_elevation[region_view.node_indices] += node_thickness

        return surface

    def ejecta_profile(self, r: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ejecta elevation profile at a given radial distance.

        Parameters
        ----------
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
        # flatten r to 1D array
        rflat = np.ravel(r)
        elevation = ejecta_functions.profile(
            rflat, self.crater.final_diameter, self.ejrim
        )
        elevation = np.array(elevation, dtype=np.float64)
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)
        return elevation

    def ejecta_distribution(
        self, r: ArrayLike, theta: ArrayLike
    ) -> NDArray[np.float64]:
        """
        Compute the ejecta thickness distribution modulated by ray patterns.

        Parameters
        ----------
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
        # flatten r and theta to 1D arrays
        thickness = self.ejecta_profile(r)
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        ray_intensity = ejecta_functions.ray_intensity(
            rflat,
            theta_flat,
            self.crater.final_diameter,
        )
        thickness = np.array(thickness, dtype=np.float64)
        ray_intensity = np.array(ray_intensity, dtype=np.float64)
        thickness *= ray_intensity
        # reshape thickness to match the shape of r and theta
        thickness = np.reshape(thickness, r.shape)
        ray_intensity = np.reshape(ray_intensity, r.shape)
        return thickness, ray_intensity

    def ray_intensity(self, r: ArrayLike, theta: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ray pattern intensity modulation at each (r, theta) pair.

        Parameters
        ----------
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
        # flatten r and theta to 1D arrays
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        intensity = ejecta_functions.ray_intensity(
            rflat,
            theta_flat,
            self.crater.final_diameter,
        )
        intensity = np.array(intensity, dtype=np.float64)
        # reshape intensity to match the shape of r and theta
        intensity = np.reshape(intensity, r.shape)
        return intensity

    def rmax(
        self,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
        crater: Crater | None = None,
    ) -> float:
        """
        Compute the maximum extent of the crater based on the minimum thickness of a feature, or the ejecta_truncation factor,
        whichever is smaller.

        Parameters
        ----------
        minimum_thickness : FloatLike
            The minimum thickness of the feature blanket in meters.
        feature : str, optional, default = "ejecta"
            The feature to compute the maximum extent. Either "crater" or "ejecta". If "crater" is chosen, the rmax is based
            on where the raised rim is smaller than minimum thickness.
        crater : Crater, optional
            The crater object to be used. If None, the current crater object is used. If passed, the current crater object is replaced.
        Returns
        -------
        float
            The maximum extent of the crater or ejecta blanket in meters.
        """

        # Compute the reference surface for the crater
        if crater is not None:
            self.crater = crater

        if feature == "ejecta":

            def _profile_invert(r):
                return self.ejecta_profile(r) - minimum_thickness
        elif feature == "crater":

            def _profile_invert(r):
                return self.crater_profile(r, np.zeros(1)) - minimum_thickness
        else:
            raise ValueError("Unknown feature type. Choose either 'crater' or 'ejecta'")

        # Get the maximum extent
        rmax = fsolve(_profile_invert, x0=self.crater.final_radius * 1.01)[0]

        if self.ejecta_truncation:
            rmax = min(rmax, self.ejecta_truncation * self.crater.final_radius)

        if feature == "crater":
            rmax = max(rmax, self.crater.final_radius)

        return float(rmax)

    @property
    def rim_height(self) -> float:
        """
        The height of the crater rim in meters.

        Returns
        -------
        float
        """
        return self._rim_height

    @property
    def rim_width(self) -> float:
        """
        The width of the crater rim in meters.

        Returns
        -------
        float
        """
        return self._rim_width

    @property
    def peak_height(self) -> float:
        """
        The height of the central peak in meters.

        Returns
        -------
        float
        """
        return self._peak_height

    @property
    def floor_diameter(self) -> float:
        """
        The diameter of the crater floor in meters.

        Returns
        -------
        float
        """
        return self._floor_diameter

    @property
    def floor_depth(self) -> float:
        """
        Return the depth of the crater floor in m

        Returns
        -------
        float
        """
        return self._floor_depth

    @property
    def ejrim(self) -> float:
        """
        The thickness of ejecta at the rim in m.

        Returns
        -------
        float
        """
        return self._ejrim

    @property
    def crater(self):
        return super().crater

    @crater.setter
    def crater(self, value):
        if value is None:
            self._crater = None
            return
        Morphology.crater.fset(self, value)
        self._set_morphology_config()

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
