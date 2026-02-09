from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from cratermaker._cratermaker import morphology_bindings
from numpy.random import Generator
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad
from scipy.optimize import root_scalar
from tqdm import tqdm

from cratermaker.components.crater import Crater, CraterFixed, CraterVariable
from cratermaker.components.morphology import Morphology, MorphologyCrater, MorphologyCraterVariable
from cratermaker.components.surface import LocalSurface, Surface
from cratermaker.constants import FloatLike
from cratermaker.utils.general_utils import format_large_units, parameter


@dataclass(frozen=True, slots=True)
class SimpleMoonCraterFixed(CraterFixed):
    rim_height: float | None = None
    rim_width: float | None = None
    floor_depth: float | None = None
    floor_diameter: float | None = None
    peak_height: float | None = None
    ejrim: float | None = None


class SimpleMoonCrater(MorphologyCrater):
    def __init__(self, crater: Crater | None = None, fixed_cls=SimpleMoonCraterFixed, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, **kwargs)
        return

    def __str__(self) -> str:
        base = super().__str__()
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
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        **kwargs: Any,
    ) -> SimpleMoonCrater:
        """
        Initialize a SimpleMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a SimpleMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        kwargs : Any
            The keyword arguments provided are passed down to :func:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.
        """
        morphology = Morphology.maker(morphology, **kwargs)
        if crater is None:
            crater = super().maker(**kwargs)
        args = {}
        diameter_m = crater.final_diameter
        diameter_km = diameter_m * 1e-3

        if crater.morphology_type in ["simple", "transitional"]:
            args["rim_height"] = 0.043 * diameter_km**1.014 * 1e3
            args["rim_width"] = 0.257 * diameter_km**1.011 * 1e3
            args["floor_depth"] = 0.224 * diameter_km**1.010 * 1e3
            args["floor_diameter"] = 0.200 * diameter_km**1.143 * 1e3
            args["peak_height"] = None
        elif crater.morphology_type in ["complex", "peakring", "multiring"]:
            args["rim_height"] = 0.236 * diameter_km**0.399 * 1e3
            args["rim_width"] = 0.467 * diameter_km**0.836 * 1e3
            args["floor_depth"] = 1.044 * diameter_km**0.301 * 1e3
            args["floor_diameter"] = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m)
            args["peak_height"] = 0.032 * diameter_km**0.900 * 1e3
        else:
            raise ValueError(f"Unknown morphology type: {crater.morphology_type}")

        args["ejrim"] = 0.14 * (diameter_m * 0.5) ** 0.74

        return cls(
            crater=crater,
            morphology=morphology,
            **args,
        )


@Morphology.register("simplemoon")
class SimpleMoon(Morphology):
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
        crater : Crater, optional
            The crater currently attached to the morphology model.
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

    def __str__(self) -> str:
        base = super().__str__()
        if self.ejecta_truncation is not None:
            base += f"\nEjecta Trunction: {self.ejecta_truncation:.2f} * crater.final_radius"
        else:
            base += "\nEjecta Truncation: Off"
        return f"{base}\nEjecta Rays: {self.dorays}"

    def emplace(self, craters: Crater, **kwargs: Any) -> None:
        """
        Convenience method to immediately emplace a crater onto the surface.

        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        crater : Crater | list[Crater] | None
            The crater or list of craters to be emplaced. If None, then a crater will be created using the provided parameters in kwargs and emplaced.
        kwargs : Any
            |kwargs|
        """
        if craters is None:
            craters = [SimpleMoonCrater.maker(morphology=self, **kwargs)]
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
                processed_craters.append(SimpleMoonCrater.maker(c, morphology=self, **kwargs))
            craters = processed_craters
        elif isinstance(craters, Crater):
            craters = [SimpleMoonCrater.maker(craters, morphology=self, **kwargs)]
        else:
            raise ValueError(
                "Invalid input for crater emplacement. Must be a Crater object, a list of Crater objects, or None with additional arguments for Crater.maker()."
            )

        super().emplace(craters)
        return

    def form_crater(self, crater: Crater | SimpleMoonCrater, **kwargs: Any) -> None:
        """
        Form a crater or list of craters on the surface.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object to be formed.
        kwargs : Any
            |kwargs|
        """
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)

        if crater.emplaceable:
            super().form_crater(crater, **kwargs)
        return

    def form_ejecta(self, crater: Crater | SimpleMoonCrater, **kwargs: Any) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        ejecta_thickness, ejecta_intensity = super().form_ejecta(crater, **kwargs)
        k_ej = self.ejecta_burial_degradation(ejecta_thickness[: crater.ejecta_region.n_face], ejecta_soften_factor=1.50)
        crater.ejecta_region.apply_diffusion(k_ej)

        k_deg = self.degradation_function(crater.diameter, fe=100) * ejecta_intensity[: crater.ejecta_region.n_face]
        crater.ejecta_region.apply_diffusion(k_deg)
        return ejecta_thickness, ejecta_intensity

    def crater_shape(self, crater: SimpleMoonCrater, region: LocalSurface, **kwargs: Any) -> NDArray[np.float64]:
        """
        Compute the crater shape based on the region view and surface.

        Parameters
        ----------
        crater : SimpleMoonCrater
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        reference_elevation = region.get_reference_surface(reference_radius=crater.final_radius)

        # Combine distances and references for nodes and faces
        distance = np.concatenate([region.face_distance, region.node_distance])

        original_elevation = np.concatenate([region.face_elevation, region.node_elevation])

        new_elevation = self.crater_profile(crater, distance, reference_elevation)
        elevation_change = new_elevation - original_elevation

        return elevation_change

    def crater_profile(self, crater: SimpleMoonCrater, r: ArrayLike, r_ref: ArrayLike | None = None) -> NDArray[np.float64]:
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
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        if r_ref is None:
            r_ref = np.zeros_like(r)

        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)

        # flatten r to 1D array
        rflat = np.ravel(r)
        r_ref_flat = np.ravel(r_ref)
        elevation = morphology_bindings.crater_profile(
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
        self, crater: SimpleMoonCrater, region: LocalSurface, **kwargs: Any
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the ejecta shape based on the region view and surface.

        Parameters
        ----------
        region : LocalSurface
            The region view of the surface mesh centered at the crat er center.
        crater : SimpleMoonCrater
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)

        distance = np.concatenate([region.face_distance, region.node_distance])
        if self.dorays:
            bearing = np.concatenate([region.face_bearing, region.node_bearing])
            thickness, intensity = self.ejecta_distribution(crater, distance, bearing)
        else:
            thickness = self.ejecta_profile(crater, distance)
            intensity = np.ones_like(thickness)

        return thickness, intensity

    def ejecta_profile(self, crater: SimpleMoonCrater, r: ArrayLike) -> NDArray[np.float64]:
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
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        if np.isscalar(r):
            r = np.array([r], dtype=np.float64)
        elif isinstance(r, (list | tuple)):
            r = np.array(r, dtype=np.float64)
        # flatten r to 1D array
        rflat = np.ravel(r)
        elevation = morphology_bindings.ejecta_profile(rflat, crater.final_diameter, crater.ejrim)
        elevation = np.array(elevation, dtype=np.float64)
        # reshape elevation to match the shape of r
        elevation = np.reshape(elevation, r.shape)
        return elevation

    def ejecta_distribution(self, crater, r: ArrayLike, theta: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        Compute the ejecta thickness distribution modulated by ray patterns.

        Parameters
        ----------
        crater : SimpleMoonCrater
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        # flatten r and theta to 1D arrays
        thickness = self.ejecta_profile(crater, r)
        rflat = np.ravel(r)
        theta_flat = np.ravel(theta)
        intensity = morphology_bindings.ray_intensity(
            radial_distance=rflat,
            initial_bearing=np.radians(theta_flat),
            crater_diameter=crater.final_diameter,
            seed=self.rng.integers(0, 2**32 - 1),
        )
        thickness = np.array(thickness, dtype=np.float64)
        intensity = np.array(intensity, dtype=np.float64)
        thickness *= intensity
        # reshape thickness to match the shape of r and theta
        return np.reshape(thickness, r.shape), np.reshape(intensity, r.shape)

    def ray_intensity(self, crater: SimpleMoonCrater, r: ArrayLike, theta: ArrayLike) -> NDArray[np.float64]:
        """
        Compute the ray pattern intensity modulation at each (r, theta) pair.

        Parameters
        ----------
        crater : SimpleMoonCrater
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)
        # flatten r and theta to 1D arrays
        rflat = np.ravel(r)
        theta_flat = np.radians(np.ravel(theta))
        intensity = morphology_bindings.ray_intensity(
            rflat,
            theta_flat,
            crater.final_diameter,
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
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021
        """
        return ejecta_soften_factor * ejecta_thickness**2

    def compute_subpixel_degradation(
        self,
        age_start: float,
        age_end: float,
        **kwargs,
    ) -> None:
        """
        Performs the subpixel degradation.

        This models the combined degradation of the part of the production population that is below the resolution of the mesh on each face. It is called between batches of craters by the `emplace` method.

        Parameters
        ----------
        age_start : float
            The age of the surface at the start of the degradation.
        age_end : float
            The age of the surface at the end of the degradation.
            |kwargs|
        """
        dc_min = 1e-8  # Minimum crater size for subpixel degradation calculation.

        if age_end >= age_start:
            raise ValueError("age_end must be less than age_start.")
        if self.production is None:
            raise RuntimeError("Production model must be set in the Morphology object to compute subpixel degradation.")

        if not hasattr(self, "_Kdiff"):
            self._Kdiff = np.zeros_like(self.surface.face_elevation)

        def _subpixel_degradation(diameter):
            fe = 100.0
            k = self.degradation_function(diameter, fe)
            n = self.production.function(
                diameter=diameter,
                age=age_start,
                age_end=age_end,
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
                self.counting.tally()

        return

    def apply_subpixel_degradation(self) -> None:
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
        if not isinstance(crater, SimpleMoonCrater):
            crater = SimpleMoonCrater.maker(crater, morphology=self)

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
        lower_limit = crater.final_radius * 1.0001
        upper_limit = self.ejecta_truncation * crater.final_radius if self.ejecta_truncation else np.pi * self.surface.target.radius

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
            ans = sol.root if sol.converged else crater.final_radius

        return float(ans)

    def degradation_function(
        self,
        final_diameter: FloatLike,
        fe: FloatLike = 100.0,
    ) -> float:
        """
        Computes the degradation function, which defines the topographic degradation that each crater contributes to the surface.

        This function implements a combination of the model by Minton et al. (2019) [#]_ for small craters and Riedel et al. (2020) [#]_ for large craters. It is currently not well-constrained, so may change in the future.

        Parameters
        ----------
        final_diameter : FloatLike
            The final diameter of the crater in meters.
        fe : FloatLike, optional
            The degradation function size factor, which is a scaling factor for the degradation function. Default is 100.0.

        Returns
        -------
        float
            The computed degradation function


        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021
        .. [#] Riedel, C., Minton, D.A., Michael, G., Orgel, C., Bogert, C.H. van der, Hiesinger, H., 2020. Degradation of Small Simple and Large Complex Lunar Craters: Not a Simple Scale Dependence. Journal of Geophysical Research: Planets 125, e2019JE006273. https://doi.org/10.1029/2019JE006273
        """

        def _kdmare(r, fe, psi):
            """
            The mare-scale degradation function from Minton et al. (2019). See eq. (32).
            """
            kv1 = 0.17
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

        return float(_kd(final_diameter / 2, fe))

    @staticmethod
    def overlap_function(crater: SimpleMoonCrater) -> tuple[set[int], set[int]]:
        """
        Get the affected node and face indices for the crater.

        Parameters
        ----------
        crater : SimpleMoonCrater
            The crater object to be used.

        Returns
        -------
        tuple[set[int], set[int]]
            The affected node and face indices for the crater.
        """
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
