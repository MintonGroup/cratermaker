from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from cratermaker.bindings import counting_bindings
from cratermaker.components.crater import Crater, CraterFixed, CraterVariable
from cratermaker.constants import PairOfFloats
from cratermaker.utils.general_utils import format_large_units

# Relative difference below which the crater and ejecta region radii are treated as the same, so the crater region is reused as
# the ejecta region. The two radii come from separate root searches that agree to within a few units in the last place.
_SAME_REGION_RADIUS_RTOL = 16 * np.finfo(np.float64).eps

if TYPE_CHECKING:
    from cratermaker.components.morphology import Morphology
    from cratermaker.components.surface import LocalSurface


from abc import abstractmethod
from collections.abc import Callable
from typing import TYPE_CHECKING

from numpy.typing import ArrayLike, NDArray
from tqdm import tqdm

from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import parameter

if TYPE_CHECKING:
    from cratermaker.components.counting import Counting
    from cratermaker.components.production import Production
    from cratermaker.components.scaling import Scaling
    from cratermaker.components.surface import LocalSurface, Surface


class MorphologyCraterVariable(CraterVariable):
    def __init__(self, morphology: Morphology | None = None, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_morphology", morphology)
        object.__setattr__(self, "_ejecta_region", None)
        object.__setattr__(self, "_crater_region", None)
        object.__setattr__(self, "_face_index", None)
        object.__setattr__(self, "_affected_face_indices", None)
        object.__setattr__(self, "_affected_node_indices", None)
        object.__setattr__(self, "_ejecta_rmax", None)
        object.__setattr__(self, "_emplaceable", None)
        object.__setattr__(self, "_Crater", None)
        object.__setattr__(self, "_frac_ejrim", None)
        frac_ejrim = kwargs.pop("frac_ejrim", None)
        if frac_ejrim is not None:
            self.frac_ejrim = frac_ejrim
        return

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

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the crater variable properties.
        """
        dict_repr = super().as_dict()
        keys = (
            "face_index",
            "ejecta_rmax",
            "emplaceable",
            "affected_face_indices",
            "affected_node_indices",
            "morphology",
            "ejecta_region",
            "crater_region",
            "emplaceable",
            "frac_ejrim",
        )

        for key in keys:
            dict_repr[key] = getattr(self, key)

        return dict_repr

    @property
    def _has_initialized_surface_data(self) -> bool:
        """
        A helper property to check if the morphology's surface data has been initialized, which is necessary for many of the variable properties to be computed.

        This is used to avoid unnecessary computations and potential errors when accessing properties that depend on surface data before it is available.
        """
        return self.morphology is not None and self.morphology.surface.uxds is not None

    @property
    def morphology(self) -> Morphology:
        """The morphology model associated with this crater variable object."""
        return self._morphology

    @property
    def emplaceable(self) -> bool | None:
        """Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted."""
        return self._emplaceable

    @property
    def face_index(self) -> int | None:
        """The index of the face on the surface mesh where the crater is centered."""
        return self._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        """The maximum radius of the ejecta region for this crater, as determined by the morphology model and is based on the distance that the ejecta thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object."""
        return self._ejecta_rmax

    @property
    def affected_face_indices(self) -> NDArray[np.int64] | None:
        """The sorted, unique face indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_face_indices

    @property
    def affected_node_indices(self) -> NDArray[np.int64] | None:
        """The sorted, unique node indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_node_indices

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model."""
        return self._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses a buffered region around the cratered rim, as determined by the morphology model."""
        return self._crater_region


def _region_indices_array(region: LocalSurface | None, attr: str, n_total: int) -> NDArray[np.int64]:
    """
    Return the global face or node indices of a region as an integer array.

    Parameters
    ----------
    region : LocalSurface or None
        The region whose indices are returned. If None, an empty array is returned.
    attr : str
        Either "face_indices" or "node_indices".
    n_total : int
        The total number of faces or nodes on the surface, used to expand a slice into explicit indices.

    Returns
    -------
    NDArray[np.int64]
        The indices covered by the region.
    """
    if region is None:
        return np.empty(0, dtype=np.int64)
    indices = getattr(region, attr)
    if isinstance(indices, slice):
        return np.arange(n_total, dtype=np.int64)[indices]
    return np.asarray(indices)


@Crater.register("morphologycrater")
class MorphologyCrater(Crater):
    def __init__(self, crater: Crater | None = None, fixed_cls=CraterFixed, variable_cls=MorphologyCraterVariable, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)

    def __str__(self) -> str:
        horiz_line = "-" * 40 + "\n"
        str_repr = super().__str__()
        str_repr += f"Ejecta region maximum radius: {format_large_units(self.ejecta_rmax, quantity='length')}\n"
        str_repr += f"Large enough to be emplaced on the grid: {self.emplaceable}\n"
        if self.emplaceable:
            str_repr += f"Face index of crater center: {self.face_index}\n"
            str_repr += horiz_line
            str_repr += f"Crater region: {self.crater_region}"
            str_repr += horiz_line
            str_repr += f"Ejecta region: {self.ejecta_region}"
            str_repr += horiz_line
        return str_repr

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        location: PairOfFloats | None = None,
        relative_location: dict | None = None,
        check_redundant_inputs: bool = True,
        **kwargs: Any,
    ) -> MorphologyCrater:
        """
        Initialize a MorphologyCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a BasicMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        location : pair of floats, optional
            The (longitude, latitude) location of the crater.
        relative_location : dict
            Set the crater's location based on the distance and bearing from a reference location. The dictionary is required to have the keys: "distance" and "bearing" with values of distance (in meters) and bearing (in degrees) from the reference point. A third key, "reference_location" should be in the form of a tuple of floats representing the longitude and latitude pair that defines the reference location. The "reference_location" is optional if the morphology.surface contains a `local_location` attribute (such as in the case of a HiResLocal surface type and its derivatives) which will be used if the key is missing or set to None, otherwise it is required. Both location and relative_location cannot be set unless check_redundant_inputs is set to False.
        check_redundant_inputs : bool, optional
            If True, check for redundant inputs such as providing both diameter and radius. Default is True.
        kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`MorphologyCrater.maker() <cratermaker.components.morphology.MorphologyCrater.maker>`.  Refer to its documentation for a detailed description of valid keyword arguments.
        """
        from cratermaker.components.morphology import Morphology

        morphology = Morphology.maker(morphology, **kwargs)
        if relative_location is not None:
            if location is not None and check_redundant_inputs:
                raise ValueError("Cannot provide both location and relative_location")
            if not isinstance(relative_location, dict) or "distance" not in relative_location or "bearing" not in relative_location:
                raise TypeError("relative_location must be a dict with keys 'distance' and 'bearing'.")
            relative_location["distances"] = relative_location.pop("distance")
            relative_location["bearings"] = relative_location.pop("bearing")
            location = morphology.surface.compute_location_from_distance_bearing(**relative_location)

        crater = super().maker(
            crater=crater, scaling=morphology.scaling, location=location, check_redundant_inputs=check_redundant_inputs, **kwargs
        )
        # Crater.maker builds an instance of this class, so only the morphology needs to be attached rather than copying it
        if isinstance(crater, cls) and isinstance(crater._var, MorphologyCraterVariable):
            object.__setattr__(crater._var, "_morphology", morphology)
            return crater
        return cls(crater=crater, morphology=morphology)

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
            ignore_keys += (
                "morphology",
                "affected_face_indices",
                "affected_node_indices",
                "ejecta_region",
                "crater_region",
            )
        dict_repr = super().as_dict(ignore_keys=ignore_keys, skip_complex_data=skip_complex_data, **kwargs)
        return dict_repr

    def remove_complex_data(self) -> None:
        """
        Remove complex data types from the crater variable properties.

        This is useful for storing a lightweight representation, as it removes complex data types that can be recomputed from the fixed properties and the morphology model when needed. The morphology model is not removed, so that the complex data can be recomputed if needed.
        """
        self._var._affected_face_indices = None
        self._var._affected_node_indices = None
        self._var._ejecta_region = None
        self._var._crater_region = None
        return

    @property
    def face_index(self) -> int | None:
        """The index of the face on the surface mesh where the crater is centered."""
        if self._var._face_index is None and self._has_initialized_surface_data:
            self._var._face_index = self.morphology.surface.find_nearest_face(self.location)
        return self._var._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        """The maximum radius of the ejecta region for this crater, as determined by the morphology model and is based on the distance that the ejecta thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object."""
        if self._var._ejecta_rmax is None and self._has_initialized_surface_data:
            self._var._ejecta_rmax = self.morphology.rmax(
                self, minimum_thickness=self.morphology.surface.smallest_length, feature="ejecta"
            )
            if self._var._crater_region is not None:
                self._var._ejecta_rmax = max(self._var._ejecta_rmax, self._var._crater_region.radius)
        return self._var._ejecta_rmax

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model.

        This is extracted from the morphology's associated Surface object based on the location of the crater and the ejecta_rmax distance. If the ejecta region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if self._var._ejecta_region is None and self.crater_region is not None:
            if self.ejecta_rmax <= self.crater_region.radius * (1.0 + _SAME_REGION_RADIUS_RTOL):
                # The ejecta region would be extracted with the same location and (to within rounding) the same radius as the crater
                # region. This is common: ejecta_rmax is never smaller than the crater region radius, and the rim and ejecta
                # profiles usually cross the threshold at the same radius, with the two root searches differing only in the last
                # bit. Reuse the crater region rather than extracting an identical region and recomputing its distances, bearings,
                # and connectivity.
                self._var._ejecta_region = self.crater_region
            else:
                self._var._ejecta_region = self.morphology.surface.extract_region(
                    location=self.location,
                    region_radius=self.ejecta_rmax,
                )
        return self._var._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses a buffered region around the cratered rim.

        This is extracted from the morphology's associated Surface object based on the location of the crater and a radius that extends by a `_RIM_BUFFER_FACTOR` constant times the crater radius. If the crater region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if self._var._crater_region is None and self._has_initialized_surface_data:
            region_radius = self.morphology.rmax(self, minimum_thickness=self.morphology.surface.smallest_length, feature="crater")
            self._var._crater_region = self.morphology.surface.extract_region(location=self.location, region_radius=region_radius)
        return self._var._crater_region

    @property
    def affected_face_indices(self) -> NDArray[np.int64] | None:
        """
        The sorted, unique face indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously. The array is the ejecta region's own index array (not a copy), so it must not be modified.
        """
        if self._var._affected_face_indices is None and self._has_initialized_surface_data and self.emplaceable:
            self._var._affected_face_indices = _region_indices_array(
                self.ejecta_region, "face_indices", self.morphology.surface.n_face
            )
        return self._var._affected_face_indices

    @property
    def affected_node_indices(self) -> NDArray[np.int64] | None:
        """
        The sorted, unique node indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously. The array is the ejecta region's own index array (not a copy), so it must not be modified.
        """
        if self._var._affected_node_indices is None and self._has_initialized_surface_data and self.emplaceable:
            self._var._affected_node_indices = _region_indices_array(
                self.ejecta_region, "node_indices", self.morphology.surface.n_node
            )
        return self._var._affected_node_indices

    def reset(self):
        """
        Resets the crater by removing saved regions and lists of affected indices.
        """
        self.remove_complex_data()
        self._var._face_index = None
        self._var._ejecta_rmax = None
        self._var._emplaceable = None
        self._var._measured_rim_height = None
        self._var._measured_floor_elevation = None
        return

    @property
    def emplaceable(self) -> bool | None:
        """Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted."""
        if self._var._emplaceable is None and self._has_initialized_surface_data:
            self._var._emplaceable = self.crater_region is not None or self.ejecta_region is not None
        return self._var._emplaceable

    @property
    def measured_rim_height(self) -> float | None:
        """The measured rim height of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_rim_height = counting_bindings.measure_rim_height(self.crater_region, self)
        return self._var._measured_rim_height

    @property
    def measured_floor_elevation(self) -> float | None:
        """The measured floor depth of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_floor_elevation = counting_bindings.measure_floor_elevation(self.crater_region, self)
        return self._var._measured_floor_elevation

    @property
    def measured_depth_to_diameter(self) -> float | None:
        """
        The measured depth to diameter ratio of the crater.

        This is computed from `measured_rim_height`-`measured_floor_elevation`
        """
        if self.crater_region is not None:
            self.crater_region._desloped_face_elevation = None
            floor_elevation = self.measured_floor_elevation
            rim_height = self.measured_rim_height
            return (rim_height - floor_elevation) / self.measured_diameter
        else:
            return None


class Morphology(ComponentBase):
    _registry: dict[str, Morphology] = {}
    """
    The base class for Morphology models.
    """

    def __init__(
        self,
        surface: Surface | str | None = None,
        production: Production | str | None = None,
        counting: Counting | str | None = None,
        scaling: Scaling | str | None = None,
        conserve_volume: bool = True,
        do_subpixel_degradation: bool = True,
        do_slope_collapse: bool = True,
        save_observed: bool | None = None,
        save_emplaced: bool | None = None,
        **kwargs: Any,
    ) -> None:
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
        production : str, Production, optional
            The name of a Production object, or an instance of Production, to be associated with the morphology model. This is used for subpixel degradation in the emplace method. It is otherwise ignored if do_subpixel_degradation is False.
        counting : str, Counting or None, optional
            The name of a Counting object, or an instance of Counting, to be associated with the morphology model. This is used to record crater counts during emplacement. If None, no counting will be performed.
        scaling : str, Scaling or None, optional
            The projectile->crater size scaling model to use from the components library. The default is "montecarlo".
        conserve_volume : bool, optional
            When True, volume will be conserved when generating craters on a Surface using a combination of adjusting input morphometric parameters and scaling ejecta thickness such that the net change of volume of material is near zero.
        do_subpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Default is True.
        do_slope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Default is True.
        save_observed : bool or None, optional
            If True, counting will be performed during emplacement if a Counting object is provided. If False, no counting will be performed even if a Counting object is provided. If None (default), the presence or absence of a Counting object will determine whether counting is performed.
        save_emplaced : bool, optional
            If True, any craters large enough to be emplaced onto the surface in each interval are stored. If False, emplaced craters will not be tracked and saved, even if a Counting object is provided. If None (default), thye presence or abscence of a Counting objecdt will determine whether emplaced craters are saved.
        **kwargs : Any
            |kwargs|

        """
        from cratermaker.components.counting import Counting
        from cratermaker.components.production import Production
        from cratermaker.components.scaling import Scaling
        from cratermaker.components.surface import Surface

        super().__init__(**kwargs)
        object.__setattr__(self, "_production", None)
        object.__setattr__(self, "_counting", None)
        object.__setattr__(self, "_scaling", None)
        object.__setattr__(self, "_conserve_volume", None)
        object.__setattr__(self, "_save_observed", None)
        object.__setattr__(self, "_save_emplaced", None)
        object.__setattr__(self, "_Crater", None)

        # Because a Surface object is associated with a Counting object, we should first check to see if we are receiving a Counting object first so that we don't end up creating a spurious Surface object that we don't want.
        if surface is None and isinstance(counting, Counting):
            self.surface = counting.surface
            self.counting = counting
        else:
            self.surface = Surface.maker(surface, **kwargs)
            if counting is not None:
                self.counting = Counting.maker(counting, surface=self.surface, **kwargs)

        if self.counting is not None and self.counting.morphology is not self:
            self.counting.morphology = self  # Associated counting and morphology with each other
        else:
            self.counting = Counting.maker(counting, surface=self.surface, morphology=self, **kwargs)

        self.scaling = Scaling.maker(scaling, **kwargs)

        self.save_observed = save_observed
        self.save_emplaced = save_emplaced

        self.conserve_volume = conserve_volume
        self.do_subpixel_degradation = do_subpixel_degradation
        self.do_slope_collapse = do_slope_collapse
        if do_subpixel_degradation:
            # Be sure to use the associated Target object from the Surface object when initializing a new Production object
            if not isinstance(production, Production):
                kwargs["target"] = self.surface.target
            self.production = Production.maker(production, **kwargs)
        return

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += f"Count observed craters: {self.save_observed}\n"
        str_repr += f"Save emplaced craters: {self.save_emplaced}\n"
        if self.counting is not None:
            str_repr += f"Counting : <{self.counting.component_name}>\n"
        str_repr += f"Do slope collapse: {self.do_slope_collapse}\n"
        str_repr += f"Do subpixel degradation: {self.do_subpixel_degradation}\n"
        if self.do_subpixel_degradation:
            str_repr += f"Production: <{self.production.component_name}>\n"
        return str_repr

    @classmethod
    def maker(
        cls,
        morphology: str | type[Morphology] | Morphology | None = None,
        surface: Surface | str | None = None,
        production: Production | str | None = None,
        counting: Counting | str | None = None,
        scaling: Scaling | str | None = None,
        do_subpixel_degradation: bool = True,
        do_slope_collapse: bool = True,
        **kwargs: Any,
    ) -> Morphology:
        """
        Initialize a Morphology model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, the default "basicmoon" is used.
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
        production : str or Production, optional
            The name of a Production object, or an instance of Production, to be associated with the morphology model. This is used for subpixel degradation in the emplace method. It is otherwise ignored.
        counting : str or Counting, optional
            The name of a Counting object, or an instance of Counting, to be associated with the morphology model. This is used to record crater counts during emplacement. If None, no counting will be performed.
        scaling : str or Scaling or None, optional
            The projectile->crater size scaling model to use from the components library. The default is "montecarlo".
        do_subpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Default is True.
        do_slope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        Morphology
            An instance of the specified Morphology model (e.g. BasicMoonMorphology).

        Raises
        ------
        KeyError
            If the specified morphology model name is not found in the registry.
        TypeError
            If the specified morphology model is not a string or a subclass of Morphology.
        """
        # Call the base class version of make and pass the morphology argument as the component argument
        if morphology is None:
            morphology = "basicmoon"
        morphology = super().maker(
            component=morphology,
            surface=surface,
            production=production,
            counting=counting,
            scaling=scaling,
            do_subpixel_degradation=do_subpixel_degradation,
            do_slope_collapse=do_slope_collapse,
            **kwargs,
        )
        return morphology

    def emplace(self, craters: Crater | list[Crater] | None = None, **kwargs: Any) -> list[Crater]:
        """
        Convenience method to immediately emplace a crater onto the surface.

        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        craters : Crater or list of Crater objects, optional
            The Crater object(s) to be emplaced. If provided, this will be used directly. Otherwise, a single crater will be generated based on the keyword arguments.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        list of Crater or None
            The crater(s) that were emplaced. If no craters were provided and a crater could not be generated based on the keyword arguments, an empty list is returned.

        Notes
        -----
        The keyword arguments provided are passed down to :py:meth:`Crater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.

        Examples
        --------
        .. code-block:: python

            from cratermaker import Morphology, Crater

            morphology = Morphology.maker()

            # Create a crater with specific diameter
            morphology.emplace(diameter=10.0e3)

            # Create a crater based on a projectile with given mass and projectile_velocity
            morphology.emplace(projectile_mass=1e15, projectile_velocity=20e3)

            # Create a crater with a specific transient diameter and location
            morphology.emplace(transient_diameter=50e3, location=(43.43, -86.92))

            # Create multiple craters
            craters = [Crater.maker(diameter=20.0e3), Crater.maker(diameter=20.0e3)]
            morphology.emplace(craters)
        """
        if craters is None:
            craters = [self.Crater.maker(**kwargs)]
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
                if isinstance(c, self._CraterType):
                    processed_craters.append(c)
                else:
                    processed_craters.append(self.Crater.maker(crater=c, **kwargs))
            craters = processed_craters
        elif isinstance(craters, self._CraterType):
            craters = [craters]
        elif isinstance(craters, Crater):
            craters = [self.Crater.maker(crater=craters, **kwargs)]
        else:
            raise ValueError(
                "Invalid input for crater emplacement. Must be a Crater object, a list of Crater objects, or None with additional arguments for Crater.maker()."
            )

        if isinstance(craters, list) and len(craters) > 0:
            if self._queue_manager is None:
                self._init_queue_manager()

            for c in tqdm(
                craters,
                total=len(craters),
                desc="Queueing craters for emplacement",
                unit="crater",
                position=0,
                leave=False,
            ):
                self._enqueue_crater(c)

        self._process_queue(**kwargs)

        return craters

    @abstractmethod
    def form_crater(self, crater: Crater, **kwargs: Any) -> None:
        """
        Form a crater and its ejecta blanket on the surface.

        Implementations alter the face and node elevations and the ejecta thickness of the surface within the crater's region,
        and record the crater with the `counting` component if there is one. They must only change the faces and nodes of the
        region returned by `overlap_function`, which is what allows non-overlapping craters to be processed as a batch.

        Morphology models implement this as a wrapper around a compiled Rust function that runs the whole formation sequence
        (see `BasicMoonMorphology.form_crater`).

        Parameters
        ----------
        crater : Crater
            The crater object to be formed.
        **kwargs : Any
            |kwargs|
        """
        ...

    def form_craters(self, craters: list[Crater], **kwargs: Any) -> None:
        """
        Form a batch of craters whose regions do not overlap.

        The emplacement queue calls this once per batch of non-overlapping craters. This default forms them one at a time with
        `form_crater`. Models override it to form the craters of a batch in parallel (see `BasicMoonMorphology.form_craters`).

        Parameters
        ----------
        craters : list of Crater
            The craters to be formed, in emplacement order. Their regions (see `overlap_function`) must not share any faces or
            nodes.
        **kwargs : Any
            |kwargs|
        """
        for crater in craters:
            self.form_crater(crater, **kwargs)

    @staticmethod
    def _scale_ejecta_rim(crater: Crater, conservation_factor: float) -> None:
        """
        Scale the ejecta fraction of the rim height of a crater by the ejecta volume conservation factor.

        The fraction is capped at 1.

        Parameters
        ----------
        crater : Crater
            The crater whose ejecta was scaled.
        conservation_factor : float
            The factor the ejecta thickness was multiplied by.
        """
        ejrim = crater.frac_ejrim * crater.rim_height
        ejrim *= conservation_factor
        crater.frac_ejrim = min(ejrim / crater.rim_height, 1.0)

    def _init_queue_manager(self) -> None:
        """
        Initialize the crater queue manager with a surface-dependent overlap function.
        """
        self._queue_manager = self.CraterQueueManager(self.overlap_function)

    def _enqueue_crater(
        self,
        crater: Crater,
        **kwarg: Any,
    ) -> None:
        """
        Add a crater to the queue for later emplacement.

        Automatically initializes the queue manager if it hasn't been set.

        Parameters
        ----------
        crater : Crater
            The crater object to enqueue.
        **kwarg : Any
            |kwargs|

        Raises
        ------
        RuntimeError
            If the queue manager must be initialized but no surface is provided.
        """
        if self._queue_manager is None:
            if self.surface is None:
                raise RuntimeError("Surface must be provided to initialize queue manager.")
            self._init_queue_manager()
        if crater.emplaceable:
            self._queue_manager.push(crater)

    def _process_queue(self, **kwargs) -> None:
        """
        Process all queued craters in the order they were added, forming non-overlapping batches and applying each to the surface.

        Raises
        ------
        RuntimeError
            If the queue manager has not been initialized.
        """
        if not hasattr(self, "_queue_manager"):
            raise RuntimeError("Queue manager has not been initialized. Call _init_queue_manager first.")

        def _batch_process(pbar=None):
            tally_cadence = 2000
            nacumulated = 0
            while not self._queue_manager.is_empty:
                batch = self._queue_manager.peek_next_batch()

                self.form_craters(batch, **kwargs)
                if pbar is not None:
                    pbar.update(len(batch))

                if self.do_subpixel_degradation and len(batch) > 1:
                    # If the craters have time values attached to them, we can perform subpixel degradation between time values
                    timevals = [crater.time for crater in batch if crater.time is not None]
                    if len(timevals) > 1:
                        maxtime = max(timevals)
                        mintime = min(timevals)
                        if maxtime > mintime:
                            self.compute_subpixel_degradation(time_start=maxtime, time_end=mintime, **kwargs)

                self._queue_manager.pop_batch(batch)
                nacumulated += len(batch)
                if self.save_observed and nacumulated >= tally_cadence:
                    measure_rim = kwargs.pop("measure_rim", False)
                    self.counting.tally(measure_rim=measure_rim, quiet=False, **kwargs)
                    nacumulated = 0

            return

        total_craters = len(self._queue_manager._queue)
        if total_craters > 10:
            with tqdm(
                total=total_craters,
                desc="Emplacing craters",
                position=0,
                leave=False,
                unit="craters",
                smoothing=10 / total_craters,
            ) as pbar:
                _batch_process(pbar)
        else:
            _batch_process()

        if self.do_subpixel_degradation:
            self.apply_subpixel_degradation(**kwargs)
        if self.save_observed:
            self.counting.tally(**kwargs)
        return

    @abstractmethod
    def rmax(
        self,
        crater: Crater,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
    ) -> float: ...

    @staticmethod
    @abstractmethod
    def overlap_function(
        crater: Crater,
    ) -> tuple[ArrayLike, ArrayLike]: ...

    @property
    def surface(self) -> Surface:
        """The Surface object associated with this morphology model."""
        return self._surface

    @surface.setter
    def surface(self, surface: Surface) -> None:
        """Set the Surface object associated with this morphology model."""
        from cratermaker.components.surface import Surface

        if isinstance(surface, Surface):
            self._surface = surface
        elif surface is None or isinstance(surface, str):
            self._surface = Surface.maker(surface)
        else:
            raise TypeError("surface must be an instance of Surface, a string, or None")
        self._queue_manager: self.CraterQueueManager | None = None

    @property
    def production(self) -> Production:
        """The Production object associated with this morphology model."""
        return self._production

    @production.setter
    def production(self, production: Production) -> None:
        from cratermaker.components.production import Production

        if isinstance(production, Production):
            self._production = production
        elif production is None or isinstance(production, str):
            self._production = Production.maker(production)
        else:
            raise TypeError("production must be an instance of Production a string or None")

    @property
    def counting(self) -> Counting:
        """The Counting object associated with this morphology model."""
        return self._counting

    @counting.setter
    def counting(self, counting: Counting) -> None:
        """Set the Counting object associated with this morphology model."""
        from cratermaker.components.counting import Counting

        if isinstance(counting, Counting):
            self._counting = counting
        elif isinstance(counting, str):
            self._counting = Counting.maker(counting, surface=self.surface, Crater=self.Crater)
        elif counting is None:
            self.save_observed = False
            self._counting = None
            return
        else:
            raise TypeError("counting must be an instance of Counting or a string")
        self.save_observed = True
        self._counting.morphology = self

    @parameter
    def conserve_volume(self) -> bool:
        """Whether to conserved volume during crater emplacement."""
        return self._conserve_volume

    @conserve_volume.setter
    def conserve_volume(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("conserve_volume must be a boolean value")
        self._conserve_volume = value

    @parameter
    def do_subpixel_degradation(self) -> bool:
        """Whether to perform subpixel degradation during crater emplacement."""
        return self._do_subpixel_degradation

    @do_subpixel_degradation.setter
    def do_subpixel_degradation(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("do_subpixel_degradation must be a boolean value")
        self._do_subpixel_degradation = value

    @parameter
    def do_slope_collapse(self) -> bool:
        """Whether to perform slope collapse during crater emplacement."""
        return self._do_slope_collapse

    @do_slope_collapse.setter
    def do_slope_collapse(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("do_slope_collapse must be a boolean value")
        self._do_slope_collapse = value

    @parameter
    def save_observed(self) -> bool:
        """Whether to perform crater counting during crater emplacement."""
        if self._save_observed is None:
            return self.counting is not None
        else:
            return self._save_observed

    @save_observed.setter
    def save_observed(self, value: bool | None) -> None:
        if not isinstance(value, (bool | None)):
            raise TypeError("save_observed must be a boolean value or None")
        if value and self.counting is None:
            raise ValueError("save_observed cannot be True with no counting object was provided or created.")
        self._save_observed = value

    @parameter
    def save_emplaced(self) -> bool:
        """Whether to save a record of emplaced craters."""
        if self._save_emplaced is None:
            return self.counting is not None
        else:
            return self._save_emplaced

    @save_emplaced.setter
    def save_emplaced(self, value: bool) -> None:
        if not isinstance(value, (bool | None)):
            raise TypeError("save_emplaced must be a boolean value or None")
        if value and self.counting is None:
            raise ValueError("save_emplaced cannot be True with no counting object was provided or created.")
        self._save_emplaced = value

    @property
    def scaling(self):
        """
        The Scaling object that defines the crater scaling relationships model. Set during initialization.
        """
        return self._scaling

    @scaling.setter
    def scaling(self, value):
        from cratermaker.components.scaling import Scaling

        if not isinstance(value, Scaling):
            raise TypeError("scaling must be of Scaling")
        self._scaling = value

    @property
    def _CraterType(self) -> type[MorphologyCrater]:
        """The Crater class associated with this morphology model."""
        return MorphologyCrater

    @property
    def Crater(self) -> type[MorphologyCrater]:
        if self._Crater is None:

            class _WrappedMorphologyCrater(self._CraterType):
                @classmethod
                def maker(cls: type[_WrappedMorphologyCrater], **kwargs):
                    kwargs["morphology"] = self
                    if "conserve_volume" not in kwargs:
                        kwargs["conserve_volume"] = self.conserve_volume
                    kwargs = {**kwargs, **vars(self.common_args)}
                    return self._CraterType.maker(**kwargs)

            self._Crater = _WrappedMorphologyCrater

        return self._Crater

    class CraterQueueManager:
        """
        A manager for craters awaiting emplacement. Craters are processed in order (FIFO) but batches of non-overlapping craters can be processed simultaneously.

        Parameters
        ----------
        overlap_fn : Callable[[Crater], tuple[ArrayLike, ArrayLike]]
            A function that takes a crater and returns a tuple of (node indices, face indices) affected. Each may be an integer
            array of global indices or a set of ints.
        """

        def __init__(self, overlap_fn: Callable[[Crater], tuple[ArrayLike, ArrayLike]]):
            self._queue: list[Crater] = []
            self._overlap_fn = overlap_fn
            # Reusable reservation masks over the global node and face indices. They are grown as needed and are always left
            # all-False between calls to peek_next_batch.
            self._reserved_nodes = np.zeros(0, dtype=bool)
            self._reserved_faces = np.zeros(0, dtype=bool)

        @staticmethod
        def _as_index_array(indices: ArrayLike) -> NDArray[np.intp]:
            """Convert the output of an overlap function (array or set of ints) to an integer index array."""
            if isinstance(indices, (set, frozenset)):
                return np.fromiter(indices, dtype=np.intp, count=len(indices))
            return np.asarray(indices, dtype=np.intp)

        @staticmethod
        def _fit_mask(mask: NDArray[np.bool_], indices: NDArray[np.intp]) -> NDArray[np.bool_]:
            """Return `mask`, grown with False entries if needed so that every index in `indices` is valid."""
            if indices.size:
                needed = int(indices.max()) + 1
                if needed > mask.size:
                    grown = np.zeros(max(needed, 2 * mask.size), dtype=bool)
                    grown[: mask.size] = mask
                    return grown
            return mask

        def push(self, crater: Crater) -> None:
            """
            Add a crater to the end of the queue.

            Parameters
            ----------
            crater : Crater
                The crater to be added to the queue.

            Raises
            ------
            TypeError
                If the provided crater is not an instance of Crater.
            """
            if not isinstance(crater, Crater):
                raise TypeError("crater must be an instance of Crater")
            self._queue.append(crater)
            return

        def peek_next_batch(self) -> list[Crater]:
            """
            Return a list of the next batch of craters that do not overlap with each other or the current active region.

            Returns
            -------
            list of Crater
                The next batch of craters that can be processed simultaneously without overlap.
            """
            batch = []
            touched_nodes = []
            touched_faces = []
            try:
                for crater in self._queue:
                    node_indices, face_indices = self._overlap_fn(crater)
                    node_indices = self._as_index_array(node_indices)
                    face_indices = self._as_index_array(face_indices)
                    self._reserved_nodes = self._fit_mask(self._reserved_nodes, node_indices)
                    self._reserved_faces = self._fit_mask(self._reserved_faces, face_indices)
                    if self._reserved_nodes[node_indices].any() or self._reserved_faces[face_indices].any():
                        break
                    batch.append(crater)
                    self._reserved_nodes[node_indices] = True
                    self._reserved_faces[face_indices] = True
                    touched_nodes.append(node_indices)
                    touched_faces.append(face_indices)
            finally:
                # Clear only the entries that were set, so the masks are all-False for the next call
                for idx in touched_nodes:
                    self._reserved_nodes[idx] = False
                for idx in touched_faces:
                    self._reserved_faces[idx] = False
            return batch

        def pop_batch(self, batch: list[Crater]) -> None:
            """
            Remove a processed batch of craters from the queue.

            Parameters
            ----------
            batch : list of Crater
                The batch of craters to be removed from the queue.

            Raises
            ------
            ValueError
                If any crater in the batch is not found in the queue.
            """
            for crater in batch:
                self._queue.remove(crater)
                crater.remove_complex_data()
            return

        @property
        def is_empty(self) -> bool:
            """Return True if the queue is empty, False otherwise."""
            return len(self._queue) == 0


import_components(__name__, __path__)
