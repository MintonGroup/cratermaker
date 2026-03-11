from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from cratermaker._cratermaker import counting_bindings

from cratermaker.components.crater import Crater, CraterFixed, CraterVariable
from cratermaker.constants import PairOfFloats
from cratermaker.utils.general_utils import format_large_units

if TYPE_CHECKING:
    from cratermaker.components.morphology import Morphology
    from cratermaker.components.surface import LocalSurface

# The factor by which the crater tagging region is extended beyond the final rim.
_RIM_BUFFER_FACTOR = 1.5


from abc import abstractmethod
from collections.abc import Callable
from math import pi
from typing import TYPE_CHECKING

from numpy.typing import NDArray
from tqdm import tqdm

from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import parameter

if TYPE_CHECKING:
    from cratermaker.components.counting import Counting
    from cratermaker.components.production import Production
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
    def affected_face_indices(self) -> set[int] | None:
        """The set of face indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        """The set of node indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_node_indices

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model."""
        return self._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses a buffered region around the cratered rim, as determined by the morphology model."""
        return self._crater_region


class MorphologyCrater(Crater):
    def __init__(self, crater: Crater | None = None, fixed_cls=CraterFixed, variable_cls=MorphologyCraterVariable, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, variable_cls=variable_cls, **kwargs)
        return

    def __str__(self) -> str:
        output = super().__str__()
        output += f"Ejecta region maximum radius: {format_large_units(self.ejecta_rmax, quantity='length')}\n"
        output += f"\nLarge enough to be emplaced on the grid: {self.emplaceable}\n"
        if self.emplaceable:
            output += f"Face index of crater center: {self.face_index}\n"
            output += f"Crater region: {self.crater_region}\n"
            output += f"Ejecta region: {self.ejecta_region}\n"
        return output

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
            location = morphology.surface.compute_location_from_distance_bearing(**relative_location)

        if crater is None:
            crater = super().maker(location=location, check_redundant_inputs=check_redundant_inputs, **kwargs)
        args = {}

        return cls(
            crater=crater,
            morphology=morphology,
            **args,
        )

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
        return self._var._ejecta_rmax

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model.

        This is extracted from the morphology's associated Surface object based on the location of the crater and the ejecta_rmax distance. If the ejecta region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if self._var._ejecta_region is None and self.crater_region is not None:
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
            self._var._crater_region = self.morphology.surface.extract_region(
                location=self.location,
                region_radius=_RIM_BUFFER_FACTOR * self.radius,
            )
        return self._var._crater_region

    @property
    def affected_face_indices(self) -> set[int] | None:
        """
        The set of face indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_face_indices is None and self._has_initialized_surface_data and self.emplaceable:
            if self.ejecta_region is not None:
                if isinstance(self.ejecta_region.face_indices, slice):
                    self._var._affected_face_indices = set(
                        np.arange(self.morphology.surface.n_face)[self.ejecta_region.face_indices]
                    )
                else:
                    self._var._affected_face_indices = set(self.ejecta_region.face_indices)
            else:
                self._var._affected_face_indices = set()

        return self._var._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        """
        The set of node indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_node_indices is None and self._has_initialized_surface_data and self.emplaceable:
            if self.ejecta_region is not None:
                if isinstance(self.ejecta_region.node_indices, slice):
                    self._var._affected_node_indices = set(
                        np.arange(self.morphology.surface.n_node)[self.ejecta_region.node_indices]
                    )
                else:
                    self._var._affected_node_indices = set(self.ejecta_region.node_indices)
            else:
                self._var._affected_node_indices = set()

        return self._var._affected_node_indices

    @property
    def emplaceable(self) -> bool | None:
        """Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted."""
        if self._var._emplaceable is None and self._has_initialized_surface_data:
            self._var._emplaceable = self.crater_region is not None
        return self._var._emplaceable

    @property
    def measured_rim_height(self) -> float | None:
        """The measured rim height of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_rim_height = counting_bindings.measure_rim_height(self.crater_region, self)
        return self._var._measured_rim_height

    @property
    def measured_floor_depth(self) -> float | None:
        """The measured floor depth of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_floor_depth = counting_bindings.measure_floor_depth(self.crater_region, self)
        return self._var._measured_floor_depth

    @property
    def measured_depth_to_diameter(self) -> float | None:
        """
        The measured depth to diameter ratio of the crater.

        This is computed from `measured_rim_height`-`measured_floor_depth`
        """
        if self.crater_region is not None:
            self.crater_region._desloped_face_elevation = None
            floor_depth = self.measured_floor_depth
            rim_height = self.measured_rim_height
            return (rim_height - floor_depth) / self.measured_diameter
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
        do_subpixel_degradation: bool = True,
        do_slope_collapse: bool = True,
        do_counting: bool | None = None,
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
        do_subpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Default is True.
        do_slope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Default is True.
        do_counting : bool or None, optional
            If True, counting will be performed during emplacement if a Counting object is provided. If False, no counting will be performed even if a Counting object is provided. If None (default), the presence or absence of a Counting object will determine whether counting is performed.
        **kwargs : Any
            |kwargs|

        """
        from cratermaker.components.counting import Counting
        from cratermaker.components.production import Production
        from cratermaker.components.surface import Surface

        super().__init__(**kwargs)
        object.__setattr__(self, "_production", None)
        object.__setattr__(self, "_counting", None)
        object.__setattr__(self, "_do_counting", None)
        object.__setattr__(self, "_excavated_volume", None)
        object.__setattr__(self, "_Crater", None)
        object.__setattr__(self, "_CraterType", MorphologyCrater)

        # Because a Surface object is associated with a Counting object, we should first check to see if we are receiving a Counting object first so that we don't end up creating a spurious Surface object that we don't want.
        if surface is None and isinstance(counting, Counting):
            self.surface = counting.surface
            self.counting = counting
        else:
            self.surface = Surface.maker(surface, **kwargs)
            if counting is not None:
                crater_cls = kwargs.pop("Crater", self.Crater)
                self.counting = Counting.maker(counting, surface=self.surface, Crater=crater_cls, **kwargs)

        if self.counting is not None and self.counting.morphology is not self:
            self.counting.morphology = self  # Associated counting and morphology with each other

        if do_counting is not None:
            if do_counting and self.counting is None:
                raise ValueError("do_counting is True but no counting object was provided or created.")
            else:
                self.do_counting = do_counting

        self.do_subpixel_degradation = do_subpixel_degradation
        self.do_slope_collapse = do_slope_collapse
        if do_subpixel_degradation:
            # Be sure to use the associated Target object from the Surface object when initializing a new Production object
            if not isinstance(production, Production):
                kwargs["target"] = self.surface.target
            self.production = Production.maker(production, **kwargs)
        return

    def __str__(self) -> str:
        base = super().__str__()
        str_repr = f"{base}\n"
        str_repr += f"\n{self.counting}\n"
        str_repr += f"Do slope collapse: {self.do_slope_collapse}\n"
        str_repr += f"Do subpixel degradation: {self.do_subpixel_degradation}\n"
        if self.do_subpixel_degradation:
            str_repr += f"\n{self.production}\n"
        return str_repr

    @classmethod
    def maker(
        cls,
        morphology: str | type[Morphology] | Morphology | None = None,
        surface: Surface | str | None = None,
        production: Production | str | None = None,
        counting: Counting | str | None = None,
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
        do_subpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Default is True.
        do_slope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Default is True.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        Morphology
            An instance of the specified Morphology model (e.g. BasicMoon).

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
        if self._queue_manager is None:
            self._init_queue_manager()

        if craters is None:
            craters = [MorphologyCrater.maker(morphology=self, **kwargs)]
        elif isinstance(craters, MorphologyCrater):
            craters = [craters]
        elif isinstance(craters, Crater):
            craters = [MorphologyCrater.maker(crater=craters, morphology=self)]

        if isinstance(craters, list) and len(craters) > 0:
            for c in craters:
                self._enqueue_crater(c)

        self._process_queue(**kwargs)

        return craters

    def form_crater(self, crater: Crater, **kwargs: Any) -> None:
        """
        Form the interior of the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        crater : Crater
            The crater object to be formed.
        **kwargs : Any
            |kwargs|
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if not isinstance(crater, MorphologyCrater):
            crater = MorphologyCrater.maker(crater=crater, morphology=self)
        if not crater.emplaceable:
            return

        crater_area = pi * crater.crater_region.radius**2

        # Check to make sure that the face at the crater location is not smaller than the crater area
        if crater_area > self.surface.face_area[crater.face_index]:
            elevation_change = self.crater_shape(crater, crater.ejecta_region)
            crater.ejecta_region.update_elevation(elevation_change)
            if self.do_slope_collapse:
                crater.ejecta_region.slope_collapse()
            self._excavated_volume = crater.ejecta_region.compute_volume(elevation_change[: crater.ejecta_region.n_face])

            # Remove any ejecta from the interior of the crater
            inner_crater_region = crater.crater_region.extract_subregion(crater.radius)
            if inner_crater_region is not None:
                inner_crater_region.add_data(
                    "ejecta_thickness",
                    long_name="ejecta thickness",
                    units="m",
                    data=0.0,
                    overwrite=True,
                )

            # Record the crater to the counting layer
            if self.do_counting:
                self.counting.add(crater, **kwargs)

            self.form_ejecta(crater, **kwargs)
        return

    def form_ejecta(self, crater: Crater, **kwargs: Any) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
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
            The computed ejecta thickness and intensity at the face and node elevations. |
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if not isinstance(crater, MorphologyCrater):
            crater = MorphologyCrater.maker(crater=crater, morphology=self)

        if not self._excavated_volume:
            return
        ejecta_thickness, ejecta_intensity = self.ejecta_shape(crater, crater.ejecta_region)
        ejecta_volume = crater.ejecta_region.compute_volume(ejecta_thickness[: crater.ejecta_region.n_face])
        conservation_factor = -self._excavated_volume / ejecta_volume
        ejecta_thickness *= conservation_factor

        crater.ejecta_region.add_data(
            "ejecta_thickness",
            long_name="ejecta thickness",
            units="m",
            data=ejecta_thickness[: crater.ejecta_region.n_face],
        )

        crater.ejecta_region.update_elevation(ejecta_thickness)

        return ejecta_thickness, ejecta_intensity

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

        from concurrent.futures import ThreadPoolExecutor, as_completed

        def _batch_process(pbar=None):
            tally_cadence = 2000
            nacumulated = 0
            while not self._queue_manager.is_empty:
                batch = self._queue_manager.peek_next_batch()

                for crater in batch:
                    self.form_crater(crater)
                    if pbar is not None:
                        pbar.update(1)

                if self.do_subpixel_degradation and len(batch) > 1:
                    # If the craters have time values attached to them, we can perform subpixel degradation between time values
                    timevals = [crater.time for crater in batch if crater.time is not None]
                    if len(timevals) > 1:
                        self.compute_subpixel_degradation(time_start=max(timevals), time_end=min(timevals), **kwargs)

                self._queue_manager.pop_batch(batch)
                nacumulated += len(batch)
                if self.do_counting and nacumulated >= tally_cadence:
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
        if self.do_counting:
            self.counting.tally(**kwargs)
        return

    @abstractmethod
    def crater_shape(
        self,
        crater: Crater,
        region: LocalSurface,
        **kwarg: Any,
    ) -> NDArray[np.float64]: ...

    @abstractmethod
    def ejecta_shape(
        self,
        crater: Crater,
        region: LocalSurface,
        **kwarg: Any,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...

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
    ) -> tuple[set[int], set[int]]: ...

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
            self.do_counting = False
            self._counting = None
            return
        else:
            raise TypeError("counting must be an instance of Counting or a string")
        self.do_counting = True
        self._counting.morphology = self

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
    def do_counting(self) -> bool:
        """Whether to perform crater counting during crater emplacement."""
        return self._do_counting

    @do_counting.setter
    def do_counting(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("do_counting must be a boolean value")
        self._do_counting = value

    @property
    def Crater(self) -> type[MorphologyCrater]:
        if self._Crater is None:

            class _WrappedMorphologyCrater(self._CraterType):
                @classmethod
                def maker(cls: type[_WrappedMorphologyCrater], **kwargs):
                    kwargs["morphology"] = self
                    return self._CraterType.maker(**kwargs)

            self._Crater = _WrappedMorphologyCrater

        return self._Crater

    class CraterQueueManager:
        """
        A manager for craters awaiting emplacement. Craters are processed in order (FIFO) but batches of non-overlapping craters can be processed simultaneously.

        Parameters
        ----------
        overlap_fn : Callable[[Crater], tuple[set[int], set[int]]]
            A function that takes a crater and returns a tuple of (node indices, face indices) affected.
        """

        def __init__(self, overlap_fn: Callable[[Crater], tuple[set[int], set[int]]]):
            self._queue: list[Crater] = []
            self._overlap_fn = overlap_fn

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
            reserved_nodes = set()
            reserved_faces = set()
            for crater in self._queue:
                node_indices, face_indices = self._overlap_fn(crater)
                if reserved_nodes.isdisjoint(node_indices) and reserved_faces.isdisjoint(face_indices):
                    batch.append(crater)
                    reserved_nodes.update(node_indices)
                    reserved_faces.update(face_indices)
                else:
                    break
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
