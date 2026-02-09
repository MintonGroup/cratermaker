from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from math import pi
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad
from tqdm import tqdm

from cratermaker.components.crater import Crater, CraterFixed, CraterVariable
from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import format_large_units, parameter

if TYPE_CHECKING:
    from cratermaker.components.counting import Counting
    from cratermaker.components.production import Production
    from cratermaker.components.surface import LocalSurface, Surface


class MorphologyCraterVariable(CraterVariable):
    def __init__(self, morphology: Morphology | None = None, **kwargs: Any) -> None:
        object.__setattr__(self, "_morphology", morphology)
        object.__setattr__(self, "_ejecta_region", None)
        object.__setattr__(self, "_crater_region", None)
        object.__setattr__(self, "_face_index", None)
        object.__setattr__(self, "_affected_face_indices", None)
        object.__setattr__(self, "_affected_node_indices", None)
        object.__setattr__(self, "_crater_rmax", None)
        object.__setattr__(self, "_ejecta_rmax", None)
        object.__setattr__(self, "_emplaceable", None)
        super().__init__(**kwargs)
        return

    @property
    def _has_initialized_surface_data(self) -> bool:
        """
        A helper property to check if the morphology's surface data has been initialized, which is necessary for many of the variable properties to be computed.

        This is used to avoid unnecessary computations and potential errors when accessing properties that depend on surface data before it is available.
        """
        return self.morphology is not None and self.morphology.surface.uxds is not None

    @property
    def morphology(self) -> Morphology:
        """
        The morphology model associated with this crater variable object.
        """
        return self._morphology

    @property
    def emplaceable(self) -> bool | None:
        """
        Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted.
        """
        return self._emplaceable

    @property
    def face_index(self) -> int | None:
        return self._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        return self._ejecta_rmax

    @property
    def crater_rmax(self) -> float | None:
        return self._crater_rmax

    @property
    def affected_face_indices(self) -> set[int] | None:
        return self._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        return self._affected_node_indices

    @property
    def ejecta_region(self) -> LocalSurface | None:
        return self._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        return self._crater_region

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the crater variable properties.
        """
        dict_repr = super().as_dict()
        keys = (
            "face_index",
            "ejecta_rmax",
            "crater_rmax",
            "emplaceable",
            "affected_face_indices",
            "affected_node_indices",
            "morphology",
            "crater_region",
            "ejecta_region",
            "emplaceable",
        )
        for key in keys:
            dict_repr[key] = getattr(self, key)

        return dict_repr


class MorphologyCrater(Crater):
    def __init__(self, crater: Crater | None = None, fixed_cls=CraterFixed, var_cls=MorphologyCraterVariable, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, var_cls=var_cls, **kwargs)
        return

    def __str__(self) -> str:
        output = super().__str__()
        output += f"Crater region maximum radius: {format_large_units(self.crater_rmax, quantity='length')}\n"
        output += f"Ejecta region maximum radius: {format_large_units(self.ejecta_rmax, quantity='length')}\n"
        output += f"\nLarge enough to be emplaced on the grid: {self.emplaceable}\n"
        if self.emplaceable:
            output += f"Face index of crater center: {self.face_index}\n"
            output += f"Crater region: {self.crater_region}\n"
            output += f"Ejecta region: {self.ejecta_region}\n"
        return output

    def as_dict(self, ignore_keys: list[str] | tuple[str] = (), skip_complex_data: bool = False, **kwargs) -> dict:
        """
        Return a dictionary representation of the crater properties, including morphology-specific variable properties.

        Parameters
        ----------
        ignore_keys : list[str] or tuple[str], optional
            A list or tuple of property names to ignore when creating the dictionary representation. Default is an empty tuple unless `skip_complex_data` is True, in which case it will be extended to include ("morphology", "affected_face_indices", "affected_node_indices", "crater_region", "ejecta_region").
        skip_complex_data : bool, optional
            If True, skip complex data types when creating the dictionary representation. This is useful when serializing the object for saving to a file, as it removes complex data types that may not be serializable. Default is False.
        """
        if skip_complex_data:
            ignore_keys += ("morphology", "affected_face_indices", "affected_node_indices", "crater_region", "ejecta_region")
        dict_repr = super().as_dict(ignore_keys=ignore_keys, skip_complex_data=skip_complex_data, **kwargs)
        return dict_repr

    def remove_complex_data(self) -> None:
        """
        Remove complex data types from the crater variable properties.

        This is useful for storing a lightweight representation, as it removes complex data types that can be recomputed from the fixed properties and the morphology model when needed.
        """
        self._var._morphology = None
        self._var._affected_face_indices = None
        self._var._affected_node_indices = None
        self._var._crater_region = None
        self._var._ejecta_region = None
        return

    @property
    def face_index(self) -> int | None:
        """
        The index of the face on the surface mesh where the crater is centered.
        """
        if self._var._face_index is None and self._has_initialized_surface_data:
            self._var._face_index = self.morphology.surface.find_nearest_face(self.location)
        return self._var._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        """
        The maximum radius of the ejecta region for this crater, as determined by the morphology model and is based on the distance that the ejecta thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object.
        """
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
        if self._var._ejecta_region is None and self._has_initialized_surface_data and self.ejecta_rmax is not None:
            self._var._ejecta_region = self.morphology.surface.extract_region(
                location=self.location,
                region_radius=self.ejecta_rmax,
            )
            if self._var._ejecta_region is None:
                self._emplaceable = False
                self._affected_face_indices = set()
                self._affected_node_indices = set()
        return self._var._ejecta_region

    @property
    def crater_rmax(self) -> float | None:
        """
        The maximum radius of the crater region for this crater, as determined by the morphology model and is based on the distance that the crater shape thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object.
        """
        if self._var._crater_rmax is None and self._has_initialized_surface_data:
            self._var._crater_rmax = self.morphology.rmax(
                self, minimum_thickness=self.morphology.surface.smallest_length, feature="crater"
            )
        return self._var._crater_rmax

    @property
    def crater_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses the crater interior, as determined by the morphology model.

        This is extracted from the morphology's associated Surface object based on the location of the crater and the crater_rmax distance. If the crater region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if (
            self._var._crater_region is None
            and self._has_initialized_surface_data
            and self.ejecta_region is not None
            and self.crater_rmax is not None
        ):
            self._var._crater_region = self.ejecta_region.extract_subregion(
                subregion_radius=self.crater_rmax,
            )
            if self._var._crater_region is not None:
                self._emplaceable = True
                region = self.ejecta_region
                if isinstance(region.node_indices, slice) or isinstance(region.face_indices, slice):
                    self._var._affected_node_indices, self._var._affected_face_indices = (
                        set(np.arange(self.morphology.surface.n_node)[region.node_indices]),
                        set(np.arange(self.morphology.surface.n_face)[region.face_indices]),
                    )
                else:
                    self._var._affected_node_indices, self._var._affected_face_indices = (
                        set(region.node_indices),
                        set(region.face_indices),
                    )
            else:
                self._var._emplaceable = False
                self._var._affected_face_indices = set()
                self._var._affected_node_indices = set()
        return self._var._crater_region

    @property
    def affected_face_indices(self) -> set[int] | None:
        """
        The set of face indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_face_indices is None and self._has_initialized_surface_data:
            _ = self.crater_region
        return self._var._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        """
        The set of node indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_node_indices is None and self._has_initialized_surface_data:
            _ = self.crater_region
        return self._var._affected_node_indices

    @property
    def emplaceable(self) -> bool | None:
        """
        Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted.
        """
        if self._var._emplaceable is None and self._has_initialized_surface_data:
            _ = self.crater_region
        return self._var._emplaceable


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
        **kwargs : Any
            |kwargs|

        """
        from cratermaker.components.counting import Counting
        from cratermaker.components.production import Production
        from cratermaker.components.surface import Surface

        super().__init__(**kwargs)
        object.__setattr__(self, "_production", None)
        object.__setattr__(self, "_counting", None)
        object.__setattr__(self, "_do_counting", False)
        object.__setattr__(self, "_excavated_volume", None)

        self.surface = Surface.maker(surface, **kwargs)
        if counting is not None:
            self.counting = Counting.maker(counting, surface=self.surface, **kwargs)
        self.do_subpixel_degradation = do_subpixel_degradation
        self.do_slope_collapse = do_slope_collapse
        if do_subpixel_degradation:
            self.production = Production.maker(production, **kwargs)
        return

    def __str__(self) -> str:
        base = super().__str__()
        return base

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
            The name of the morphology model to use, or an instance of Morphology. If None, the default "simplemoon" is used.
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
            An instance of the specified Morphology model (e.g. SimpleMoon).

        Raises
        ------
        KeyError
            If the specified morphology model name is not found in the registry.
        TypeError
            If the specified morphology model is not a string or a subclass of Morphology.
        """
        # Call the base class version of make and pass the morphology argument as the component argument
        if morphology is None:
            morphology = "simplemoon"
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

    def emplace(self, craters: Crater | list[Crater], **kwargs: Any) -> None:
        """
        Convenience method to immediately emplace a crater onto the surface.

        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        crater : Crater or list[Crater]
            The crater to be emplaced.
        **kwargs : Any
            |kwargs|
        """
        if self._queue_manager is None:
            self._init_queue_manager()

        if isinstance(craters, list) and len(craters) > 0:
            for c in craters:
                self._enqueue_crater(c)
        elif isinstance(craters, Crater):
            self._enqueue_crater(craters)
        self._process_queue()

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

        crater_area = pi * crater.crater_rmax**2

        # Check to make sure that the face at the crater location is not smaller than the crater area
        if crater_area > self.surface.face_area[crater.face_index]:
            elevation_change = self.crater_shape(crater, crater.crater_region)
            crater.crater_region.update_elevation(elevation_change)
            if self.do_slope_collapse:
                crater.crater_region.slope_collapse()
            self._excavated_volume = crater.crater_region.compute_volume(elevation_change[: crater.crater_region.n_face])

            # Remove any ejecta from the surface
            inner_crater_region = crater.crater_region.extract_subregion(crater.final_radius)
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
                self.counting.add(crater, count_region=crater.ejecta_region, **kwargs)

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
        self._queue_manager = CraterQueueManager(self.overlap_function)

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

        self._queue_manager.push(crater)

    def _process_queue(self) -> None:
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
            tally_cadence = 10000
            nacumulated = 0
            while not self._queue_manager.is_empty():
                batch = self._queue_manager.peek_next_batch()

                with ThreadPoolExecutor() as executor:
                    # executor.map(process, batch)
                    futures = [executor.submit(self.form_crater, crater) for crater in batch]
                    for future in as_completed(futures):
                        try:
                            future.result()
                        except Exception as e:
                            raise RuntimeError("Error processing crater in batch\n") from e
                        else:
                            if pbar is not None:
                                pbar.update(1)

                if self.do_subpixel_degradation and len(batch) > 1:
                    # If the craters have age values attached to them, we can perform subpixel degradation between time values
                    agevals = [crater.age for crater in batch if crater.age is not None]
                    if len(agevals) > 1:
                        self.compute_subpixel_degradation(age_start=max(agevals), age_end=min(agevals))

                self._queue_manager.pop_batch(batch)
                nacumulated += len(batch)
                if self.do_counting and nacumulated >= tally_cadence:
                    self.counting.tally(quiet=False)
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
            self.apply_subpixel_degradation()
        if self.do_counting:
            self.counting.tally()
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
        """
        The surface object associated with this morphology model.
        """
        return self._surface

    @surface.setter
    def surface(self, surface: Surface) -> None:
        """
        Set the surface object associated with this morphology model.
        """
        from cratermaker.components.surface import Surface

        if isinstance(surface, Surface):
            self._surface = surface
        elif surface is None or isinstance(surface, str):
            self._surface = Surface.maker(surface)
        else:
            raise TypeError("surface must be an instance of Surface, a string, or None")
        self._queue_manager: CraterQueueManager | None = None

    @property
    def production(self) -> Production:
        """
        The production object associated with this morphology model.
        """
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
        """
        The counting object associated with this morphology model.
        """
        return self._counting

    @counting.setter
    def counting(self, counting: Counting) -> None:
        """
        Set the Counting object associated with this morphology model.
        """
        from cratermaker.components.counting import Counting

        if isinstance(counting, Counting):
            self._counting = counting
        elif isinstance(counting, str):
            self._counting = Counting.maker(counting, surface=self.surface)
        elif counting is None:
            self.do_counting = False
            self._counting = None
            return
        else:
            raise TypeError("counting must be an instance of Counting or a string")
        self.do_counting = True

    @parameter
    def do_subpixel_degradation(self) -> bool:
        """
        Whether to perform subpixel degradation during crater emplacement.
        """
        return self._do_subpixel_degradation

    @do_subpixel_degradation.setter
    def do_subpixel_degradation(self, value: bool) -> None:
        """
        Set whether to perform subpixel degradation during crater emplacement.
        """
        if not isinstance(value, bool):
            raise TypeError("do_subpixel_degradation must be a boolean value")
        self._do_subpixel_degradation = value

    @parameter
    def do_slope_collapse(self) -> bool:
        """
        Whether to perform slope collapse during crater emplacement.
        """
        return self._do_slope_collapse

    @do_slope_collapse.setter
    def do_slope_collapse(self, value: bool) -> None:
        """
        Set whether to perform slope collapse during crater emplacement.
        """
        if not isinstance(value, bool):
            raise TypeError("do_slope_collapse must be a boolean value")
        self._do_slope_collapse = value

    @parameter
    def do_counting(self) -> bool:
        """
        Whether to perform crater counting during crater emplacement.
        """
        return self._do_counting

    @do_counting.setter
    def do_counting(self, value: bool) -> None:
        """
        Set whether to perform crater counting during crater emplacement.
        """
        if not isinstance(value, bool):
            raise TypeError("do_counting must be a boolean value")
        self._do_counting = value


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
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._queue.append(crater)

    def peek_next_batch(self) -> list[Crater]:
        """
        Return a list of the next batch of craters that do not overlap with each other or the current active region.
        """
        batch = []
        reserved_nodes = set()  # set(self._active_nodes)
        reserved_faces = set()  # set(self._active_faces)
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
        """
        for crater in batch:
            self._queue.remove(crater)
            crater.remove_complex_data()
        return

    # def clear_active(self) -> None:
    #     """
    #     Clear the active region set after batch processing is complete.
    #     """
    #     self._active_nodes.clear()
    #     self._active_faces.clear()
    #     return

    def is_empty(self) -> bool:
        return len(self._queue) == 0


import_components(__name__, __path__)
