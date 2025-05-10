from __future__ import annotations

from abc import abstractmethod
from math import pi
from typing import TYPE_CHECKING, Any, Callable

from numpy.typing import NDArray
from tqdm import tqdm

from cratermaker.constants import FloatLike
from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import Surface, SurfaceView


class Morphology(ComponentBase):
    def __init__(self, **kwargs: Any) -> None:
        """
        Initialize the Morphology class.

        Parameters
        ----------
        **kwargs : Any
            Additional keyword arguments.

        Raises
        -------

        TypeError
            If the crater is not an instance of Crater.

        """
        super().__init__(**kwargs)
        self._queue_manager: CraterQueueManager | None = None

    def __str__(self) -> str:
        base = super().__str__()
        return base

    @classmethod
    def maker(
        cls,
        morphology: str | type[Morphology] | Morphology | None = None,
        **kwargs: Any,
    ) -> Morphology:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, the default "simplemoon" is used.
        **kwargs : Any
            Additional keyword arguments that are required for the specific morphology model being created.

        Returns
        -------
        component
            An instance of the specified component model.

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
        morphology = super().maker(component=morphology, **kwargs)
        return morphology

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
        **kwargs : Any
            Additional keyword arguments.
        """
        if self._queue_manager is None:
            self.init_queue_manager(surface)
        self.enqueue_crater(crater)
        self.process_queue(surface)

    def form_crater(self, crater: Crater, surface: Surface, **kwargs: Any) -> Surface:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        crater : Crater
            The crater object to be formed.
        surface : Surface
            The surface to be altered. This will be modified in place.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions.

        Returns
        -------
        None
        """
        from cratermaker.components.surface import Surface

        if not isinstance(surface, Surface):
            raise TypeError("surface must be an instance of Surface")
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        # Find the node and face center of the crater
        self.node_index, self.face_index = surface.find_nearest_index(crater.location)

        # Test if the ejecta is big enough to modify the surface
        ejecta_rmax = self.rmax(
            crater, minimum_thickness=surface.smallest_length, feature="ejecta"
        )
        ejecta_region_view = surface.extract_region(crater.location, ejecta_rmax)
        ejecta_area = pi * ejecta_rmax**2
        if (
            ejecta_region_view is None
            or ejecta_area < surface.face_areas[self.face_index]
        ):  # The crater is too small to change the surface
            return

        crater_rmax = self.rmax(
            crater, minimum_thickness=surface.smallest_length, feature="crater"
        )
        crater_region_view = surface.extract_region(crater.location, crater_rmax)
        if (
            crater_region_view is not None
        ):  # The crater is big enough to affect the surface
            crater_area = pi * crater_rmax**2

            # Check to make sure that the face at the crater location is not smaller than the crater area
            if crater_area > surface.face_areas[self.face_index]:
                # Form the crater shape
                surface = self.crater_shape(crater, crater_region_view, surface)

        # Now form the ejecta blanket
        surface = self.ejecta_shape(crater, ejecta_region_view, surface)
        return

    def affected_indices(
        self, crater: Crater, surface: Surface
    ) -> tuple[set[int], set[int]]:
        """
        Determine the set of node and face indices affected by a crater's ejecta blanket.

        Parameters
        ----------
        crater : Crater
            The crater object whose effect region is to be computed.
        surface : Surface
            The surface object on which the crater would be emplaced.

        Returns
        -------
        affected_nodes : set of int
            The set of node indices affected by the crater.
        affected_faces : set of int
            The set of face indices affected by the crater.
        """
        rmax = self.rmax(
            crater, minimum_thickness=surface.smallest_length, feature="ejecta"
        )
        region_view = surface.extract_region(crater.location, rmax)
        if region_view is None:
            return set(), set()
        return set(region_view.node_indices), set(region_view.face_indices)

    def init_queue_manager(self, surface: Surface) -> None:
        """
        Initialize the crater queue manager with a surface-dependent overlap function.

        Parameters
        ----------
        surface : Surface
            The surface object used to determine which nodes and faces each crater affects.
        """

        def overlap_fn(crater: Crater) -> tuple[set[int], set[int]]:
            return self.affected_indices(crater, surface)

        self._queue_manager = CraterQueueManager(overlap_fn)

    def enqueue_crater(
        self,
        crater: Crater | None = None,
        surface: "Surface" | None = None,
        **kwarg: Any,
    ) -> None:
        """
        Add a crater to the queue for later emplacement. Automatically initializes
        the queue manager if it hasn't been set.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to enqueue. If None, one is created from keyword args.
        surface : Surface, optional
            The surface used to determine overlap regions if initialization is needed.
        **kwarg : Any
            Additional keyword arguments for crater construction.

        Raises
        ------
        RuntimeError
            If the queue manager must be initialized but no surface is provided.
        """
        if self._queue_manager is None:
            if surface is None:
                raise RuntimeError(
                    "Surface must be provided to initialize queue manager."
                )
            self.init_queue_manager(surface)

        if crater is None:
            crater = Crater.maker(**kwarg)
        self._queue_manager.push(crater)

    def process_queue(self, surface: Surface) -> None:
        """
        Process all queued craters in the order they were added, forming non-overlapping
        batches and applying each to the surface.

        Parameters
        ----------
        surface : Surface
            The surface object to be modified in place.

        Raises
        ------
        RuntimeError
            If the queue manager has not been initialized.
        """
        if not hasattr(self, "_queue_manager"):
            raise RuntimeError(
                "Queue manager has not been initialized. Call init_queue_manager first."
            )

        def _batch_process(pbar=None):
            while not self._queue_manager.is_empty():
                batch = self._queue_manager.peek_next_batch()
                for crater in batch:
                    self.form_crater(crater, surface)
                    if pbar is not None:
                        pbar.update(1)
                self._queue_manager.pop_batch(batch)
                self._queue_manager.clear_active()
            return

        total_craters = len(self._queue_manager._queue)
        if total_craters > 10:
            with tqdm(total=total_craters, desc="Processing craters") as pbar:
                _batch_process(pbar)
        else:
            _batch_process()
        return

    @abstractmethod
    def crater_shape(
        self,
        crater: Crater,
        region_view: SurfaceView | NDArray,
        surface: Surface | NDArray,
    ) -> Surface | NDArray: ...

    @abstractmethod
    def ejecta_shape(
        self,
        crater: Crater,
        region_view: SurfaceView | NDArray,
        surface: Surface | NDArray,
    ) -> Surface | NDArray: ...

    @abstractmethod
    def rmax(
        self,
        crater: Crater,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
    ) -> float: ...


class CraterQueueManager:
    """
    A manager for craters awaiting emplacement. Craters are processed in order (FIFO),
    but batches of non-overlapping craters can be processed simultaneously.

    Parameters
    ----------
    overlap_fn : Callable[[Crater], tuple[set[int], set[int]]]
        A function that takes a crater and returns a tuple of (node indices, face indices) affected.
    """

    def __init__(self, overlap_fn: Callable[[Crater], tuple[set[int], set[int]]]):
        self._queue: list[Crater] = []
        self._active_entities: set[int] = set()
        self._overlap_fn = overlap_fn

    def push(self, crater: Crater) -> None:
        self._queue.append(crater)

    def peek_next_batch(self) -> list["Crater"]:
        """
        Return a list of the next batch of craters that do not overlap with each other
        or the current active region.
        """
        batch = []
        reserved = set(self._active_entities)
        for crater in self._queue:
            node_indices, face_indices = self._overlap_fn(crater)
            combined = node_indices | face_indices
            if reserved.isdisjoint(combined):
                batch.append(crater)
                reserved.update(combined)
            else:
                break
        return batch

    def pop_batch(self, batch: list["Crater"]) -> None:
        """
        Remove a processed batch of craters from the queue.
        """
        for crater in batch:
            self._queue.remove(crater)
            node_indices, face_indices = self._overlap_fn(crater)
            self._active_entities.update(node_indices | face_indices)

    def clear_active(self) -> None:
        """
        Clear the active region set after batch processing is complete.
        """
        self._active_entities.clear()

    def is_empty(self) -> bool:
        return len(self._queue) == 0


import_components(__name__, __path__, ignore_private=True)
