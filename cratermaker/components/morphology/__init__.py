from __future__ import annotations

from abc import abstractmethod
from math import pi
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from cratermaker.constants import FloatLike
from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import Surface, SurfaceView


class Morphology(ComponentBase):
    def __init__(self, surface: Surface | str | None = None, **kwargs: Any) -> None:
        """
        Initialize the Morphology class.

        Parameters
        ----------
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
        **kwargs : Any
            Additional keyword arguments.

        Raises
        -------

        TypeError
            If the crater is not an instance of Crater.

        """
        from cratermaker.components.surface import Surface

        super().__init__(**kwargs)
        self._surface = Surface.maker(surface, **kwargs)
        self._queue_manager: CraterQueueManager | None = None

    def __str__(self) -> str:
        base = super().__str__()
        return base

    @classmethod
    def maker(
        cls,
        morphology: str | type[Morphology] | Morphology | None = None,
        surface: Surface | str | None = None,
        **kwargs: Any,
    ) -> Morphology:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, the default "simplemoon" is used.
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
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
        morphology = super().maker(component=morphology, surface=surface, **kwargs)
        return morphology

    def emplace(self, crater: Crater | list[Crater], **kwargs: Any) -> None:
        """
        Convenience method to immediately emplace a crater onto the surface.
        Initializes and uses the queue system behind the scenes.

        Parameters
        ----------
        crater : Crater or list[Crater]
            The crater to be emplaced.
        **kwargs : Any
            Additional keyword arguments.
        """
        if self._queue_manager is None:
            self._init_queue_manager()

        if isinstance(crater, list) and len(crater) > 0:
            for c in crater:
                self._enqueue_crater(c)
        elif isinstance(crater, Crater):
            self._enqueue_crater(crater)
        self._process_queue()

    def form_crater(self, crater: Crater, **kwargs: Any) -> None:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        crater : Crater
            The crater object to be formed.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions.

        Returns
        -------
        None
        """

        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        # Find the node and face center of the crater
        self.face_index, self.node_index = self.surface.find_nearest_index(
            crater.location
        )

        # Test if the ejecta is big enough to modify the surface
        ejecta_rmax = self.rmax(
            crater, minimum_thickness=self.surface.smallest_length, feature="ejecta"
        )
        ejecta_region_view = self.surface.extract_region(crater.location, ejecta_rmax)
        ejecta_area = pi * ejecta_rmax**2
        if (
            ejecta_region_view is None
            or ejecta_area < self.surface.face_areas[self.face_index]
        ):  # The crater is too small to change the surface
            return

        crater_rmax = self.rmax(
            crater, minimum_thickness=self.surface.smallest_length, feature="crater"
        )
        crater_region_view = self.surface.extract_region(crater.location, crater_rmax)
        crater_volume = None
        if (
            crater_region_view is not None
        ):  # The crater is big enough to affect the surface
            crater_area = pi * crater_rmax**2

            # Check to make sure that the face at the crater location is not smaller than the crater area
            if crater_area > self.surface.face_areas[self.face_index]:
                # Form the crater shape
                elevation_change = self.crater_shape(crater, crater_region_view)
                crater_region_view.update_elevation(elevation_change)
                crater_volume = crater_region_view.compute_volume(
                    elevation_change[: crater_region_view.n_face]
                )

                # Remove any ejecta from the surface
                inner_crater_region = self.surface.extract_region(
                    crater.location, crater.final_radius
                )
                if inner_crater_region is not None:
                    inner_crater_region.add_data(
                        "ejecta_thickness",
                        long_name="ejecta thickness",
                        units="m",
                        data=0.0,
                        overwrite=True,
                    )

        # Now form the ejecta blanket
        ejecta_thickness, ejecta_intensity = self.ejecta_shape(
            crater, ejecta_region_view
        )

        if crater_volume:
            ejecta_volume = ejecta_region_view.compute_volume(
                ejecta_thickness[: ejecta_region_view.n_face]
            )
            conservation_factor = -crater_volume / ejecta_volume
            ejecta_thickness *= conservation_factor

        ejecta_region_view.add_data(
            "ejecta_thickness",
            long_name="ejecta thickness",
            units="m",
            data=ejecta_thickness[: ejecta_region_view.n_face],
        )

        ejecta_region_view.update_elevation(ejecta_thickness)

        self.degradation_function(
            crater, ejecta_region_view, ejecta_thickness, ejecta_intensity
        )

        ejecta_region_view.slope_collapse()

        return

    def _affected_indices(self, crater: Crater) -> tuple[set[int], set[int]]:
        """
        Determine the set of node and face indices affected by a crater's ejecta blanket.

        Parameters
        ----------
        crater : Crater
            The crater object whose effect region is to be computed.

        Returns
        -------
        affected_nodes : set of int
            The set of node indices affected by the crater.
        affected_faces : set of int
            The set of face indices affected by the crater.
        """
        rmax = self.rmax(
            crater, minimum_thickness=self.surface.smallest_length, feature="ejecta"
        )
        region_view = self.surface.extract_region(crater.location, rmax)
        if region_view is None:
            return set(), set()
        return set(region_view.node_indices), set(region_view.face_indices)

    def _init_queue_manager(self) -> None:
        """
        Initialize the crater queue manager with a surface-dependent overlap function.
        """

        def overlap_fn(crater: Crater) -> tuple[set[int], set[int]]:
            return self._affected_indices(crater)

        self._queue_manager = CraterQueueManager(overlap_fn)

    def _enqueue_crater(
        self,
        crater: Crater | None = None,
        **kwarg: Any,
    ) -> None:
        """
        Add a crater to the queue for later emplacement. Automatically initializes
        the queue manager if it hasn't been set.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to enqueue. If None, one is created from keyword args.
        **kwarg : Any
            Additional keyword arguments for crater construction.

        Raises
        ------
        RuntimeError
            If the queue manager must be initialized but no surface is provided.
        """
        if self._queue_manager is None:
            if self.surface is None:
                raise RuntimeError(
                    "Surface must be provided to initialize queue manager."
                )
            self._init_queue_manager()

        if crater is None:
            crater = Crater.maker(**kwarg)
        self._queue_manager.push(crater)

    def _process_queue(self) -> None:
        """
        Process all queued craters in the order they were added, forming non-overlapping
        batches and applying each to the surface.

        Raises
        ------
        RuntimeError
            If the queue manager has not been initialized.
        """
        if not hasattr(self, "_queue_manager"):
            raise RuntimeError(
                "Queue manager has not been initialized. Call _init_queue_manager first."
            )

        # import threading
        from concurrent.futures import ThreadPoolExecutor

        def _batch_process(pbar=None):
            while not self._queue_manager.is_empty():
                batch = self._queue_manager.peek_next_batch()

                def process(crater):
                    try:
                        self.form_crater(crater)
                        if pbar is not None:
                            pbar.update(1)
                    except Exception as e:
                        print(f"Exception during form_crater: {e}")

                # max_workers=1 because something needs access to HDF files (probably grid.nc) that is not thread safe
                with ThreadPoolExecutor(max_workers=1) as executor:
                    executor.map(process, batch)

                self._queue_manager.pop_batch(batch)
                self._queue_manager.clear_active()
            return

        total_craters = len(self._queue_manager._queue)
        if total_craters > 10:
            with tqdm(
                total=total_craters,
                desc="Processing craters",
                position=0,
                leave=False,
                unit="craters",
            ) as pbar:
                _batch_process(pbar)
        else:
            _batch_process()
        return

    @abstractmethod
    def degradation_function(self) -> None: ...

    @abstractmethod
    def crater_shape(
        self,
        crater: Crater,
        region_view: SurfaceView,
        **kwarg: Any,
    ) -> NDArray[np.float64]: ...

    @abstractmethod
    def ejecta_shape(
        self,
        crater: Crater,
        region_view: SurfaceView,
        **kwarg: Any,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...

    @abstractmethod
    def rmax(
        self,
        crater: Crater,
        minimum_thickness: FloatLike,
        feature: str = "ejecta",
    ) -> float: ...

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

        if not isinstance(surface, (Surface, str)):
            raise TypeError("surface must be an instance of Surface or a string")
        self._surface = Surface.maker(surface)
        self._queue_manager: CraterQueueManager | None = None


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
        self._active_nodes: set[int] = set()
        self._active_faces: set[int] = set()
        self._overlap_fn = overlap_fn

    def push(self, crater: Crater) -> None:
        self._queue.append(crater)

    def peek_next_batch(self) -> list["Crater"]:
        """
        Return a list of the next batch of craters that do not overlap with each other
        or the current active region.
        """
        batch = []
        reserved_nodes = set(self._active_nodes)
        reserved_faces = set(self._active_faces)
        for crater in self._queue:
            node_indices, face_indices = self._overlap_fn(crater)
            if reserved_nodes.isdisjoint(node_indices) and reserved_faces.isdisjoint(
                face_indices
            ):
                batch.append(crater)
                reserved_nodes.update(node_indices)
                reserved_faces.update(face_indices)
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
            self._active_nodes.update(node_indices)
            self._active_faces.update(face_indices)

    def clear_active(self) -> None:
        """
        Clear the active region set after batch processing is complete.
        """
        self._active_nodes.clear()
        self._active_faces.clear()

    def is_empty(self) -> bool:
        return len(self._queue) == 0


import_components(__name__, __path__, ignore_private=True)
