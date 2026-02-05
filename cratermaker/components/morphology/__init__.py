from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from dataclasses import asdict, dataclass
from math import pi
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad
from tqdm import tqdm

from cratermaker.components.crater import Crater
from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import format_large_units, parameter

if TYPE_CHECKING:
    from cratermaker.components.counting import Counting
    from cratermaker.components.production import Production
    from cratermaker.components.surface import LocalSurface, Surface


@dataclass(frozen=True, slots=True)
class MorphologyCrater(Crater):
    face_index: int | None = None
    ejecta_region: LocalSurface | None = None
    crater_region: LocalSurface | None = None
    crater_rmax: float | None = None
    ejecta_rmax: float | None = None
    affected_face_indices: set[int] | None = None
    affected_node_indices: set[int] | None = None

    def __str__(self) -> str:
        base = super().__str__()

        return (
            f"{base}\n"
            f"Face index of crater center: {self.face_index}\n"
            f"Maximum affected radius: {format_large_units(self.rmax, quantity='length')}\n"
            f"Crater region maximum radius: {format_large_units(self.crater_rmax, quantity='length')}\n"
            f"Crater region: {self.crater_region}\n"
            f"Ejecta region maximum radius: {format_large_units(self.ejecta_rmax, quantity='length')}\n"
            f"Ejecta region: {self.ejecta_region}\n"
        )

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        **kwargs: Any,
    ) -> MorphologyCrater:
        """
        Initialize a MorphologyCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters. It is a stub that is intended to be extended by specific morphology models.

        Parameters
        ----------
        crater : Crater | None
            The crater object to be converted into a MorphologyCrater. If None, a new Crater object is created using the provided keyword arguments.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        kwargs : Any
            The keyword arguments provided are passed down to :func:`cratermaker.crater.Crater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.
        """
        if crater is None:
            crater = super(MorphologyCrater, cls).maker(**kwargs)
        base_fields = asdict(crater)
        return cls(
            **base_fields,
        )


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

        if crater.crater_region is not None:  # The crater is big enough to affect the surface
            crater_area = pi * crater.crater_rmax**2

            # Check to make sure that the face at the crater location is not smaller than the crater area
            if crater_area > self.surface.face_area[crater.face_index]:
                # Form the crater shape
                # Because the crater_region is locked in a mutable dataclass, we can't use it directly to update the elevation, so we have to rebuild it
                crater_region = self._get_mutable_region(crater, self.surface, "crater")
                elevation_change = self.crater_shape(crater, crater_region)
                crater_region.update_elevation(elevation_change)
                if self.do_slope_collapse:
                    crater_region.slope_collapse()
                self._excavated_volume = crater_region.compute_volume(elevation_change[: crater_region.n_face])

                # Remove any ejecta from the surface
                inner_crater_region = crater_region.extract_subregion(crater.final_radius)
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
                    ejecta_region = self._get_mutable_region(crater, self.surface, "ejecta")
                    self.counting.add(crater, count_region=ejecta_region, **kwargs)

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
            |kwargs

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
        ejecta_region = self._get_mutable_region(crater, self.surface, "ejecta")
        ejecta_thickness, ejecta_intensity = self.ejecta_shape(crater, ejecta_region)
        ejecta_volume = ejecta_region.compute_volume(ejecta_thickness[: ejecta_region.n_face])
        conservation_factor = -self._excavated_volume / ejecta_volume
        ejecta_thickness *= conservation_factor

        ejecta_region.add_data(
            "ejecta_thickness",
            long_name="ejecta thickness",
            units="m",
            data=ejecta_thickness[: ejecta_region.n_face],
        )

        ejecta_region.update_elevation(ejecta_thickness)

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

        from concurrent.futures import ThreadPoolExecutor

        def _batch_process(pbar=None):
            tally_cadence = 10000
            nacumulated = 0
            while not self._queue_manager.is_empty():
                batch = self._queue_manager.peek_next_batch()

                def process(crater):
                    try:
                        self.form_crater(crater)
                        if pbar is not None:
                            pbar.update(1)
                    except Exception as e:
                        raise RuntimeError(f"Error processing crater {crater}: {e}") from e

                # max_workers=1 because something needs access to HDF files (probably grid.nc) that is not thread safe
                with ThreadPoolExecutor(max_workers=1) as executor:
                    executor.map(process, batch)

                if self.do_subpixel_degradation and len(batch) > 1:
                    # If the craters have age values attached to them, we can perform subpixel degradation between time values
                    agevals = [crater.age for crater in batch if crater.age is not None]
                    if len(agevals) > 1:
                        self.compute_subpixel_degradation(age_start=max(agevals), age_end=min(agevals))

                self._queue_manager.pop_batch(batch)
                self._queue_manager.clear_active()
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

    @staticmethod
    def _get_mutable_region(crater: MorphologyCrater, surface: Surface, feature: Literal["crater", "ejecta"]) -> LocalSurface:
        from cratermaker.components.surface import LocalSurface

        if feature == "crater":
            return LocalSurface(
                surface=surface,
                face_indices=crater.crater_region.face_indices,
                node_indices=crater.crater_region.node_indices,
                edge_indices=crater.crater_region.edge_indices,
                location=crater.crater_region.location,
                region_radius=crater.crater_region.region_radius,
            )
        elif feature == "ejecta":
            return LocalSurface(
                surface=surface,
                face_indices=crater.ejecta_region.face_indices,
                node_indices=crater.ejecta_region.node_indices,
                edge_indices=crater.ejecta_region.edge_indices,
                location=crater.ejecta_region.location,
                region_radius=crater.ejecta_region.region_radius,
            )
        else:
            raise ValueError("feature must be either 'crater' or 'ejecta'")

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
        self._active_nodes: set[int] = set()
        self._active_faces: set[int] = set()
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
        reserved_nodes = set(self._active_nodes)
        reserved_faces = set(self._active_faces)
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


import_components(__name__, __path__)
