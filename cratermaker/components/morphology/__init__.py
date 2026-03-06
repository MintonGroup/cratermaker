from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from math import pi
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from cratermaker.components.crater import Crater
from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import parameter

if TYPE_CHECKING:
    from cratermaker.components.counting import Counting
    from cratermaker.components.crater.morphologycrater import MorphologyCrater
    from cratermaker.components.production import Production
    from cratermaker.components.surface import LocalSurface, Surface


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

        # Because a Surface object is associated with a Counting object, we should first check to see if we are receiving a Counting object first so that we don't end up creating a spurious Surface object that we don't want.
        if surface is None and isinstance(counting, Counting):
            self.surface = counting.surface
            self.counting = counting
        else:
            self.surface = Surface.maker(surface, **kwargs)
            if counting is not None:
                crater_cls = kwargs.pop("Crater", self.Crater)
                self.counting = Counting.maker(counting, surface=self.surface, Crater=crater_cls, **kwargs)

        if self.counting is not None:
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
        from cratermaker.components.crater.morphologycrater import MorphologyCrater

        if self._queue_manager is None:
            self._init_queue_manager()

        if craters is None:
            craters = [MorphologyCrater.maker(**kwargs)]
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
        from cratermaker.components.crater.morphologycrater import MorphologyCrater

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
        from cratermaker.components.crater.morphologycrater import MorphologyCrater

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
        self._queue_manager: CraterQueueManager | None = None

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
            # Make sure the associated Crater class is correct for this Morphology class
            if not issubclass(counting.Crater, self.Crater):
                counting.Crater = self.Crater
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
    def Crater(self) -> type[Crater]:
        """The Crater class used for this counting component, which is determined by the morphology component."""
        from cratermaker.components.crater.morphologycrater import MorphologyCrater

        return MorphologyCrater


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
