from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from math import pi
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad
from tqdm import tqdm

from cratermaker.components.production import Production
from cratermaker.constants import FloatLike
from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import parameter

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface


class Morphology(ComponentBase):
    def __init__(
        self,
        surface: Surface | str | None = None,
        production: Production | str | None = None,
        dosubpixel_degradation: bool = False,
        doslope_collapse: bool = True,
        **kwargs: Any,
    ) -> None:
        """
        Initialize the Morphology class.

        Parameters
        ----------
        surface : str or Surface, optional
            The name of a Surface object, or an instance of Surface, to be associated the morphology model.
        production : str or Production, optional
            The name of a Production object, or an instance of Production, to be associated with the morphology model. This is used for subpixel degradation in the emplace method. It is otherwise ignored.
        dosubpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Defaults to True.
        doslope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Defaults to True.
        **kwargs : Any

        """
        from cratermaker.components.surface import Surface

        super().__init__(**kwargs)
        object.__setattr__(self, "_production", None)
        self._surface = Surface.maker(surface, **kwargs)
        self._queue_manager: CraterQueueManager | None = None
        if production is not None:
            self._production = Production.maker(production, **kwargs)
        self.dosubpixel_degradation = dosubpixel_degradation
        self.doslope_collapse = doslope_collapse
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
        dosubpixel_degradation: bool = False,
        doslope_collapse: bool = True,
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
        dosubpixel_degradation : bool, optional
            If True, subpixel degradation will be performed during the emplacement of craters. Defaults to True.
        doslope_collapse : bool, optional
            If True, slope collapse will be performed during the emplacement of craters. Defaults to True.
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
        morphology = super().maker(
            component=morphology,
            surface=surface,
            production=production,
            dosubpixel_degradation=dosubpixel_degradation,
            doslope_collapse=doslope_collapse,
            **kwargs,
        )
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
        Form the interior of the crater by altering the elevation variable of the surface mesh.

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
        self._face_index = self.surface.find_nearest_face(crater.location)

        # Test if the ejecta is big enough to modify the surface
        ejecta_rmax = self.rmax(crater, minimum_thickness=self.surface.smallest_length, feature="ejecta")
        ejecta_region = self.surface.extract_region(crater.location, ejecta_rmax)
        ejecta_area = pi * ejecta_rmax**2
        if (
            ejecta_region is None or ejecta_area < self.surface.face_area[self.face_index]
        ):  # The crater is too small to change the surface
            return

        crater_rmax = self.rmax(crater, minimum_thickness=self.surface.smallest_length, feature="crater")
        crater_region = ejecta_region.extract_subregion(crater_rmax)
        crater_volume = None
        if crater_region is not None:  # The crater is big enough to affect the surface
            crater_area = pi * crater_rmax**2

            # Check to make sure that the face at the crater location is not smaller than the crater area
            if crater_area > self.surface.face_area[self.face_index]:
                # Form the crater shape
                elevation_change = self.crater_shape(crater, crater_region)
                crater_region.update_elevation(elevation_change)
                if self.doslope_collapse:
                    crater_region.slope_collapse()
                crater_volume = crater_region.compute_volume(elevation_change[: crater_region.n_face])

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

        # Now form the ejecta blanket
        ejecta_thickness, ejecta_intensity = self.ejecta_shape(crater, ejecta_region)

        if crater_volume:
            ejecta_volume = ejecta_region.compute_volume(ejecta_thickness[: ejecta_region.n_face])
            conservation_factor = -crater_volume / ejecta_volume
            ejecta_thickness *= conservation_factor

        ejecta_region.add_data(
            "ejecta_thickness",
            long_name="ejecta thickness",
            units="m",
            data=ejecta_thickness[: ejecta_region.n_face],
        )

        ejecta_region.update_elevation(ejecta_thickness)

        k_ej = self.ejecta_burial_degradation(ejecta_thickness[: ejecta_region.n_face], ejecta_soften_factor=1.50)
        ejecta_region.apply_diffusion(k_ej)

        k_deg = self.degradation_function(crater.final_diameter, fe=100) * ejecta_intensity[: ejecta_region.n_face]
        ejecta_region.apply_diffusion(k_deg)

        return

    def compute_subpixel_degradation(
        self,
        age_start: float,
        age_end: float,
        production: Production | str | None = None,
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
        production : Production or str, optional
            The production object containing the production population that will contribute to the degradation. This is only necessary if the Morphology model was not initialized with a production object. If not provided, the production object associated with this morphology model will be used. If passed, this will override the current production object.
        **kwargs : Any
            Additional keyword arguments for the degradation function.
        """
        dc_min = 1e-8  # Minimum crater size for subpixel degradation calculation.

        if age_end >= age_start:
            raise ValueError("age_end must be less than age_start.")

        production = Production.maker(production, **kwargs) if production is not None else self.production

        if not hasattr(self, "_Kdiff"):
            self._Kdiff = np.zeros_like(self.surface.face_elevation)

        def _subpixel_degradation(final_diameter):
            fe = 100.0
            k = self.degradation_function(final_diameter, fe)
            n = production.function(
                diameter=final_diameter,
                age=age_start,
                age_end=age_end,
                validate_inputs=False,
            ).item()
            degradation_region_area = np.pi * (final_diameter / 2) * fe
            return k * n * degradation_region_area

        for face_indices, dc_max in zip(self.surface.face_bin_indices, self.surface.face_bin_max_sizes, strict=False):
            delta_kdiff, _ = quad(_subpixel_degradation, dc_min, dc_max)
            self._Kdiff[face_indices] += delta_kdiff

        # If any Kdiff values reaches a threshold where a meaningful amount of diffusion will occur on the surface, then go ahead and apply it.
        # Otherwise, degradation will continue to accumulate until the next batch of craters is processed.
        if np.any(self._Kdiff / self.surface.face_area > 1):
            self.apply_subpixel_degradation()

        return

    def apply_subpixel_degradation(self) -> None:
        """
        Apply subpixel degradation to the surface using the current Kdiff values.

        This method is called after all craters have been processed and is used to
        apply the accumulated degradation effects.
        """
        if not hasattr(self, "_Kdiff"):
            raise RuntimeError("Kdiff has not been initialized. Call compute_subpixel_degradation first.")

        self.surface.apply_diffusion(self._Kdiff)
        self._Kdiff = np.zeros_like(self.surface.face_elevation)
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
        rmax = self.rmax(crater, minimum_thickness=self.surface.smallest_length, feature="ejecta")
        region = self.surface.extract_region(crater.location, rmax)
        if region is None:
            return set(), set()
        if isinstance(region.node_indices, slice) or isinstance(region.face_indices, slice):
            return set(np.arange(self.surface.n_node)[region.node_indices]), set(
                np.arange(self.surface.n_face)[region.face_indices]
            )
        return set(region.node_indices), set(region.face_indices)

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
        Add a crater to the queue for later emplacement.

        Automatically initializes the queue manager if it hasn't been set.

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
                raise RuntimeError("Surface must be provided to initialize queue manager.")
            self._init_queue_manager()

        if crater is None:
            crater = Crater.maker(**kwarg)
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

                if self.dosubpixel_degradation and len(batch) > 1:
                    # If the craters have age values attached to them, we can perform subpixel degradation between time values
                    agevals = [crater.age for crater in batch if crater.age is not None]
                    if len(agevals) > 1:
                        self.compute_subpixel_degradation(age_start=max(agevals), age_end=min(agevals))

                self._queue_manager.pop_batch(batch)
                self._queue_manager.clear_active()

            if self.dosubpixel_degradation:
                self.apply_subpixel_degradation()
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

    @abstractmethod
    def degradation_function(self, final_diameter: FloatLike, fe: FloatLike) -> float: ...

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

        if not isinstance(surface, (Surface | str)):
            raise TypeError("surface must be an instance of Surface or a string")
        self._surface = Surface.maker(surface)
        self._queue_manager: CraterQueueManager | None = None

    @property
    def production(self) -> Production:
        """
        The production object associated with this morphology model.
        """
        return self._production

    @production.setter
    def production(self, production: Production) -> None:
        """
        Set the production object associated with this morphology model.
        """
        from cratermaker.components.production import Production

        if not isinstance(production, (Production | str)):
            raise TypeError("production must be an instance of Production or a string")
        self._production = Production.maker(production)

    @parameter
    def dosubpixel_degradation(self) -> bool:
        """
        Whether to perform subpixel degradation during crater emplacement.
        """
        return self._dosubpixel_degradation

    @dosubpixel_degradation.setter
    def dosubpixel_degradation(self, value: bool) -> None:
        """
        Set whether to perform subpixel degradation during crater emplacement.
        """
        if not isinstance(value, bool):
            raise TypeError("dosubpixel_degradation must be a boolean value")
        self._dosubpixel_degradation = value

    @parameter
    def doslope_collapse(self) -> bool:
        """
        Whether to perform slope collapse during crater emplacement.
        """
        return self._doslope_collapse

    @doslope_collapse.setter
    def doslope_collapse(self, value: bool) -> None:
        """
        Set whether to perform slope collapse during crater emplacement.
        """
        if not isinstance(value, bool):
            raise TypeError("doslope_collapse must be a boolean value")
        self._doslope_collapse = value

    @property
    def face_index(self):
        """
        The index of the face closest to the crater location.

        Returns
        -------
        int
        """
        return self._face_index


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


import_components(__name__, __path__, ignore_private=True)
