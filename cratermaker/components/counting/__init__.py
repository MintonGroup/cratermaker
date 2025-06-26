from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any

import numpy as np
import uxarray as uxr
from numpy.random import Generator

from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_NAME = "crater_id"
_TALLY_LONG_NAME = "Unique crater identification number"
_N_LAYER = 8
_RIM_BUFFER_FACTOR = 1.2


class Counting(ComponentBase):
    _registry: dict[str, Counting] = {}

    """
    Base class for all crater counting models. It defines the interface for tallying the observable craters on a surface.

    Parameters
    ----------
    surface : Surface | LocalSurface
        The surface or local surface view to be counted. 
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        surface: Surface | LocalSurface,
        **kwargs: Any,
    ):
        self.surface = surface
        rng = kwargs.pop("rng", surface.rng)
        super().__init__(rng=rng, **kwargs)
        self.reset()

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
        **kwargs: Any,
    ) -> Counting:
        """
        Initialize a crater counting model based on the provided name or class.

        Parameters
        ----------
        counting : str, Counting or None, default=None
            The name of the counting model to initialize. If None, the default model is used.
        surface : Surface | LocalSurface
            The surface or local surface view to be counted.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        Counting
            An instance of the specified counting model.

        Raises
        ------
        KeyError
            If the specified counting model name is not found in the registry.
        TypeError
            If the specified counting model is not a string or a subclass of Scaling.
        """
        if counting is None:
            counting = "minton2019"

        if surface is None:
            raise ValueError("surface must be provided")

        counting = super().maker(
            component=counting,
            surface=surface,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nSurface: {self.surface}\n"

    def reset(self):
        """
        Remove all craters count records from the surface.
        """
        dims = ("n_face", "layer")
        data = np.zeros((self.surface.n_face, self.n_layer), dtype=np.uint32)
        uxda = uxr.UxDataArray(
            data=data,
            dims=dims,
            attrs={"long_name": _TALLY_LONG_NAME},
            name=_TALLY_NAME,
            uxgrid=self.surface.uxgrid,
        )

        self.surface._uxds[_TALLY_NAME] = uxda
        self._observed = {}
        return

    def add(self, crater: Crater):
        """
        Add a crater to the surface.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        self.observed[crater.id] = crater
        # Tag a region just outside crater rim with the id
        crater_region = self.surface.extract_region(
            location=crater.location, region_radius=_RIM_BUFFER_FACTOR * crater.final_radius
        )
        if crater_region:
            insert_layer = -1
            for i in reversed(range(self.n_layer)):
                if np.any(self.surface.uxds[_TALLY_NAME].isel(layer=i).data[crater_region.face_indices] > 0):
                    # Gather the unique id values for the current layer
                    unique_ids = np.unique(self.surface.uxds[_TALLY_NAME].data[crater_region.face_indices, i])
                    removes = [
                        id for id, v in self.observed.items() if v.id in unique_ids and v.final_diameter < crater.final_diameter
                    ]

                    # For every id that appears in the removes list, set it to 0 in the data array
                    if removes:
                        data = self.surface.uxds[_TALLY_NAME].data[crater_region.face_indices, :]
                        for remove in removes:
                            data[data == remove] = 0
                        self.surface.uxds[_TALLY_NAME].data[crater_region.face_indices, :] = data
                if insert_layer == -1 and np.all(
                    self.surface.uxds[_TALLY_NAME].isel(layer=i).data[crater_region.face_indices] == 0
                ):
                    insert_layer = i
            if insert_layer == -1:
                raise ValueError("Crater counting layers are full")
            data = self.surface.uxds[_TALLY_NAME].data[crater_region.face_indices, :]
            data[:, insert_layer] = crater.id
            self.surface.uxds[_TALLY_NAME].data[crater_region.face_indices, :] = data

        return

    @property
    def surface(self):
        """
        Surface mesh data for the simulation. Set during initialization.
        """
        return self._surface

    @surface.setter
    def surface(self, value):
        from cratermaker.components.surface import LocalSurface, Surface

        if not isinstance(value, (Surface | LocalSurface)):
            raise TypeError("surface must be an instance of Surface or LocalSurface")
        self._surface = value

    @property
    def n_layer(self) -> int:
        """
        Number of layers in the counting model.
        """
        return _N_LAYER

    @property
    def observed(self) -> dict[int, Crater]:
        """
        List of observed craters on the surface.

        Returns
        -------
            A dict with the crater id as the key and the Crater object as the values
        """
        return self._observed


import_components(__name__, __path__, ignore_private=True)
