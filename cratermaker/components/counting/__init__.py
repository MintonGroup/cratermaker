from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.random import Generator

from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_NAME = "crater_id"
_TALLY_LONG_NAME = "Unique crater identification number"


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
        self.surface.add_data(name=_TALLY_NAME, long_name=_TALLY_LONG_NAME, data=0, dtype=np.int64, overwrite=True)

    def add(self, crater: Crater):
        """
        Add a crater to the surface.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        # Tag a region twice the size of the crater rim with the id
        crater_region = self.surface.extract_region(location=crater.location, region_radius=crater.final_diameter)
        if crater_region:
            crater_region.add_data(name=_TALLY_NAME, long_name=_TALLY_LONG_NAME, data=crater._id, overwrite=True, dtype=np.int64)
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


import_components(__name__, __path__, ignore_private=True)
