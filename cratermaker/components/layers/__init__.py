from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING, Any

import awkward as ak
import numpy as np
import ragged
import uproot

from cratermaker.components.surface import Surface
from cratermaker.core.base import ComponentBase, import_components


class Layers(ComponentBase):
    """
    Base class for layer models. It defines the interface for tallying the observable craters on a surface.
    """

    _registry: dict[str, Layers] = {}

    def __init__(
        self,
        surface: Surface,
        reset: bool = True,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        surface : Surface
            The surface to associate with this Layers object.
        reset : bool, optional
            Flag to indicate whether to reset the layers and delete any old output files. Default is True.
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.components.surface import Surface

        super().__init__(reset=reset, **kwargs)

        return

    @classmethod
    def maker(
        cls,
        layers: str | Layers | None = None,
        surface: Surface | None = None,
        reset: bool = True,
        **kwargs: Any,
    ) -> Layers:
        """
        Initialize a Layers model with the given name or instance.

        Parameters
        ----------
        layers : str, Layers or None, default=None
            The name of the layers model to initialize. If None, the default model is used.
        surface : Surface | None
            The surface to associate with this Layers object. If None, the surface will be set later and must be set before using the Layers object.
        reset : bool, optional
            Flag to indicate whether to reset the layers and delete any old output files. Default is True
        **kwargs : Any
            |kwargs|

        Returns
        -------
        Layers
            An instance of the specified layers component.

        Raises
        ------
        KeyError
            If the specified layers component name is not found in the registry.
        TypeError
            If the specified layers component type is not a string or a subclass of Layers.
        """
        layers = super().maker(
            component=layers,
            surface=surface,
            reset=reset,
            **kwargs,
        )

        return layers

    def __str__(self) -> str:
        str_repr = super().__str__()
        return str_repr

    def reset(self, **kwargs: Any) -> None:
        """
        Reset the layers layer structure.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|
        """
        super().reset(**kwargs)
        return

    @property
    def surface(self):
        """The Surface object associated with this Layers object."""
        return self._surface

    @surface.setter
    def surface(self, value):
        from cratermaker.components.surface import LocalSurface, Surface

        if not isinstance(value, (Surface | LocalSurface)):
            raise TypeError("surface must be an instance of Surface or LocalSurface")
        self._surface = value


import_components(__name__, __path__)
