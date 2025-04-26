from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any
from cratermaker.core.crater import Crater
from cratermaker.core.surface import Surface
from cratermaker.utils.general_utils import parameter
from cratermaker.core.base import CratermakerBase
from cratermaker.utils.component_loader import import_components

class Morphology(CratermakerBase, ABC):
    """
    Abstract base class for morphology models. A morphology model defines how the surface mesh is modified by a given crater.
    """
    _registry: dict[str, type[Morphology]] = {}

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_crater" , None)

    @classmethod
    def make(cls, morphology : str | type[Morphology] | Morphology | None = None, **kwargs: Any) -> type[Morphology]:
        """
        Initialize the morphology model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, defaults to "simplemoon".
        kwargs : Any
            Additional keyword arguments to pass to the morphology model constructor.

        Returns
        -------
        Morphology
            An instance of the specified morphology model.

        Raises
        ------
        KeyError
            If the specified morphology model name is not found in the registry.
        TypeError
            If the specified morphology model is not a string or a subclass of Morphology.
        """

        if morphology is None:
            morphology = "simplemoon"
        if isinstance(morphology, str):
            if morphology not in cls.available():
                raise KeyError(f"Unknown morphology model: {morphology}. Available models: {cls.available()}")
            return cls._registry[morphology](**kwargs)
        elif isinstance(morphology, type) and issubclass(morphology, Morphology):
            return morphology(**kwargs)
        elif isinstance(morphology, Morphology):
            return morphology
        else:
            raise TypeError(f"morphology must be a string or a subclass of Morphology, not {type(morphology)}")

    @abstractmethod
    def form_crater(self, 
                    surf: Surface,
                    crater: Crater | None = None,
                    **kwargs) -> None: ...    

    @parameter
    def name(self):
        """
        The registered name of this scaling model set by the @register_scaling_model decorator.
        """ 
        return self._name
    
    @property
    def crater(self):
        """
        The crater to be created.
        
        Returns
        -------
        Crater
        """ 
        if self._crater is None:
            raise RuntimeError("No crater has been added to the morphology model yet.")
        return self._crater
    
    @crater.setter
    def crater(self, value):
        if value is not None and not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value

    @classmethod
    def register(cls, name: str):
        """
        Class decorator to register a morphology model component under the given key.
        """
        def decorator(subcls):
            subcls._name = name 
            subcls._registry[name] = subcls
            return subcls
        return decorator

    @classmethod 
    def available(cls) -> list[str]:
        """Return list of all registered catalogue names."""
        return list(cls._registry.keys())

import_components(__name__, __path__, ignore_private=True)

