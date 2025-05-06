from __future__ import annotations
from abc import abstractmethod
from typing import Any
from cratermaker.core.crater import Crater
from cratermaker.components.surface import Surface
from cratermaker.utils.component_utils import ComponentBase, import_components

class Morphology(ComponentBase):
    def __init__(self, crater : Crater | None = None, **kwargs: Any) -> None:
        """
        Initialize the Morphology class.

        Parameters
        ----------
        crater : Crater, optional
            The crater currently attached to the morphology model.
        **kwargs : Any
            Additional keyword arguments.

        Raises
        -------

        TypeError
            If the crater is not an instance of Crater.

        """
        super().__init__(**kwargs)
        object.__setattr__(self, "_crater" , crater)

    def __repr__(self) -> str:
        base = super().__repr__()
        if self.crater is None:
            return base
        return (
            f"{base}\n"
            f"{self.crater}"
        )

    @classmethod
    def maker(cls, 
             morphology: str | type[Morphology] | Morphology| None = None, 
             crater : Crater | None = None,
             **kwargs: Any) -> Morphology:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, the default "simplemoon" is used.
        crater : Crater, optional
            The crater currently attached to the morphology model.
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
        morphology = super().maker(component=morphology, crater=crater, **kwargs)
        return morphology

    @abstractmethod
    def form_crater(self, 
                    surface: Surface,
                    crater: Crater | None = None,
                    **kwargs) -> None: ...    
    
    @property
    def crater(self):
        """
        The crater to be created.
        
        Returns
        -------
        Crater
        """ 
        return self._crater
    
    @crater.setter
    def crater(self, value):
        if value is not None and not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value


import_components(__name__, __path__, ignore_private=True)

