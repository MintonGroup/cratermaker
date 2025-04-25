import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
from cratermaker.core.crater import Crater
from cratermaker.core.surface import Surface
from cratermaker.utils.general_utils import parameter
from cratermaker.core.base import CratermakerBase

class MorphologyModel(CratermakerBase, ABC):
    """
    Abstract base class for morphology models. A morphology model defines how the surface mesh is modified by a given crater.
    """
    def __init__(self, **kwargs: Any) -> None:
        """
        Initialize the morphology model with given parameters.
        """
        super().__init__(**kwargs)
        object.__setattr__(self, "_crater" , None)

    @abstractmethod
    def form_crater(self, 
                    surf: Surface,
                    crater: Crater | None = None,
                    **kwargs) -> None: ...    

    @parameter
    def model(self):
        """
        The registered name of this scaling model set by the @register_scaling_model decorator.
        """ 
        return self._model
    
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

    
_registry: dict[str, MorphologyModel] = {}

def register_morphology_model(name: str):
    """
    Class decorator to register a morphology model component under the given key.
    """
    def decorator(cls):
        cls._model = name 
        _registry[name] = cls
        return cls
    return decorator

def available_morphology_models() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_morphology_model(name: str):
    """Return the component instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_morphology_model) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")


def _init_morphology(morphology : str | MorphologyModel | None = None, 
                     **kwargs: Any) -> MorphologyModel:
    """
    Initialize the morphology model with the given name or instance.

    Parameters
    ----------
    morphology : str or MorphologyModel or None
        The name of the morphology model to use, or an instance of MorphologyModel. If None, defaults to "simplemoon".
    kwargs : Any
        Additional keyword arguments to pass to the morphology model constructor.

    Returns
    -------
    MorphologyModel
        An instance of the specified morphology model.

    Raises
    ------
    KeyError
        If the specified morphology model name is not found in the registry.
    TypeError
        If the specified morphology model is not a string or a subclass of MorphologyModel.
    """

    if morphology is None:
        morphology = "simplemoon"
    if isinstance(morphology, str):
        if morphology not in available_morphology_models():
            raise KeyError(f"Unknown morphology model: {morphology}. Available models: {available_morphology_models()}")
        return get_morphology_model(morphology)(**kwargs)
    elif issubclass(morphology, MorphologyModel):
        return morphology(**kwargs)
    elif isinstance(morphology, MorphologyModel):
        return morphology
    else:
        raise TypeError(f"morphology must be a string or a subclass of MorphologyModel, not {type(morphology)}")