import pkgutil
import importlib
from abc import ABC, abstractmethod
from cratermaker.utils.custom_types import FloatLike, PairOfFloats
from collections.abc import Sequence
from numpy.typing import ArrayLike
from typing import Any, Union
import numpy as np

class ProductionModel(ABC):
    @abstractmethod
    def function(self,
            diameter: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
            age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
            age_end: FloatLike | Sequence[FloatLike] | ArrayLike | None = None,
            **kwargs: Any,
            ) -> Union[FloatLike, ArrayLike]: ...
    @abstractmethod
    def function_inverse(self,
             diameter: FloatLike | Sequence[FloatLike] | ArrayLike,
             cumulative_number_density: FloatLike | Sequence[FloatLike] | ArrayLike,
             **kwargs: Any,
             ) -> Union[FloatLike, ArrayLike]: ...
    @abstractmethod
    def sample(self,
               age: FloatLike | None = None,
               age_end: FloatLike | None = None,
               diameter_number: PairOfFloats | None = None,
               diameter_number_end: PairOfFloats | None = None,
               diameter_range: PairOfFloats | None = None,
               area: FloatLike | None = None, 
               return_age: bool = True
               ) -> np.ndarray: ...
    @abstractmethod
    def chronology(self,
             age: FloatLike | Sequence[FloatLike] | ArrayLike = 1.0,
             check_valid_time: bool=True
             ) -> Union[FloatLike, ArrayLike]: ...
    @abstractmethod
    def to_config(self) -> dict: ...

    def __init__(self):
        object.__setattr__(self, "_user_defined", set())
        self._user_defined.add("model")

    def to_config(self) -> dict:
        """
        Only include those parameters the user actually set.
        """
        return {name: getattr(self, name) for name in self._user_defined}

    @property
    def model(self):
        """
        The registered name of this scaling model set by the @register_scaling_model decorator.
        """ 
        return self._model

_registry: dict[str, ProductionModel] = {}

def register_production_model(name: str):
    """
    Class decorator to register a production model plugin under the given key.
    """
    def decorator(cls):
        cls._model = name 
        _registry[name] = cls
        return cls
    return decorator

def available_production_models() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_production_model(name: str):
    """Return the plugin instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_production_model) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")