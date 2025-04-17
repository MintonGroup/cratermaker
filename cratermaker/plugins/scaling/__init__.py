import pkgutil
import importlib
from abc import ABC, abstractmethod

class ScalingModel(ABC):
    @abstractmethod
    def projectile_to_crater(self, projectile): ...
    @abstractmethod
    def crater_to_projectile(self, crater): ...
    @abstractmethod
    def to_config(self) -> dict: ...

    def __init__(self):
        object.__setattr__(self, "_user_defined", set())
        self._user_defined.add("scaling_model")

    @property
    def scaling_model(self):
        """
        The registered name of this scaling model set by the @register_scaling_model decorator.
        """ 
        return self._scaling_model

_registry: dict[str, ScalingModel] = {}

def register_scaling_model(name: str):
    """
    Class decorator to register an impactor->crater size scaling plugin under the given key.
    """
    def decorator(cls):
        cls._scaling_model = name 
        _registry[name] = cls
        return cls
    return decorator

def available_scaling_models() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_scaling_model(name: str):
    """Return the plugin instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_scaling_model) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")