import pkgutil
import importlib
from abc import ABC, abstractmethod

class ScalingModel(ABC):
    @abstractmethod
    def get_morphology_type(self, final_diameter): ...
    @abstractmethod
    def projectile_to_crater(self, projectile): ...
    @abstractmethod
    def crater_to_projectile(self, crater): ...

_registry: dict[str, ScalingModel] = {}

def register_scaling_model(name: str):
    """
    Class decorator to register an impactor->crater size scaling plugin under the given key.
    """
    def decorator(cls):
        instance = cls()
        _registry[name] = instance
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