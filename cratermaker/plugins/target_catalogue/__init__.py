import pkgutil
import importlib
from abc import ABC, abstractmethod

class TargetCataloguePlugin(ABC):
    """
    Abstract base class that all target‑catalogue plugins must implement.
    """
    @abstractmethod
    def get_targets(self) -> dict[str, dict]:
        """
        Return a mapping from body name to its property dict, e.g.:
           { "Earth": { "radius": 6371000.0, ... }, ... }
        """
        ...

_registry: dict[str, TargetCataloguePlugin] = {}

def register_target_catalogue(name: str):
    """
    Class decorator to register a target‑catalogue plugin under the given key.
    """
    def decorator(cls):
        instance = cls()
        _registry[name] = cls
        return cls
    return decorator

def available_target_catalogues() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_target_catalogue(name: str):
    """Return the plugin instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_target_catalogue) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")