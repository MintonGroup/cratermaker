import pkgutil
import importlib
from abc import ABC, abstractmethod

class MaterialCataloguePlugin(ABC):
    """
    Abstract base class that all material‑catalogue plugins must implement.
    """
    @abstractmethod
    def get_materials(self) -> dict[str, dict]:
        """
        Return a mapping from body name to its property dict, e.g.:
           { "Earth": { "radius": 6371000.0, ... }, ... }
        """
        ...

_registry: dict[str, MaterialCataloguePlugin] = {}

def register_material_catalogue(name: str):
    """
    Class decorator to register a material‑catalogue plugin under the given key.
    """
    def decorator(cls):
        instance = cls()
        _registry[name] = instance
        return cls
    return decorator

def available_material_catalogues() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_material_catalogue(name: str):
    """Return the plugin instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_material_catalogue) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")