import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.general_utils import _to_config, parameter
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils import montecarlo as mc
from cratermaker.core import Target

class ImpactorModel(ABC):
    def __init__(self, 
                 target: Target | str = "Moon",
                 rng: Generator | None = None,
                 **kwargs):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_rng", None)
        object.__setattr__(self, "_density", None)
        if isinstance(target, str):
            try:
                target = Target(target,**kwargs)
            except:
                raise ValueError(f"Invalid target name {target}")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target or a valid name of a target body")
        self.rng = rng


    def to_config(self, **kwargs: Any) -> dict:
        return _to_config(self)

    @parameter
    def model(self):
        """
        The registered name of this impactor model set by the @register_impactor_model decorator.
        """ 
        return self._model
    
    @property
    def target(self):
        """
        The target body for the impact.
        
        Returns
        -------
        Target
        """ 
        return self._target
    
    @target.setter
    def target(self, value):
        if value is None:
            self._target = Target(name="Moon")
            return
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value
        return 


    @property
    def density(self):
        """
        Volumetric density of the projectile in kg/m^3.
        
        Returns
        -------
        float 
        """
        return self._density
    
    @density.setter
    def density(self, value):
        if value is not None:
            if not isinstance(value, FloatLike): 
                raise TypeError("density must be a numeric value or None")
            if value < 0:
                raise ValueError("density must be a positive number")
            self._density = float(value)
        else:
            if self.target.surface_density is not None:
                self._density = self.target.surface_density

    @property
    def angle(self):
        """
        The impact angle in degrees.
        
        Returns
        -------
        float 
        """
        return mc.get_random_impact_angle(self._angle)


    @property
    def rng(self):
        """
        A random number generator instance.
        
        Returns
        -------
        Generator
        """ 
        return self._rng
    
    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()       


_registry: dict[str, ImpactorModel] = {}

def register_impactor_model(name: str):
    """
    Class decorator to register an impactor->crater size impactor component under the given key.
    """
    def decorator(cls):
        cls._model = name 
        cls._user_defined = set()
        cls._user_defined.add("model")
        _registry[name] = cls
        return cls
    return decorator

def available_impactor_models() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_impactor_model(name: str):
    """Return the component instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_impactor_model) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")
