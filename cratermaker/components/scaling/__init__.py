import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.general_utils import _to_config, parameter
from cratermaker.utils.custom_types import FloatLike
from cratermaker.core.target import Target
from cratermaker.components.impactor import ImpactorModel, get_impactor_model

class ScalingModel(ABC):
    """
    This is the abstract base class for all scaling models. It defines the interface for converting between projectile and crater diameters.
    """
    def __init__(self, 
                 target: Target | str = "Moon",
                 impactor: ImpactorModel | str = "asteroids",
                 rng : Generator | None = None,
                 **kwargs):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_target_density", None)
        object.__setattr__(self, "_impactor", None)
        if isinstance(target, str):
            try:
                self.target = Target(target,**kwargs)
            except:
                raise ValueError(f"Invalid target name {target}")
        else:
            self.target = Target

    @abstractmethod
    def projectile_to_transient(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_projectile(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_final(self, transient_diameter: FloatLike) -> tuple[np.float64, str]: ...
    @abstractmethod
    def final_to_transient(self, final_diameter: FloatLike, morphology_type: str | None = None, **kwargs) -> np.float64: ...


    @parameter
    def model(self):
        """
        The registered name of this scaling model set by the @register_scaling_model decorator.
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
    def target_density(self):
        """
        Volumentric density of material in kg/m^3.
        
        Returns
        -------
        float 
        """
        return self._target_density
    
    @target_density.setter
    def target_density(self, value):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("target_density must be a numeric value or None")
            if value < 0:
                raise ValueError("target_density must be a positive number")
            self._target_density = float(value)

    @property
    def impactor(self):
        """
        The impactor model for the impact.
        
        Returns
        -------
        ImpactorModel
        """ 
        return self._impactor
    
    @impactor.setter
    def impactor(self, value):
        if value is None:
            value = get_impactor_model("asteroids")
        elif isinstance(value, str):
            try:
                value = get_impactor_model(value)
            except KeyError:
                raise ValueError(f"Invalid impactor name {value}")
        elif not isinstance(value, ImpactorModel):
            raise TypeError("impactor must be an instance of ImpactorModel")
        self._impactor = value
        return


_registry: dict[str, ScalingModel] = {}

def register_scaling_model(name: str):
    """
    Class decorator to register an impactor->crater size scaling component under the given key.
    """
    def decorator(cls):
        cls._model = name 
        cls._user_defined = set()
        cls._user_defined.add("model")
        _registry[name] = cls
        return cls
    return decorator

def available_scaling_models() -> list[str]:
    """Return list of all registered catalogue names."""
    return list(_registry.keys())

def get_scaling_model(name: str):
    """Return the component instance for the given name (KeyError if not found)."""
    return _registry[name]

# This loop will import every .py in this folder, causing those modules
# (which use @register_scaling_model) to run and register themselves.
package_dir = __path__[0]
for finder, module_name, is_pkg in pkgutil.iter_modules([package_dir]):
    importlib.import_module(f"{__name__}.{module_name}")
