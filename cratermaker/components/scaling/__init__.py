import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.general_utils import  parameter
from cratermaker.utils.custom_types import FloatLike
from cratermaker.core.target import Target, _init_target
from cratermaker.components.impactor import ImpactorModel, _init_impactor
from cratermaker.core.base import CratermakerBase
class ScalingModel(CratermakerBase, ABC):
    """
    This is the abstract base class for all scaling models. It defines the interface for converting between projectile and crater diameters.

    Parameters
    ----------
    target : Target | str, default="Moon"
        The target body for the impact. Can be a Target object or a string representing the target name.
    impactor : ImpactorModel | str, default="asteroids"
        The impactor model for the impact. Can be an ImpactorModel object or a string representing the impactor name.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    """
    def __init__(self, 
                 target: Target | str = "Moon",
                 impactor: ImpactorModel | str = "asteroids",
                 rng : Generator | None = None,
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs: Any):
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_target_density", None)
        object.__setattr__(self, "_impactor", None)
        # combine the kwargs with the common_args, giving common_args priority
        kwargs = {**kwargs, **vars(self.common_args)}
        self._target = _init_target(target, **kwargs)
        self._impactor = _init_impactor(impactor, **kwargs)

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
        self._target = _init_target(value, **vars(self.common_args))
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
        self._impactor = _init_impactor(value, **vars(self.common_args))
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


def _init_scaling(scaling: str | ScalingModel | None = None, 
                  **kwargs: Any) -> ScalingModel:
    """
    Initialize a scaling model based on the provided name or class.
    Parameters
    ----------
    scaling : str, ScalingModel, or None, default=None
        The name of the scaling model to initialize. If None, the default model is used.
    kwargs : Any
        Additional keyword arguments to pass to the scaling model constructor.

    Returns
    -------
    ScalingModel
        An instance of the specified scaling model.
    Raises
    ------
    KeyError
        If the specified scaling model name is not found in the registry.
    TypeError
        If the specified scaling model is not a string or a subclass of ScalingModel.
    """

    if scaling is None:
        scaling = "simplemoon"
    if isinstance(scaling, str):
        if scaling not in available_scaling_models():
            raise KeyError(f"Unknown scaling model: {scaling}. Available models: {available_scaling_models()}")
        return get_scaling_model(scaling)(**kwargs)
    elif isinstance(scaling, type) and issubclass(scaling, ScalingModel):
        return scaling(**kwargs)
    elif isinstance(scaling, ScalingModel):
        return scaling
    else:
        raise TypeError(f"scaling must be a string or a subclass of ScalingModel, not {type(scaling)}")