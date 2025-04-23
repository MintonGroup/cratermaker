import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.general_utils import _to_config, parameter
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils import montecarlo as mc

class ImpactorModel(ABC):
    """
    This is the abstract base class for all impactor models. It defines the interface for generating impactor velocities, angles, and densities for a given target body.
    """
    def __init__(self, 
                 target_name : str = "Moon",
                 mean_velocity : FloatLike = 22100.0, 
                 density : FloatLike = 2500.0,
                 sample_velocities : bool = True,
                 sample_angles : bool = True,
                 sample_directions : bool = True,
                 rng: Generator | None = None,
                 **kwargs):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_target_name", None)
        object.__setattr__(self, "_model", None)
        object.__setattr__(self, "_sample_angles", None)
        object.__setattr__(self, "_sample_velocities", None)
        object.__setattr__(self, "_sample_directions", None)
        object.__setattr__(self, "_mean_velocity", None)
        object.__setattr__(self, "_density", None)
        object.__setattr__(self, "_rng", None)
        self.target_name = target_name
        self.sample_velocities = sample_velocities
        self.sample_angles = sample_angles
        self.sample_directions = sample_directions
        self.mean_velocity = mean_velocity
        self.density = density
        self.rng = rng


    def to_config(self, **kwargs: Any) -> dict:
        return _to_config(self)
    
    def get_projectile(self, **kwargs: Any) -> dict:
        """
        Returns a dictionary of projectile properties that can be passed as arguments to the Crater class.
        
        Parameters
        ----------
        **kwargs : Any
            Additional keyword arguments to be passed to internal functions.
        
        Returns
        -------
        dict
            A dictionary containing the impactor properties.
        """
        return {
            "projectile_velocity": self.velocity,
            "projectile_angle": self.angle,
            "projectile_density": self.density,
            "projectile_direction": self.direction, 
        }

    @parameter
    def model(self):
        """
        The registered name of this impactor model set by the @register_impactor_model decorator.
        """ 
        return self._model
    
    @parameter
    def target_name(self):
        """
        The name of the target body.
        
        Returns
        -------
        str
        """ 
        return self._target_name
    
    @target_name.setter
    def target_name(self, value):
        if not isinstance(value, str):
            raise TypeError("target_name must be a string")
        self._target_name = value
        return 

    @parameter
    def sample_angles(self):
        """
        Flag that determines whether to sample impact angles from a distribution. If set to False, impact angles will be set to 90 degrees (vertical impact). 
        
        Returns
        -------
        bool
        """
        return self._sample_angles
    
    @sample_angles.setter
    def sample_angles(self, value):
        if not isinstance(value, bool):
            raise TypeError("sample_angles must be a boolean value")
        self._sample_angles = value
        return
    
    @parameter
    def sample_velocities(self):
        """
        Flag that determines whether to sample impact velocities from a distribution. If set to False, impact velocities will be set to the mean velocity.
        
        Returns
        -------
        bool
        """
        return self._sample_velocities
    
    @sample_velocities.setter
    def sample_velocities(self, value):
        if not isinstance(value, bool):
            raise TypeError("sample_velocities must be a boolean value")
        self._sample_velocities = value
        return
    
    @parameter
    def sample_directions(self):
        """
        Flag that determines whether to sample impact directions from a distribution. If set to False, impact directions will be set to the mean velocity.
        
        Returns
        -------
        bool
        """
        return self._sample_directions
    
    @sample_directions.setter
    def sample_directions(self, value):
        if not isinstance(value, bool):
            raise TypeError("sample_directions must be a boolean value")
        self._sample_directions = value
        return
    
    @parameter
    def mean_velocity(self):
        """
        The mean velocity of the projectile in m/s.
        
        Returns
        -------
        float 
        """
        return self._mean_velocity
    
    @mean_velocity.setter
    def mean_velocity(self, value):
        """
        Sets the mean velocity of the projectile. If none is provided, it will look at the target name, and if it matches a known body, we will draw from a distribution using the predefined mean velocity value. Otherwise, it will raise an error.
        
        Parameters
        ----------
        value : float | None
            The mean velocity in m/s or None to use a default value based on the target name.
        
        Raises
        ------
        ValueError
            If the provided value is negative or None.
        TypeError
            If the provided value is not a numeric value or None.
        """
        if isinstance(value, FloatLike):
            if value < 0:
                raise ValueError("mean_velocity must be a positive number")
            self._mean_velocity = float(value)
        else: 
            raise TypeError("mean_velocity must be a numeric value") 
        return

    @property 
    def angle(self):
        """
        The impact angle in degrees.
        
        Returns
        -------
        float 
        """
        if self.sample_angles:
            return mc.get_random_impact_angle(rng=self.rng)
        else:
            return 90.0

    @property
    def direction(self):
        """
        The impact direction in degrees.
        
        Returns
        -------
        float 
        """
        if self.sample_directions:
            return mc.get_random_impact_direction(rng=self.rng)
        else:
            return 0.0

    @property
    def velocity(self):
        """
        The impact velocity in m/s.
        
        Returns
        -------
        float 
        """
        if self.sample_velocities:
            return mc.get_random_velocity(self.mean_velocity, rng=self.rng)
        else:
            return self.mean_velocity

    @parameter
    def density(self):
        """
        The density of the impactor in kg/m^3.
        
        Returns
        -------
        float 
        """
        return self._density
    
    @density.setter
    def density(self, value):
        """
        Sets the density of the impactor in kg/m^3. 
        
        Parameters
        ----------
        value : float | None
            The density in kg/m^3 or None to use a default value.
        
        Raises
        ------
        ValueError
            If the provided value is negative or None.
        TypeError
            If the provided value is not a numeric value or None.
        """
        if not isinstance(value, FloatLike):
            raise TypeError("density must be a numeric value")
        if value < 0:
            raise ValueError("density must be a positive number")
        self._density = float(value)

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
