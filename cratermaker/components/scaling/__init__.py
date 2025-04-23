import pkgutil
import importlib
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from cratermaker.utils.general_utils import _to_config, parameter
from cratermaker.utils.custom_types import FloatLike
from cratermaker.core import Target

class ScalingModel(ABC):
    def __init__(self, 
                 target: Target | str = "Moon",
                 **kwargs):
        object.__setattr__(self, "_user_defined", set())
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_projectile_density", None)
        object.__setattr__(self, "_projectile_velocity", None)
        object.__setattr__(self, "_projectile_angle", None)
        if isinstance(target, str):
            try:
                target = Target(target,**kwargs)
            except:
                raise ValueError(f"Invalid target name {target}")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target or a valid name of a target body")

        self.projectile_density = projectile_density
        self.projectile_vertical_velocity = projectile_vertical_velocity

    @abstractmethod
    def projectile_to_transient(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_projectile(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_final(self, transient_diameter: FloatLike) -> tuple[np.float64, str]: ...
    @abstractmethod
    def final_to_transient(self, final_diameter: FloatLike, morphology_type: str | None = None, **kwargs) -> np.float64: ...

    def to_config(self, **kwargs: Any) -> dict:
        return _to_config(self)

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
    def projectile_density(self):
        """
        Volumetric density of the projectile in kg/m^3.
        
        Returns
        -------
        float 
        """
        if self._projectile_density is None:
            self.projectile_density = None
        return self._projectile_density
    
    @projectile_density.setter
    def projectile_density(self, value):
        if value is not None:
            if not isinstance(value, FloatLike): 
                raise TypeError("projectile_density must be a numeric value or None")
            if value < 0:
                raise ValueError("projectile_density must be a positive number")
            self._projectile_density = float(value)
        else:
            if self.target_density is not None:
                self._projectile_density = self.target_density

    @property
    def projectile_mean_velocity(self):
        """
        The mean velocity of the projectile in m/s.
        
        Returns
        -------
        float 
        """
        predefined_models = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'MBA']
        predefined_velocities = [41100.0, 29100.0, 24600.0, 22100.0, 10700.0, 5300.0]
        predefined = dict(zip(predefined_models, predefined_velocities))
        if self.target.name in predefined_models:
            pmv = float(predefined[self.target.name])
        elif self.target.name in ["Ceres", "Vesta"]:
            pmv = float(predefined["MBA"])
        else:
            pmv = None
        return pmv
    

    @property
    def projectile_vertical_velocity(self):
        """Get the impact velocity model name."""
        if self._projectile_vertical_velocity is None:
            self.projectile_vertical_velocity = None
        return self._projectile_vertical_velocity

    @projectile_vertical_velocity.setter
    def projectile_vertical_velocity(self, value):
        """"
        Sets the vertical component of the impact velocity. If none is provided, it will look at the target name, and if it matches a known body, we will draw from a distribution using the predefined mean velocityvalue. Otherwise, it will raise an error.
            
        Notes
        ----- 
        Mean velocities for terrestrial planets and the Moon are based on analysis of simulations of main-belt derived asteroids from Minton & Malhotra (2010) [1]_  and Yue et al. (2013) [2]_. Mean velocities for the asteroids are from Bottke et al. (1994) [3]_.
        
        References
        ----------
        .. [1] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744-757. https://doi.org/10.1016/j.icarus.2009.12.008
        .. [2] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
        .. [3] Bottke, W.F., Nolan, M.C., Greenberg, R., Kolvoord, R.A., 1994. Velocity distributions among colliding asteroids. Icarus 107, 255-268. https://doi.org/10.1006/icar.1994.1021
        """

        if value is None: 
            pmv = self.projectile_mean_velocity
            vencounter_mean = np.sqrt(pmv**2 - self.target.escape_velocity**2)
            vencounter = mc.get_random_velocity(vencounter_mean, rng=self.rng)
            pv = np.sqrt(vencounter**2 + self.target.escape_velocity**2)
            pang = mc.get_random_impact_angle(rng=self.rng)
            self._projectile_vertical_velocity = pv * np.sin(np.deg2rad(pang))
        elif isinstance(value, (int, float)):
            if value < 0:
                raise ValueError("projectile_vertical_velocity must be a positive number")
            self._projectile_vertical_velocity = float(value)
        else: 
            raise TypeError("projectile_vertical_velocity must be a numeric value or None") 

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
