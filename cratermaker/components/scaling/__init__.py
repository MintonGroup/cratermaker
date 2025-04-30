from __future__ import annotations
from abc import  abstractmethod
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.custom_types import FloatLike
from cratermaker.core.target import Target
from cratermaker.components.projectile import Projectile
from cratermaker.utils.component_utils import ComponentBase, import_components

class Scaling(ComponentBase):
    _registry: dict[str, Scaling] = {}

    """
    This is the abstract base class for all scaling models. It defines the interface for converting between projectile and crater diameters.

    Parameters
    ----------
    target : Target | str, default="Moon"
        The target body for the impact. Can be a Target object or a string representing the target name.
    projectile : Projectile | str, default="asteroids"
        The projectile model for the impact. Can be an Projectile object or a string representing the projectile name.
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
                 target: Target | str | None = None,
                 projectile: Projectile | str | None = None,
                 rng : Generator | None = None,
                 rng_seed: int | None = None,
                 rng_state: dict | None = None,
                 **kwargs: Any):
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_target_density", None)
        object.__setattr__(self, "_projectile", None)
        # combine the kwargs with the common_args, giving common_args priority
        kwargs = {**kwargs, **vars(self.common_args)}
        self._target = Target.maker(target, **kwargs)
        self._projectile = Projectile.maker(projectile, **kwargs)

    @abstractmethod
    def projectile_to_transient(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_projectile(self, **kwargs: Any) -> np.float64: ...
    @abstractmethod
    def transient_to_final(self, transient_diameter: FloatLike) -> tuple[np.float64, str]: ...
    @abstractmethod
    def final_to_transient(self, final_diameter: FloatLike, morphology_type: str | None = None, **kwargs) -> np.float64: ...

    @classmethod
    def maker(cls,
             scaling: str | Scaling | None = None, 
             target: Target | str | None = None,
             projectile: Projectile | str | None = None,
             rng : Generator | None = None,
             rng_seed: int | None = None,
             rng_state: dict | None = None,
             **kwargs: Any) -> Scaling:
        """
        Initialize a scaling model based on the provided name or class.

        Parameters
        ----------
        scaling : str, Scaling, or None, default=None
            The name of the scaling model to initialize. If None, the default model is used.
        target : Target | str, default="Moon"
            The target body for the impact. Can be a Target object or a string representing the target name.
        projectile : Projectile | str, default="asteroids"
            The projectile model for the impact. Can be an Projectile object or a string representing the projectile name.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        Scaling
            An instance of the specified scaling model.

        Raises
        ------
        KeyError
            If the specified scaling model name is not found in the registry.
        TypeError
            If the specified scaling model is not a string or a subclass of Scaling.
        """

        if scaling is None:
            scaling = "richardson2009"
        return super().maker(component=scaling, target=target, projectile=projectile, rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    def __repr__(self) -> str:
        return (
            f"<Scaling Model: {self.model}\n"
            f"Target: {self.target}\n"
            f"Projectile: {self.projectile}>"
        )

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
        self._target = Target.maker(value, **vars(self.common_args))
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
    def projectile(self):
        """
        The projectile model for the impact.
        
        Returns
        -------
        Projectile
        """ 
        return self._projectile
    
    @projectile.setter
    def projectile(self, value):
        self._projectile = Projectile.maker(value, **vars(self.common_args))
        return
    
    @property
    def model(self):
        """
        The name of the scaling model.
        
        Returns
        -------
        str
        """ 
        return self._component_name

import_components(__name__, __path__, ignore_private=True)

