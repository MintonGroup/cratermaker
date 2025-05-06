from __future__ import annotations
from typing import Any
import numpy as np
import math
from numpy.random import Generator
from cratermaker.utils.general_utils import parameter, format_large_units, _set_properties
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils import montecarlo_utils as mc
from cratermaker.components.target import Target

class Projectile(ComponentBase):
    _registry: dict[str, Projectile] = {}
    def __init__(self, 
                 target : Target | str | None = None,
                 mean_velocity : FloatLike | None = None,
                 density : FloatLike | None = None,
                 sample : bool | None = None,
                 angle: FloatLike | None = None,
                 velocity : FloatLike | None = None,
                 direction : FloatLike | None = None,
                 rng: Generator | None = None,
                 rng_seed : int | None = None,
                 rng_state : dict | None = None, 
                 **kwargs):
        """
        This is the abstract base class for all projectile models. It defines the interface for generating projectile velocities, angles, and densities for a given target body.

        Parameters
        ----------
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        mean_velocity : float
            The mean velocity of the projectile in m/s. 
        density : float
            The density of the projectile in kg/m^3. 
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.
        angle : float
            The impact angle in degrees. 
        velocity : float | None
            The impact velocity in m/s. If None, the velocity will be sampled from a distribution.
        direction : float | None
            The impact direction in degrees. If None, the direction will be sampled from a distribution.
        rng : numpy.random.Generator | None
            A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
        rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
            The rng_rng_seed for the RNG. If None, a new RNG is created.
        rng_state : dict, optional
            The state of the random number generator. If None, a new state is created.
        **kwargs : Any
            Additional keyword arguments.
        """
        super().__init__(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_sample", sample)
        object.__setattr__(self, "_mean_velocity", mean_velocity)
        object.__setattr__(self, "_density", density)
        object.__setattr__(self, "_velocity", velocity)
        object.__setattr__(self, "_direction", direction)
        object.__setattr__(self, "_angle", angle)
        if self.sample is None:
            self.sample = True
        self.target = Target.maker(target, **kwargs)

    def __repr__(self) -> str:
        base = super().__repr__()
        mean_velocity = format_large_units(self.mean_velocity, quantity="velocity")
        return (
            f"{base}\n"
            f"Target: {self.target.name}\n"
            f"Density: {self.density:.1f} kg/m³\n"
            f"Mean Velocity: {mean_velocity}\n"
            f"Sample from distributions: {self.sample}"
        )


    @classmethod
    def maker(cls,
             projectile: Projectile | str | None = None, 
             target : Target | str | None = None,
             mean_velocity : FloatLike | None = None,
             density : FloatLike | None = None,
             sample : bool = True,
             angle: FloatLike | None = None,
             velocity : FloatLike | None = None,
             direction : FloatLike | None = None,
             rng: Generator | None = None,
             rng_seed : int | None = None,
             rng_state : dict | None = None, 
             **kwargs: Any) -> Projectile:
        """
        Initialize an projectile model based on the provided name or class.
        
        Parameters
        ----------
        projectile : Projectile or str
            The projectile model to initialize. Can be a class or a string representing the model name.
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        mean_velocity : float
            The mean velocity of the projectile in m/s. 
        density : float
            The density of the projectile in kg/m^3. 
        sample : bool
            Flag that determines whether to sample impact velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.
        angle : float
            The impact angle in degrees. 
        velocity : float | None
            The impact velocity in m/s. If None, the velocity will be sampled from a distribution.
        direction : float | None
            The impact direction in degrees. If None, the direction will be sampled from a distribution.
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
        Projectile
            The initialized projectile model.
        
        Raises
        ------
        KeyError
            If the specified projectile model is not found in the registry.
        TypeError
            If the specified projectile model is not a string or a subclass of Projectile.
        """
        target = Target.maker(target, **kwargs)
        if projectile is None:
            target_name = target.name.capitalize()
            if target_name in ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Ceres', 'Vesta']:
                projectile = "asteroids"
            else:
                projectile = "comets"

        projectile = super().maker(component=projectile,
                            target=target,
                            mean_velocity=mean_velocity,
                            density=density,
                            sample=sample,
                            angle=angle,
                            velocity=velocity,
                            direction=direction,
                            rng=rng,
                            rng_seed=rng_seed,
                            rng_state=rng_state, 
                            **kwargs) 
        projectile.new_projectile(**kwargs)
        return projectile

        
    def new_projectile(self, **kwargs: Any) -> dict:
        """
        Updates the values of the velocities and angles and returns them as a dictionary of projectile properties that can be passed as arguments to the Crater class.
        
        Parameters
        ----------
        **kwargs : Any
            Additional keyword arguments to be passed to internal functions.
        
        Returns
        -------
        dict
            A dictionary containing the projectile properties.
        """

        if self.sample:
            self._velocity = float(mc.get_random_velocity(self.mean_velocity, rng=self.rng)[0])
        elif self.velocity is None:
            self._velocity = float(self.mean_velocity)

        if self.sample:
            self._angle = float(mc.get_random_impact_angle(rng=self.rng)[0])
        elif self.angle is None:
            self._angle = float(90.0)

        if self.sample:
            self._direction = float(mc.get_random_impact_direction(rng=self.rng)[0])
        elif self._direction is None:
            self._direction = float(0.0)

        return {
            "projectile_velocity": self.velocity,
            "projectile_angle": self.angle,
            "projectile_density": self.density,
            "projectile_direction": self.direction, 
        }

    @parameter
    def target_name(self):
        """
        The name of the target body.
        
        Returns
        -------
        str
        """ 
        return self._target.name
    
    @property
    def target(self):
        """
        The target object for the projectile model.
        
        Returns
        -------
        Target
        """
        return self._target
    
    @target.setter
    def target(self, value):
        self._target = Target.maker(value)

    @parameter
    def sample(self):
        """
        Flag that determines whether to sample velocities, angles, and directions from distributions. If set to False, impact velocities will be set to the mean velocity, impact angles will be set to 90 degrees (vertical impact), and directions will be 0.
        
        Returns
        -------
        bool
        """
        return self._sample
    
    @sample.setter
    def sample(self, value):
        if not isinstance(value, bool):
            raise TypeError("sample must be a boolean value")
        self._sample = value
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
        if isinstance(value, np.ndarray):
            value = value.item()
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
        return self._angle
    
    @angle.setter
    def angle(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("angle must be a numeric value")
        if value < 0:
            raise ValueError("angle must be a positive number")
        self._angle = float(value)

    @property
    def direction(self):
        """
        The impact direction in degrees.
        
        Returns
        -------
        float 
        """
        return self._direction
    
    @direction.setter
    def direction(self, value):
        if isinstance(value, np.ndarray):
            value = value.item()
        if not isinstance(value, FloatLike):
            raise TypeError("direction must be a numeric value")
        if value < 0:
            raise ValueError("direction must be a positive number")
        self._direction = float(value)

    @property
    def velocity(self):
        """
        The impact velocity in m/s.
        
        Returns
        -------
        float 
        """
        return self._velocity
    
    @velocity.setter
    def velocity(self, value):
        if value is not None:
            if not isinstance(value, FloatLike):
                raise TypeError("velocity must be a numeric value")
            if value < 0:
                raise ValueError("velocity must be a positive number")
            self._velocity = float(value)
        else:
            self._velocity = None
        return

    @parameter
    def density(self):
        """
        The density of the projectile in kg/m^3.
        
        Returns
        -------
        float 
        """
        return self._density
    
    @density.setter
    def density(self, value):
        """
        Sets the density of the projectile in kg/m^3. 
        
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
    def vertical_velocity(self):
        """
        The vertical component of the impact velocity in m/s.
        
        Returns
        -------
        float 
        """
        return self.velocity * math.sin(math.radians(self.angle))
    
    @property
    def population(self):
        """
        The name of the population of the projectile model.
        
        Returns
        -------
        int
        """
        return self._component_name

import_components(__name__, __path__, ignore_private=True)

