from __future__ import annotations
from typing import Any
import numpy as np
from numpy.random import Generator
from cratermaker.utils.general_utils import parameter
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.custom_types import FloatLike
from cratermaker.utils import montecarlo as mc
from cratermaker.core.target import Target

class Impactor(ComponentBase):
    _registry: dict[str, Impactor] = {}
    def __init__(self, 
                 target : Target | str | None = None,
                 mean_velocity : FloatLike = 22100.0, 
                 density : FloatLike = 2250.0,
                 sample_velocities : bool = True,
                 sample_angles : bool = True,
                 sample_directions : bool = True,
                 angle: FloatLike = 90.0,
                 velocity : FloatLike | None = None,
                 direction : FloatLike | None = None,
                 rng: Generator | None = None,
                 rng_seed : int | None = None,
                 rng_state : dict | None = None, 
                 **kwargs):
        """
        This is the abstract base class for all impactor models. It defines the interface for generating impactor velocities, angles, and densities for a given target body.

        Parameters
        ----------
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        mean_velocity : float
            The mean velocity of the projectile in m/s. Default is 22100.0 m/s.
        density : float
            The density of the impactor in kg/m^3. Default is 2250.0 kg/m^3.
        sample_velocities : bool
            Flag that determines whether to sample impact velocities from a distribution. If set to False, impact velocities will be set to the mean velocity.
        sample_angles : bool
            Flag that determines whether to sample impact angles from a distribution. If set to False, impact angles will be set to 90 degrees (vertical impact).
        sample_directions : bool
            Flag that determines whether to sample impact directions from a distribution. If set to False, impact directions will be set to the mean velocity.
        angle : float
            The impact angle in degrees. Default is 90.0 degrees.
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
        object.__setattr__(self, "_sample_angles", None)
        object.__setattr__(self, "_sample_velocities", None)
        object.__setattr__(self, "_sample_directions", None)
        object.__setattr__(self, "_mean_velocity", mean_velocity)
        object.__setattr__(self, "_density", density)
        object.__setattr__(self, "_velocity", velocity)
        object.__setattr__(self, "_direction", direction)
        object.__setattr__(self, "_angle", angle)

        self.target = Target.maker(target, **kwargs)
        self.sample_velocities = sample_velocities
        self.sample_angles = sample_angles
        self.sample_directions = sample_directions


    @classmethod
    def maker(cls,
             impactor: Impactor | str | None = None, 
             target : Target | str | None = None,
             mean_velocity : FloatLike = 22100.0, 
             density : FloatLike = 2250.0,
             sample_velocities : bool = True,
             sample_angles : bool = True,
             sample_directions : bool = True,
             angle: FloatLike = 90.0,
             velocity : FloatLike | None = None,
             direction : FloatLike | None = None,
             rng: Generator | None = None,
             rng_seed : int | None = None,
             rng_state : dict | None = None, 
             **kwargs: Any) -> Impactor:
        """
        Initialize an impactor model based on the provided name or class.
        
        Parameters
        ----------
        impactor : Impactor or str
            The impactor model to initialize. Can be a class or a string representing the model name.
        target : Target or str.
            The name of the target body for the impact. Default is "Moon"
        mean_velocity : float
            The mean velocity of the projectile in m/s. Default is 22100.0 m/s.
        density : float
            The density of the impactor in kg/m^3. Default is 2250.0 kg/m^3.
        sample_velocities : bool
            Flag that determines whether to sample impact velocities from a distribution. If set to False, impact velocities will be set to the mean velocity.
        sample_angles : bool
            Flag that determines whether to sample impact angles from a distribution. If set to False, impact angles will be set to 90 degrees (vertical impact).
        sample_directions : bool
            Flag that determines whether to sample impact directions from a distribution. If set to False, impact directions will be set to the mean velocity.
        angle : float
            The impact angle in degrees. Default is 90.0 degrees.
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
        Impactor
            The initialized impactor model.
        
        Raises
        ------
        KeyError
            If the specified impactor model is not found in the registry.
        TypeError
            If the specified impactor model is not a string or a subclass of Impactor.
        """
        target = Target.maker(target, **kwargs)
        if impactor is None:
            target_name = target.name.capitalize()
            if target_name in ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Ceres', 'Vesta']:
                impactor = "asteroids"
            else:
                impactor = "comets"

        return super().maker(component=impactor,
                            target=target,
                            mean_velocity=mean_velocity,
                            density=density,
                            sample_velocities=sample_velocities,
                            sample_angles=sample_angles,
                            sample_directions=sample_directions,
                            angle=angle,
                            velocity=velocity,
                            direction=direction,
                            rng=rng,
                            rng_seed=rng_seed,
                            rng_state=rng_state, 
                            **kwargs) 

        
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
            A dictionary containing the impactor properties.
        """

        if self.sample_velocities:
            self._velocity = mc.get_random_velocity(self.mean_velocity, rng=self.rng)

        if self.sample_angles:
            self.angle = mc.get_random_impact_angle(rng=self.rng)

        if self.sample_directions:
            self._direction = mc.get_random_impact_direction(rng=self.rng)

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
        The target object for the impactor model.
        
        Returns
        -------
        Target
        """
        return self._target
    
    @target.setter
    def target(self, value):
        self._target = Target.maker(value)

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
    def vertical_velocity(self):
        """
        The vertical component of the impact velocity in m/s.
        
        Returns
        -------
        float 
        """
        return self.velocity * np.sin(np.radians(self.angle))
    

import_components(__name__, __path__, ignore_private=True)

