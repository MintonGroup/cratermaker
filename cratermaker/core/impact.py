import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike
from typing import Type, Any
from .target import Target
from .scale import Scale 
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..utils.custom_types import FloatLike, PairOfFloats
from .morphology import Morphology
from .production import Production
from abc import ABC, abstractmethod

class Impact(ABC):
    """
    Abstract base class representing an impact event in the simulation. 

    This class provides a generic framework for impact-related entities, defining common attributes and methods such as diameter,
    radius, and location. It is not intended to be instantiated directly but serves as a base for specific impact event classes like 
    Crater and Projectile.


    Parameters
    ----------
    diameter : FloatLike, optional
        The diameter of the impact event in meters.
    radius : FloatLike, optional
        The radius of the impact event in meters.
    location : PairOfFloats, optional
        The location (latitude, longitude) of the impact event.
    target : Target, optional
        The target body of the impact event.
    scale_cls : Type[Scale], optional
        The class used for scaling calculations.
    rng : Generator, optional
        Random number generator instance.

    Attributes
    ----------
    _diameter : FloatLike
        Diameter of the impact event.
    _radius : FloatLike
        Radius of the impact event.
    _location : PairOfFloats
        Location of the impact event.

    Abstract Methods
    ----------------
    diameter
        Property that must be implemented to handle the diameter of the impact event.
    radius
        Property that must be implemented to handle the radius of the impact event.

    Notes
    -----
    Subclasses must provide implementations for all abstract methods and properties. 
    They may also extend the __init__ method to handle additional parameters specific 
    to the subclass.

    Examples
    --------
    Subclassing `Impact` to create a `Crater` class:

    >>> class Crater(Impact):
    ...     def __init__(self, diameter, ...):
    ...         super().__init__(diameter=diameter, ...)
    ...     # Implement abstract methods and properties here

    See Also
    --------
    Crater, Projectile : Subclasses of `Impact`.
    """
    
    def __init__(self, 
                 diameter: FloatLike = None,
                 radius: FloatLike = None,
                 location: PairOfFloats = None,
                 target: Target = None,
                 scale_cls: Type[Scale] = None,
                 rng: Generator = None,
                 **kwargs: Any):
        
        # Evaluate and check diameter/radius inputs
        self._target = None
        self._rng = None
        self._scale_cls = None
        self._scale = None
        self._location = None
        self._diameter = None
        self._radius = None
        self._face_index = None
        self._node_index = None
        
        self.target = target
        self.rng = rng 
        self.scale_cls = scale_cls 
        self.scale = self.scale_cls(target=self.target, rng=self.rng)
        self.location = location 
        self.diameter = diameter
        self.radius = radius
        
        values_set = sum(x is not None for x in [diameter, radius])
        if values_set > 1:
            raise ValueError("Only one of diameter or radius may be set") 

    @property
    @abstractmethod
    def diameter(self):
        """
        The diameter of the the impact object in m.
        
        Returns
        -------
        np.float64
        """         
        pass

    @property
    @abstractmethod
    def radius(self):
        """
        The radius of the the impact object in m.
        
        Returns
        -------
        np.float64
        """             
        pass

    @property
    def location(self):
        """
        The longitude and latitude of the impact in degrees. 
       
        Returns
        -------
        np.ndarray
        """
        return self._location

    @location.setter
    def location(self, value):
        if value is None:
            self._location = mc.get_random_location(rng=self.rng)
        else:    
            self._location = validate_and_convert_location(value)
            
            
    @property
    def face_index(self):
        """
        The index of the face closest to the impact location.
        
        Returns
        -------
        int
        """
        return self._face_index
    
    @face_index.setter
    def face_index(self, value):
        if value is None:
            self._face_index = None
        else:
            self._face_index = int(value)
        return
    
    @property
    def node_index(self):
        """
        The index of the node closest to the impact location.
        
        Returns
        -------
        int
        """
        return self._node_index
    
    @node_index.setter
    def node_index(self, value):
        if value is None:
            self._node_index = None
        else:
            self._node_index = int(value)
        return

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
    def scale_cls(self):
        """
        The class to use for computing crater scaling relationships.
        
        Returns
        -------
        Type[Scale]
        """ 
        return self._scale_cls

    @scale_cls.setter
    def scale_cls(self, cls):
        if cls is None:
            self._scale_cls = Scale
            return
        if not issubclass(cls, Scale):
            raise ValueError("The class must be a subclass of scale")
        self._scale_cls = cls
        return
        
    @property
    def scale(self):
        """
        An object that defines parameters relavant to scaling the crater size and projectile->crater scaling.
        
        Returns
        -------
        Scale
        """ 
        return self._scale
    
    @scale.setter
    def scale(self, value):
        if value is not None and not isinstance(value, Scale):
            raise TypeError("scale must be an instance of Scale")
        self._scale = value
        return    

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


class Crater(Impact):
    """
    Represents a crater formed by an impact in the simulation.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features.

    """

    def __init__(self, 
                transient_diameter: FloatLike = None,
                transient_radius: FloatLike = None,
                morphology_cls: Type[Morphology] = None, 
                **kwargs: Any):
        """
        Constructor for the Crater class.
        
        Parameters
        ----------
        transient_diameter : FloatLike
            The diameter of the transient crater in m.
        transient_radius : FloatLIke
            The radius of the transient crater in m.    
        morphology_cls : Type[Morphology]
            The class to use for computing the crater morphology.
        **kwargs : Any 
            Additional keyword arguments to pass to the scale and morphology classes.
        """ 
        self._transient_diameter = None
        self._transient_radius = None
        super().__init__(**kwargs)
        
        #  Evaluate and check diameter/radius values 
        values_set = sum(x is not None for x in [self.diameter, transient_diameter, transient_radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius, transient_diameter, transient_radius may be set")
        elif values_set == 0:
            raise ValueError("A crater must include one of diameter, radius, transient_diameter, or transient_radius!")                     
        
        # Now call the setters to ensure that all related calculations and checks are performed
        self.transient_diameter = transient_diameter
        self.transient_radius = transient_radius    
        
        self.morphology_cls = morphology_cls
        self.morphology = self.morphology_cls(crater=self, target=self.target, rng=self.rng)
            
        return

    
    def __repr__(self):
        return (f"Crater(diameter={self.diameter}, radius={self.radius}, "
                f"transient_diameter={self.transient_diameter}, transient_radius={self.transient_radius}, "
                f"morphology_type={self.morphology_type} "
                f"location={self.location}")
        
    @property
    def diameter(self):
        """
        The diameter of the crater rim in m. Setting diameter automatically sets radius, transient_diameter, and transient_radius.
        
        Returns
        -------
        np.float64
        """
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of crater rim must be finite and positive!")
            self._diameter = np.float64(value)
            self._radius = np.float64(value) / 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2
    
    @property
    def radius(self):
        """
        The radius of the crater rim in m. Setting radius automatically sets diameter, transient_diameter, and transient_radius.
        
        Returns
        -------
        np.float64
        """        
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of crater rim must be finite and positive!")            
            self._radius = np.float64(value)
            self._diameter = np.float64(value) * 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2    
            
    @property
    def transient_diameter(self):
        """
        The diameter of the transient crater in m. Setting transint diameter automatically sets transient radius, diameter, and radius.
        
        Returns
        -------
        np.float64
        """        
        return self._transient_diameter

    @transient_diameter.setter
    def transient_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of transient crater must be finite and positive!")
            self._transient_diameter = np.float64(value)
            self._transient_radius = np.float64(value) / 2
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
        return
    
    @property
    def transient_radius(self):
        """
        The radius the transient crater in m. Setting transint radius automatically sets transient diameter, diameter, and radius.
        
        Returns
        -------
        np.float64
        """            
        return self._transient_radius

    @transient_radius.setter
    def transient_radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of transient crater must be finite and positive!")            
            self._transient_radius = np.float64(value)
            self._transient_diameter = np.float64(value) * 2
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
    

    @property
    def morphology_cls(self):
        """
        The class to use for computing the crater morphology.
        
        Returns
        -------
        Type[Morphology]
        """ 
        return self._morphology_cls

    @morphology_cls.setter
    def morphology_cls(self, cls):
        if cls is None:
            self._morphology_cls = Morphology
        elif not issubclass(cls, Morphology):
            raise ValueError("The class must be a subclass of Morphology")
        else:
            self._morphology_cls = cls
        
    @property
    def morphology(self):
        """
        An object that defines parameters relavant to describing the morphology of the crater.
        
        Returns
        -------
        Morphology
        """ 
        return self._morphology
    
    @morphology.setter
    def morphology(self, value):
        if not isinstance(value, Morphology):
            raise TypeError("morphology must be an instance of Morphology")
        self._morphology = value
        return 


class Projectile(Impact):
    """
    Represents the projectile body in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    """
    
    def __init__(self, 
                mass: FloatLike = None,
                density: FloatLike = None,
                velocity: FloatLike = None,
                vertical_velocity: FloatLike = None,
                angle: FloatLike = None,
                **kwargs):
        """
        
        Constructor for the Projectile class.
        
        Parameters
        ----------
        mass : float
            The mass of the projectile in kg.
        density : float
            The mass density of the projectile in kg/m**3.            
        velocity : float
            The velocity of the projectile upon impact, in m/s.
        vertical_velocity : float
            The vertical component of the projectile velocity upon impact, in m/s.
        angle : float
            The angle of impact, in degrees.
        **kwargs : Any
        
        Notes
        -----
        The Projectile class is initialized with at least one of diameter, radius, or mass, but only two of mass, density, or
        radius or two of mass, density, and diameter, and only two of velocity, vertical_velocity, or angle.  If density cannot be
        determined from the inputs, the density of the target body is used as a default. If the velocity cannot be determined from 
        the inputs, the mean impact velocity of the target is used as an input to the get_random_velocity function. If the angle 
        cannot be determined from the inputs, a random angle is generated using the get_random_impact_angle function.
        """ 

        # ensure that all related calculations and checks are performed
        self._density = density
        self._mass = mass
        
        self._velocity = None
        self._vertical_velocity = None
        self._angle = None
        
        super().__init__(**kwargs)
        values_set = sum(x is not None for x in [self.diameter, self.radius, mass])
        # Evaluate and check mass/density/radius inputs 
        if values_set == 0:
            raise ValueError("A projectile must include one of diameter, radius, or mass!")
        values_set = sum(x is not None for x in [mass, density, self.radius])
        if values_set > 2:
            raise ValueError("Only two of mass, density, and radius may be set")
        values_set = sum(x is not None for x in [mass, density, self.diameter])
        if values_set > 2:
            raise ValueError("Only two of mass, density, and diameter may be set")        
        
        if self.mass is None or self.radius is None: # Default to target density if we are given no way to figure it out
            self.density = self.target.material.density 
        
        # Evaluate and check velocity/impact angle inputs
        values_set = sum(x is not None for x in [velocity, vertical_velocity, angle])
        if values_set > 2:
            raise ValueError("Only two of velocity, vertical_velocity, angle may be set")        
        self._velocity = velocity    
        self._angle = angle    
        self._vertical_velocity = vertical_velocity
        self._initialize_velocities(self.target,self.rng)
        
        return

    def __repr__(self):
        return (f"Projectile(diameter={self.diameter} mm, radius={self.radius} m, "
                f"mass={self.mass} kg, density={self.density} kg/m^3, "
                f"velocity={self.velocity} m/s, angle={self.angle} deg, "
                f"vertical_velocity={self.vertical_velocity} m/s, "
                f"lon: {self.location[0]}, lat {self.location[1]}")


    @property
    def diameter(self):
        """
        The diameter of the projectile in m. Setting the diameter automatically sets the radius. Mass and density are updated accordingly.
        
        Returns
        -------
        np.float64 
        """       
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of projectile must be finite and positive!")
            self._diameter = np.float64(value)
            self._radius = np.float64(value) / 2
            self._update_mass()
        return

    @property
    def radius(self):
        """
        The radius of the projectile in m. Setting the radius automatically sets the diameter. Mass and density are updated accordingly.
        
        Returns
        -------
        np.float64 
        """
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of projectile must be finite and positive!")            
            self._radius = np.float64(value)
            self._diameter = np.float64(value) * 2
            self._update_mass()
        return

    @property
    def mass(self):
        """
        The mass of the projectile in kg. Setting the mass automatically updates the diameter and radius from a given density.
        
        Returns
        -------
        np.float64 
        """
        return self._mass

    @mass.setter
    def mass(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Mass of projectile must be finite and positive!")            
            self._mass = np.float64(value)
            self._update_volume_based_properties()
        return

    @property
    def density(self):
        """
        The mass density of the projectile in kg/m**3. Setting the density automatically updates the mass from a given radius.
        
        Returns
        -------
        np.float64 
        """
        return self._density

    @density.setter
    def density(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Density of projectile must be finite and positive!")                 
            self._density = np.float64(value)
            self._update_mass()

    def _update_mass(self):
        if self._density is not None and self._radius is not None:
            self._mass = 4.0/3.0 * np.pi * self._radius**3 * self._density

    def _update_volume_based_properties(self):
        if self._mass is not None and self._density is not None:
            volume = self._mass / self._density
            self._radius = ((3.0 * volume) / (4.0 * np.pi))**(1.0/3.0)
            self._diameter = self._radius * 2

    @property
    def velocity(self):
        """
        The velocity of the projectile upon impact, in m/s. Setting the velocity automatically updates the vertical velocity from a given angle.
        
        Returns
        -------
        np.float64 
        """
        if self._velocity is None and self._vertical_velocity is not None and self._angle is not None:
            self._velocity = self._vertical_velocity / np.sin(np.deg2rad(self._angle))
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Velocity of projectile must be finite and positive!")        
        self._velocity = np.float64(value)
            
        if value is not None and self._angle is not None:
            self._vertical_velocity = np.float64(value) * np.sin(np.deg2rad(self._angle))

    @property
    def vertical_velocity(self):
        """
        The vertical component of the projectile velocity upon impact, in m/s. Setting the vertical velocity automatically updates the velocity from a given angle. 
        The vertical velocity must be less than or equal to the velocity.
        
        Returns
        -------
        np.float64 
        """
        if self._vertical_velocity is None and self._velocity is not None and self._angle is not None:
            self._vertical_velocity = self._velocity * np.sin(np.deg2rad(self._angle))
        return self._vertical_velocity

    @vertical_velocity.setter
    def vertical_velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Vertical velocity of projectile must be finite and positive!")        
        self._vertical_velocity = np.float64(value)
        
        # Update velocity only if angle is already set
        if self._angle is not None and value is not None:
            try:
                self._velocity = np.float64(value) / np.sin(np.deg2rad(self._angle))
            except ValueError:
                raise ValueError(f"Invalid vertical velocity value {value} for a given angle value {self._angle}!")

    @property
    def angle(self):
        """
        The angle of impact, in degrees. Setting the angle automatically updates the vertical velocity from a given velocity.
        
        Returns
        -------
        np.float64
        """
        return self._angle

    @angle.setter
    def angle(self, value):
        if value is not None:
            if value < 0.0 or value > 90.0:
                raise ValueError("Impact angle of projectile must be between 0 and 90 degrees")
            self._angle = np.float64(value)
            # Update vertical_velocity only if velocity is already set
            if self._velocity is not None:
                self._vertical_velocity = self._velocity * np.sin(np.deg2rad(self._angle))

    def _initialize_velocities(self, target: Target, rng: Generator | None = None):
        if self._velocity is None:
            if target.mean_impact_velocity is None:
                raise ValueError("No mean impact velocity is defined for target {target.name}. Projectiles will not be computed.")
        
            # Check if the impact velocity is greater than the escape velocity. If so, be sure to add the escape velocity to the mean impact velocity in quadrature
            # Otherwise, just use the mean velocity value
            if target.mean_impact_velocity > target.escape_velocity:
                vencounter_mean = np.sqrt(target.mean_impact_velocity**2 - target.escape_velocity**2)
                vencounter = mc.get_random_velocity(vencounter_mean, rng=rng)
                self.velocity = np.sqrt(vencounter**2 + target.escape_velocity**2)
            else: 
                # Be sure not to generate velocities above escape in this scenario
                while True:  
                    self.velocity = mc.get_random_velocity(target.mean_impact_velocity, rng=rng)
                    if self.velocity < target.escape_velocity:
                        break

        if self._angle is None:
            if rng:
                self.angle = mc.get_random_impact_angle(rng=rng)
            else:
                self.angle = mc.get_random_impact_angle()
        return

