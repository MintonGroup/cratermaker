import numpy as np
from numpy.random import Generator
from typing import Type
from .target import Target
from ..utils.general_utils import validate_and_convert_location, float_like
from ..utils import montecarlo as mc
from .scale import Scale 
from .morphology import Morphology
class Crater:
    """
    Represents a crater formed by an impact in the simulation.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features.

    Attributes
    ----------
    diameter : float
        The diameter of the crater rim in m.
    radius : float
        The radius of the crater rim in m.
    transient_diameter : float
        The diameter of the transient crater in m.
    transient_radius : float
        The radius of the transient crater in m.    
    location : (2,) float
        The lat. and lon. of the impact point in degrees
    """

    def __init__(self, 
                diameter: float_like = None,
                radius: float_like = None,
                transient_diameter: float_like = None,
                transient_radius: float_like = None,
                location: np.ndarray = None,
                target: Target = None,
                scale: Scale = None,
                morphology: Morphology = None, 
                rng: Generator = None):
        
        if target is None:
            target = Target(name="Moon")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")
        
        if scale is None:
            self.scale = Scale(target,rng) 
        elif isinstance(scale, Scale):
            self.scale = scale
        else:
            raise TypeError("scale must be an instance of Scale") 
        
        if rng is None:
            self.rng = np.random.default_rng()
        elif isinstance(rng, Generator):
            self.rng = rng
        else:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        
        #  Evaluate and check diameter/radius values 
        values_set = sum(x is not None for x in [diameter, radius, transient_diameter, transient_radius])
        if values_set > 1:
            raise ValueError("Only one of diameter, radius, transient_diameter, transient_radius may be set")
        elif values_set == 0:
            raise ValueError("A crater must include one of diameter, radius, transient_diameter, or transient_radius!")
        
        # Now call the setters to ensure that all related calculations and checks are performed
        self._diameter = diameter
        self._radius = radius
        self._transient_diameter = transient_diameter
        self._transient_radius = transient_radius    
        self._location = location
        
        if diameter is not None:
            self.diameter = diameter
        elif radius is not None:
            self.radius = radius
        elif transient_diameter is not None:
            self.transient_diameter = transient_diameter
        elif transient_radius is not None:  
            self.transient_radius = transient_radius

        if location is not None:
            self.location = location              
        self._initialize_location(self.rng) 
        
        if morphology is None:
            self.morphology = Morphology(self,target,self.rng) 
        elif isinstance(morphology, Morphology):
            self.morphology = morphology
        else:
            raise TypeError("morphology must be an instance of Morphology")
            
        return

    
    def __repr__(self):
        return (f"Crater(diameter={self.diameter}, radius={self.radius}, "
                f"transient_diameter={self.transient_diameter}, transient_radius={self.transient_radius}, "
                f"morphology_type={self.morphology_type} "
                f"location={self.location}")
    
    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of crater rim must be finite and positive!")
            self._diameter = value
            self._radius = value / 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2
        return
    
    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of crater rim must be finite and positive!")            
            self._radius = value
            self._diameter = value * 2
            self._transient_diameter, self.morphology_type = self.scale.final_to_transient(value)
            self._transient_radius = self._transient_diameter / 2
            
    @property
    def transient_diameter(self):
        return self._transient_diameter

    @transient_diameter.setter
    def transient_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of transient crater must be finite and positive!")
            self._transient_diameter = value
            self._transient_radius = value / 2
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
        return
    
    @property
    def transient_radius(self):
        return self._transient_radius

    @transient_radius.setter
    def transient_radius(self, value):
        self._transient_radius = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of transient crater must be finite and positive!")            
            self._transient_diameter = value * 2            
            self._diameter,self.morphology_type = self.scale.transient_to_final(value)
            self._radius = self._diameter / 2
    
    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, value):
        self._location = value
    
    def _initialize_location(self, rng):
        if self._location is None:
            self.location = mc.get_random_location(rng=rng)
        else:    
            self.location = validate_and_convert_location(self.location)


class Projectile:
    """
    Represents the self.in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    Attributes
    ----------
    diameter : float
        The diameter of the projectile in m.
    radius : float
        The radius of the projectile in m.
    density : float
        The mass density of the projectile in kg/m**3.
    mass : float
        The mass of the projectile in kg.
    velocity : float
        The velocity of the projectile upon impact, in m/s.
    vertical_velocity : float
        The vertical component of the projectile velocity upon impact, in m/s.
    angle : float
        The angle of impact, in degrees.
    location : (2,) float
        The lat. and lon. of the impact point in degrees
    """
    
    def __init__(self, 
                diameter: float_like = None,
                radius: float_like = None,
                density: float_like = None,
                mass: float_like = None,
                velocity: float_like = None,
                angle: float_like = None,
                vertical_velocity: float_like = None,
                location: np.ndarray = None,
                target: Target = None, 
                scale = None,
                rng: Generator = None):
        from .scale import Scale 
        if target is None:
            target = Target(name="Moon")
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target")        
        if rng and not isinstance(rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        if scale is None:
            self.scale = Scale(target, rng) 
        elif isinstance(scale, Scale):
            self.scale = scale
        else:
            raise TypeError("scale must be an instance of Scale") 

        
        # Evaluate and check diameter/radius inputs
        values_set = sum(x is not None for x in [diameter, radius])
        if values_set > 1:
            raise ValueError("Only one of diameter or radius may be set")
        
        values_set = sum(x is not None for x in [diameter, radius, mass])
        if values_set == 0:
            raise ValueError("A projectile must include one of diameter, radius, or mass!")
        
        # Evaluate and check mass/density/radius inputs 
        values_set = sum(x is not None for x in [mass, density, radius])
        if values_set > 2:
            raise ValueError("Only two of mass, density, and radius may be set")
        
        # Evaluate and check mass/density/diameter inputs 
        values_set = sum(x is not None for x in [mass, density, diameter])
        if values_set > 2:
            raise ValueError("Only two of mass, density, and diameter may be set")
        
        # Evaluate and check velocity/impact angle inputs
        values_set = sum(x is not None for x in [velocity, vertical_velocity, angle])
        if values_set > 2:
            raise ValueError("Only two of velocity, vertical_velocity, angle may be set")
        
        # Now call the setters to ensure that all related calculations and checks are performed
        self._diameter = diameter
        self._radius = radius
        self._density = density
        self._mass = mass    
        self._velocity = velocity    
        self._angle = angle    
        self._vertical_velocity = vertical_velocity
        self._location = location
        
        if diameter is not None:
            self.diameter = diameter
        elif radius is not None:
            self.radius = radius

        if mass is not None:
            self.mass = mass
       
        if density is not None:
            self.density = density 
        else:
            if self.mass is None or self.radius is None: # Default to target density if we are given no way to figure it out
                self.density = target.material.density 
                
        if location is not None:
            self.location = location              
        
        if velocity is not None:
            self.velocity = velocity
        if angle is not None:
            self.angle = angle
        if vertical_velocity is not None:
            self.vertical_velocity = vertical_velocity 
            
        self._initialize_velocities(target,rng)
        self._initialize_location(rng)
        
        self.scale = scale 
        
        return

    def __repr__(self):
        return (f"Projectile(diameter={self.diameter} mm, radius={self.radius} m, "
                f"mass={self.mass} kg, density={self.density} kg/m^3, "
                f"velocity={self.velocity} m/s, angle={self.angle} deg, "
                f"vertical_velocity={self.vertical_velocity} m/s, "
                f"lon: {self.location[0]}, lat {self.location[1]}")


    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of projectile must be finite and positive!")
            self._radius = value / 2
            self._update_mass()
        return

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of projectile must be finite and positive!")            
            self._diameter = value * 2
            self._update_mass()
        return

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        if value is not None and value <= 0.0:
            raise ValueError("Mass of projectile must be finite and positive!")            
        self._update_volume_based_properties()
        return

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value
        if value is not None and value <= 0.0:
            raise ValueError("Density of projectile must be finite and positive!")                 
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
        if self._velocity is None and self._vertical_velocity is not None and self._angle is not None:
            self._velocity = self._vertical_velocity / np.sin(np.deg2rad(self._angle))
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Velocity of projectile must be finite and positive!")        
        self._velocity = value
            
        if value is not None and self._angle is not None:
            self._vertical_velocity = value * np.sin(np.deg2rad(self._angle))

    @property
    def vertical_velocity(self):
        if self._vertical_velocity is None and self._velocity is not None and self._angle is not None:
            self._vertical_velocity = self._velocity * np.sin(np.deg2rad(self._angle))
        return self._vertical_velocity

    @vertical_velocity.setter
    def vertical_velocity(self, value):
        if value is not None and value <= 0.0:
            raise ValueError("Vertical velocity of projectile must be finite and positive!")        
        self._vertical_velocity = value
        
        # Update velocity only if angle is already set
        if self._angle is not None and value is not None:
            try:
                self._velocity = value / np.sin(np.deg2rad(self._angle))
            except ValueError:
                raise ValueError(f"Invalid vertical velocity value {value} for a given angle value {self._angle}!")

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, value):
        if value is not None:
            if value < 0.0 or value > 90.0:
                raise ValueError("Impact angle of projectile must be between 0 and 90 degrees")
            self._angle = value
            # Update vertical_velocity only if velocity is already set
            if self._velocity is not None:
                self._vertical_velocity = self._velocity * np.sin(np.deg2rad(self._angle))

    def _initialize_velocities(self, target: Target, rng: Generator | None = None):
        if self._velocity is None:
            vencounter_mean = np.sqrt(target.mean_impact_velocity**2 - target.escape_velocity**2)
            vencounter = mc.get_random_velocity(vencounter_mean, rng=rng)
            self.velocity = np.sqrt(vencounter**2 + target.escape_velocity**2)

        if self._angle is None:
            if rng:
                self.angle = mc.get_random_impact_angle(rng=rng)
            else:
                self.angle = mc.get_random_impact_angle()
        return

    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, value):
        self._location = value
    
    
    def _initialize_location(self, rng):
        if self._location is None:
            self.location = mc.get_random_location(rng=rng)
        else:    
            self.location = validate_and_convert_location(self.location)

