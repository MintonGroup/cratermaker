import numpy as np
from numpy.random import Generator
from typing import Type, Any
from .target import Target
from .surface import Surface
from ..utils.general_utils import validate_and_convert_location
from ..utils import montecarlo as mc
from ..utils.custom_types import FloatLike, PairOfFloats
from ..components.scaling import ScalingModel, get_scaling_model
from ..components.morphology import MorphologyModel, get_morphology_model

class Crater:
    """
    Represents a crater formed by an impact in the simulation, including both the target and projectile properties.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features, as well as the properties of the projectile that formed it.

    Parameters
    ----------
    final_diameter : FloatLike, optional
        The final diameter of the crater rim in meters.
    final_radius : FloatLike, optional
        The final radius of the crater rim in meters.
    transient_diameter : FloatLike, optional
        The diameter of the transient crater in m.
    transient_radius : FloatLike, optional
        The radius of the transient crater in m.
    projectile_diameter : float, optional
        The diameter of the projectile in m.
    projectile_radius : float, optional
        The radius of the projectile in m.
    projectile_mass : float, optional
        The mass of the projectile in kg.
    projectile_density : float, optional
        The mass density of the projectile in kg/m**3. If not set, the density of the target is used.
    projectile_mean_velocity : float, optional
        The mean impact velocity of the projectile population in m/s. 
    projectile_velocity : float, optional
        The velocity of the projectile upon impact, in m/s.
    projectile_vertical_velocity : float, optional
        The vertical component of the projectile velocity upon impact, in m/s.
    projectile_angle : float, optional
        The angle of impact, in degrees.
    projectile_direction: float, optional
        The direction of impact in degrees relative to north (the bearing on the surface)
    location : PairOfFloats, optional
        The location (latitude, longitude) of the impact event.
    age : FloatLike, optional
        The age of the impact event in My before present.
    target : Target, optional
        The target body of the impact event.
    scale : ScalingModel, optional
        The model used to compute impactor to crater scaling.
    morphology_cls : Type[MorphologyModel], optional
        The class to use for computing the crater morphology.
    rng : Generator, optional
        Random number generator instance.
    **kwargs : Any 
        Additional keyword arguments to pass to the scale and morphology classes.
    """
    __slots__ = (
        "_target",
        "_rng",
        "_scale",
        "_morphology_cls",
        "_morphology_type",
        "_morphology",
        "_location",
        "_age",
        "_final_diameter",
        "_transient_diameter",
        "_face_index",
        "_node_index",
        "_simple_enlargement_factor",
        "_final_exp",
        "_morphology",
        "_projectile_diameter",
        "_projectile_mass",
        "_projectile_density",
        "_projectile_velocity",
        "_projectile_vertical_velocity",
        "_projectile_mean_velocity",
        "_projectile_angle",
        "_projectile_direction")

    def __init__(self,
                 final_diameter: FloatLike = None,
                 final_radius: FloatLike = None,
                 transient_diameter: FloatLike = None,
                 transient_radius: FloatLike = None,
                 projectile_diameter: FloatLike = None,
                 projectile_radius: FloatLike = None,
                 projectile_mass: FloatLike = None,
                 projectile_density: FloatLike = None,
                 projectile_mean_velocity: FloatLike = None,
                 projectile_velocity: FloatLike = None,
                 projectile_vertical_velocity: FloatLike = None,
                 projectile_angle: FloatLike = None,
                 projectile_direction: FloatLike = None,
                 location: PairOfFloats = None,
                 age: FloatLike = None,
                 target: Target = None,
                 scale: ScalingModel | None = None,
                 morphology_cls: Type[MorphologyModel] = None,
                 rng: Generator = None,
                 **kwargs: Any):
        
        self.rng = rng
        self.target = target
        self.scale = scale
        self.morphology_cls = morphology_cls
        self.location = location
        self.age = age

        # validate and set projectile properties
        if projectile_mean_velocity is not None:
            if projectile_velocity is not None or projectile_vertical_velocity is not None:
                raise ValueError("projectile_mean_velocity cannot be used with projectile_velocity or projectile_vertical_velocity")
        else:
            v_values_set = sum(x is not None for x in [projectile_velocity, projectile_vertical_velocity, projectile_angle])
            if v_values_set > 2:
                raise ValueError("Only two of projectile_velocity, projectile_vertical_velocity, projectile_angle may be set")

        self.projectile_mean_velocity = projectile_mean_velocity
        self.projectile_velocity = projectile_velocity
        self.projectile_vertical_velocity = projectile_vertical_velocity
        self.projectile_angle = projectile_angle 
        self.projectile_direction = projectile_direction
        self._initialize_projectile_velocities()

        # Validate and set crater size properties
        values_set = sum(x is not None for x in [final_diameter, final_radius, transient_diameter, transient_radius, projectile_radius, projectile_diameter, projectile_mass])
        if values_set > 1:
            raise ValueError("Only one of final_diameter, final_radius, transient_diameter, transient_radius, projectile_diameter, projectile_radius, or projectile_mass may be set")
        elif values_set == 0:
            raise ValueError("A crater must include one of final_diameter, final_radius, transient_diameter, transient_radius, projectile_diameter, projectile_radius, or projectile_mass.")

        if projectile_density is None:
            self._projectile_density = self.scale.target_density
        else:
            self._projectile_density = projectile_density        

        self.projectile_mass = projectile_mass
        self.final_diameter = final_diameter
        self.final_radius = final_radius
        self.transient_diameter = transient_diameter
        self.transient_radius = transient_radius

        self.morphology = self.morphology_cls(crater=self, target=self.target, rng=self.rng, **kwargs)

    @property
    def final_diameter(self):
        """
        The diameter of the crater rim in m. Setting diameter automatically sets radius, transient_diameter, and transient_radius.
        """
        return self._final_diameter

    @final_diameter.setter
    def final_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of crater rim must be finite and positive!")
            self._final_diameter = np.float64(value)
            self._final_radius = np.float64(value) / 2
            self._transient_diameter, self._morphology_type = self.scale.final_to_transient(value)
            self._projectile_diameter = self.scale.transient_to_projectile(self._transient_diameter)
            self._update_projectile_mass()

    @property
    def final_radius(self):
        """
        The final radius of the crater rim in m. Setting radius automatically sets diameter, transient_diameter, and transient_radius.
        """
        return self._final_diameter / 2

    @final_radius.setter
    def final_radius(self, value):
        if value is not None:
            self.final_diameter = np.float64(value) * 2

    @property
    def location(self):
        """
        The longitude and latitude of the impact in degrees.
        """
        return self._location

    @location.setter
    def location(self, value):
        if value is None:
            self._location = mc.get_random_location(rng=self.rng)
        else:
            if len(value) == 1:
                value = value.item()
            self._location = validate_and_convert_location(value)

    @property
    def age(self):
        """
        The age of the impact in My before present.
        """
        return self._age

    @age.setter
    def age(self, value):
        if value is None:
            self._age = None
        else:
            self._age = np.float64(value)
        return

    @property
    def face_index(self):
        """The index of the face closest to the impact location."""
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
        """The index of the node closest to the impact location."""
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
        """The target body for the impact."""
        return self._target

    @target.setter
    def target(self, value):
        if value is None:
            self._target = Target(name="Moon")
        elif isinstance(value, str):
            try:
                self._target = Target(name=value)
            except:
                raise ValueError(f"Target '{value}' not found. Please provide a valid target name.")
        elif isinstance(value, Target):
            self._target = value
        else:
            raise TypeError("target must be an instance of Target, a string, or None")
        return

    @property
    def scale(self):
        """An object that defines parameters relevant to scaling the crater size and projectile->crater scaling."""
        return self._scale

    @scale.setter
    def scale(self, value):
        if value is None:
            self._scale = get_scaling_model("richardson2009")(target=self.target, rng=self.rng)
        elif isinstance(value, ScalingModel):
            self._scale = value
        else:
            raise TypeError("scale must be an instance of ScalingModel")
        return

    @property
    def rng(self):
        """A random number generator instance."""
        return self._rng

    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        self._rng = value or np.random.default_rng()

    @property
    def transient_diameter(self):
        """
        The diameter of the transient crater in m. Setting transient diameter automatically sets transient radius, diameter, and radius.
        """
        return self._transient_diameter

    @transient_diameter.setter
    def transient_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of transient crater must be finite and positive!")
            self._transient_diameter = np.float64(value)
            self._final_diameter, self._morphology_type = self.scale.transient_to_final(value)
            self._projectile_diameter = self.scale.transient_to_projectile(self._transient_diameter)
            self._update_projectile_mass()
        return

    @property
    def transient_radius(self):
        """
        The radius the transient crater in m. Setting transient radius automatically sets transient diameter, diameter, and radius.
        """
        return self._transient_diameter / 2

    @transient_radius.setter
    def transient_radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of transient crater must be finite and positive!")
            self.transient_diameter = np.float64(value) * 2

    def emplace(self, surf: Surface):
        """
        Emplace the crater on the surface.
        """
        self.node_index, self.face_index = surf.find_nearest_index(self.location)
        self.morphology.form_crater(surf)
        return

    @property
    def morphology_cls(self):
        """The class to use for computing the crater morphology."""
        return self._morphology_cls

    @morphology_cls.setter
    def morphology_cls(self, cls: Type[MorphologyModel] = None):
        if cls is None:
            self._morphology_cls = get_morphology_model("simplemoon")
        elif not issubclass(cls, MorphologyModel):
            raise ValueError("The class must be a subclass of Morphology")
        else:
            self._morphology_cls = cls

    @property
    def morphology(self):
        """An object that defines parameters relevant to describing the morphology of the crater."""
        return self._morphology

    @morphology.setter
    def morphology(self, value):
        if not isinstance(value, MorphologyModel):
            raise TypeError("morphology must be an instance of Morphology")
        self._morphology = value
        return

    @property
    def morphology_type(self):
        """The type of morphology to use for the crater."""
        return self._morphology_type

    @property
    def projectile_mass(self):
        """
        The mass of the projectile in kg. Setting the mass automatically updates the diameter and radius from a given density.
        """
        return self._projectile_mass

    @projectile_mass.setter
    def projectile_mass(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Mass of projectile must be finite and positive!")
            self._projectile_mass = np.float64(value)
            self._update_projectile_volume_based_properties()
        return

    @property
    def projectile_density(self):
        """
        The mass density of the projectile in kg/m**3. Setting the density automatically updates the mass from a given radius.
        """
        return self._projectile_density

    @projectile_density.setter
    def projectile_density(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Density of projectile must be finite and positive!")
            self._projectile_density = np.float64(value)
            self._update_projectile_mass()

    @property
    def projectile_radius(self):
        """
        The radius of the projectile in m.
        """
        # If mass and density are set, compute radius if not present
        if self._projectile_diameter is not None:
            return self._projectile_diameter / 2
        else:
            return None

    @projectile_radius.setter
    def projectile_radius(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Radius of projectile must be finite and positive!")
            self.projectile_diameter = 2 * np.float64(value)

    @property
    def projectile_diameter(self):
        """
        The diameter of the projectile in m.
        """
        return self._projectile_diameter

    @projectile_diameter.setter
    def projectile_diameter(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Diameter of projectile must be finite and positive!")
            self._projectile_diameter = np.float64(value)
            self._projectile_radius = np.float64(value) / 2
            self._update_projectile_mass()
            self._transient_diameter = self.scale.projectile_to_transient(self._projectile_diameter)
            self._transient_radius = self._transient_diameter / 2
            self._final_diameter, self._morphology_type = self.scale.transient_to_final(self._transient_diameter)


    def _update_projectile_mass(self):
        if self._projectile_density is not None and self.projectile_radius is not None:
            self._projectile_mass = 4.0 / 3.0 * np.pi * self._projectile_radius ** 3 * self._projectile_density

    def _update_projectile_volume_based_properties(self):
        if self._projectile_mass is not None and self._projectile_density is not None:
            volume = self._projectile_mass / self._projectile_density
            self.projectile_diameter = 2 * ((3.0 * volume) / (4.0 * np.pi)) ** (1.0 / 3.0)

    @property
    def projectile_mean_velocity(self):
        """
        The mean impact velocity of the projectile population in m/s.
        """
        return self._projectile_mean_velocity

    @projectile_mean_velocity.setter
    def projectile_mean_velocity(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Mean velocity of projectile must be finite and positive!")
            self._projectile_velocity = None
            self._projectile_vertical_velocity = None
            self._projectile_mean_velocity = np.float64(value)
            self._initialize_projectile_velocities()
        return

    @property
    def projectile_velocity(self):
        """
        The velocity of the projectile upon impact, in m/s.
        """
        if self._projectile_velocity is None and self._projectile_vertical_velocity is not None and self._projectile_angle is not None:
            self._projectile_velocity = np.float64(self._projectile_vertical_velocity / np.sin(np.deg2rad(self._projectile_angle)))
        return self._projectile_velocity

    @projectile_velocity.setter
    def projectile_velocity(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Velocity of projectile must be finite and positive!")
            self._projectile_velocity = np.float64(value)
            if value is not None and self._projectile_angle is not None:
                self._projectile_vertical_velocity = np.float64(value) * np.sin(np.deg2rad(self._projectile_angle))
            self._initialize_projectile_velocities()

    @property
    def projectile_vertical_velocity(self):
        """
        The vertical component of the projectile velocity upon impact, in m/s.
        """
        if self._projectile_vertical_velocity is None and self._projectile_velocity is not None and self._projectile_angle is not None:
            self._projectile_vertical_velocity = np.float64(self._projectile_velocity * np.sin(np.deg2rad(self._projectile_angle)))
        return self._projectile_vertical_velocity

    @projectile_vertical_velocity.setter
    def projectile_vertical_velocity(self, value):
        if value is not None:
            if value <= 0.0:
                raise ValueError("Vertical velocity of projectile must be finite and positive!")
            self._projectile_vertical_velocity = np.float64(value)
            if self._projectile_angle is not None and value is not None:
                try:
                    self._projectile_velocity = np.float64(value) / np.sin(np.deg2rad(self._projectile_angle))
                except ValueError:
                    raise ValueError(f"Invalid vertical velocity value {value} for a given angle value {self._projectile_angle}!")
            self._initialize_projectile_velocities()

    @property
    def projectile_angle(self):
        """
        The angle of impact, in degrees.
        """
        return self._projectile_angle

    @projectile_angle.setter
    def projectile_angle(self, value):
        if value is not None:
            if value < 0.0 or value > 90.0:
                raise ValueError("Impact angle of projectile must be between 0 and 90 degrees")
            self._projectile_angle = np.float64(value)
            if self._projectile_velocity is not None:
                self._projectile_vertical_velocity = self._projectile_velocity * np.sin(np.deg2rad(self._projectile_angle))
            self._initialize_projectile_velocities()

    def _initialize_projectile_velocities(self):
        """
        Initialize the impact velocity and angle of the projectile.
        """
        if self._projectile_velocity is None:
            if self._projectile_mean_velocity is None:
                raise ValueError("No projectile_mean_velocity is defined. Projectiles will not be computed.")
            if self._projectile_mean_velocity > self.target.escape_velocity: # Primary impact
                vencounter_mean = np.sqrt(self._projectile_mean_velocity ** 2 - self.target.escape_velocity ** 2)
                vencounter = mc.get_random_velocity(vencounter_mean, rng=self.rng)
                self._projectile_velocity = np.sqrt(vencounter ** 2 + self.target.escape_velocity ** 2, dtype=np.float64)
            else: # Secondary impact
                while True:
                    self._projectile_velocity = mc.get_random_velocity(self._projectile_mean_velocity, rng=self.rng)
                    if self._projectile_velocity < self.target.escape_velocity:
                        break
        if self._projectile_angle is None:
            self._projectile_angle = mc.get_random_impact_angle(rng=self.rng)
        if self._projectile_direction is None:
            self._projectile_direction = np.float64(self.rng.uniform(0.0, 360.0))
        return

    @property
    def projectile_direction(self):
        """
        The direction of impact in degrees relative to north (the bearing on the surface)
        """
        return self._projectile_direction

    @projectile_direction.setter
    def projectile_direction(self, value):
        if value is not None:
            if value < 0.0 or value >= 360.0:
                raise ValueError("Direction of impact must be between 0 and 360 degrees")
            self._projectile_direction = np.float64(value)
        return

    def __repr__(self):
        return (f"Crater(final_diameter={self.final_diameter} m, " 
                f"transient_diameter={self.transient_diameter} m, "
                f"morphology_type={self.morphology_type} "
                f"projectile_diameter={self.projectile_diameter} m, "
                f"projectile_mass={self.projectile_mass} kg, projectile_density={self.projectile_density} kg/m^3, "
                f"projectile_velocity={self.projectile_velocity} m/s, projectile_angle={self.projectile_angle} deg, "
                f"projectile_vertical_velocity={self.projectile_vertical_velocity} m/s, projectile_direction={self.projectile_direction} deg, "
                f"lon: {self.location[0]}, lat {self.location[1]} age={self.age} My)")