from . import util
from . import montecarlo as mc
from  ._bind import _BodyBind
from dataclasses import dataclass, field
import numpy as np
from numpy.random import default_rng
import jigsawpy
import os
import trimesh
import json
from typing import Union, Tuple, List


@dataclass    
class Projectile:
    """
    Represents the projectile in the crater simulation.

    This class defines the properties of the impacting object, such as its size,
    velocity, material, and angle of impact.

    Attributes
    ----------
    velocity : float
        The velocity of the projectile upon impact, in m/s.
    angle : float
        The angle of impact, in degrees.
    material : Material
        The material composition of the projectile. 
    """
    diameter: float = None              # The diameter of the projectile (m)
    velocity: float = None              # The velocity of the projectile (m/s)
    density: float  = None              # The mass density of the projectile (kg/m**3)
    location: (float,float) = None      # Tuple that specifies a location of the impact onto the target surface: (lat,lon)? (theta,phi)? some other measure of location?
    angle: float = None          # The impact angle of the projectile (deg)
        

@dataclass
class Crater:
    """
    Represents a crater formed by an impact in the simulation.

    This class models the crater resulting from an impact, including its size,
    shape, depth, and other morphological features.

    Attributes
    ----------
    TBD
    """
    location: np.ndarray = field(default=None)
    diameter: np.float64 = field(default=None)
    radius: np.float64 = field(default=None)
    transient_diameter: np.float64 = field(default=None)
    transient_radius: np.float64 = field(default=None)
    
    def __post_init__(self):
        values_set = sum(x is not None for x in [self.diameter, self.radius, 
                                                 self.transient_diameter, self.transient_radius])

        if values_set > 1:
            raise ValueError("Only one of diameter, radius, transient_diameter, transient_radius may be set")

        if self.diameter is not None:
            self.radius = self.diameter / 2
        elif self.radius is not None:
            self.diameter = self.radius * 2
        elif self.transient_diameter is not None:
            self.transient_radius = self.transient_diameter / 2
        elif self.transient_radius is not None:
            self.transient_diameter = self.transient_radius * 2

        if self.location is not None:
            if not isinstance(self.location, np.ndarray):
                self.location = np.array(self.location, dtype=np.float64)
            if self.location.shape != (2,):
                raise ValueError("location must be a 2-element array")
        
        return
   

@dataclass
class Material:
    """
    Represents the material properties relevant to the crater simulation.

    This class defines various physical properties of the material involved in the cratering process.
    

    Attributes
    ----------
    name : str
        The name of the material. If the material is matched to one that is present in the catalogue, the rest of the properties will be retrieved for it unless specified. If the name is not known from the catalogue, then all other properties must be supplied and in order to build a custom material.
    Ybar : float
        The strength of the material, typically defined in Pa. 
    other_properties : dict
        Other relevant properties of the material.

    Methods
    -------
    set_properties(name, **kwargs):
        Add a custom property to the material.

    """

    # Define all valid properties for the Target object
    name: str = None
    K1: float = None
    mu: float = None
    Ybar: float = None
    density: float = None 

    config_ignore = ['catalogue']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        # Define some default crater scaling relationship terms (see Richardson 2009, Table 1, and Kraus et al. 2011 for Ice) 
        material_properties = [
            "name",       "K1",     "mu",   "Ybar",     "density" 
        ]
        material_values = [
            ("Water",     2.30,     0.55,   0.0,        1000.0),
            ("Sand",      0.24,     0.41,   0.0,        1750.0),
            ("Dry Soil",  0.24,     0.41,   0.18e6,     1500.0),
            ("Wet Soil",  0.20,     0.55,   1.14e6,     2000.0),
            ("Soft Rock", 0.20,     0.55,   7.60e6,     2250.0),
            ("Hard Rock", 0.20,     0.55,   18.0e6,     2500.0),
            ("Ice",       15.625,   0.48,   0.0,        900.0), 
        ]        
        
        self.catalogue = util._create_catalogue(material_properties, material_values)
        
        # Set properties for the Material object based on the catalogue value)
        if self.name:
            self.set_properties(catalogue=self.catalogue, key=self.name)
        else:
            raise ValueError('No material defined!')    
        
        return    
    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util._set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """         
        util._set_properties(self,**kwargs)
        return

@dataclass
class Target:
    """
    Represents the target body in the crater simulation.

    This class encapsulates the properties of the target that is impacted, including
    its material composition, size, and other relevant physical characteristics.

    Attributes
    ----------
    material : Material
        The material composition of the target.
    size : float
        The size of the target, in relevant units. 
    """
       
    # Set up instance variables
    name: str = None
    radius: float = None
    gravity: float = None
    material_name: str = None
    material: Material = None
    mean_impact_velocity: float = None
    transition_scale_type: str = "silicate" # Options are silicate and ice
    
    config_ignore = ['catalogue','material']  # Instance variables to ignore when saving to file
    def __post_init__(self):
        # Define some built-in catalogue values for known solar system targets of interest
        gEarth = 9.80665 # 1 g in SI units
        
        body_properties = [
            "name",    "radius",   "gravity",      "material_name", "mean_impact_velocity"
        ]
        body_values = [
            ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft Rock", 41100.0),
            ("Venus",   6051.84e3, 0.905 * gEarth, "Hard Rock", 29100.0),
            ("Earth",   6371.01e3, 1.0   * gEarth, "Wet Soil" , 24600.0),
            ("Moon",    1737.53e3, 0.1657* gEarth, "Soft Rock", 22100.0),
            ("Mars",    3389.92e3, 0.379 * gEarth, "Soft Rock", 10700.0),
            ("Ceres",   469.7e3,   0.029 * gEarth, "Ice"      , 5300.0),
            ("Vesta",   262.7e3,   0.025 * gEarth, "Soft Rock", 5300.0),
        ]      
        # Mean velocities for terrestrial planets based on analysis of simulations from Minton & Malhotra (2010) of main belt-derived asteroid
        # Mean velocities for the asteroids are from Bottke et al. (1994)
        
        self.catalogue = util._create_catalogue(body_properties, body_values)
        
        # Set properties for the Target object based on the arguments passed to the function
        if self.name:
            self.material = "TEMP" 
            self.set_properties(catalogue=self.catalogue, key=self.name)
            self.material = Material(name=self.material_name)
        else: 
            raise ValueError('No target defined!')    
        
        return
   
    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util._set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """         
        util._set_properties(self,**kwargs)
        return
    
class Simulation():
    """
    This is a class that defines the basic Cratermaker body object. 
    """
    def __init__(self, target_name="Moon", material_name=None, **kwargs):
        if material_name:
            material = Material(name=material_name)
            self.target = Target(name=target_name, material=material, **kwargs) 
        else: 
            self.target = Target(name=target_name, **kwargs)
       
        # Set some default values for the simulation parameters
        self.pix = kwargs.get('pix', self.target.radius / 1e3)
        self.time_function = kwargs.get('time_function', None)
        self.tstart = kwargs.get('tstart', 0.0)  # Simulation start time (in y)
        self.tstop = kwargs.get('tstop', 4.31e9)    # Simulation stop time (in y)
        
        # Set some default values for the production function population
        self.impactor_sfd  = kwargs.get('impactor_sfd', None)
        self.impactor_velocity = kwargs.get('impactor_velocity', None)
        
        # Set the random number generator seed
        self.seed = kwargs.get('seed', None) 
        self.rng = default_rng(seed=self.seed)
        
        self.cachedir = os.path.join(os.getcwd(),'.cache')
        self.mesh_file = kwargs.get('mesh_file', os.path.join(self.cachedir,"target_mesh.glb") )
        if not os.path.exists(self.cachedir):
            os.mkdir(self.cachedir)
       
        self.mesh = None 
        # Check if a mesh exists, and if so, load it up
        if os.path.exists(self.mesh_file):
            self.load_target_mesh()
        else:
            self.make_target_mesh()
    
   
    @property
    def crater(self):
        return self._crater

    @crater.setter
    def crater(self, value):
        self._crater = value
        if hasattr(value, 'diameter') and value.diameter is not None:
            self.final_to_transient()
            self.crater.transient_radius = self.crater.transient_diameter / 2
        elif hasattr(value, 'transient_diameter') and value.transient_diameter is not None:
            self.transient_to_final()
            self.crater.radius = self.crater.diameter / 2
                     
    def populate(self):
        """
        Populate the surface with craters
        """
        
        if not self.mesh:
            print("Generating new mesh. Please be patient.")
            self.make_target_mesh
            
        
    config_ignore = ['target', 'projectile', 'crater']  # Instance variables to ignore when saving to file
    def to_json(self, filename):
        
        # Get the simulation configuration into the correct structure
        material_config = util._to_config(self.target.material)
        target_config = {**util._to_config(self.target), 'material' : material_config}
        sim_config = {**util._to_config(self),'target' : target_config} 
        
        # Write the combined configuration to a JSON file
        with open(filename, 'w') as f:
            json.dump(sim_config, f, indent=4)
        
    @property
    def elevation(self):
        return self._body.fobj.elevation
    
    @property
    def name(self):
        return self._body.fobj.name

    def get_elevation(self):
        return self._body.get_elevation()

    
    def set_elevation(self, elevation_array):
        self._body.set_elevation(elevation_array)


    def make_target_mesh(self):
        """
        Generate a tessellated mesh of a sphere using the jigsawpy library and convert it to GLB format.

        This function sets up Jigsaw mesh files and a mesh body, defines a basic sphere using an ellipsoid mesh model,
        sets mesh options, generates a tessellated mesh, and converts the mesh to a trimesh object. The generated mesh 
        is then saved in GLB format.

        Returns
        -------
        None
            The function does not return a value but updates the `self.mesh` attribute with the generated trimesh object.
        """
        
        # This will get updated evantually after we're done testing
        # Set up jigsaw objects
        opts = jigsawpy.jigsaw_jig_t()
        geom = jigsawpy.jigsaw_msh_t()
        mesh = jigsawpy.jigsaw_msh_t()
        
        # Set up Jigsaw mesh files and mesh body if a file doesn't already exist
        opts.mesh_file = self.mesh_file.replace(".glb",".msh")
        opts.jcfg_file = self.mesh_file.replace(".glb",".jig")
        opts.geom_file = self.mesh_file.replace(".glb","_geom.msh")
        
        # Define basic sphere using the built-in ellipsoid mesh model
        geom.mshID = "ellipsoid-mesh"
        geom.radii = np.full(3, self.target.radius, dtype=geom.REALS_t)
        jigsawpy.savemsh(opts.geom_file,geom)

        # Set mesh options
        #opts.numthread = int(os.getenv('OMP_NUM_THREADS'))
        opts.hfun_scal = "absolute"  #scaling type for mesh-size function, either "relative" or "absolute." For "relative" mesh-size values as percentages of the (mean) length of the axis-aligned bounding-box (AABB) associated with the geometry. "absolute" interprets mesh-size values as absolute measures.
        opts.hfun_hmax = np.sqrt(2.0) * self.pix       # mesh-size function value. The "pix" value is adjusted to keep it scaled  to roughly the same as the older CTEM pix
        opts.mesh_dims = +2          # Number of mesh dimensions (2 for body mesh, 3 for volume)
        opts.optm_qlim = +9.5E-01    # threshold on mesh cost function above which gradient-based optimisation is attempted.
        opts.optm_iter = +64         # max. number of mesh optimisation iterations. 
        opts.optm_qtol = +1.0E-05    # tolerance on mesh cost function for convergence. Iteration on a given node is terminated if adjacent element cost-functions are improved by less than QTOL.
        opts.mesh_kern = "delfront"  # meshing kernel, choice of the standard Delaunay-refinement algorithm ('delaunay') or the Frontal-Delaunay method ('delfront').
        # Generate tesselated mesh
        jigsawpy.cmd.jigsaw(opts, mesh)
        
        # Convert the mesh to trimesh
        # Ensure the vertex and face data is in numpy array format
        vertices = np.array([vert[0] for vert in mesh.vert3])
        faces = np.array([face[0] for face in mesh.tria3])

        # Create a trimesh object
        self.mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
       
        # Save generated mesh 
        export = self.mesh.export(self.mesh_file)
        
        return
    
    
    def load_target_mesh(self):
        """
        Load a target mesh from a file into the `self.mesh` attribute.

        This function uses the trimesh library to read a mesh file, which is expected to be in GLB format. The function
        handles the peculiarity of trimesh reading the mesh as a Scene instead of a Trimesh object and extracts the 
        actual mesh from the scene.

        Returns
        -------
        None
            The function does not return a value but updates the `self.mesh` attribute with the loaded trimesh object.
        """
    
        # This is not well documented, but trimesh reads the mesh in as a Scene instead of a Trimesh, so we need to extract the actual mesh
        scene = trimesh.load_mesh(self.mesh_file)
        
        # Assuming that the only geometry in the file is the body mesh, this should work
        self.mesh = next(iter(scene.geometry.values()))
        
        return
    
    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `util._set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """        
        util._set_properties(self,**kwargs)
        return 

    
    def get_simple_to_complex_transition_factors(self):
        simple_enlargement_factor = 0.84 # See Melosh (1989) pg. 129 just following eq. 8.2.1
        transition_exp_mean = 0.15 # See Croft (1985) eq. 9
        final_exp_mean = 0.85      # See Croft (1985) eq. 9
        # Standard deviations for the exponents
        transition_exp_std = 0.04
        final_exp_std = 0.04
       
        # The transition values come from CTEM and are a synthesis of Pike (1980), Croft (1985), Schenk et al. (2004).
        if self.target.transition_scale_type == "silicate":
            simple_complex_exp = -1.0303 
            simple_complex_mean = 2*16533.8 
            simple_complex_std = 0.05
        elif self.target.transition_scale_type == "ice":
            simple_complex_exp = -1.22486
            simple_complex_mean = 2*3081.39
            simple_complex_std = 0.05
           
        # Draw from a normal distribution for each exponent
        transition_exp = self.rng.normal(transition_exp_mean, transition_exp_std)
        final_exp = self.rng.normal(final_exp_mean, final_exp_std)
        simple_complex_fac = simple_complex_mean * np.exp(self.rng.normal(scale=simple_complex_std))
        transition_diameter = simple_complex_fac * self.target.gravity**simple_complex_exp
        
        return transition_diameter, simple_enlargement_factor, transition_exp, final_exp 
     

    def final_to_transient(self):
        transition_diameter, simple_enlargement_factor, transition_exp, final_exp = self.get_simple_to_complex_transition_factors() 
        if self.crater.diameter < transition_diameter:
            self.crater.transient_diameter = simple_enlargement_factor * self.crater.diameter  # Simple crater scaling
        else:
            self.crater.transient_diameter = 1e3 * (transition_diameter * 1e-3)**transition_exp * (self.crater.diameter * 1e-3)**final_exp # Complex crater scaling (in m)
        self.crater.transient_diameter = np.float64(self.crater.transient_diameter)
        
        return 

    def transient_to_final(self):
        transition_diameter, simple_enlargement_factor, transition_exp, final_exp = self.get_simple_to_complex_transition_factors()
        
        final_diameter_simple = np.float64(self.crater.transient_diameter / simple_enlargement_factor)
        final_diameter_complex = np.float64(1e3 * ((1e-3 * self.crater.transient_diameter) / (transition_diameter * 1e-3)**transition_exp)**(1.0 / final_exp))
       
        if final_diameter_simple < transition_diameter and final_diameter_complex < transition_diameter: # Unambiguosly simple
            self.crater.diameter = final_diameter_simple
            self.crater.morphology_type = "simple"
        elif final_diameter_simple > transition_diameter and final_diameter_complex > transition_diameter: # Unambiguously complex
            self.crater.diameter = final_diameter_complex
            self.crater.morphology_type = "complex"
        else: # Could be either complex or simple. We'll just draw which one from a hat weighted in the direction of whichever size is closest to the transition
            is_simple = self.rng.random() < np.abs(final_diameter_complex - transition_diameter) / np.abs(final_diameter_simple - final_diameter_complex)
            self.crater.morphology_type = self.rng.choice(["simple", "complex"])
            if is_simple:
                self.crater.diameter = final_diameter_simple
                self.crater.morphology_type = "simple"
            else:
                self.crater.diameter = final_diameter_complex
                self.crater.morphology_type = "complex"
            
        return    


 