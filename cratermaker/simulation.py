from ._bind import _BodyBind
import numpy as np
import jigsawpy
import os
import trimesh
from numpy.random import default_rng


class Target:
    # Define some built-in catalogue values for known solar system targets of interest
    gEarth = 9.80665 # 1 g in SI units
   
    body_properties = [
        "name",    "radius",   "gravity",      "material"
    ]
    body_values = [
        ("Mercury", 2440.0e3,  0.377 * gEarth, "Soft rock"),
        ("Venus",   6051.84e3, 0.905 * gEarth, "Hard rock"),
        ("Earth",   6371.01e3, 1.0   * gEarth, "Wet soil"),
        ("Moon",    1737.53e3, 0.1657* gEarth, "Soft rock"),
        ("Mars",    3389.92e3, 0.379 * gEarth, "Soft rock"),
        ("Ceres",   469.7e3,   0.29  * gEarth, "Ice"     ),
        ("Vesta",   262.7e3,   0.25  * gEarth, "Soft rock"),
    ]
   
    # Define some default crater scaling relationship terms (see Richardson 2009, Table 1) 
    material_properties = [
        "name",       "K1",     "mu",   "Ybar",     "density" 
    ]
    material_values = [
        ("Water",     2.30,     0.55,   0.0,        1000.0),
        ("Sand",      0.24,     0.41,   0.0,        1750.0),
        ("Dry Soil",  0.24,     0.41,   0.18,       1500.0),
        ("Wet Soil",  0.20,     0.55,   1.14,       2000.0),
        ("Soft Rock", 0.20,     0.55,   7.60,       2250.0),
        ("Hard Rock", 0.20,     0.55,   18.0,       2500.0),
        ("Ice",       2.30,     0.39,   0.0,        900.0), # TODO: Update these based on Kraus, Senft, and Stewart (2011) 
    ]

    def __init__(self, **kwargs):
        # Define all valid properties for the Target object
        self.name = None
        self.radius = None
        self.gravity = None 
        self.material = 'Soft Rock'
        self.K1 = None
        self.mu = None
        self.Ybar = None
        self.crustal_density = None
        
        # Set properties for the Target object based on the arguments passed to the function
        self.set_properties(**kwargs)
        
        return
    
    @staticmethod
    def create_body_catalogue(self):
        # Create the catalogue dictionary using the class variables
        body_catalogue = {
            body[0]: dict(zip(self.body_properties, body))
            for body in self.body_values
        }

        # Remove the 'name' key from each dictionary in the catalogue
        for body_name in list(body_catalogue):
            del body_catalogue[body_name]['name']

        return body_catalogue
    
    
    @staticmethod
    def create_material_catalogue(self):
        # Create the catalogue dictionary using the class variables
        material_catalogue = {
            material[0]: dict(zip(self.material_properties, material))
            for material in self.material_values
        }

        # Remove the 'name' key from each dictionary in the catalogue
        for material_name in list(material_catalogue):
            del material_catalogue[material_name]['name']

        return material_catalogue 
    
   
    def set_properties(self, **kwargs):
        
        # If the "name" variable is passed, look the name up in the catalogue and set properties to those from the catalogue
        self.name = kwargs.get('name', None)
        if self.name is not None:
            self.name = self.name.title()
        
        body_properties = self.body_catalogue.get(self.name)
        if body_properties: # The name was found in the catalogue
            print(f"{self.name} was found in the catalogue of known bodies.")
            self.radius = body_properties.get('radius', self.radius)
            self.gravity = body_properties.get('gravity', self.gravity)
            self.material = body_properties.get('material', self.material)
        elif self.name is not None:
            print(f"{self.name} was not found in the catalogue. Setting custom properties for this target body")
        else:
            raise ValueError("A name must be supplied to generate a custom target body!")
        
        # Override body properties with input arguments if provided 
        self.radius = kwargs.get('radius', self.radius)  
        self.gravity = kwargs.get('gravity', self.gravity)
        self.material = kwargs.get('material',self.material)
        
        # Look up material in the catalogue
        if self.material is not None:
            self.material = self.material.title()
        material_properties = self.material_catalogue.get(self.material)
        if material_properties: # The material was found in the catalogue
            print(f"{self.material} was found in the catalogue of known materials.")
            self.K1 = material_properties.get('K1', self.K1)
            self.mu = material_properties.get('mu', self.mu)
            self.Ybar = material_properties.get('Ybar', self.Ybar)
            self.crustal_density = material_properties.get('density', self.crustal_density)
        elif self.material is not None:
            print(f"{self.material} was not found in the catalogue. Setting custom properties for this material")
        else:
            raise ValueError("A valid material must must be supplied, or a name must be supplied to generate a custom target material!")
           
        # Override material properties with input arguments if provided
        self.K1 = kwargs.get('K1', self.K1) 
        self.mu = kwargs.get('mu', self.mu) 
        self.Ybar = kwargs.get('Ybar', self.Ybar) 
        self.crustal_density = kwargs.get('crustal_density', self.crustal_density)

        # Check for any unset properties
        for property_name, value in self.__dict__.items():
            if value is None:
                raise ValueError(f"The property {property_name} has not been set!")
        return
        
 
Target.body_catalogue = Target.create_body_catalogue(Target)
Target.material_catalogue = Target.create_material_catalogue(Target)
class Projectile:
    def __init__(self):
        # Initialize the projectile's attributes
        self.sfd   = None
        self.radius = None
        self.diameter = None
        self.velocity = None
        self.sin_impact_angle = None
        self.vertical_velocity = None
        
        return

class Crater:
    def __init__(self):
        # Initialize the crater's attributes
        self.diameter = None
        self.radius = None
        self.morphotype = None
        return

class Simulation(object):
    """
    This is a class that defines the basic Cratermaker body object. 
    """
    def __init__(self, name="Moon", **kwargs):
        kwargs['name'] = name
        self.target = Target(**kwargs)
        # Set default configuration options
        # TODO: Initialize with configure options as arguments or read in configuration file
        self.config = {
            "body" : "Moon", # Which body to simulation (Options are "Moon","Custom" )
            "pix"  : 6.16e3, # Approximate cell global cell size to make a body with 1e6 faces 
            "cachedir" : os.path.join(os.getcwd(),".cache"), # Directory location of output files
            "meshname" : "body", # Name of files generated during mesh process
            "seed" : 5421109845, # RNG seed
        }
        
        # Determine the radius of the target body
        #if self.config['body'] != 'Custom':
        #    self.config['body_radius'] = body_radius_values[self.config['body']]
           
        # if not os.path.exists(self.config['cachedir']):
        #     os.mkdir(self.config['cachedir'])
            
        # if os.path.exists(os.path.join(self.config['cachedir'],f"{self.config['meshname']}.glb")):
        #     self.load_body_mesh()
        # else:
        #     self.make_body_mesh()
               
        # self.rng = default_rng(seed=self.config['seed'])
        # # self._body = _BodyBind(gridshape)     
        # # lat = np.linspace(-90, 90, gridshape[0])
        # # lon = np.linspace(0, 360, gridshape[1])
        # # self.ds = xr.Dataset(
        # #     {'elevation': (['lat', 'lon'], self.get_elevation())},
        # #                 coords={'lat': lat, 'lon': lon}
        # )
    
    def set_target(self, body_name):
        self.target = Target(body_name)
        
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

    def make_body_mesh(self):
        '''
            This will use jigsawpy to tesselate a sphere to make our intial mesh. This will then be converted to a GLB format
        '''
        
        # This will get updated evantually after we're done testing
        # Set up jigsaw objects
        opts = jigsawpy.jigsaw_jig_t()
        geom = jigsawpy.jigsaw_msh_t()
        mesh = jigsawpy.jigsaw_msh_t()
        
        # Set up Jigsaw mesh files and mesh body if a file doesn't already exist
        meshname=self.config['meshname']
        cachedir=self.config['cachedir']  
        opts.jcfg_file = os.path.join(cachedir,f"{meshname}.jig")
        opts.mesh_file = os.path.join(cachedir,f"{meshname}.msh")
        opts.geom_file = os.path.join(cachedir,f"{meshname}_base.msh")
        
        # Define basic sphere using the built-in ellipsoid mesh model
        geom.mshID = "ellipsoid-mesh"
        geom.radii = np.full(3, self.config['body_radius'], dtype=geom.REALS_t)
        jigsawpy.savemsh(opts.geom_file,geom)

        # Set mesh options
        opts.hfun_scal = "absolute"  #scaling type for mesh-size function, either "relative" or "absolute." For "relative" mesh-size values as percentages of the (mean) length of the axis-aligned bounding-box (AABB) associated with the geometry. "absolute" interprets mesh-size values as absolute measures.
        opts.hfun_hmax = np.sqrt(2.0) * self.config['pix']       # mesh-size function value. The "pix" value is adjusted to keep it scaled  to roughly the same as the older CTEM pix
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
        export = self.mesh.export(os.path.join(cachedir,f"{meshname}.glb"))
        
        return
    
    
    def load_body_mesh(self):
        
        # Set up Jigsaw mesh files and mesh body if a file doesn't already exist
        meshname=self.config['meshname']
        cachedir=self.config['cachedir']  
          
        # This is not well documented, but trimesh reads the mesh in as a Scene instead of a Trimesh, so we need to extract the actual mesh
        scene = trimesh.load_mesh(os.path.join(cachedir,f"{meshname}.glb"))
        
        # Assuming that the only geometry in the file is the body mesh, this should work
        self.mesh = next(iter(scene.geometry.values()))
        
        return
       