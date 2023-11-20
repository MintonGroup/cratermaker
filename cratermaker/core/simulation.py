from numpy.random import default_rng
import os
import json
from .target import Target
from .material import Material
from .projectile import Projectile
from .crater import Crater 
from .mesh import Mesh, load_target_mesh
from ..utils import general_utils 
from ..models.crater_scaling import get_simple_to_complex_transition_factors
from ..utils.montecarlo import get_random_location


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
        mesh_file = kwargs.get('mesh_file', os.path.join(self.cachedir,"target_mesh.glb") )
        if not os.path.exists(self.cachedir):
            os.mkdir(self.cachedir)
        if os.path.exists(mesh_file):
            self.mesh = load_target_mesh(mesh_file)
        else:
            self.mesh = Mesh(mesh_file,self.target,self.pix) 
        self._crater = None
        self._projectile = None

        return
    
    @property
    def crater(self):
        return self._crater
    
    @crater.setter
    def crater(self, value):
        if isinstance(value, Crater):
            self._crater = value
        else:
            self._crater = Crater(target=self.target, rng=self.rng, **value)
            
    def add_crater(self, **kwargs):
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        crater = Crater(**kwargs)
        if crater.diameter is not None:
            crater.radius = crater.diameter / 2
            crater.final_to_transient(self.target,self.rng)
        elif crater.radius is not None:
            crater.diameter = crater.radius * 2
            crater.final_to_transient(self.target,self.rng)
        elif crater.transient_diameter is not None:
            crater.transient_radius = crater.transient_diameter / 2
            crater.transient_to_final(self.target,self.rng)
        elif crater.transient_radius is not None:
            crater.transient_diameter = crater.transient_radius * 2
            crater.transient_to_final(self.target,self.rng)

        # if crater.morphology_type is None:
        #     crater.set_morphology_type(self.target,self.rng)
        if crater.location is None:
            crater.location = get_random_location(rng=self.rng)
            
        self.crater = crater
        
        return
        
            
    config_ignore = ['target', 'projectile', 'crater']  # Instance variables to ignore when saving to file
    def to_json(self, filename):
        
        # Get the simulation configuration into the correct structure
        material_config = general_utils.to_config(self.target.material)
        target_config = {**general_utils.to_config(self.target), 'material' : material_config}
        sim_config = {**general_utils.to_config(self),'target' : target_config} 
        
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

    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `utils.general_utils.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """        
        general_utils.set_properties(self,**kwargs)
        return 
