import numpy as np
from numpy.random import default_rng
import os
import json
from .target import Target
from .material import Material
from .projectile import Projectile
from .crater import Crater 
from .mesh import Mesh, load_target_mesh
from ..utils import general_utils  as gu
from . import montecarlo as mc
from ..models import craterscaling as cs 


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
            
    def add_crater(self, **kwargs):
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        crater = Crater(**kwargs)
        
        if crater.diameter is not None:
            crater.radius = crater.diameter / 2
        elif crater.radius is not None:
            crater.diameter = crater.radius * 2
            
        if crater.transient_diameter is not None:
            crater.transient_radius = crater.transient_diameter / 2
        elif crater.transient_radius is not None:
            crater.transient_diameter = crater.transient_radius * 2

        if crater.diameter is not None and crater.transient_diameter is None:
            crater.transient_diameter, crater.morphology_type = cs.final_to_transient(crater.diameter,self.target,self.rng)
            crater.transient_radius = crater.transient_diameter / 2            
        elif crater.transient_diameter is not None and crater.diameter is None:
            crater.diameter, crater.morphology_type = cs.transient_to_final(crater.transient_diameter,self.target,self.rng)
            crater.radius = crater.diameter / 2            
            
        if crater.morphology_type is None:
            crater.set_morphology_type(self.target,self.rng)
        if crater.location is None:
            crater.location = mc.get_random_location(rng=self.rng)
            
        self.crater = crater
        
        return
    

    def add_projectile(self, **kwargs):
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        projectile = Projectile(**kwargs)
        
        if projectile.location is None:
            projectile.location = mc.get_random_location(rng=self.rng)
            
        if projectile.angle is None:
            projectile.angle = mc.get_random_impact_angle(rng=self.rng)
        else: 
            projectile.angle = np.deg2rad(projectile.angle)
        
        if projectile.velocity is None:
            vencounter_mean = np.sqrt(self.target.mean_impact_velocity**2 - self.target.escape_velocity**2)
            vencounter = mc.get_random_velocity(vencounter_mean)
            projectile.velocity = np.sqrt(vencounter**2 + self.target.escape_velocity**2)
                       
        if projectile.velocity is None and projectile.vertical_velocity is not None:
            projectile.velocity = projectile.vertical_velocity / np.sin(projectile.angle)
        elif projectile.vertical_velocity is None and projectile.velocity is not None: 
            projectile.vertical_velocity = projectile.velocity * np.sin(projectile.angle)
        elif projectile.vertical_velocity is not None and projectile.velocity is not None:
            projectile.angle = np.arcsin(projectile.vertical_velocity / projectile.velocity)
            
        if projectile.density is None:
            if projectile.mass is not None and projectile.radius is not None:
                volume = 4.0/3.0 * np.pi * projectile.radius**3
                projectile.density = projectile.mass / volume    
            else: # Default to target density if we are given no way to figure it out
                projectile.density = self.target.material.density 
                
        if projectile.mass is None:
            projectile.mass = 4.0/3.0 * np.pi * projectile.density * projectile.radius**3
        elif projectile.radius is None:
            projectile.radius = (3.0 * projectile.mass / (4.0 * np.pi * projectile.density))**(1.0/3.0)
            projectile.diameter = projectile.radius * 2         
            
        # Now add the crater that this projectile would make
        transient_diameter, _ = cs.projectile_to_transient(projectile, self.target, self.rng) 
        self.add_crater(transient_diameter=transient_diameter,location=projectile.location)
        self.projectile = projectile
        
        return
        
            
    config_ignore = ['target', 'projectile', 'crater']  # Instance variables to ignore when saving to file
    def to_json(self, filename):
        
        # Get the simulation configuration into the correct structure
        material_config = gu.to_config(self.target.material)
        target_config = {**gu.to_config(self.target), 'material' : material_config}
        sim_config = {**gu.to_config(self),'target' : target_config} 
        
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
        property setting is handled by the `utils.gu.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """        
        gu.set_properties(self,**kwargs)
        return 
