import numpy as np
from numpy.random import default_rng
import os
import json
from .target import Target
from .material import Material
from . import projectile
from . import crater 
from .mesh import Mesh, load_target_mesh
from ..utils.general_utils import to_config, set_properties


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

        return
    

    @property
    def projectile(self):
        return self._projectile

    @projectile.setter
    def projectile(self, value):
        self._projectile = value
        if hasattr(value, 'diameter') and value.diameter is not None:
            self.projectile.radius = self.projectile.diameter / 2
        elif hasattr(value, 'radius') and value.radius is not None:
            self.projectile.diameter = self.projectile.radius * 2
                     
   
    @property
    def crater(self):
        return self._crater

    @crater.setter
    def crater(self, value):
        self._crater = value
        if hasattr(value, 'diameter') and value.diameter is not None:
            self.crater.radius = self.crater.diameter / 2
            self.final_to_transient()
            self.crater.transient_radius = self.crater.transient_diameter / 2
        elif hasattr(value, 'radius') and value.radius is not None:
            self.crater.diameter = self.crater.radius * 2
            self.final_to_transient()
            self.crater.transient_radius = self.crater.transient_diameter / 2    
        elif hasattr(value, 'transient_diameter') and value.transient_diameter is not None:
            self.crater.transient_radius = self.crater.transient_diameter / 2
            self.transient_to_final()
            self.crater.radius = self.crater.diameter / 2
        elif hasattr(value, 'transient_radius') and value.transient_radius is not None:
            self.crater.transient_diameter = self.crater.transient_radius * 2
            self.transient_to_final()
            self.crater.radius = self.crater.diameter / 2    
            

    config_ignore = ['target', 'projectile', 'crater']  # Instance variables to ignore when saving to file
    def to_json(self, filename):
        
        # Get the simulation configuration into the correct structure
        material_config = to_config(self.target.material)
        target_config = {**to_config(self.target), 'material' : material_config}
        sim_config = {**to_config(self),'target' : target_config} 
        
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
        property setting is handled by the `util.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        Returns
        -------
        None
            The function does not return a value.
        """        
        set_properties(self,**kwargs)
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


 