import numpy as np
from numpy.random import default_rng
import json
from dacite import from_dict
import uxarray as uxr
from typing import Any
import os
from pathlib import Path
from .target import Target, Material
from .crater import Crater, Projectile
from .surface import Surface, initialize_surface
from ..utils import general_utils as gu
from ..utils.general_utils import float_like

class Simulation():
    """
    This class orchestrates the processes involved in running a crater simulation.

    Attributes
    ----------
    pix : float
        Pixel resolution for the mesh.
    target : Target
        The target body for the impact simulation.
    rng : numpy.random.Generator
        Random number generator instance.

    Methods
    -------
    set_properties(**kwargs):
        Set properties of the current object based on the provided keyword arguments.
    to_json(filename):
        Export the current simulation configuration to a JSON file.
    generate_crater(**kwargs):
        Create a new Crater object and its corresponding Projectile.
    generate_projectile(**kwargs):
        Create a new Projectile object and its corresponding Crater.
    emplace_crater(from_projectile=False, **kwargs):
        Emplace a crater in the simulation, optionally based on a projectile.
    """

    def __init__(self, 
                target_name: str="Moon",
                material_name: str | None = None,
                pix: float_like | None = None,
                reset_surface: bool = True,
                *args: Any,
                **kwargs: Any):
        """
        Initialize the Simulation object.

        Parameters
        ----------
        target_name : str, optional
            Name of the target body for the simulation, default is "Moon".
        material_name : str, optional
            Name of the material for the target body, default is None.
        reset_surface : bool, optional
            Flag to reset the surface elevation, default is True.
        pix : float, optional
            Pixel resolution for the mesh, default is None.
        **kwargs : Any
            Additional keyword arguments.
        """
            
        if material_name:
            material = Material(name=material_name)
            self.target = from_dict(data_class=Target,data=dict({"name":target_name,"material":material}, **kwargs)) 
        else: 
            self.target = from_dict(data_class=Target,data=dict({"name":target_name}, **kwargs))
            
        if pix is not None:
            self.pix = np.float64(pix)
        else:    
            self.pix = np.sqrt(4 * np.pi * self.target.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid

        self.surf = initialize_surface(pix=self.pix, target=self.target, reset_surface=reset_surface, *args, **kwargs)
        
        # Set some default values for the simulation parameters
        self.time_function = kwargs.get('time_function', None)
        self.tstart = kwargs.get('tstart', 0.0)  # Simulation start time (in y)
        self.tstop = kwargs.get('tstop', 4.31e9)    # Simulation stop time (in y)
        
        # Set some default values for the production function population
        self.impactor_sfd  = kwargs.get('impactor_sfd', None)
        self.impactor_velocity = kwargs.get('impactor_velocity', None)
        
        # Set the random number generator seed
        self.seed = kwargs.get('seed', None) 
        self.rng = default_rng(seed=self.seed)
        
        self._crater = None
        self._projectile = None

        return

    
    def set_properties(self, **kwargs):
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `utils.gu.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        """        
        gu.set_properties(self,**kwargs)
        return 

    
    def to_json(self, filename):
        """
        Export the current simulation configuration to a JSON file.

        Parameters
        ----------
        filename : str
            The file path where the JSON configuration will be saved.
        """        
        #TODO: Re-do this once the dust settles a bit
        # Get the simulation configuration into the correct structure
        material_config = gu.to_config(self.target.material)
        target_config = {**gu.to_config(self.target), 'material' : material_config}
        sim_config = {**gu.to_config(self),'target' : target_config} 
        
        # Write the combined configuration to a JSON file
        with open(filename, 'w') as f:
            json.dump(sim_config, f, indent=4)
            
        return
    

    def generate_crater(self, **kwargs):
        """
        Create a new Crater object and its corresponding Projectile.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments for initializing the Crater object.

        Returns
        -------
        (Crater, Projectile)
            A tuple containing the newly created Crater and Projectile objects.
        """        
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        crater = Crater(target=self.target, rng=self.rng, **kwargs)
        projectile = crater.scale.crater_to_projectile(crater)
        
        return crater, projectile
    
    
    def generate_projectile(self, **kwargs):
        """
        Create a new Projectile object and its corresponding Crater.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments for initializing the Projectile object.

        Returns
        -------
        (Projectile, Crater)
            A tuple containing the newly created Projectile and Crater objects.
        """
        projectile = Projectile(target=self.target, rng=self.rng, **kwargs)
        crater = projectile.scale.projectile_to_crater(projectile)
        
        return projectile, crater
   
   
    def emplace_crater(self, from_projectile=False, **kwargs):
        """
        Emplace a crater in the simulation, optionally based on a projectile.

        Parameters
        ----------
        from_projectile : bool, optional
            Flag to create a crater based on a projectile, default is False.
        **kwargs : dict
            Keyword arguments for initializing the Crater or Projectile object.
        """        
        if from_projectile:
            self.projectile, self.crater = self.generate_projectile(**kwargs)
        else:
            self.crater, self.projectile = self.generate_crater(**kwargs)
        self.surf['crater_distance'] = self.surf.get_cell_distance(self.crater.location, self.target.radius)
        self.surf['crater_bearing'] = self.surf.get_cell_initial_bearing(self.crater.location)
        
        self.crater.average_surface_normal_vector = self.surf.get_average_surface(self.crater.location, self.crater.radius)
        
        return  

    # The following are placeholders for if/when we need to pass data back and forth to the Fortran library     
    # @property
    # def elevation(self):
    #     return self._body.fobj.elevation
    
    # @property
    # def name(self):
    #     return self._body.fobj.name


    # def get_elevation(self):
    #     return self._body.get_elevation()

    
    # def set_elevation(self, elevation_array):
    #     self._body.set_elevation(elevation_array)
