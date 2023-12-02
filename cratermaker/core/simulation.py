import numpy as np
from numpy.random import default_rng
import json
from dacite import from_dict
import xarray as xr
from typing import Any
import os
from pathlib import Path
from .target import Target, Material
from .crater import Crater, Projectile, CraterScaling
from ..utils import general_utils as gu
from ..utils.general_utils import float_like
from ..utils import mesh_tools
from ..utils import montecarlo as mc

class Simulation():
    """
    This is a class that defines the Simulation object in Cratermaker, which orchestrates the processes involved in running a simulation. 
    
 
    Attributes
    ----------
    mesh_temp_dir : str, Path, or os.PathLike or None
        Directory path for caching data.
    mesh_file : str, Path, or os.PathLike or None
        File path for the dataset file describing the target mesh.
    dem_file : str, Path, or os.PathLike or None
        File path for the dataset file describing the target elevation data. 
    ds : xarray.Dataset
        xarray Dataset representing the surface mesh and associated data.    
    pix : float_like or None
        Pixel resolution for the mesh.
        
        
    """
    def __init__(self, 
                target_name: str="Moon",
                material_name: str | None = None,
                make_new_mesh:  bool | None = None,
                reset_surface: bool = True,
                mesh_temp_dir: str | Path | os.PathLike = ".mesh",
                mesh_file: str | Path | os.PathLike = "surface_mesh.nc",
                dem_file: str | Path | os.PathLike = "surface_dem.nc",
                pix: float_like | None = None,
                mesh_data: xr.Dataset | None=None,
                dem_data: xr.Dataset | None=None,
                **kwargs: Any):
      
        self.mesh_temp_dir = Path(mesh_temp_dir)
        self.mesh_file = Path(mesh_file)
        if not self.mesh_file.is_absolute():
            self.mesh_file = Path.cwd() / self.mesh_file
            
        self.dem_file = Path(dem_file)
        if not self.dem_file.is_absolute():
            self.dem_file = Path.cwd() / self.dem_file
            
        if material_name:
            material = Material(name=material_name)
            self.target = from_dict(data_class=Target,data=dict({"name":target_name,"material":material}, **kwargs)) 
        else: 
            self.target = from_dict(data_class=Target,data=dict({"name":target_name}, **kwargs))
            
        if pix is not None:
            self.pix = np.float64(pix)
        else:    
            self.pix = np.sqrt(4 * np.pi * self.target.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid
                
        if not os.path.exists(self.mesh_temp_dir):
            os.mkdir(self.mesh_temp_dir)

        # Generate a new surface if either it is explicitly requested via parameter or a data file doesn't yet exist 
        if make_new_mesh is None:
            make_new_mesh = not os.path.exists(self.mesh_file)
        
        if make_new_mesh:
            self.mesh = mesh_tools.generate(self.target.radius,self.pix,self.mesh_file,self.mesh_temp_dir)
            reset_surface = True
        else:
            self.mesh = xr.open_dataset(self.mesh_file)
            
        if reset_surface:
            self.dem = mesh_tools.set_elevation(nCells=self.mesh.nCells.size)
            
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

        Returns
        -------
        None
            The function does not return a value.
        """        
        gu.set_properties(self,**kwargs)
        return 

    
    def to_json(self, filename):
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
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        crater = Crater(self.target, self.rng, **kwargs)
        
        self.crater = crater
        
        return
    
    
    def generate_projectile(self, **kwargs):
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        projectile = Projectile(self.target, self.rng, **kwargs)

        # Now add the crater that this projectile would make
        #transient_diameter, _ = .projectile_to_transient(projectile, self.target, self.rng) 
        self.generate_crater(transient_diameter=transient_diameter,location=projectile.location)
        self.projectile = projectile
        
        return
   
   
    def emplace_crater(self, from_projectile=False, **kwargs):
        if from_projectile:
            self.generate_projectile(**kwargs)
        else:
            self.generate_crater(**kwargs)
        self.dem['crater_distance'] = mesh_tools.get_cell_distance(self.mesh, self.crater.location, self.target.radius)
        self.dem['crater_bearing'] = mesh_tools.get_cell_initial_bearing(self.mesh, self.crater.location)
        
        self.crater.average_surface_normal_vector = mesh_tools.get_average_surface(self.mesh, self.dem, self.crater.location, self.crater.radius)
        
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
