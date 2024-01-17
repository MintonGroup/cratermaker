import numpy as np
from numpy.random import default_rng, Generator
import xarray as xr
import json
import os
import tempfile
from typing import Any, Tuple, Type
from .target import Target, Material
from .impact import Crater, Projectile
from .surface import Surface, initialize_surface, save_surface, elevation_to_cartesian
from .scale import Scale
from .morphology import Morphology
from .production import Production, NeukumProduction
from ..utils.general_utils import to_config, set_properties
from ..utils.custom_types import FloatLike, PairOfFloats
from mpas_tools.viz.paraview_extractor import extract_vtk
from ..cython.perlin import apply_noise
import warnings

_POISSON_BATCH_SIZE = 1000

class Simulation:
    """
    This class orchestrates the processes involved in running a crater simulation.

    """
    def __init__(self, 
                 target: str | Target = "Moon",
                 material: str | Material | None = None,
                 pix: FloatLike | None = None,
                 reset_surface: bool = True,
                 simdir: os.PathLike | None = None, 
                 scale_cls: Type[Scale] | None = None,
                 morphology_cls: Type[Morphology] | None = None,
                 production_cls: Type[Production] | None = None,
                 *args: Any,
                 **kwargs: Any):
        """
        Initialize the Simulation object.

        Parameters
        ----------
        target: str or Target, optional, default "Moon"
            Name of the target body or Target object for the simulation, default is "Moon".
        material : str or Material, optional
            Name of the material or Material object for the target body, if None is passed, the default material for the target body 
            is used.
        pix : float, optional
            Pixel resolution for the mesh, default is None.
        reset_surface : bool, optional
            Flag to reset the surface elevation, default is True.
        simdir: PathLike, optional
            Path to the simulation directory, default is current working directory.
        scale_cls : Type[Scale], optional
            The Scale class that defines the crater scaling law. If none provided, then the default will be 
        morphology_cls : Type[Morphology], optional
            The Morphology class that defines the model used to describe the morphology of the crater. If none provided, then the 
            default will be based on the default morphology model.
        production_cls: Type[Production], optional
            The Production class that defines the production function used to populate the surface with craters. If none provided, 
            then the default will be based on the target body, with the NeukumProduction crater-based scaling law used if the target 
            body is the Moon or Mars, the NeukumProduction projectile-based scaling law if the target body is Mercury, Venus, or 
            Earth, and a simple power law model otherwise.
        **kwargs : Any
            Additional keyword arguments that can be passed to any of the method of the class, such as arguments to set the scale, 
            morphology, or production function constructors.
        """
     
        self.simdir = simdir
        
        # # Set some default values for the simulation parameters
        # self.tstart = kwargs.get('tstart', 0.0)  # Simulation start time (in y)
        # self.tstop = kwargs.get('tstop', 4.31e9)    # Simulation stop time (in y)
        
        # Set the random number generator seed
        self.seed = kwargs.get('seed', None) 
        self.rng = kwargs.get('rng', default_rng(seed=self.seed))
        if not isinstance(self.rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        
        self._crater = None
        self._projectile = None        
        
        # First we need to establish the production function. This will allow us to compute the mean impact velocity, which is needed
        # in order to instantiate the target body.
        #  
        # Check to see if the impact velocity model is set with an argument. If not, set it based on a combination of target body 
        # and production function 
        impact_velocity_model = kwargs.get('impact_velocity_model', None)
        
        if not target:
            target_name = "Moon"
        elif isinstance(target, str):
            target_name = target
        elif isinstance(target, Target):
            target_name = target.name
        
        # Set the production function model for this simulation 
        impact_velocity_model = None
        if production_cls is None:
            if target_name in ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars']:
                production_cls = NeukumProduction
                if target_name in ['Moon', 'Mars']:
                    self.production = production_cls(model=target_name, rng=self.rng, **kwargs) 
                else:
                    self.production = production_cls(model='projectile', rng=self.rng, **kwargs)
                impact_velocity_model = target_name + "_MBA"
            else:
                self.production = Production(rng=self.rng, **kwargs)
        elif issubclass(production_cls, Production):
            self.production = production_cls(rng=self.rng, **kwargs)
        else:
            raise TypeError("production must be a subclass of Production")
       
        impact_velocity_model = kwargs.get('impact_velocity_model', impact_velocity_model) 
        if impact_velocity_model is None:
            if target_name in ['Ceres', 'Vesta']:
                impact_velocity_model = "MBA_MBA" 
            
        # Allow an argument to override the value from the production function 
        mean_impact_velocity = kwargs.get("mean_impact_velocity",
                                          self.production.set_mean_impact_velocity(impact_velocity_model=impact_velocity_model)
                                          )
        
        if material:
            if isinstance(material, str):
                try:
                    material = Material(material, **kwargs)
                except:
                    raise ValueError(f"Invalid material name {material}")
            elif not isinstance(material, Material):
                raise TypeError("materiat must be an instance of Material or a valid name of a material")        
            
        if not target:
            if material:
                target = Target("Moon", material=material, mean_impact_velocity=mean_impact_velocity, **kwargs)
            else:
                target = Target("Moon", mean_impact_velocity=mean_impact_velocity, **kwargs)
        elif isinstance(target, str):
            try:
                if material:
                    target = Target(target, material=material, mean_impact_velocity=mean_impact_velocity, **kwargs)
                else:
                    target = Target(target, mean_impact_velocity=mean_impact_velocity, **kwargs)
            except:
                raise ValueError(f"Invalid target name {target}")
        elif isinstance(target, Target):
            if material:
                target.material = material
        elif not isinstance(target, Target):
            raise TypeError("target must be an instance of Target or a valid name of a target body")

        self.target = target

        if pix is not None:
            self.pix = np.float64(pix)
        else:    
            self.pix = np.sqrt(4 * np.pi * self.target.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid

        self.initialize_surface(pix=self.pix, target=self.target, reset_surface=reset_surface, simdir=simdir, *args, **kwargs)
        
        # Set the scaling law model for this simulation 
        if scale_cls is None:
            self.scale_cls = Scale
        elif issubclass(scale_cls, Scale):
            self.scale_cls = scale_cls
        else:
            raise TypeError("scale must be a subclass of Scale") 
      
        # Set the morphology model for this simulation 
        if morphology_cls is None:
            self.morphology_cls = Morphology 
        elif issubclass(morphology_cls, Morphology):
            self.morphology_cls = morphology_cls
        else:
            raise TypeError("morphology must be a subclass of Morphology")
      
        self._craterlist = []
       
        # Find the minimum and maximum possible crater diameter (will be roughly the area of the minimum face size)
        self.smallest_crater = np.sqrt(self.surf['face_areas'].min().item() / np.pi) * 2
        self.largest_crater = self.target.radius * 2
        self.surface_area = 4 * np.pi * self.target.radius**2

        return

    
    def set_properties(self, 
                       **kwargs: Any
                      ) -> None:
        """
        Set properties of the current object based on the provided keyword arguments.

        This function is a utility to update the properties of the current object. The actual implementation of the 
        property setting is handled by the `utils.set_properties` method.

        Parameters
        ----------
        **kwargs : dict
            A dictionary of keyword arguments that represent the properties to be set on the current object.

        """        
        set_properties(self,**kwargs)
        return 

    
    def to_json(self, 
                filename: os.PathLike,
                ) -> str:
        """
        Export the current simulation configuration to a JSON file.

        Parameters
        ----------
        filename : str
            The file path where the JSON configuration will be saved.
        """        
        #TODO: Re-do this once the dust settles a bit
        # Get the simulation configuration into the correct structure
        material_config = to_config(self.target.material)
        target_config = {**to_config(self.target), 'material' : material_config}
        sim_config = {**to_config(self),'target' : target_config} 
        
        # Write the combined configuration to a JSON file
        with open(filename, 'w') as f:
            json.dump(sim_config, f, indent=4)
            
        return
    

    def initialize_surface(self, 
                           *args: Any, 
                           **kwargs: Any
                          ) -> None:
        """
        Initialize the surface mesh.

        Parameters
        ----------
        *args : dict
            Variable length argument list to pass to initialize_surface.
        **kwargs : dict
            Keyword arguments for initializing the surface mesh.
        """        
        self.surf = initialize_surface(*args, **kwargs)
        return
   
    
    def generate_crater(self, 
                        **kwargs: Any
                       ) -> Tuple[Crater, Projectile]:
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
        crater = Crater(target=self.target, morphology_cls=self.morphology_cls, scale_cls=self.scale_cls, rng=self.rng, **kwargs)
        if self.target.mean_impact_velocity is None:
            projectile = None
        else:
            projectile = crater.scale.crater_to_projectile(crater)
        
        return crater, projectile
    
    
    def generate_projectile(self, 
                            **kwargs: Any
                           ) -> Tuple[Projectile, Crater]:
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
        if self.target.mean_impact_velocity is None:
            raise RuntimeError("The mean impact velocity is not set for this simulation")
        projectile = Projectile(target=self.target, rng=self.rng, scale_cls=self.scale_cls, **kwargs)
        crater = projectile.scale.projectile_to_crater(projectile, morphology_cls=self.morphology_cls)
        
        return projectile, crater
   
   
    def emplace_crater(self, 
                       from_projectile: bool=False, 
                       **kwargs: Any
                      ) -> None:
        """
        Emplace a crater in the simulation, optionally based on a projectile.

        Parameters
        ----------
        from_projectile : bool, optional
            Flag to create a crater based on a projectile, default is False.
        **kwargs : Any
            Keyword arguments for initializing the Crater or Projectile object.
        """        
        if from_projectile:
            self.projectile, self.crater = self.generate_projectile(**kwargs)
        else:
            self.crater, self.projectile = self.generate_crater(**kwargs)
        self.surf['node_crater_distance'], self.surf['face_crater_distance'] = self.surf.get_distance(self.crater.location)
        self.surf['node_crater_bearing'], self.surf['face_crater_bearing'] = self.surf.get_initial_bearing(self.crater.location)
        self.crater.node_index, self.crater.face_index = self.surf.find_nearest_index(self.crater.location)
       
        self.crater.morphology.form_crater(self.surf)
        
        return  


    def populate(self, 
                 age: FloatLike | None = None,
                 reference_age: FloatLike | None = None,
                 cumulative_number_at_diameter: PairOfFloats | None = None,
                 reference_cumulative_number_at_diameter: PairOfFloats | None = None,
                 **kwargs: Any,
                ) -> None:
        """
        Populate the surface with craters over a specified interval using the current production function.
        
        
        Parameters
        ----------
        age : FloatLike or ArrayLike, optional
            Age in the past relative to the reference age to compute cumulative SFD in units of My. 
        reference_age, FloatLike or ArrayLike, optional
            The reference used when computing age in My. The default is 0 (present day). If a non-zero value is passed, the `age` is 
            interpreted as a delta on the reference age. So for instance, if `age=500` and `reference_age=3500`, then this means 
            "4.0 Gy to 3.5 Gy ago". 
        cumulative_number_at_diameter : PairOfFloats, optional
            A pair of cumulative number and diameter values, in the form of a (N, D). If provided, the function will convert this value
            to a corresponding age and use the production function for a given age.
        reference_cumulative_number_at_diameter : PairOfFloats, optional
            A pair of cumulative number and diameter values, in the form of a (N, D). If provided, the function will convert this
            value to a corresponding reference age and use the production function for a given age.        
        """
        
        if not hasattr(self, 'production'):
            raise RuntimeError("No production function defined for this simulation")
        elif not hasattr(self.production, 'generator_type'):
            raise RuntimeError("The production function is not properly defined. Missing 'generator_type' attribute")
        elif self.production.generator_type not in ['crater', 'projectile']:
            raise RuntimeError(f"Invalid production function type {self.production.generator_type}")
        
        impacts_this_interval = self.production.sample(age=age, 
                                                       reference_age=reference_age, 
                                                       cumulative_number_at_diameter=cumulative_number_at_diameter, 
                                                       reference_cumulative_number_at_diameter=reference_cumulative_number_at_diameter, 
                                                       diameter_range=(self.smallest_crater, self.largest_crater),
                                                       area=self.surface_area, 
                                                       **kwargs)
        for diameter in impacts_this_interval:
            print(f"Emplacing crater of diameter {diameter*1e-3:.2f} km")
            self.emplace_crater(diameter=diameter)
        return 
         
    def save(self, 
             *args: Any, 
             **kwargs: Any
            ) -> None:
        """
        Save the current simulation state to a file.
        """
        save_surface(self.surf, *args, **kwargs)
        
        return
    
    
    def export_vtk(self, 
                   out_dir: os.PathLike | None = None,
                   *args, **kwargs
                   ) -> None:
        """
        Export the surface mesh to a VTK file.

        Parameters
        ----------
        out_dir : str, Default "vtk_files"
            Directory to store the VTK files.
        """
        
        if out_dir is None:
            out_dir = os.path.join(self.simdir, "vtk_files")
            
        # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
        # The API change does not affect the functionality of the code, so we can safely ignore the warning
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", FutureWarning)
            warnings.simplefilter("ignore", DeprecationWarning) # Ignores a warning issued in bar.py
            ignore_time = "time" not in self.surf.dims
        
            # Save the surface data to a combined netCDF file
            with tempfile.TemporaryDirectory() as temp_dir:
                self.save(combine_data_files=True, out_dir=temp_dir)
                
                # Use elevation data to modify the mesh for visualization purposes
                grid = xr.open_dataset(self.grid_file)
            
                vert_vars = ['xVertex', 'yVertex', 'zVertex']
                
                ds_new = elevation_to_cartesian(grid[vert_vars], self.surf['node_elevation'])
                for var in vert_vars:
                    grid[var] = ds_new[var]
                    
                face_vars = ['xCell', 'yCell', 'zCell']
                ds_new = elevation_to_cartesian(grid[face_vars], self.surf['face_elevation'])
                for var in face_vars:
                    grid[var] = ds_new[var]
                
                grid.to_netcdf(os.path.join(temp_dir, "surface_mesh.nc"))
                
                # Combine the grid and data into one file
                try:
                    extract_vtk(
                        filename_pattern=os.path.join(temp_dir, "*.nc"),
                        mesh_filename=os.path.join(temp_dir, "surface_mesh.nc"),
                        variable_list=['allOnVertices', 'allOnCells'] , 
                        dimension_list=['maxEdges=','vertexDegree='], 
                        combine=True,
                        include_mesh_vars=True,
                        out_dir=out_dir, 
                        ignore_time=ignore_time, 
                        *args, **kwargs)
                except:
                    raise RuntimeError("Error in xtract_vtk. Cannot export VTK files")
        
        return
    

    def apply_noise(self, 
                    model="turbulence",
                    noise_width=1000e3,
                    noise_height=20e3,
                    to_nodes=True,
                    to_faces=True,
                    **kwargs,
                    ) -> None:
        """
        Applies a specified noise model to the simulation's surface elevation.

        This method adjusts the surface elevation of the simulation based on the chosen noise model. It supports various noise models, each with its own set of default and customizable parameters. The method also ensures that the applied noise is volume-conserving.

        Parameters
        ----------
        model : str, optional
            The noise model to apply. Supported models include 'turbulence', 'billowed', 'plaw', 'ridged', 'swiss', and 'jordan'. The default is 'turbulence'.
        noise_width : float, optional
            The width scale of the noise in meters. The default is 1000e3 (1000 km).
        noise_height : float, optional
            The height scale of the noise in meters. The default is 20e3 (20 km).
        to_nodes : bool, optional
            Flag to apply the noise to the nodes of the mesh. The default is True.
        to_faces : bool, optional
            Flag to apply the noise to the faces of the mesh. The default is False.
        **kwargs :
            Additional keyword arguments specific to the noise model. Common parameters include 'num_octaves' and 'anchor'. Model-specific parameters like 'freq', 'pers', 'slope', 'lacunarity', 'gain', etc., can also be set.

        Returns
        -------
        None
            This method modifies the simulation's surface elevation in-place and does not return a value.

        Notes
        -----
        - The noise is scaled to be volume-conserving, ensuring the mean of the noise is zero.
        - The method internally calculates normalized coordinates based on the target radius and scales the noise appropriately.
        - Default values for noise parameters are set based on the chosen model.
        - For details on thes noise models, see https://www.decarpentier.nl/scape-procedural-basics

        Examples
        --------
        Apply default turbulence noise:

        >>> sim = cratermaker.Simulation()
        >>> sim.apply_noise()

        Apply ridged noise with custom Parameters

        >>> sim.apply_noise(model="ridged", noise_width=500e3, num_octaves=10, freq=1.5)
        """      
        scale = self.target.radius / noise_width
        num_octaves = kwargs.pop("num_octaves", 12)
        anchor = kwargs.pop("anchor", self.rng.uniform(0.0,scale, size=(num_octaves, 3))) 
        
        # Set reasonable default values for the different models
        if model == "turbulence" or model == "billowed" or model == "plaw" or model == "ridged":
            kwargs.setdefault("noise_height", noise_height)
            kwargs.setdefault("freq", 2.00)
            kwargs.setdefault("pers", 0.5)
        if model == "plaw":
            kwargs.setdefault("slope", 2.0)
        if model == "swiss" or model == "jordan":
            kwargs.setdefault("lacunarity", 2.00)
            kwargs.setdefault("gain", 0.5)
            kwargs.setdefault("warp", 0.35)
        if model == "jordan":
            kwargs.setdefault("gain0", 0.8)
            kwargs.setdefault("warp0", 0.4)
            kwargs.setdefault("damp0", 1.0)
            kwargs.setdefault("damp", 0.8)
            kwargs.setdefault("damp_scale", 1.0) 
            
        if "noise_height" in kwargs:
            kwargs["noise_height"] = kwargs["noise_height"] / self.target.radius
            
        def _noisemaker(vars):
            ds_norm = self.surf.uxgrid._ds[vars] * scale / self.target.radius
            x = ds_norm[vars[0]].values
            y = ds_norm[vars[1]].values
            z = ds_norm[vars[2]].values
            noise = apply_noise(model, x, y, z, num_octaves, anchor, **kwargs)
        
            # Make sure the noise is volume-conserving (i.e., the mean is zero)
            # TODO: Take into account the nodes are not uniformly distributed on the sphere
            noise = noise - np.mean(noise)
            return noise                
        
        if to_nodes:
            vars = ['node_x', 'node_y', 'node_z']
            node_noise = _noisemaker(vars)
        else:
            node_noise = None
        if to_faces:
            vars = ['face_x', 'face_y', 'face_z']
            face_noise = _noisemaker(vars)
        else:
            face_noise = None
            
        if model =="swiss" or model == "jordan":
            if node_noise is not None:
                self.surf['node_elevation'] += node_noise * noise_height
            if face_noise is not None:
                self.surf['face_elevation'] += face_noise * noise_height
        else:
            if node_noise is not None:
                self.surf['node_elevation'] += node_noise * self.target.radius 
            if face_noise is not None:
                self.surf['face_elevation'] += face_noise * self.target.radius
        
        return
    
     
    def set_elevation(self, 
                      *args: Any, 
                      **kwargs: Any
                     ) -> None:
        """
        Set the elevation on the surface. Delegates to the Surface object.

        Parameters
        ----------
        *args: Variable length argument list to pass to self.surf.set_elevation.
        **kwargs: Arbitrary keyword arguments to pass to self.surf.set_elevation.
        """
        return self.surf.set_elevation(*args, **kwargs)   
    

    @property
    def simdir(self):
        """
        Directory where the simulation data is stored. Set during initialization.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value):
        if value is None:
            self._simdir = os.getcwd() 
        elif not os.path.isabs(value):
            self._simdir = os.path.abspath(value)
        else:
            self._simdir = value
            
        if not os.path.exists(self._simdir):
            os.makedirs(self._simdir)

    @property
    def seed(self):
        """
        Seed for the random number generator. Set during initialization.
        """
        return self._seed

    @seed.setter
    def seed(self, value):
        self._seed = value

    @property
    def rng(self):
        """
        Random number generator instance. Set during initialization.
        """
        return self._rng

    @rng.setter
    def rng(self, value):
        if not isinstance(value, Generator) and value is not None:
            raise TypeError("rng must be a numpy.random.Generator instance or None")
        self._rng = value

    @property
    def target(self):
        """
        The target body for the impact simulation. Set during initialization.
        """
        return self._target

    @target.setter
    def target(self, value):
        if not isinstance(value, Target):
            raise TypeError("target must be an instance of Target")
        self._target = value

    @property
    def pix(self):
        """
        Pixel resolution for the mesh. Set during initialization.
        """
        return self._pix

    @pix.setter
    def pix(self, value):
        if value <= 0:
            raise ValueError("pix must be greater than zero")
        self._pix = value

    @property
    def surf(self):
        """
        Surface mesh data for the simulation. Set during initialization.
        """
        return self._surf
    
    @surf.setter
    def surf(self, value):
        if not isinstance(value, Surface):
            raise TypeError("surf must be an instance of Surface")
        self._surf = value

    @property
    def production(self):
        """
        The Production class instance used for crater production. Set during initialization.
        """
        return self._production

    @production.setter
    def production(self, value):
        if not issubclass(value.__class__, Production):
            raise TypeError("production must be a subclass of Production")
        self._production = value

    @property
    def scale_cls(self):
        """
        The Scale class that defines the crater scaling law. Set during initialization.
        """
        return self._scale_cls

    @scale_cls.setter
    def scale_cls(self, value):
        if not issubclass(value, Scale):
            raise TypeError("scale_cls must be a subclass of Scale")
        self._scale_cls = value

    @property
    def morphology_cls(self):
        """
        The Morphology class that defines the crater morphology model. Set during initialization.
        """
        return self._morphology_cls

    @morphology_cls.setter
    def morphology_cls(self, value):
        if not issubclass(value, Morphology):
            raise TypeError("morphology_cls must be a subclass of Morphology")
        self._morphology_cls = value

    @property
    def crater(self):
        """
        The current Crater object in the simulation. Set during runtime.
        """
        return self._crater

    @crater.setter
    def crater(self, value):
        if not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value

    @property
    def projectile(self):
        """
        The current Projectile object in the simulation. Set during runtime.
        """
        return self._projectile

    @projectile.setter
    def projectile(self, value):
        if value is not None and not isinstance(value, Projectile):
            raise TypeError("projectile must be an instance of Projectile")
        self._projectile = value

    @property
    def data_dir(self):
        """
        Directory where the data files are stored. Dynamically set based on `surf` attribute.
        """
        return self.surf.data_dir

    @property
    def grid_file(self):
        """
        File path of the grid file. Dynamically set based on `surf` attribute.
        """
        return self.surf.grid_file


    @property
    def n_node(self):
        """
        Number of nodes in the simulation mesh. Dynamically set based on `surf` attribute.
        """
        return self.surf.uxgrid.n_node
   
    @property
    def n_node(self):
        """
        Number of nodes in the simulation mesh. Dynamically set based on `surf` attribute.
        """
        return self.surf.uxgrid.n_node

    @property
    def n_face(self):
        """
        Number of faces in the simulation mesh. Dynamically set based on `surf` attribute.
        """
        return self.surf.uxgrid.n_face
