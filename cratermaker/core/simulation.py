import numpy as np
from numpy.random import default_rng, Generator
import xarray as xr
import json
import os
import shutil
from glob import glob
import tempfile
from typing import Any, Tuple, Type
from .target import Target, Material
from .impact import Crater, Projectile
from .surface import Surface, save, initialize_surface, elevation_to_cartesian
from .scale import Scale
from .morphology import Morphology
from .production import Production, NeukumProduction
from ..utils.general_utils import to_config, set_properties
from ..utils.custom_types import FloatLike, PairOfFloats
from mpas_tools.viz.paraview_extractor import extract_vtk
from ..cython.perlin import apply_noise
import warnings
from tqdm import tqdm

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
        
        # Set the random number generator seed
        self.seed = kwargs.get('seed', None) 
        self.rng = kwargs.get('rng', default_rng(seed=self.seed))
        if not isinstance(self.rng, Generator):
            raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
        
        self._crater = None
        self._projectile = None
        self._interval_number = 0
        self._elapsed_time = 0.0
        self._current_age = 0.0
        self._elapsed_n1 = 0.0
        
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
        **kwargs : Any
            Keyword arguments for initializing the :class:`Crater` object. Refer to 
            its documentation for details on valid keyword arguments.

        Returns
        -------
        tuple of (Crater, Projectile)
            A tuple containing the newly created Crater and Projectile objects.
            
        Notes
        -----
        The keyword arguments provided are passed to the constructor of 
        :class:`Crater`. Additionally, these arguments are used in 
        :meth:`crater.scale.crater_to_projectile` method. Refer to the 
        documentation of these classes and methods for a detailed description 
        of valid keyword arguments.

        Examples
        --------
        .. code-block:: python
        
            # Create a crater and projectile pair with a specific diameter
            crater, projectile = sim.generate_crater(diameter=1000.0)
    
            # Create a crater with a specific transient diameter and location (but ignore the projectile)
            crater, _ = sim.generate_crater(transient_diameter=5e3, location=(43.43, -86.92))
        """       
        # Create a new Crater object with the passed arguments and set it as the crater of this simulation
        crater = Crater(target=self.target, morphology_cls=self.morphology_cls, scale_cls=self.scale_cls, rng=self.rng, **kwargs)
        if self.target.mean_impact_velocity is None:
            projectile = None
        else:
            projectile = crater.scale.crater_to_projectile(crater,**kwargs)
        
        return crater, projectile
    
    
    def generate_projectile(self, 
                            **kwargs: Any
                           ) -> Tuple[Projectile, Crater]:
        """
        Create a new Projectile object and its corresponding Crater.

        Parameters
        ----------
        **kwargs : Any
            Keyword arguments for initializing the :class:`Projectile` object. Refer to its documentation 
            for details on valid keyword arguments.

        Returns
        -------
        (Projectile, Crater)
            A tuple containing the newly created Projectile and Crater objects.
            
        Notes
        -----
        The keyword arguments provided are passed to the constructor of 
        :class:`Projectile`. Refer to its documentation for a detailed description 
        of valid keyword arguments.    
        
        Examples
        --------
        .. code-block:: python
        
            # Create a projectile and crater pair with a specific diameter
            projectile, crater = sim.generate_projectile(diameter=1000.0)

            # Create a projectile and crater with a specific mass and velocity and impact angle
            projectile, crater = sim.generate_projectile(mass=1e14, velocity=20e3, angle=10.0)       
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

        This method orchestrates the creation and placement of a crater in the
        simulation. It can create a crater directly or based on the characteristics
        of a projectile.

        Parameters
        ----------
        from_projectile : bool, optional
            Flag to create a crater based on a projectile, default is False.
        **kwargs : Any
            Keyword arguments for initializing the :class:`Crater` or 
            :class:`Projectile` object. Refer to the documentation of these 
            classes for details on valid keyword arguments.

        Notes
        -----
        The keyword arguments provided are passed down to :meth:`generate_crater` 
        or :meth:`generate_projectile`, and subsequently to the constructors of 
        :class:`Crater`, :class:`Projectile`, or :class:`Impact`. Refer to the 
        documentation of these classes for a detailed description of valid 
        keyword arguments.

        Examples
        --------
        .. code-block:: python        
        
            # Create a crater with specific diameter
            sim.emplace_crater(diameter=1000.0)

            # Create a crater based on a projectile with given mass and velocity
            sim.emplace_crater(from_projectile=True, mass=1e14, velocity=20e3)
            
            # Create a crater with a specific transient diameter and location
            sim.emplace_crater(transient_diameter=5e3, location=(43.43, -86.92))
        """ 
        if from_projectile:
            self.projectile, self.crater = self.generate_projectile(**kwargs)
        else:
            self.crater, self.projectile = self.generate_crater(**kwargs)
       
        self.crater.morphology.form_crater(self.surf,**kwargs)
        
        return  


    def populate(self, 
                 age: FloatLike | None = None,
                 age_end: FloatLike | None = None,
                 diameter_number: PairOfFloats | None = None,
                 diameter_number_end: PairOfFloats | None = None,
                 **kwargs: Any,
                ) -> None:
        """
        Populate the surface with craters over a specified interval using the current production function.
        
        Parameters
        ----------
        age : FloatLike, optional
            Age in the past in units of My relative to the present, which is used compute the cumulative SFD.
        age_end, FloatLike, optional
            The ending age in units of My relative to the present, which is used to compute the cumulative SFD. The default is 0 (present day).
        diameter_number : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function will convert this value
            to a corresponding age and use the production function for a given age.
        diameter_number_end : PairOfFloats, optional
            A pair of diameter and cumulative number values, in the form of a (D, N). If provided, the function will convert this
            value to a corresponding reference age and use the production function for a given age.        
        """

        if not hasattr(self, 'production'):
            raise RuntimeError("No production function defined for this simulation")
        elif not hasattr(self.production, 'generator_type'):
            raise RuntimeError("The production function is not properly defined. Missing 'generator_type' attribute")
        elif self.production.generator_type not in ['crater', 'projectile']:
            raise RuntimeError(f"Invalid production function type {self.production.generator_type}")
        
        impacts_this_interval, impact_ages = self.production.sample(age=age, 
                                                                    age_end=age_end, 
                                                                    diameter_number=diameter_number, 
                                                                    diameter_number_end=diameter_number_end, 
                                                                    diameter_range=(self.smallest_crater, self.largest_crater),
                                                                    area=self.surf.area.item(), 
                                                                    **kwargs)
        
        if impacts_this_interval is not None:
            if impacts_this_interval.size == 1:
                self.emplace_crater(diameter=impacts_this_interval, age=impact_ages)
            else:
                for i, diameter in tqdm(enumerate(impacts_this_interval), total=len(impacts_this_interval)):
                    self.emplace_crater(diameter=diameter, age=impact_ages[i])
        return 
    
    
    def run(self,
            age: FloatLike | None = None,
            age_end: FloatLike | None = None,
            age_interval: FloatLike | None = None,
            diameter_number: PairOfFloats | None = None,
            diameter_number_end: PairOfFloats | None = None,
            diameter_number_interval: FloatLike | None = None,
            ninterval: int | None = None,
            **kwargs: Any,
           ) -> None:
        """
        Run the simulation over a specified interval using the current production function.

        Parameters
        ----------
        age : FloatLike, optional
            Start age in My relative to the present for the simulation, used to compute the starting point of the production function.
            Default is None, which requires `diameter_number` to be set.
        age_end : FloatLike, optional
            End age in My relative to the present for the simulation, used to compute the ending point of the production function.
            Default is 0 (present day) if not provided.
        age_interval : FloatLike, optional
            Interval in My for outputting intermediate results. If not provided, calculated as `age - age_end` / `ninterval` if `ninterval` is provided, otherwise set to the total simulation duration.
        diameter_number : PairOfFloats, optional
            Starting cumulative number and diameter pair (D, N) to define the simulation start in terms of crater frequency. 
            Default is None, which requires `age` to be set.
        diameter_number_end : PairOfFloats, optional
            Ending cumulative number and diameter pair (D, N) to define the simulation end in terms of crater frequency.
            Default is the present-day values if not provided.
        diameter_number_interval : PairOfFloats, optional
            Interval for outputting results in terms of cumulative number and diameter (D, N). Calculated based on `ninterval` if provided.
        ninterval : int, optional
            Number of intervals for outputting results. This has a special use case where one can specify age-based inputs but output
            in equal cumulative number intervals and vice versa.  

        Notes
        -----
        This function allows defining the simulation parameters either in terms of age or crater frequency (cumulative number). The 
        arguments `age`, `age_end`, `age_interval` are mutually exclusive with `diameter_number`, `diameter_number_end`, `diameter_number_interval`.
        
        The input diameter values used in diameter_number, diameter_number_end, and diameter_number_interval need not be the same.
        However, the production function will be used to convert all numbers to a common diameter value on output.

        Examples
        --------
        .. code-block:: python
        
            # Create a simulation object with default parameters (Moon, NeukumProduction, etc.)
            sim = cratermaker.Simulation()
            
            # Run the simulation for 3.8 billion years, saving the results every 100 million years
            sim.run(age=3.8e3, age_interval=100.0)
        
            # Run the simulation to create 80 craters larger than 300 km in diameter and output the results after approximately 100 
            # craters larger than 50km in diameter have formed
            sim.run(diameter_number=(300e3, 80), diameter_number_interval=(50e3, 100))
        
            # Run the simulation for 3.8 billion years, saving 100 intervals with equal cumulative number intervals
            sim.run(age=3.8e3, ninterval=100)
        
            # Run the simulation to create 80 craters larger than 300 km and output 100 intervals with equal age intervals
            sim.run(diameter_number=(300e3, 80), ninterval=100)
        
            # Run the simulation from 3.8 billion years to 3.0 billion years, saving the results every 100 million years
            sim.run(age=3.8e3, age_end=3.0e3, age_interval=100.0)

        """
        arguments = locals().copy()
        arguments.pop('self')
        arguments = self._validate_run_args(**arguments) 
        age = arguments.pop('age', None)
        age_end = arguments.pop('age_end', None)
        age_interval = arguments.pop('age_interval', None)
        diameter_number = arguments.pop('diameter_number', None)
        diameter_number_end = arguments.pop('diameter_number_end', None)
        diameter_number_interval = arguments.pop('diameter_number_interval', None)
        ninterval = arguments.pop('ninterval', None)
        is_age_interval = arguments.pop('is_age_interval', None)
        
        
        self.current_age = age
        self.elapsed_time = 0.0
        self.elapsed_n1 = 0.0
        for i in tqdm(range(ninterval+1), total=ninterval+1):
            self.interval_number = i
            if i > 0: # This allows us to save the initial state of the simulation
                if is_age_interval:
                    current_age = age - (i-1) * age_interval
                    current_age_end = age - i * age_interval
                    self.populate(age=current_age, age_end=current_age_end)
                else:
                    current_diameter_number = (diameter_number[0], diameter_number[1] - (i-1) * diameter_number_interval[1])
                    current_diameter_number_end = (diameter_number[0], diameter_number[1] - i * diameter_number_interval[1])
                    self.populate(diameter_number=current_diameter_number, diameter_number_end=current_diameter_number_end)
                    current_age = self.production.function_inverse(*current_diameter_number)
                    current_age_end = self.production.function_inverse(*current_diameter_number_end)
                    age_interval = current_age - current_age_end 
                self.elapsed_time += age_interval
                self.elapsed_n1 += (self.production.function(diameter=1000.0, age=current_age)  - self.production.function(diameter=1000.0, age = current_age_end))
                self.current_age = current_age_end
                
            self.save()
            
        self.export_vtk()
        return

    def _validate_run_args(self,**kwargs) -> dict:
        """
        Validate all the input arguments to the sample method. This function will raise a ValueError if any of the arguments are invalid.
        It will also convert age arguments to diameter_number and vice versa.
        
        Parameters
        ----------
        age : FloatLike, optional
            Start age in My relative to the present for the simulation, used to compute the starting point of the production function.
            Default is None, which requires `diameter_number` to be set.
        age_end : FloatLike, optional
            End age in My relative to the present for the simulation, used to compute the ending point of the production function.
            Default is 0 (present day) if not provided.
        age_interval : FloatLike, optional
            Interval in My for outputting intermediate results. If not provided, calculated as `age - age_end` / `ninterval` if `ninterval` is provided, otherwise set to the total simulation duration.
        diameter_number : PairOfFloats, optional
            Starting cumulative number and diameter pair (D, N) to define the simulation start in terms of crater frequency. 
            Default is None, which requires `age` to be set.
        diameter_number_end : PairOfFloats, optional
            Ending cumulative number and diameter pair (D, N) to define the simulation end in terms of crater frequency.
            Default is the present-day values if not provided.
        diameter_number_interval : PairOfFloats, optional
            Interval for outputting results in terms of cumulative number and diameter (D, N). Calculated based on `ninterval` if provided.
        ninterval : int, optional
            Number of intervals for outputting results. This has a special use case where one can specify age-based inputs but output
            in equal cumulative number intervals and vice versa.  
                
        Returns
        -------
        A dict containing all arguments listed in Parameters above, as well as `is_age_interval`, which is a boolean flag indicating
        whether or not the simulation is being run in equal age intervals or equal number intervals.
            
        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - Neither the age nore the diameter_number argument is provided.
            - Both the age and diameter_number arguments are provided.
            - Both the age_end and diameter_number_end arguments are provided.
            - The age argument is provided but is not a scalar.
            - The age_end argument is provided but is not a scalar. 
            - The age_interval is provided but is not a positive scalar
            - The age_interval provided is negative, or is greater than age - age_end
            - The diameter_number argument is not a pair of values, or any of them are less than 0
            - The diameter_number_end argument is not a pair of values, or any of them are less than 0
            - The diameter_number_interval argument is not a pair of values, or any of them are less than 0
            - The age_interval and diameter_number_interval arguments are both provided.
            - The diameter_number_interval provided is negative, or is greater than diameter_number - diameter_number_end
            - The ninterval is provided but is not an integer or is less than 1.
            - The ninterval is provided and either age_interval or diameter_number_interval is also provided
        """         
        # Determine whether we are going to do equal time intervals or equal number intervals
        age = kwargs.get('age', None)
        age_interval = kwargs.get('age_interval', None)
        diameter_number_interval = kwargs.get('diameter_number_interval', None)
        ninterval = kwargs.pop('ninterval', None)
        if age_interval is not None and diameter_number_interval is not None:
            raise ValueError("Cannot specify both ninterval and age_interval or diameter_number_interval")
        if ninterval is not None:
            if not isinstance(ninterval, int):
                raise TypeError("ninterval must be an integer")
            if ninterval < 1:
                raise ValueError("ninterval must be greater than zero") 
            if age_interval is not None or diameter_number_interval is not None:
                raise ValueError("Cannot specify both ninterval and age_interval or diameter_number_interval")
            
        is_age_interval = age_interval is not None or (ninterval is not None and diameter_number_interval is None)
        
        # Validate arguments using the production function validator first 
        if 'diameter_range' not in kwargs:
            kwargs['diameter_range'] = (self.smallest_crater, self.largest_crater)
        if 'area' not in kwargs:
            kwargs['area'] = self.surf.area.item()
        kwargs  = self.production._validate_sample_args(**kwargs)
       
        if is_age_interval:
            age = kwargs.get('age', None)
            if age is None:
                raise ValueError("Something went wrong! age should be set by self.production_validate_sample_args")
            age_end = kwargs.get('age_end', None)
            if age_end is None:
                raise ValueError("Something went wrong! age_end should be set by self.production_validate_sample_args")
            if age_interval is None:
                if ninterval is None:
                    ninterval = 1
                age_interval = (age - kwargs['age_end']) / ninterval
            else:
                if age_interval > age - age_end:
                    raise ValueError("age_interval must be less than age - age_end")
                elif age_interval <= 0:
                    raise ValueError("age_interval must be greater than zero")
                ninterval = int(np.ceil((age - age_end) / age_interval)) 
            
            kwargs['age_interval'] = age_interval
            kwargs['ninterval'] = ninterval
        else:
            diameter_number = kwargs.get('diameter_number', None)
            if diameter_number is None:
                raise ValueError("Something went wrong! diameter_number should be set by self.production_validate_sample_args")
            diameter_number_end = kwargs.get('diameter_number_end', None)
            if diameter_number_end is None:
                raise ValueError("Something went wrong! diameter_number_end should be set by self.production_validate_sample_args")
            if diameter_number_interval is None:
                if ninterval is None:
                    ninterval = 1
                diameter_number_interval = (diameter_number[0], (diameter_number[1] - diameter_number_end[1]) / ninterval)
            else:
                if len(diameter_number_interval) != 2:
                    raise ValueError("The 'diameter_number_interval' must be a pair of values in the form (D,N)")
                # Check to be sure that the diameter in the diameter_number_interval is the same as diameter_number. 
                # If not, we need to adjust the end diameter value it so that they match 
                diameter_number_interval = self.production._validate_csfd(*diameter_number_interval)
                if diameter_number_interval[0] != diameter_number[0]:
                    area = kwargs.get('area', None)
                    if area is None:
                        raise ValueError("Something went wrong! area should be set by self.production_validate_sample_args")
                    diameter_number_density_interval = (diameter_number_interval[0], diameter_number_interval[1] / area)
                    age_val = self.production.function_inverse(*diameter_number_density_interval)
                    diameter_number_density_interval = (diameter_number[0], self.production.function(diameter=diameter_number[0], age=age_val))
                    diameter_number_interval = (diameter_number[0], diameter_number_density_interval[1] * area)
                    
                if diameter_number_interval[1] >= diameter_number[1] - diameter_number_end[1]:
                    raise ValueError("diameter_number_interval must be less than diameter_number - diameter_number_end")
                if diameter_number_interval[1] <= 0:
                    raise ValueError("diameter_number_interval must be greater than zero")
                ninterval = int(np.ceil((diameter_number[1] - diameter_number_end[1]) / diameter_number_interval[1]))
                kwargs['diameter_number_interval'] = diameter_number_interval
                kwargs['ninterval'] = ninterval
                
        kwargs['is_age_interval'] = is_age_interval
        
        # Remove unecessary arguments that came out of the production._validate_sample_args method        
        kwargs.pop('diameter_range')
        kwargs.pop('area')
        kwargs.pop('return_age')
         
        return kwargs
    
    def save(self, **kwargs: Any) -> None:
        """
        Save the current simulation state to a file.
        """

        time_variables = {
            "current_age": self.current_age,
            "elapsed_time": self.elapsed_time,
            "elapsed_n1": self.elapsed_n1
            }
         
        save(self.surf, interval_number=self.interval_number, time_variables=time_variables, **kwargs)
        
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
        
        self.save()  
        if out_dir is None:
            out_dir = os.path.join(self.simdir, "vtk_files")
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(out_dir)
            
        data_file_list = glob(os.path.join(self.surf.data_dir, "*_*.nc"))
        data_file_list.append(self.surf.time_file)
        
        # This will suppress the warning issued by xarray starting in version 2023.12.0 about the change in the API regarding .dims
        # The API change does not affect the functionality of the code, so we can safely ignore the warning
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore", FutureWarning)
            warnings.simplefilter("ignore", DeprecationWarning) # Ignores a warning issued in bar.py
        
            # Save the surface data to a combined netCDF file
            with tempfile.TemporaryDirectory() as temp_dir:
                filename_pattern = ""
                for f in data_file_list:
                    filename_pattern += f"{f};"
                try:
                    extract_vtk(
                        filename_pattern=filename_pattern,
                        mesh_filename=self.surf.grid_file,
                        variable_list=['allOnVertices', 'allOnCells'], 
                        dimension_list=['maxEdges=','vertexDegree=','maxEdges2=','TWO='], 
                        combine=False,
                        include_mesh_vars=True,
                        ignore_time=False,
                        time="0:",
                        xtime='Time',
                        out_dir=out_dir)
                            
                except:
                    raise RuntimeError("Error in extract_vtk. Cannot export VTK files")
        
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
        .. code-block:: python
        
            # Apply default turbulence noise:
            sim = cratermaker.Simulation()
            sim.apply_noise()

            #Apply ridged noise with custom Parameters
            sim.apply_noise(model="ridged", noise_width=500e3, num_octaves=10, freq=1.5)
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

    @property
    def interval_number(self):
        """
        The index of the current time step. 
        """
        return self._interval_number

    @interval_number.setter
    def interval_number(self, value):
        if not isinstance(value, int):
            raise TypeError("interval_number must be an integer")
        if value < 0:
            raise ValueError("interval_number must be greater than or equal to zero")
        
        self._interval_number = value    
        
    @property
    def elapsed_time(self):
        """
        The elasped time in My since the start of the simulation.
        """
        return self._elapsed_time
    
    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = np.float64(value)
        
        
    @property
    def current_age(self):
        """
        The age of the current time step in My relative to the present from the chronology of the production function.
        """
        return self._current_age
    
    @current_age.setter
    def current_age(self, value):
        self._current_age = np.float64(value)
        
    @property
    def elapsed_n1(self):
        """
        The elapsed number of craters larger than 1 km in diameter.
        """
        return self._elapsed_n1
    
    @elapsed_n1.setter
    def elapsed_n1(self, value):
        self._elapsed_n1 = np.float64(value)
        