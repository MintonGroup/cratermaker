import numpy as np
from numpy.random import Generator
import xarray as xr
import os
from pathlib import Path
import shutil
from tqdm import tqdm
from glob import glob
from collections.abc import Sequence
from typing import Any
from numpy.typing import ArrayLike
import warnings
import yaml
from .target import Target
from .crater import Crater
from .surface import Surface, _save_surface
from ..utils.general_utils import _convert_for_yaml, _set_properties, _to_config, parameter
from ..utils.custom_types import FloatLike, PairOfFloats
from ..components.scaling import ScalingModel, get_scaling_model, available_scaling_models
from ..components.production import ProductionModel, get_production_model, available_production_models
from ..components.morphology import MorphologyModel, get_morphology_model, available_morphology_models

class Simulation:
    """
    This class orchestrates the processes involved in running a crater simulation.

    """
    def __init__(self, *, # Enforce keyword-only arguments
                 target_name: str = "Moon",
                 seed: int | None = None,
                 scaling_model: str = "richardson2009",
                 production_model: str | None = None,
                 morphology_model: str = "simplemoon",
                 reset_surface: bool = True,
                 simdir: os.PathLike = Path.cwd(), 
                 **kwargs: Any):
        """
        Initialize the Simulation object.

        Parameters
        ----------
        target_name: str, optional, default "Moon"
            Name of the target body or Target object for the simulation, default is "Moon".
        seed : int, optional
            Seed for the random number generator, default is None which will generate a random seed.
        scaling_model : str, optional
            The name of the impactor->crater size scaling model to use from the components library. The default is "richardson2009".
        production_model: str, optional
            The name of the production function model to use from the components library that defines the production function used to populate the surface with craters. If none provided, 
            then the default will be based on the target body, with the NeukumProduction crater-based scaling law used if the target 
            body is the Moon or Mars, the NeukumProduction projectile-based scaling law if the target body is Mercury, Venus, or 
            Earth, and a simple power law model otherwise.
        morphology_model : str, optional
            The name of the component model used to describe the morphology of the crater. If none provided, then the default will "simplemoon", which is similar to the one used by CTEM.
        reset_surface : bool, optional
            Flag to reset the surface elevation, default is True. If False, the surface will be preserved and the simulation will
            continue from the last saved state.
        simdir: PathLike, optional
            Path to the simulation directory, default is current working directory.
        **kwargs : Any
            Additional keyword arguments that can be passed to other cratermaker components, such as arguments to set the surface, scale, 
            morphology, or production function constructors. Refer to the documentation of each component module for details.
            
        See Also
        --------
        cratermaker.core.surface.Surface.initialize : Parameters for initializing a surface mesh.
        """
        object.__setattr__(self, "_target_name", None)
        object.__setattr__(self, "_target", None)
        object.__setattr__(self, "_production_model", None)
        object.__setattr__(self, "_production", None)
        object.__setattr__(self, "_scaling_model", None)
        object.__setattr__(self, "_scale", None)
        object.__setattr__(self, "_morphology_model", None)
        object.__setattr__(self, "_morphology", None)
        object.__setattr__(self, "_craterlist", None)
        object.__setattr__(self, "_crater", None)
        object.__setattr__(self, "_interval_number", None)
        object.__setattr__(self, "_elapsed_time", None)
        object.__setattr__(self, "_current_age", None)
        object.__setattr__(self, "_elapsed_n1", None)
        object.__setattr__(self, "_smallest_crater", None)
        object.__setattr__(self, "_smallest_projectile", None)
        object.__setattr__(self, "_largest_crater", None)
        object.__setattr__(self, "_largest_projectile", None)
        object.__setattr__(self, "_simdir", None)
        object.__setattr__(self, "_rng", None)
        object.__setattr__(self, "_surf", None)
        object.__setattr__(self, "_seed", None)
        object.__setattr__(self, "_reset_surface", None)

        self.simdir = simdir
        
        if self.config_file.exists():
            filename = self.config_file
        else:
            filename = None
        matched, unmatched = _set_properties(self, target_name=target_name, seed=seed, scaling_model=scaling_model, production_model=production_model, morphology_model=morphology_model, reset_surface=reset_surface, filename=filename)
        production_model_parameters = unmatched.pop("production_model_parameters", {})
        scaling_model_parameters = unmatched.pop("scaling_model_parameters", {})
        surface_parameters = unmatched.pop("surface_parameters", {})
        morphology_model_parameters = unmatched.pop("morphology_model_parameters", {})
        target_parameters = unmatched.pop("target_parameters", {})
        if "name" in target_parameters:
            self.target_name = target_parameters.pop("name")

        kwargs.update(unmatched)

        self.rng = np.random.default_rng(self.seed)

        try:
            self.target = Target(name=self.target_name, **target_parameters, **kwargs)
        except:
            raise ValueError(f"Error initializing target {self.target_name}")
         
        # Set the production function model for this simulation 
        if production_model is None:
            if self.target_name in ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars']:
                self.production_model = "neukum"
                if self.target_name in ['Moon', 'Mars']:
                    self.production = get_production_model(self.production_model)(version=self.target_name, rng=self.rng, **kwargs)
                else:
                    self.production = get_production_model(self.production_model)(version="projectile", rng=self.rng, **kwargs)
            else:
                self.production = get_production_model("powerlaw")(rng=self.rng, **kwargs)
        else:
            self.production_model = production_model
        try:
            self.production = get_production_model(self.production_model)(rng=self.rng, **production_model_parameters, **kwargs)
        except:
            raise ValueError(f"Error initializing production model {self.production_model}")
            
        # Set the scaling law model for this simulation 
        self.scaling_model = scaling_model
        try:
            self.scale = get_scaling_model(self.scaling_model)(target=self.target, rng=self.rng, **scaling_model_parameters, **kwargs)
        except:
            raise ValueError(f"Error initializing scaling model {self.scaling_model}")
      
        self.morphology_model = morphology_model
        try:
            self.morphology = get_morphology_model(self.morphology_model)(rng=self.rng, **morphology_model_parameters, **kwargs)
        except:
            raise ValueError(f"Error initializing {morphology_model}")
        
        grid_type = kwargs.get('grid_type', None)
        if grid_type is not None and grid_type == 'hires local':
            if 'superdomain_scale_factor' not in kwargs:
                # Determine the scale factor for the superdomain based on the smallest crater whose ejecta can reach the edge of the 
                # superdomain. This will be used to set the superdomain scale factor. TODO: Streamline this a bit
                for d in np.logspace(np.log10(self.target.radius*2), np.log10(self.target.radius / 1e6), 1000):
                    crater, _ = self.generate_crater(diameter=d, angle=90.0, projectile_velocity=self.scale.projectile_mean_velocity*10)
                    rmax = crater.morphology.compute_rmax(minimum_thickness=1e-3) 
                    if rmax < self.target.radius * 2 * np.pi:
                        superdomain_scale_factor = rmax / crater.final_radius
                        break
                kwargs['superdomain_scale_factor'] = superdomain_scale_factor
        self.surf = Surface.initialize(target=self.target, reset_surface=self.reset_surface, simdir=self.simdir, rng=self.rng, **surface_parameters, **kwargs)

        self._craterlist = []
        self._crater = None
        self._projectile = None
        self._interval_number = 0
        self._elapsed_time = 0.0
        self._current_age = 0.0
        self._elapsed_n1 = 0.0
        self._smallest_crater = 0.0 # The smallest crater will be determined by the smallest face area
        self._smallest_projectile = 0.0 # The smallest crater will be determined by the smallest face area
        self._largest_crater = np.inf # The largest crater will be determined by the target body radius
        self._largest_projectile = np.inf # The largest projectile will be determined by the target body radius

        return

    
    def get_smallest_diameter(self, 
                              face_areas: ArrayLike | None = None, 
                              from_projectile: bool = False) -> float:
        """
        Get the smallest possible crater or projectile be formed on the surface.
        """
        if face_areas is None:
            face_areas = self.surf.face_areas
        else:
            face_areas = np.asarray(face_areas)
        smallest_crater = np.sqrt(face_areas.min().item() / np.pi) * 2        
        if from_projectile:
            crater = self.generate_crater(final_diameter=smallest_crater, angle=90.0, projectile_velocity=self.scale.projectile_mean_velocity*10)
            return crater.projectile_diameter 
        else:
            return smallest_crater 
        
    
    def get_largest_diameter(self,
                             from_projectile: bool = False) -> float:
        """
        Get the largest possible crater or projectile that can be formed on the surface.
        """
        largest_crater = self.target.radius * 2
        if from_projectile:
            crater = self.generate_crater(final_diameter=largest_crater, angle=1.0, projectile_velocity=self.scale.projectile_mean_velocity/10.0)
            return crater.projectile_diameter
        else:
            return largest_crater


    def generate_crater(self, 
                        **kwargs: Any
                       ) -> Crater:
        """
        Create a new Crater object 

        Parameters
        ----------
        **kwargs : Any
            Keyword arguments for initializing the :class:`Crater` object. Refer to 
            its documentation for details on valid keyword arguments.

        Returns
        -------
        Crater
            The newly created Crater object
            
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
        
        if "projectile_velocity" not in kwargs and "projectile_mean_velocity" not in kwargs and "projectile_vertical_velocity" not in kwargs:
            pmv = self.scale.projectile_mean_velocity
            if pmv is None:
                raise RuntimeError("No projectile_velocity value is set for this projectile")
            kwargs['projectile_mean_velocity'] = pmv
        crater = Crater.from_args(target=self.target, scale=self.scale, rng=self.rng, **kwargs)
        
        return crater
    
    
    def emplace_crater(self, **kwargs: Any
                      ) -> None:
        """
        Emplace a crater in the simulation, optionally based on a projectile.

        This method orchestrates the creation and placement of a crater in the
        simulation. It can create a crater directly or based on the characteristics
        of a projectile.

        Parameters
        ----------
        **kwargs : Any
            Keyword arguments for initializing the :class:`Crater.
            Refer to the documentation of this class for details on valid keyword arguments.

        Notes
        -----
        The keyword arguments provided are passed down to :meth:`generate_crater` 
        or :meth:`generate_projectile`, and subsequently to the constructor of 
        :class:`Crater` . Refer to its documentation for a detailed description of valid 
        keyword arguments.

        Examples
        --------
        .. code-block:: python        
        
            # Create a crater with specific diameter
            sim.emplace_crater(final_diameter=1000.0)

            # Create a crater based on a projectile with given mass and projectile_velocity
            sim.emplace_crater(projectile_mass=1e14, projectile_velocity=20e3)
            
            # Create a crater with a specific transient diameter and location
            sim.emplace_crater(transient_diameter=5e3, location=(43.43, -86.92))
        """ 
        self.crater = self.generate_crater(**kwargs)
        self.morphology.form_crater(self.surf,self.crater,**kwargs)
        
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
       
        from_projectile = self.production.generator_type == 'projectile'
        if from_projectile:
            diam_key = 'projectile_diameter'
        else:
            diam_key = 'final_diameter'
        Dmax = self.get_largest_diameter(from_projectile=from_projectile)
        Dmin = self.get_smallest_diameter(from_projectile=from_projectile)
        
        # Loop over each face in the mesh to build up a population of craters in this interval. This is done because faces may
        # not all have the same surface area, the range of crater sizes that can be formed on each face may be different.
        impact_diameters = []
        impact_ages = []
        impact_locations = []
        face_areas = self.surf.face_areas
        min_area = face_areas.min()
        n_face = face_areas.size
        surface_area = self.surf.area.item() 
         
        # Group surfaces into bins based on their area. All bins within a factor of 2 in surface area are grouped together.
        max_bin_index = np.ceil(np.log2(face_areas.max() / min_area)).astype(int)
        bins = {i: [] for i in range(max_bin_index + 1)}
        
        for face_index, area in enumerate(face_areas):
            bin_index = np.floor(np.log2(area / min_area)).astype(int)
            bins[bin_index].append(face_index)        
            
        # Process each bin
        for bin_index, face_indices in bins.items():
            if not face_indices:
                continue  # Skip empty bins
            face_indices = np.array(face_indices)

            bin_areas = face_areas[face_indices]
            total_bin_area = bin_areas.sum()
            area_ratio = total_bin_area / surface_area 
             
            Dmin = self.get_smallest_diameter(bin_areas, from_projectile=from_projectile)
            if diameter_number is not None:
                diameter_number_local = (diameter_number[0], diameter_number[1] * area_ratio)
            else:
                diameter_number_local = None
            if diameter_number_end is not None:
                diameter_number_end_local = (diameter_number_end[0], diameter_number_end[1] * area_ratio)
            else:
                diameter_number_end_local = None
        
            diameters, ages= self.production.sample(age=age, 
                                                    age_end=age_end, 
                                                    diameter_number=diameter_number_local, 
                                                    diameter_number_end=diameter_number_end_local, 
                                                    diameter_range=(Dmin, Dmax),
                                                    area=total_bin_area,
                                                    **kwargs)
            if diameters.size > 0:
                impact_diameters.extend(diameters.tolist())
                impact_ages.extend(ages.tolist())
                
                # Get the probability of impact onto any particular face then get the locations of the impacts
                p = bin_areas / total_bin_area
                locations = []
                generated = self.rng.choice(face_indices, size=diameters.shape)
                for i in range(len(diameters)): 
                    face_index = generated[i]
                    locations = self.surf.get_random_location_on_face(face_index)
                    impact_locations.append(locations) 
            
        if len(impact_diameters) > 0: 
            if len(impact_diameters) == 1:
                diam_arg = {diam_key: impact_diameters[0]}
                self.emplace_crater(age=impact_ages[0], location=impact_locations[0], **diam_arg)
            else:
                # Sort the ages, diameters, and locations so that they are in order of decreasing age
                sort_indices = np.argsort(impact_ages)[::-1]
                impact_diameters = np.asarray(impact_diameters)[sort_indices]
                impact_ages = np.asarray(impact_ages)[sort_indices] 
                impact_locations = np.array(impact_locations, dtype=[('lon', 'float64'), ('lat', 'float64')])[sort_indices]
                
                for i, diameter in tqdm(enumerate(impact_diameters), total=len(impact_diameters)):
                    location = impact_locations[i][0], impact_locations[i][1]
                    diam_arg = {diam_key: diameter}
                    self.emplace_crater(age=impact_ages[i], location=location, from_projectile=from_projectile, **diam_arg)
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
            kwargs['diameter_range'] = (self.get_smallest_diameter(), self.get_largest_diameter())
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


    def to_config(self, **kwargs: Any) -> dict:
        sim_config = _to_config(self)
        sim_config['target_parameters'] = self.target.to_config()
        sim_config['scaling_model_parameters'] = self.scale.to_config()
        sim_config['production_model_parameters'] = self.production.to_config()
        sim_config['surface_parameters'] = self.surf.to_config()
        # Write the combined configuration to a YAML file
        with open(self.config_file, 'w') as f:
            yaml.safe_dump(sim_config, f, indent=4)

        return _to_config(self)

    def save(self, **kwargs: Any) -> None:
                
        """
        Save the current simulation state to a file.
        """

        time_variables = {
            "current_age": self.current_age,
            "elapsed_time": self.elapsed_time,
            "elapsed_n1": self.elapsed_n1
            }
         
        _save_surface(self.surf, interval_number=self.interval_number, time_variables=time_variables, **kwargs)

        self.to_config(**kwargs)
        
        return
    
    
    def export_vtk(self, 
                   out_dir: os.PathLike | None = None,
                   *args, **kwargs
                   ) -> None:
        """
        Export the surface mesh to a VTK file.

        Parameters
        ----------
        out_dir : str, Default "vtk_files" in the simulation directory
            Directory to store the VTK files.
        """
        from vtk import vtkUnstructuredGrid, vtkPoints, VTK_POLYGON, vtkWarpScalar, vtkXMLPolyDataWriter
        from vtkmodules.util.numpy_support import numpy_to_vtk
        from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
        from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
        
        self.save()  
        if out_dir is None:
            out_dir = os.path.join(self.simdir, "vtk_files")
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(out_dir)
            
        data_file_list = glob(os.path.join(self.surf.data_dir, "*.nc"))
        if self.surf.grid_file in data_file_list:
            data_file_list.remove(self.surf.grid_file)
       
        # Convert uxarray grid arrays to regular numpy arrays for vtk processing 
        n_node = self.surf.data.uxgrid.n_node
        n_face = self.surf.data.uxgrid.n_face
        node_x = self.surf.data.uxgrid.node_x.values * self.target.radius
        node_y = self.surf.data.uxgrid.node_y.values * self.target.radius
        node_z = self.surf.data.uxgrid.node_z.values * self.target.radius
        n_nodes_per_face = self.surf.data.uxgrid.n_nodes_per_face.values
        face_node_connectivity = self.surf.data.uxgrid.face_node_connectivity.values
        
        vtk_data = vtkUnstructuredGrid()
        nodes = vtkPoints()
        for i in range(n_node):
            nodes.InsertNextPoint(node_x[i], node_y[i], node_z[i])
        vtk_data.SetPoints(nodes)
        vtk_data.Allocate(n_face)
        for i,n in enumerate(n_nodes_per_face):
            point_ids=face_node_connectivity[i][0:n]
            vtk_data.InsertNextCell(VTK_POLYGON, n, point_ids) 
       
        warp = vtkWarpScalar()
        warp.SetInputArrayToProcess(0, 0, 0,
                            vtkUnstructuredGrid.FIELD_ASSOCIATION_POINTS,
                            "node_elevation")                
            
        writer = vtkXMLPolyDataWriter()
        writer.SetDataModeToBinary()
        writer.SetCompressorTypeToZLib()
        print("Exporting VTK files...")
        
        with xr.open_mfdataset(data_file_list) as ds:
            # Warp the surface based on node_elevation data
            for i in tqdm(range(len(ds.time))):
                
                ids = ds.isel(time=i)
                current_grid = vtkUnstructuredGrid()
                current_grid.DeepCopy(vtk_data) 
                 
                for v in ds.variables:
                    array = numpy_to_vtk(ids[v].values, deep=True)
                    array.SetName(v)
                    n = ids[v].size
                    if 'n_face' in ids[v].dims:
                        current_grid.GetCellData().AddArray(array)
                    elif 'n_node' in ids[v].dims:
                        current_grid.GetPointData().AddArray(array)
                        if v == 'node_elevation':
                            current_grid.GetPointData().SetActiveScalars(v) 
                    elif n == 1:
                        current_grid.GetFieldData().AddArray(array)
                        
                geomFilter = vtkGeometryFilter()
                geomFilter.SetInputData(current_grid)
                geomFilter.Update()    
                polyData = geomFilter.GetOutput()    

                normalsFilter = vtkPolyDataNormals()
                normalsFilter.SetInputData(polyData)
                normalsFilter.ComputeCellNormalsOn()
                normalsFilter.ConsistencyOn()           # Tries to make normals consistent across shared edges
                normalsFilter.AutoOrientNormalsOn()     # Attempt to orient normals consistently outward/inward
                normalsFilter.SplittingOff()   
                normalsFilter.Update()
                polyDataWithNormals = normalsFilter.GetOutput()        
        
                warp.SetInputData(polyDataWithNormals)
                warp.Update()
                warped_output = warp.GetOutput()
                output_filename = os.path.join(out_dir, f"surf{i:06d}.vtp")
                writer.SetFileName(output_filename)
                writer.SetInputData(warped_output) 
                writer.Write()               
        
        return


    def make_circle_file(self,
                         diameters: FloatLike | Sequence[FloatLike] | ArrayLike ,
                         longitudes: FloatLike | ArrayLike,
                         latitudes: FloatLike | ArrayLike,
                         output_filename: os.PathLike | None = None,
                         *args, **kwargs
                        ) -> None:
        """
        Plot circles of diameter D centered at the given location.
    
        Parameters
        ----------
        diameters : FloatLike or ArrayLike
            Diameters of the circles in m.
        longitudes : FloatLike or ArrayLike of Floats
            Longitudes of the circle centers in degrees.
        latitudes : FloatLike or ArrayLike of Floats
            Latitudes of the circle centers in degrees.
        out_dir : str, Default "vtk_files" in the simulation directory
            Directory to store the VTK files.
        """ 
        import vtk 

        if output_filename is None:
            output_filename = os.path.join(self.simdir, "circles.vtp") 
        else:
            output_filename = os.path.join(self.simdir, output_filename)
        
        diameters = np.atleast_1d(diameters)
        longitudes = np.atleast_1d(longitudes)
        latitudes = np.atleast_1d(latitudes)
        # Check for length consistency
        if len(diameters) != len(longitudes) or len(diameters) != len(latitudes):
                raise ValueError("The diameters, latitudes, and longitudes, arguments must have the same length")
            
        # Validate non-negative values
        if np.any(diameters < 0):
            raise ValueError("All values in 'diameters' must be non-negative")
        
        sphere_radius = self.target.radius 
        def create_circle(lon, lat, circle_radius, num_points=360):
            """
            Create a circle on the sphere's surface with a given radius and center.
            
            Parameters
            ----------
            lon : float
                Longitude of the circle's center in degrees.
            lat : float
                Latitude of the circle's center in degrees.
            circle_radius : float
                Radius of the circle in meters.
            num_points : int, optional
                Number of points to use to approximate the circle. The default is 360. 
            """
            # Create an array of angle steps for the circle
            radians = np.linspace(0, 2 * np.pi, num_points)
            # Convert latitude and longitude to radians
            lat_rad = np.deg2rad(lat)
            lon_rad = np.deg2rad(lon)

            # Calculate the Cartesian coordinates for the circle's center
            center_x = sphere_radius * np.cos(lat_rad) * np.cos(lon_rad)
            center_y = sphere_radius * np.cos(lat_rad) * np.sin(lon_rad)
            center_z = sphere_radius * np.sin(lat_rad)

            # Calculate the vectors for the local east and north directions on the sphere's surface
            east = np.array([-np.sin(lon_rad), np.cos(lon_rad), 0])
            north = np.array([-np.cos(lon_rad)*np.sin(lat_rad), -np.sin(lon_rad)*np.sin(lat_rad), np.cos(lat_rad)])
            
            # Initialize arrays to hold the circle points
            x = np.zeros_like(radians)
            y = np.zeros_like(radians)
            z = np.zeros_like(radians)

            # Calculate the points around the circle
            for i in range(num_points):
                x[i] = center_x + circle_radius * np.cos(radians[i]) * east[0] + circle_radius * np.sin(radians[i]) * north[0]
                y[i] = center_y + circle_radius * np.cos(radians[i]) * east[1] + circle_radius * np.sin(radians[i]) * north[1]
                z[i] = center_z + circle_radius * np.cos(radians[i]) * east[2] + circle_radius * np.sin(radians[i]) * north[2]

            return x, y, z 

        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        point_id = 0  # Keep track of the point ID across all circles
        
        for lon, lat, diameter in zip(longitudes, latitudes, diameters):
            circle_radius = diameter / 2
            x, y, z = create_circle(lon, lat, circle_radius)
            
            for i in range(len(x)):
                points.InsertNextPoint(x[i], y[i], z[i])
            
            polyline = vtk.vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(len(x))
            for i in range(len(x)):
                polyline.GetPointIds().SetId(i, point_id)
                point_id += 1
            
            lines.InsertNextCell(polyline)
        
        # Create a polydata object and add points and lines to it
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        
        # Write the polydata to a VTK file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_filename)
        writer.SetInputData(polydata)
        
        # Optional: set the data mode to binary to save disk space
        writer.SetDataModeToBinary()
        writer.Write()
        
        return    
    

    def apply_noise(self, 
                    model="turbulence",
                    noise_width=1000e3,
                    noise_height=1e3,
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
            raise NotImplementedError
            # with warnings.catch_warnings():
            #     warnings.simplefilter("ignore", Warning)           
            #     x = self.surf.uxgrid[vars[0]].values * scale 
            #     y = self.surf.uxgrid[vars[1]].values * scale 
            #     z = self.surf.uxgrid[vars[2]].values * scale
            #     noise = realistic.apply_noise(model, x, y, z, num_octaves, anchor, **kwargs)
        
            # # Make sure the noise is volume-conserving (i.e., the mean is zero)
            # # TODO: Take into account the nodes are not uniformly distributed on the sphere
            # noise = noise - np.mean(noise)
            # return noise                
        
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

    @parameter
    def production_model(self):
        """
        The production model used for crater production. Set during initialization.
        """
        return self._production_model

    @production_model.setter
    def production_model(self, value):
        if not isinstance(value, str):
            raise TypeError("production_model must be a string")
        if value.lower() not in available_production_models():
            raise ValueError(f"production_model must be one of {available_production_models()}")
        self._production_model = value.lower()

    @parameter
    def target_name(self):
        """
        The name of the target body for the impact simulation. Set during initialization.
        """
        return self._target_name
    
    @target_name.setter
    def target_name(self, value):
        if not isinstance(value, str):
            raise TypeError("target_name must be a string")
        self._target_name = value.title()

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

    @parameter
    def seed(self):
        """
        Seed for the random number generator. Set during initialization.
        """
        return self._seed

    @seed.setter
    def seed(self, value):
        if value is None:
            value = np.random.default_rng().integers(0, 2**32)
        elif not isinstance(value, int):
            raise TypeError("seed must be an integer or None")
        self._seed = value

    @property
    def rng(self):
        """
        Random number generator instance. Set during initialization.
        """
        return self._rng

    @rng.setter
    def rng(self, value):
        if value is None:
            value = np.random.default_rng(seed=self.seed) 
        elif not isinstance(value, Generator): 
            raise TypeError("rng must be a numpy.random.Generator instance or None")
        self._rng = value

    @parameter
    def reset_surface(self):
        """
        Flag to reset the surface before running the simulation. Set during initialization.
        """
        return self._reset_surface

    @reset_surface.setter
    def reset_surface(self, value):
        if not isinstance(value, bool):
            raise TypeError("reset_surface must be a boolean")
        self._reset_surface = value
    
    @parameter
    def simdir(self):
        """
        Directory where the simulation data is stored. Set during initialization.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value):
        if value is None:
            self._simdir = Path.cwd()
        else:
            p = Path(value)
            if not p.is_absolute():
                p = Path.cwd() / p
            self._simdir = p

        self._simdir.mkdir(parents=True, exist_ok=True)

    @property
    def production(self):
        """
        The Production class instance used for crater production. Set during initialization.
        """
        return self._production

    @production.setter
    def production(self, value):
        if not issubclass(value.__class__, ProductionModel):
            raise TypeError("production must be a subclass of Production")
        self._production = value

    @parameter
    def scaling_model(self):
        """
        The scaling model used for crater scaling. Set during initialization.
        """
        return self._scaling_model

    @scaling_model.setter
    def scaling_model(self, value):
        if not isinstance(value, str):
            raise TypeError("scaling_model must be a string")
        if value.lower() not in available_scaling_models():
            raise ValueError(f"scaling_model must be one of {available_scaling_models()}")
        self._scaling_model = value.lower()

    @property
    def scale(self):
        """
        The ScalingModel object that defines the crater scaling relationships model. Set during initialization.
        """
        return self._scale

    @scale.setter
    def scale(self, value):
        if not isinstance(value, ScalingModel):
            raise TypeError("scale must be of ScalingModel type")
        self._scale = value

    @property
    def morphology(self):
        """
        The crater morphology model. Set during initialization.
        """
        return self._morphology

    @morphology.setter
    def morphology(self, value):
        if not isinstance(value, MorphologyModel):
            raise TypeError("morpholog must be of MorphologyModel type")
        self._morphology = value

    @parameter
    def morphology_model(self):
        """
        The name of the morphology model to load from the components library.
        """
        return self._morphology_model

    @morphology_model.setter
    def morphology_model(self, value):
        if not isinstance(value, str):
            raise TypeError("morphology_model must be a string")
        if value.lower() not in available_morphology_models():
            raise ValueError(f"morphology_model must be one of {available_morphology_models()}")
        self._morphology_model = value.lower()

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

    @parameter
    def target_name(self):
        """
        The name of the target body for the impact simulation. Set during initialization.
        """
        return self._target_name

    @target_name.setter
    def target_name(self, value):
        if not isinstance(value, str):
            raise TypeError("target_name must be a string")
        self._target_name = value.title()

    @parameter
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
        
    @parameter
    def elapsed_time(self):
        """
        The elasped time in My since the start of the simulation.
        """
        return self._elapsed_time
    
    @elapsed_time.setter
    def elapsed_time(self, value):
        self._elapsed_time = float(value)
        
    @parameter
    def current_age(self):
        """
        The age of the current time step in My relative to the present from the chronology of the production function.
        """
        return self._current_age
    
    @current_age.setter
    def current_age(self, value):
        self._current_age = float(value)
        
    @parameter
    def elapsed_n1(self):
        """
        The elapsed number of craters larger than 1 km in diameter.
        """
        return self._elapsed_n1
    
    @elapsed_n1.setter
    def elapsed_n1(self, value):
        self._elapsed_n1 = float(value)
       
    @parameter
    def smallest_crater(self):
        """
        The smallest crater diameter in meters. Set during initialization.
        """
        return self._smallest_crater
    
    @smallest_crater.setter
    def smallest_crater(self, value):
        if value is None:
            self._smallest_crater = 0.0
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("smallest_crater must be a scalar value")
        elif value < 0:
            raise ValueError("smallest_crater must be greater than or equal to zero")
        elif self._largest_crater is not None and value > self._largest_crater:
            raise ValueError("smallest_crater must be less than or equal to largest_crater")
        self._smallest_crater = float(value)
        
    @parameter
    def largest_crater(self):
        """
        The largest crater diameter in meters. Set during initialization.
        """
        return self._largest_crater
    
    @largest_crater.setter
    def largest_crater(self, value):
        if value is None or np.isinf(value):
            self._largest_crater = np.inf
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("largest_crater must be a scalar value")
        elif value <= 0:
            raise ValueError("largest_crater must be greater than zero")
        elif self._smallest_crater is not None and value < self._smallest_crater:
            raise ValueError("largest_crater must be greater than or equal to smallest_crater")
        self._largest_crater = float(value)
       
    @parameter
    def smallest_projectile(self):
        """
        The smallest projectile diameter in meters. Set during initialization.
        """
        return self._smallest_projectile
    
    @smallest_projectile.setter
    def smallest_projectile(self, value):
        if value is None:
            self._smallest_projectile = 0.0
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("smallest_projectile must be a scalar value")
        elif value < 0:
            raise ValueError("smallest_projectile must be greater or equal to zero")
        elif self._largest_projectile is not None and value > self._largest_projectile:
            raise ValueError("smallest_projectile must be less than or equal to largest_projectile")  
        self._smallest_projectile = float(value)
        
    @parameter
    def largest_projectile(self):
        """
        The largest projectile diameter in meters. Set during initialization.
        """
        return self._largest_projectile
    
    @largest_projectile.setter
    def largest_projectile(self, value):
        if value is None or np.isinf(value):
            self._largest_projectile = np.inf
            return
        elif not isinstance(value, FloatLike):
            raise TypeError("largest_projectile must be a scalar value")
        elif value <= 0:
            raise ValueError("largest_projectile must be greater than zero")
        elif self._smallest_projectile is not None and value < self._smallest_projectile:
            raise ValueError("largest_projectile must be greater than or equal to smallest_projectile")
        self._largest_projectile = float(value)

    @property
    def name(self):
        """
        The name of the simulation. 
        """
        return "Cratermaker Simulation object"
    
    @property
    def config_file(self):
        """
        The path to the configuration file for the simulation.
        """
        return self.simdir / "cratermaker.yaml" 

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        # Avoid recursive calls during initialization or early access
        if hasattr(self, "to_config") and callable(getattr(self, "to_config", None)):
            if _convert_for_yaml(value) is not None:
                try:
                    self.to_config()
                except Exception:
                    pass