.. currentmodule:: cratermaker

.. _api:

#########################
Cratermaker API reference
#########################

This section of the documentation provides a detailed reference for the core functionality and modular components of the Cratermaker project.

.. _api-Simulation:


Simulation
==========

The Simulation class is the main class for the cratermaker project. It is used to create a simulation of a crater on a given target 
body. The Simulation class is used to generate craters of a given size and morphology based on the production function, morphology
function, and crater scaling relationship model. The surface of the target body is represented by a Surface attribute called
`surface`, which contains a UxDataset object called `surface.uxds`. This is an unstructured grid dataset that contains data for the target body surface.

Creating a Simulation
---------------------

.. autosummary::
    :toctree: generated/

    Simulation

Methods
-------

.. autosummary::
    :toctree: generated/

    Simulation.run
    Simulation.populate
    Simulation.emplace_crater
    Simulation.generate_crater
    Simulation.save
    Simulation.export
    Simulation.set_elevation
    Simulation.to_config
    Simulation.get_smallest_diameter
    Simulation.get_largest_diameter


Attributes
----------

.. autosummary::
    :toctree: generated/

    Simulation.surface
    Simulation.production
    Simulation.scaling
    Simulation.crater
    Simulation.projectile
    Simulation.target
    Simulation.morphology
    Simulation.crater
    Simulation.data_dir
    Simulation.n_face
    Simulation.n_node
    Simulation.rng
    Simulation.rng_seed
    Simulation.simdir
    Simulation.interval_number
    Simulation.elapsed_time
    Simulation.current_age
    Simulation.elapsed_n1
    Simulation.smallest_crater
    Simulation.largest_crater
    Simulation.smallest_projectile
    Simulation.largest_projectile


.. _api-Crater:

Crater
======

The ``Crater`` class represents a single crater in the simulation. It is used to model the crater resulting from an impact, including its size, shape, depth, and other morphological features. It also defines the properties of the projectile, such as its size, velocity, material, and angle of impact.


Methods
-------

.. autosummary::
    :toctree: generated/

    Crater.maker

Attributes
----------

.. autosummary::
    :toctree: generated/

    Crater.final_diameter
    Crater.final_radius
    Crater.transient_diameter
    Crater.transient_radius
    Crater.projectile_diameter
    Crater.projectile_radius
    Crater.projectile_density
    Crater.projectile_mass
    Crater.projectile_velocity
    Crater.projectile_vertical_velocity
    Crater.projectile_angle
    Crater.projectile_direction
    Crater.location
    Crater.age

.. _api-Target:

Target
======

The Target class represents the target body in a crater simulation. It encapsulates properties of the target, such as its material composition, size, and other physical characteristics.

Creating a Target
-----------------

.. autosummary::
    :toctree: generated/

    Target.maker

Methods
-------

.. autosummary::
    :toctree: generated/

    Target.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    Target.name
    Target.material_name
    Target.catalogue
    Target.diameter
    Target.escape_velocity
    Target.gravity
    Target.radius
    Target.transition_scale_type

.. _api-Surface:

.. currentmodule:: cratermaker.components.surface


Surface
=======

The Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

Creating a surface
------------------

.. autosummary::
    :toctree: generated/

    Surface.maker

Methods
-------

.. autosummary::
    :toctree: generated/

    Surface.available
    Surface.load_from_files
    Surface.save_to_files
    Surface.generate_grid
    Surface.regrid_if_needed
    Surface.reset
    Surface.generate_data
    Surface.set_elevation
    Surface.calculate_haversine_distance
    Surface.calculate_initial_bearing
    Surface.find_nearest_index
    Surface.get_reference_surface
    Surface.get_distance
    Surface.get_initial_bearing
    Surface.elevation_to_cartesian
    Surface.extract_region
    Surface.get_random_location_on_face
    Surface.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    Surface.uxds
    Surface.uxgrid
    Surface.node_tree
    Surface.face_tree
    Surface.data_dir
    Surface.grid_file
    Surface.smallest_length
    Surface.area
    Surface.target
    Surface.rng
    Surface.rng_seed
    Surface.simdir

.. _api-IcosphereSurface:

.. currentmodule:: cratermaker.components.surface.icosphere

Icosphere grid
==============

.. autosummary::
    :toctree: generated/

    IcosphereSurface.maker

Methods
--------

.. autosummary::
    :toctree: generated/

    IcosphereSurface.maker
    IcosphereSurface.available
    IcosphereSurface.load_from_files
    IcosphereSurface.save_to_files
    IcosphereSurface.generate_grid
    IcosphereSurface.regrid_if_needed
    IcosphereSurface.reset
    IcosphereSurface.generate_data
    IcosphereSurface.set_elevation
    IcosphereSurface.calculate_haversine_distance
    IcosphereSurface.calculate_initial_bearing
    IcosphereSurface.find_nearest_index
    IcosphereSurface.get_reference_surface
    IcosphereSurface.get_distance
    IcosphereSurface.get_initial_bearing
    IcosphereSurface.elevation_to_cartesian
    IcosphereSurface.extract_region
    IcosphereSurface.get_random_location_on_face
    IcosphereSurface.to_config
    IcosphereSurface.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    IcosphereSurface.uxds
    IcosphereSurface.uxgrid
    IcosphereSurface.node_tree
    IcosphereSurface.face_tree
    IcosphereSurface.data_dir
    IcosphereSurface.grid_file
    IcosphereSurface.smallest_length
    IcosphereSurface.area
    IcosphereSurface.target
    IcosphereSurface.rng
    IcosphereSurface.rng_seed
    IcosphereSurface.simdir
    IcosphereSurface.name
    IcosphereSurface.radius
    IcosphereSurface.gridlevel

Usage example
-------------

.. code-block:: python

    from cratermaker import Surface
    surface = Surface.maker("icosphere",gridlevel=7)

.. _api-ArbitraryResolutionSurface:

.. currentmodule:: cratermaker.components.surface.arbitrary_resolution

Arbitrary resolution grid
=========================

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionSurface.maker

Methods
--------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionSurface.maker
    ArbitraryResolutionSurface.available
    ArbitraryResolutionSurface.load_from_files
    ArbitraryResolutionSurface.save_to_files
    ArbitraryResolutionSurface.generate_grid
    ArbitraryResolutionSurface.regrid_if_needed
    ArbitraryResolutionSurface.reset
    ArbitraryResolutionSurface.generate_data
    ArbitraryResolutionSurface.set_elevation
    ArbitraryResolutionSurface.calculate_haversine_distance
    ArbitraryResolutionSurface.calculate_initial_bearing
    ArbitraryResolutionSurface.find_nearest_index
    ArbitraryResolutionSurface.get_reference_surface
    ArbitraryResolutionSurface.get_distance
    ArbitraryResolutionSurface.get_initial_bearing
    ArbitraryResolutionSurface.elevation_to_cartesian
    ArbitraryResolutionSurface.extract_region
    ArbitraryResolutionSurface.get_random_location_on_face
    ArbitraryResolutionSurface.to_config
    ArbitraryResolutionSurface.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionSurface.uxds
    ArbitraryResolutionSurface.uxgrid
    ArbitraryResolutionSurface.node_tree
    ArbitraryResolutionSurface.face_tree
    ArbitraryResolutionSurface.data_dir
    ArbitraryResolutionSurface.grid_file
    ArbitraryResolutionSurface.smallest_length
    ArbitraryResolutionSurface.area
    ArbitraryResolutionSurface.target
    ArbitraryResolutionSurface.rng
    ArbitraryResolutionSurface.rng_seed
    ArbitraryResolutionSurface.simdir
    ArbitraryResolutionSurface.pix
    ArbitraryResolutionSurface.radius

Usage example
-------------

.. code-block:: python

    from cratermaker import Surface
    surface = Surface.maker("arbitrary_resolution",pix=100)

.. _api-HiResLocalSurface:

.. currentmodule:: cratermaker.components.surface.hireslocal

Hi-res local grid
=================

.. autosummary::
    :toctree: generated/

    HiResLocalSurface.maker

Methods
--------

.. autosummary::
    :toctree: generated/

    HiResLocalSurface.maker
    HiResLocalSurface.available
    HiResLocalSurface.load_from_files
    HiResLocalSurface.save_to_files
    HiResLocalSurface.generate_grid
    HiResLocalSurface.regrid_if_needed
    HiResLocalSurface.reset
    HiResLocalSurface.generate_data
    HiResLocalSurface.set_elevation
    HiResLocalSurface.calculate_haversine_distance
    HiResLocalSurface.calculate_initial_bearing
    HiResLocalSurface.find_nearest_index
    HiResLocalSurface.get_reference_surface
    HiResLocalSurface.get_distance
    HiResLocalSurface.get_initial_bearing
    HiResLocalSurface.elevation_to_cartesian
    HiResLocalSurface.extract_region
    HiResLocalSurface.get_random_location_on_face
    HiResLocalSurface.to_config
    HiResLocalSurface.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    HiResLocalSurface.uxds
    HiResLocalSurface.uxgrid
    HiResLocalSurface.node_tree
    HiResLocalSurface.face_tree
    HiResLocalSurface.data_dir
    HiResLocalSurface.grid_file
    HiResLocalSurface.smallest_length
    HiResLocalSurface.area
    HiResLocalSurface.target
    HiResLocalSurface.rng
    HiResLocalSurface.rng_seed
    HiResLocalSurface.simdir
    HiResLocalSurface.pix
    HiResLocalSurface.radius
    HiResLocalSurface.local_radius
    HiResLocalSurface.local_location
    HiResLocalSurface.superdomain_scale_factor

Usage example
-------------

.. code-block:: python

    from cratermaker import Surface
    surface = Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9))


.. _api-Production:

.. currentmodule:: cratermaker.components.production

Production
==========

.. autosummary::
    :toctree: generated/

    Production.maker


.. _api-NeukumProduction:

.. currentmodule:: cratermaker.components.production.neukum

NeukumProduction
================

.. autosummary::
    :toctree: generated/

    NeukumProduction.maker
    NeukumProduction.sample
    NeukumProduction.function
    NeukumProduction.chronology
    NeukumProduction.size_frequency_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    NeukumProduction.sfd_coef
    NeukumProduction.sfd_range
    NeukumProduction.valid_time
    NeukumProduction.tau
    NeukumProduction.Cexp
    NeukumProduction.Clin
    NeukumProduction.generator_type

Usage example
-------------

.. code-block:: python

    from cratermaker import Production
    production = Production.maker("neukum", version="projectile")


.. _api-PowerLawProduction:

.. currentmodule:: cratermaker.components.production.powerlaw

Power law production function
=============================

.. autosummary::
    :toctree: generated/

    PowerLawProduction.maker
    PowerLawProduction.sample
    PowerLawProduction.function
    PowerLawProduction.chronology

Attributes
----------

.. autosummary::
    :toctree: generated/

    PowerLawProduction.N1_coef
    PowerLawProduction.slope

Usage example
-------------

.. code-block:: python

    from cratermaker import Production
    production = Production.maker("powerlaw", slope=-4.0, N1_coef=1.0e-6)


.. _api-Scaling:

.. currentmodule:: cratermaker.components.scaling

Scaling
=======

The `Scaling` class is an abstract base class for crater scaling relationships. 
Use the `Scaling.maker` method to create a specific scaling model. 
Available models: see `Scaling.available()`.


Methods
-----------

.. autosummary::
    :toctree: generated/

    Scaling.maker
    Scaling.available

Attributes
----------

.. autosummary::
    :toctree: generated/

    Scaling.target
    Scaling.projectile


.. _api-MonteCarloScaling:

.. currentmodule:: cratermaker.components.scaling.montecarlo

Monte Carlo scaling model
=========================

This implements the scaling laws similar to those implemented in CTEM. However, unlike in CTEM, we apply monte carlo methods to the scaling laws to account for the uncertainty in the scaling laws.


Methods
--------

.. autosummary::
    :toctree: generated/

    MonteCarloScaling.final_to_transient
    MonteCarloScaling.transient_to_final
    MonteCarloScaling.projectile_to_transient
    MonteCarloScaling.transient_to_projectile
    MonteCarloScaling.get_morphology_type

Attributes
----------

.. autosummary::
    :toctree: generated/

    MonteCarloScaling.target
    MonteCarloScaling.projectile
    MonteCarloScaling.material_catalogue
    MonteCarloScaling.material_name
    MonteCarloScaling.K1
    MonteCarloScaling.mu
    MonteCarloScaling.Ybar
    MonteCarloScaling.catalogue_key
    MonteCarloScaling.transition_diameter
    MonteCarloScaling.transition_nominal
    MonteCarloScaling.complex_enlargement_factor
    MonteCarloScaling.simple_enlargement_factor
    MonteCarloScaling.final_exp

Usage example
-------------

.. code-block:: python

    from cratermaker import Scaling
    scaling_model = Scaling.maker("montecarlo", target="Mars", projectile="asteroids")

.. _api-CTEMScaling:

.. currentmodule:: cratermaker.components.scaling.ctem


CTEM Scaling model
===================

This implements the scaling laws similar to those implemented in CTEM. It is identical to the default scaling model, but it does not apply monte carlo methods to the scaling laws. 


Methods
--------

.. autosummary::
    :toctree: generated/

    CTEMScaling.final_to_transient
    CTEMScaling.transient_to_final
    CTEMScaling.projectile_to_transient
    CTEMScaling.transient_to_projectile
    CTEMScaling.get_morphology_type

Attributes
----------

.. autosummary::
    :toctree: generated/

    CTEMScaling.target
    CTEMScaling.projectile
    CTEMScaling.material_catalogue
    CTEMScaling.material_name
    CTEMScaling.K1
    CTEMScaling.mu
    CTEMScaling.Ybar
    CTEMScaling.catalogue_key
    CTEMScaling.transition_diameter
    CTEMScaling.transition_nominal
    CTEMScaling.complex_enlargement_factor
    CTEMScaling.simple_enlargement_factor
    CTEMScaling.final_exp

Usage example
-------------

.. code-block:: python

    from cratermaker import Scaling
    scaling_model = Scaling.maker("ctem", target="Mars", projectile="asteroids")


.. currentmodule:: cratermaker.components.morphology

.. _api-Morphology:

Morphology
==========

The Morphology class is an operations class for computing the morphology of a crater based on its size and target properties. It encapsulates the logic for altering the topography of the surface based on the crater properties.


.. autosummary::
    :toctree: generated/

    Morphology.maker
    Morphology.available

Attributes
----------

.. autosummary::
    :toctree: generated/

    Morphology.crater


.. currentmodule:: cratermaker.components.morphology.simplemoon

SimpleMoon
==========

.. autosummary::
    :toctree: generated/

    SimpleMoon.maker
    SimpleMoon.crater_profile
    SimpleMoon.ejecta_profile
    SimpleMoon.ejecta_distribution
    SimpleMoon.form_crater
    SimpleMoon.form_ejecta
    SimpleMoon.ray_intensity
    SimpleMoon.maker

Attributes
----------

.. autosummary::
    :toctree: generated/

    SimpleMoon.crater
    SimpleMoon.dorays
    SimpleMoon.ejrim
    SimpleMoon.ejecta_truncation
    SimpleMoon.floordepth
    SimpleMoon.floor_diameter
    SimpleMoon.rimheight
    SimpleMoon.rimwidth
    SimpleMoon.peakheight
    SimpleMoon.rng

Usage example
-------------

.. code-block:: python

    from cratermaker import Morphology
    morphology = Morphology.maker("simplemoon")


.. currentmodule:: cratermaker.components.projectile

.. _api-Projectile:

Projectile
==========

The Projectile class is an operations class defining the interface for generating projectile velocities, angles, and densities for a given target body.

.. autosummary::
    :toctree: generated/

    Projectile.maker
    Projectile.available

Attributes
----------

.. autosummary::
    :toctree: generated/

    Projectile.target_name
    Projectile.target
    Projectile.sample_angles
    Projectile.sample_velocities
    Projectile.sample_directions
    Projectile.mean_velocity
    Projectile.angle
    Projectile.direction
    Projectile.velocity
    Projectile.density
    Projectile.vertical_velocity


.. currentmodule:: cratermaker.components.projectile.asteroids

Asteroid projectiles
==================== 

.. autosummary::
    :toctree: generated/

    AsteroidProjectiles.maker
    AsteroidProjectiles.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    AsteroidProjectiles.target_name
    AsteroidProjectiles.target
    AsteroidProjectiles.sample_angles
    AsteroidProjectiles.sample_velocities
    AsteroidProjectiles.sample_directions
    AsteroidProjectiles.mean_velocity
    AsteroidProjectiles.angle
    AsteroidProjectiles.direction
    AsteroidProjectiles.velocity
    AsteroidProjectiles.density
    AsteroidProjectiles.vertical_velocity

Usage example
-------------

.. code-block:: python

    from cratermaker import Projectile
    asteroids = Projectile.maker("asteroids")


.. currentmodule:: cratermaker.components.projectile.comets

Comet projectiles
=================

.. autosummary::
    :toctree: generated/

    CometProjectiles.maker
    CometProjectiles.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    CometProjectiles.target_name
    CometProjectiles.target
    CometProjectiles.sample_angles
    CometProjectiles.sample_velocities
    CometProjectiles.sample_directions
    CometProjectiles.mean_velocity
    CometProjectiles.angle
    CometProjectiles.direction
    CometProjectiles.velocity
    CometProjectiles.density
    CometProjectiles.vertical_velocity

Usage example
-------------

.. code-block:: python

    from cratermaker import Projectile
    comets = Projectile.maker("comets")


.. _api-Utility:

.. currentmodule:: cratermaker.utils

Utility functions
=================

Monte Carlo
-----------

.. autosummary::
   :toctree: generated/

    montecarlo.get_random_location
    montecarlo.get_random_impact_angle
    montecarlo.get_random_velocity
    montecarlo.get_random_size
    montecarlo.bounded_norm

General utilities
-----------------

.. autosummary::
    :toctree: generated/

    general_utils.Parameter
    general_utils.validate_and_normalize_location
    general_utils.normalize_coords
    general_utils.R_to_CSFD

Export functions
----------------

.. autosummary::
    :toctree: generated/

    export.to_vtk
    export.make_circle_file

Custom type definitions
-----------------------

.. autosummary::
    :toctree: generated/

    custom_types

.. _api-ComponentAPI:


Component API utilities
-----------------------

.. autosummary::
    :toctree: generated/

    component_utils.ComponentBase
    component_utils.ComponentBase.maker
    component_utils.ComponentBase.name
    component_utils.ComponentBase.register
    component_utils.ComponentBase.available
    component_utils.import_components
