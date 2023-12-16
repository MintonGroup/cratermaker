.. currentmodule:: cratermaker

.. _api:

#########################
Cratermaker API reference
#########################

This section of the documentation provides a detailed reference for the Production classes in the Cratermaker project.

Simulation
==========

The Simulation class is the main class for the cratermaker project. It is used to create a simulation of a crater on a given target 
body. The Simulation class is used to generate craters of a given size and morphology based on the production function, morphology
function, and crater scaling relationship model. The surface of the target body is represented by a Surface attribute called
`surf`, which is derived from a UxDataset object. This is an unstructured grid dataset that contains data for the target body surface.

Creating a Simulation
---------------------

.. autosummary::
    :toctree: generated/

    Simulation

Attributes
----------

.. autosummary::
    :toctree: generated/

    Simulation.crater
    Simulation.data_dir
    Simulation.elevation_file
    Simulation.grid_file
    Simulation.morphology_cls
    Simulation.n_face
    Simulation.n_node
    Simulation.pix
    Simulation.production
    Simulation.projectile
    Simulation.rng
    Simulation.scale_cls
    Simulation.seed
    Simulation.simdir
    Simulation.surf
    Simulation.target

Methods
-------

.. autosummary::
    :toctree: generated/

    Simulation.apply_noise
    Simulation.emplace_crater
    Simulation.export_vtk
    Simulation.generate_crater
    Simulation.generate_projectile
    Simulation.initialize_surface
    Simulation.populate
    Simulation.save
    Simulation.set_elevation
    Simulation.set_properties
    Simulation.to_json

Surface
=======

The Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

Creating a surface
------------------

.. autosummary::
    :toctree: generated/

    Surface
    initialize_surface

Attributes
----------

.. autosummary::
    :toctree: generated/

    Surface.data_dir
    Surface.elevation_file
    Surface.grid_file
    Surface.grid_temp_dir
    Surface.pix
    Surface.target_radius

Methods
-------

.. autosummary::
    :toctree: generated/

    Surface.calculate_haversine_distance
    Surface.calculate_initial_bearing
    Surface.find_nearest_face_index
    Surface.find_nearest_node_index
    Surface.get_average_surface
    Surface.get_face_distance
    Surface.get_face_initial_bearing
    Surface.get_node_distance
    Surface.get_node_initial_bearing
    Surface.set_elevation

Production
==========

The Production class serves as a base class for computing the production function for craters and impactors.

Creating a Production Function
------------------------------

.. autosummary::
    :toctree: generated/

    Production

Attributes
----------

.. autosummary::
    :toctree: generated/

    Production.generator_type
    Production.impact_velocity_model
    Production.model
    Production.N1_coef
    Production.slope
    Production.rng
    Production.valid_models
    Production.valid_generator_types
    Production.valid_time

Methods
-------

.. autosummary::
    :toctree: generated/

    Production.set_model_parameters
    Production.set_mean_impact_velocity
    Production.function
    Production.function_inverse
    Production.sample

NeukumProduction
================

The NeukumProduction class extends the Production class to implement the Neukum production function for the Moon and Mars.

Creating a Neukum Production Function
-------------------------------------

.. autosummary::
    :toctree: generated/

    NeukumProduction

Attributes
----------

.. autosummary::
    :toctree: generated/

    NeukumProduction.sfd_coef
    NeukumProduction.sfd_range
    NeukumProduction.tau
    NeukumProduction.Cexp
    NeukumProduction.Clin

Methods
-------

.. autosummary::
    :toctree: generated/

    NeukumProduction.set_model_parameters
    NeukumProduction.function

.. currentmodule:: cratermaker

Target
======

The Target class represents the target body in a crater simulation. It encapsulates properties of the target, such as its material composition, size, and other physical characteristics.

Creating a Target
-----------------

.. autosummary::
    :toctree: generated/

    Target

Attributes
----------

.. autosummary::
    :toctree: generated/

    Target.catalogue
    Target.diameter
    Target.escape_velocity
    Target.gravity
    Target.material
    Target.material_name
    Target.mean_impact_velocity
    Target.name
    Target.radius
    Target.transition_scale_type

Methods
-------

.. autosummary::
    :toctree: generated/

    Target.set_properties


Material
========

The Material class represents a material that can be used in the simulation. The properties defined in a Material object are used in crater scaling calculations.

Creating a matererial
---------------------

.. autosummary::
    :toctree: generated/

    Material

Attributes
----------

.. autosummary::
    :toctree: generated/

    Material.name
    Material.catalogue
    Material.density
    Material.K1
    Material.mu
    Material.Ybar

Crater
======

The Crater class represents a single crater in the simulation. It is used to model the crater resulting from an impact, including its size, shape, depth, and other morphological features.

Creating a Crater
-----------------

.. autosummary::
    :toctree: generated/

    Crater


Attributes
----------

.. autosummary::
    :toctree: generated/

    Crater.diameter
    Crater.location
    Crater.morphology
    Crater.morphology_cls
    Crater.radius
    Crater.rng
    Crater.scale
    Crater.scale_cls
    Crater.target
    Crater.transient_diameter
    Crater.transient_radius



Projectile
==========

The Projectile class represents a single projectile in the simulation. It defines the properties of the impacting object, such as its size, velocity, material, and angle of impact.

Creating a Projectile
---------------------

.. autosummary::
    :toctree: generated/

    Projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    Projectile.angle
    Projectile.density
    Projectile.diameter
    Projectile.location
    Projectile.mass
    Projectile.radius
    Projectile.rng
    Projectile.scale
    Projectile.scale_cls
    Projectile.target
    Projectile.velocity
    Projectile.vertical_velocity


Scale
=====

Scale is an operations class for computing the scaling relationships between impactors and craters.  This class encapsulates the logic for converting between projectile properties and crater properties, as well as determining crater morphology based on size and target properties.


Morphology
==========

The Morphology class represents the dimensions of different morphological characteristics of a crater.


Utility functions
=================

Monte Carlo
-----------

.. autosummary::
   :toctree: generated/

    cratermaker.utils.montecarlo.get_random_location
    cratermaker.utils.montecarlo.get_random_impact_angle
    cratermaker.utils.montecarlo.get_random_velocity
    cratermaker.utils.montecarlo.get_random_size
    cratermaker.utils.montecarlo.bounded_norm

General utilities
-----------------

.. autosummary::
    :toctree: generated/

    cratermaker.utils.general_utils.to_config
    cratermaker.utils.general_utils.set_properties
    cratermaker.utils.general_utils.check_properties
    cratermaker.utils.general_utils.create_catalogue
    cratermaker.utils.general_utils.validate_and_convert_location
    cratermaker.utils.general_utils.normalize_coords
    cratermaker.utils.general_utils.R_to_CSFD

Perlin noise
------------

.. autosummary::
    :toctree: generated/

    cratermaker.perlin.apply_noise


Custom type definitions
-----------------------

.. autosummary::
    :toctree: generated/

    cratermaker.utils.custom_types

