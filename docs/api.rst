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
`surf`, which contains a UxDataset object called `surf.uxds`. This is an unstructured grid dataset that contains data for the target body surface.

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
    Simulation.export_vtk
    Simulation.set_elevation
    Simulation.to_config
    Simulation.make_circle_file
    Simulation.get_smallest_diameter
    Simulation.get_largest_diameter


Attributes
----------

.. autosummary::
    :toctree: generated/

    Simulation.surf
    Simulation.production
    Simulation.scaling
    Simulation.crater
    Simulation.impactor
    Simulation.target
    Simulation.morphology
    Simulation.crater
    Simulation.grid
    Simulation.data_dir
    Simulation.grid_file
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


.. _api-Surface:

Surface
=======

The Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

Creating a surface
------------------

.. autosummary::
    :toctree: generated/

    Surface.make

Methods
-------

.. autosummary::
    :toctree: generated/

    Surface.calculate_haversine_distance
    Surface.calculate_initial_bearing
    Surface.find_nearest_index
    Surface.get_reference_surface
    Surface.get_distance
    Surface.get_initial_bearing
    Surface.generate_data
    Surface.set_elevation
    Surface.elevation_to_cartesian
    Surface.extract_region
    Surface.get_random_location_on_face
    Surface.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    Surface.uxds
    Surface.grid
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

.. _api-Target:

Target
======

The Target class represents the target body in a crater simulation. It encapsulates properties of the target, such as its material composition, size, and other physical characteristics.

Creating a Target
-----------------

.. autosummary::
    :toctree: generated/

    Target.make

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

.. _api-Crater:

Crater
======

The ``Crater`` class represents a single crater in the simulation. It is used to model the crater resulting from an impact, including its size, shape, depth, and other morphological features. It also defines the properties of the projectile, such as its size, velocity, material, and angle of impact.


Methods
-------

.. autosummary::
    :toctree: generated/

    Crater.make

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

.. _api-Grid:

.. currentmodule:: cratermaker.components.grid

Grid
====

.. autosummary::
    :toctree: generated/

    Grid.make



Icosphere grid
==============

.. autosummary::
    :toctree: generated/

    icosphere.IcosphereGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    icosphere.IcosphereGrid.available
    icosphere.IcosphereGrid.generate_face_distribution
    icosphere.IcosphereGrid.generate_grid
    icosphere.IcosphereGrid.check_if_regrid
    icosphere.IcosphereGrid.create_grid
    icosphere.IcosphereGrid.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    icosphere.IcosphereGrid.name
    icosphere.IcosphereGrid.radius
    icosphere.IcosphereGrid.uxgrid
    icosphere.IcosphereGrid.file
    icosphere.IcosphereGrid.regrid
    icosphere.IcosphereGrid.gridlevel

.. currentmodule:: cratermaker.components.grid

Arbitrary resolution grid
=========================

.. autosummary::
    :toctree: generated/

    arbitrary_resolution.ArbitraryResolutionGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    arbitrary_resolution.ArbitraryResolutionGrid.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    arbitrary_resolution.ArbitraryResolutionGrid.pix
    arbitrary_resolution.ArbitraryResolutionGrid.radius


.. currentmodule:: cratermaker.components.grid

Hi-res local grid
=================

.. autosummary::
    :toctree: generated/

    hireslocal.HiResLocalGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    hireslocal.HiResLocalGrid.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    hireslocal.HiResLocalGrid.pix
    hireslocal.HiResLocalGrid.radius
    hireslocal.HiResLocalGrid.local_radius
    hireslocal.HiResLocalGrid.local_location
    hireslocal.HiResLocalGrid.superdomain_scale_factor

.. _api-Production:

.. currentmodule:: cratermaker.components.production

Production
==========

.. autosummary::
    :toctree: generated/

    Production.make


.. _api-NeukumProduction:

.. currentmodule:: cratermaker.components.production

NeukumProduction
================

The NeukumProduction class extends the Production class to implement the Neukum production function for the Moon and Mars.

.. autosummary::
    :toctree: generated/

    neukum.NeukumProduction.make

Methods
-------

.. autosummary::
    :toctree: generated/

    neukum.NeukumProduction.sample
    neukum.NeukumProduction.function
    neukum.NeukumProduction.chronology
    neukum.NeukumProduction.size_frequency_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    neukum.NeukumProduction.sfd_coef
    neukum.NeukumProduction.sfd_range
    neukum.NeukumProduction.valid_time
    neukum.NeukumProduction.tau
    neukum.NeukumProduction.Cexp
    neukum.NeukumProduction.Clin
    neukum.NeukumProduction.generator_type



.. _api-PowerLawProduction:

.. currentmodule:: cratermaker.components.production

Power law production function
=============================

The PowerLawProduction class extends the Production class to implement a power law production function for the Moon and Mars.

.. autosummary::
    :toctree: generated/

    powerlaw.PowerLawProduction.make

Methods
-------

.. autosummary::
    :toctree: generated/

    powerlaw.PowerLawProduction.sample
    powerlaw.PowerLawProduction.function
    powerlaw.PowerLawProduction.chronology


Attributes
----------

.. autosummary::
    :toctree: generated/

    powerlaw.PowerLawProduction.N1_coef
    powerlaw.PowerLawProduction.slope

.. _api-Scaling:

.. currentmodule:: cratermaker.components.scaling

Scaling
=======

The Scaling class is an operations class for computing the scaling relationships between impactors and craters. It encapsulates the logic for converting between projectile properties and crater properties, as well as determining crater morphology based on size and target properties.

.. autosummary::
    :toctree: generated/

    Scaling.make


Attributes
----------

.. autosummary::
    :toctree: generated/

    Scaling.target
    Scaling.impactor
    Scaling.target_density


.. currentmodule:: cratermaker.components.scaling

Richardson 2009 scaling
========================

.. autosummary::
    :toctree: generated/

    richardson2009.Richardson2009.make

Methods
-------

    richardson2009.Richardson2009.final_to_transient
    richardson2009.Richardson2009.transient_to_final
    richardson2009.Richardson2009.projectile_to_transient
    richardson2009.Richardson2009.transient_to_projectile
    richardson2009.Richardson2009.get_morphology_type

Attributes
----------

    richardson2009.Richardson2009.target
    richardson2009.Richardson2009.impactor
    richardson2009.Richardson2009.target_density
    richardson2009.Richardson2009.material_catalogue
    richardson2009.Richardson2009.material_name
    richardson2009.Richardson2009.K1
    richardson2009.Richardson2009.mu
    richardson2009.Richardson2009.Ybar
    richardson2009.Richardson2009.catalogue_key
    richardson2009.Richardson2009.transition_diameter
    richardson2009.Richardson2009.transition_nominal
    richardson2009.Richardson2009.complex_enlargment_factor
    richardson2009.Richardson2009.simple_enlargment_factor
    richardson2009.Richardson2009.final_exp


.. currentmodule:: cratermaker.components.morphology

.. _api-Morphology:

Morphology
==========

The Morphology class is an operations class for computing the morphology of a crater based on its size and target properties. It encapsulates the logic for altering the topography of the surface based on the crater properties.


.. autosummary::
    :toctree: generated/

    Morphology.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    Morphology.crater


.. currentmodule:: cratermaker.components.morphology

SimpleMoon
==========

.. autosummary::
    :toctree: generated/

    simplemoon.SimpleMoon.make



Methods
-------

.. autosummary::
    :toctree: generated/

    simplemoon.SimpleMoon.crater_profile
    simplemoon.SimpleMoon.ejecta_profile
    simplemoon.SimpleMoon.ejecta_distribution
    simplemoon.SimpleMoon.form_crater
    simplemoon.SimpleMoon.form_ejecta
    simplemoon.SimpleMoon.ray_intensity
    simplemoon.SimpleMoon.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    simplemoon.SimpleMoon.crater
    simplemoon.SimpleMoon.dorays
    simplemoon.SimpleMoon.ejrim
    simplemoon.SimpleMoon.ejecta_truncation
    simplemoon.SimpleMoon.floordepth
    simplemoon.SimpleMoon.floor_diameter
    simplemoon.SimpleMoon.rimheight
    simplemoon.SimpleMoon.rimwidth
    simplemoon.SimpleMoon.peakheight
    simplemoon.SimpleMoon.rng

.. currentmodule:: cratermaker.components.impactor

.. _api-Impactor:

Impactor
==========

The Impactor class is an operations class defining the interface for generating impactor velocities, angles, and densities for a given target body.

.. autosummary::
    :toctree: generated/

    Impactor.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    Impactor.target_name
    Impactor.target
    Impactor.sample_angles
    Impactor.sample_velocities
    Impactor.sample_directions
    Impactor.mean_velocity
    Impactor.angle
    Impactor.direction
    Impactor.velocity
    Impactor.density
    Impactor.vertical_velocity


.. currentmodule:: cratermaker.components.impactor

Asteroid impactors
================== 

.. autosummary::
    :toctree: generated/

    asteroids.AsteroidImpactors.make


Methods
-------

.. autosummary::
    :toctree: generated/

    asteroids.AsteroidImpactors.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    asteroids.AsteroidImpactors.target_name
    asteroids.AsteroidImpactors.target
    asteroids.AsteroidImpactors.sample_angles
    asteroids.AsteroidImpactors.sample_velocities
    asteroids.AsteroidImpactors.sample_directions
    asteroids.AsteroidImpactors.mean_velocity
    asteroids.AsteroidImpactors.angle
    asteroids.AsteroidImpactors.direction
    asteroids.AsteroidImpactors.velocity
    asteroids.AsteroidImpactors.density
    asteroids.AsteroidImpactors.vertical_velocity


.. currentmodule:: cratermaker.components.impactor

Comet impactors
===============

.. autosummary::
    :toctree: generated/

    comets.CometImpactors.make


Methods
-------

.. autosummary::
    :toctree: generated/

    comets.CometImpactors.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    comets.CometImpactors.target_name
    comets.CometImpactors.target
    comets.CometImpactors.sample_angles
    comets.CometImpactors.sample_velocities
    comets.CometImpactors.sample_directions
    comets.CometImpactors.mean_velocity
    comets.CometImpactors.angle
    comets.CometImpactors.direction
    comets.CometImpactors.velocity
    comets.CometImpactors.density
    comets.CometImpactors.vertical_velocity


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
    general_utils.validate_and_convert_location
    general_utils.normalize_coords
    general_utils.R_to_CSFD


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
    component_utils.ComponentBase.make
    component_utils.ComponentBase.name
    component_utils.ComponentBase.register
    component_utils.ComponentBase.available
    component_utils.import_components

