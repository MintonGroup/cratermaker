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

    Simulation.make

Methods
-------

.. autosummary::
    :toctree: generated/

    Simulation.run
    Simulation.populate
    Simulation.emplace_crater
    Simulation.generate_crater
    Simulation.generate_projectile
    Simulation.apply_noise
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
    Simulation.ejecta_truncation
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


Components
==========

.. currentmodule:: cratermaker.components.grid

.. _api-Grid:

Generating grids
----------------

.. autosummary::
    :toctree: generated/


Generating a uniform Icosphere grid
-----------------------------------

.. autosummary::
    :toctree: generated/

    IcosphereGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    IcosphereGrid.available
    IcosphereGrid.generate_face_distribution
    IcosphereGrid.generate_grid
    IcosphereGrid.check_if_regrid
    IcosphereGrid.create_grid
    IcosphereGrid.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    IcosphereGrid.name
    IcosphereGrid.pix
    IcosphereGrid.radius
    IcosphereGrid.uxgrid
    IcosphereGrid.file
    IcosphereGrid.regrid
    IcosphereGrid.gridlevel

Generating a uniform arbitrary resolution grid
----------------------------------------------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionGrid.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionGrid.pix
    ArbitraryResolutionGrid.radius


Generating a non-uniform grid with a high resolution local region
-----------------------------------------------------------------

.. autosummary::
    :toctree: generated/

    HiResLocalGrid.make

Methods
-------

.. autosummary::
    :toctree: generated/

    HiResLocalGrid.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    HiResLocalGrid.pix
    HiResLocalGrid.radius
    HiResLocalGrid.local_radius
    HiResLocalGrid.local_location
    HiResLocalGrid.superdomain_scale_factor

.. _api-Production:


The Production class serves as a base class for computing the production function for craters and impactors.

Creating a Production Function
------------------------------

.. autosummary::
    :toctree: generated/

    Production.make


.. _api-NeukumProduction:

.. currentmodule:: cratermaker.components.production

NeukumProduction
================

The NeukumProduction class extends the Production class to implement the Neukum production function for the Moon and Mars.

Creating a Neukum Production Function
-------------------------------------

.. autosummary::
    :toctree: generated/

    NeukumProduction.make

Methods
-------

.. autosummary::
    :toctree: generated/

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

.. currentmodule:: cratermaker.components.production


.. _api-Scaling:

Scaling
=======

The Scaling class is an operations class for computing the scaling relationships between impactors and craters. It encapsulates the logic for converting between projectile properties and crater properties, as well as determining crater morphology based on size and target properties.

Creating Scaling
--------------

.. autosummary::
    :toctree: generated/

    Scaling.make

Methods
-------

.. autosummary::
    :toctree: generated/

    Scaling.get_morphology_type
    Scaling.final_to_transient
    Scaling.transient_to_final
    Scaling.projectile_to_crater
    Scaling.crater_to_projectile
    Scaling.projectile_to_transient
    Scaling.transient_to_projectile
    Scaling.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    Scaling.target
    Scaling.material
    Scaling.material_name
    Scaling.rng
    Scaling.transition_diameter
    Scaling.transition_nominal
    Scaling.simple_enlargement_factor
    Scaling.complex_enlargement_factor
    Scaling.final_exp

.. currentmodule:: cratermaker.components.morphology

.. _api-Morphology:

Morphology
==========

The Morphology class is an operations class for computing the morphology of a crater based on its size and target properties. It encapsulates the logic for altering the topography of the surface based on the crater properties.

Creating Morphology
-------------------

.. autosummary::
    :toctree: generated/

    Morphology

Methods
-------

.. autosummary::
    :toctree: generated/

    Morphology.compute_rmax
    Morphology.crater_profile
    Morphology.ejecta_profile
    Morphology.ejecta_distribution
    Morphology.form_crater
    Morphology.form_ejecta
    Morphology.form_secondaries
    Morphology.ray_intensity
    Morphology.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    Morphology.crater
    Morphology.final_diameter
    Morphology.dorays
    Morphology.ejrim
    Morphology.ejecta_truncation
    Morphology.floordepth
    Morphology.floor_diameter
    Morphology.morphology_type
    Morphology.radius
    Morphology.rimheight
    Morphology.rimwidth
    Morphology.peakheight
    Morphology.target
    Morphology.rng

.. currentmodule:: cratermaker.components.impactor

.. _api-Morphology:

Impactor
==========

The Impactor class is an operations class defining the interface for generating impactor velocities, angles, and densities for a given target body.

Creating asteroids
-------------------

.. autosummary::
    :toctree: generated/

    Impactor

Methods
-------

.. autosummary::
    :toctree: generated/

    Impactor.compute_rmax
    Impactor.crater_profile
    Impactor.ejecta_profile
    Impactor.ejecta_distribution
    Impactor.form_crater
    Impactor.form_ejecta
    Impactor.form_secondaries
    Impactor.ray_intensity
    Impactor.make

Attributes
----------

.. autosummary::
    :toctree: generated/

    Impactor.crater
    Impactor.final_diameter
    Impactor.dorays
    Impactor.ejrim
    Impactor.ejecta_truncation
    Impactor.floordepth
    Impactor.floor_diameter
    Impactor.morphology_type
    Impactor.radius
    Impactor.rimheight
    Impactor.rimwidth
    Impactor.peakheight
    Impactor.target
    Impactor.rng


.. _api-Utility:

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

    cratermaker.utils.general_utils.Parameter
    cratermaker.utils.general_utils.parameter
    cratermaker.utils.general_utils.validate_and_convert_location
    cratermaker.utils.general_utils.normalize_coords
    cratermaker.utils.general_utils.R_to_CSFD


Custom type definitions
-----------------------

.. autosummary::
    :toctree: generated/

    cratermaker.utils.custom_types

.. _api-Fortran:


Component API utilities
-----------------------

.. autosummary::
    :toctree: generated/

    cratermaker.utils.component_utils.ComponentBase
    cratermaker.utils.component_utils.ComponentBase.make
    cratermaker.utils.component_utils.ComponentBase.name
    cratermaker.utils.component_utils.ComponentBase.register
    cratermaker.utils.component_utils.ComponentBase.available
    cratermaker.utils.component_utils.import_components

