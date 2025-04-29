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

.. currentmodule:: cratermaker.components.grid.icosphere

Grid
====

.. autosummary::
    :toctree: generated/

    IcosphereGrid.make
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
    IcosphereGrid.radius
    IcosphereGrid.uxgrid
    IcosphereGrid.file
    IcosphereGrid.regrid
    IcosphereGrid.gridlevel

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import IcosphereGrid
    grid = IcosphereGrid.make(3, radius=1737.4)


.. currentmodule:: cratermaker.components.grid.arbitrary_resolution

Arbitrary resolution grid
=========================

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionGrid.make
    ArbitraryResolutionGrid.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionGrid.pix
    ArbitraryResolutionGrid.radius

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import ArbitraryResolutionGrid
    grid = ArbitraryResolutionGrid.make(pix=100, radius=1737.4)


.. currentmodule:: cratermaker.components.grid.hireslocal

Hi-res local grid
=================

.. autosummary::
    :toctree: generated/

    HiResLocalGrid.make
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

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import HiResLocalGrid
    grid = HiResLocalGrid.make(pix=50, radius=1737.4, local_radius=100)


.. _api-Production:

.. currentmodule:: cratermaker.components.production

Production
==========

.. autosummary::
    :toctree: generated/

    Production.make

Usage example
-------------

.. code-block:: python

    from cratermaker import Production
    production = Production.make("neukum")


.. _api-NeukumProduction:

.. currentmodule:: cratermaker.components.production.neukum

NeukumProduction
================

.. autosummary::
    :toctree: generated/

    NeukumProduction.make
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

    from cratermaker.components.production.neukum import NeukumProduction
    neukum = NeukumProduction.make()


.. _api-PowerLawProduction:

.. currentmodule:: cratermaker.components.production.powerlaw

Power law production function
=============================

.. autosummary::
    :toctree: generated/

    PowerLawProduction.make
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

    from cratermaker.components.production.powerlaw import PowerLawProduction
    powerlaw = PowerLawProduction.make()


.. _api-Scaling:

.. currentmodule:: cratermaker.components.scaling

Scaling
=======

The `Scaling` class is an abstract base class for crater scaling relationships. 
Use the `Scaling.make` method to create a specific scaling model. 
Available models: see `Scaling.available()`.

Example:

.. code-block:: python

    from cratermaker import Scaling
    scaling_model = Scaling.make("richardson2009")

Main Method
-----------

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

Usage example
-------------

.. code-block:: python

    from cratermaker.components.scaling import Scaling
    scaling_model = Scaling.make("richardson2009")


.. _api-Richardson2009:

.. currentmodule:: cratermaker.components.scaling.richardson2009

Richardson 2009 Scaling Model
=============================

The `Richardson2009` model implements crater scaling based on Richardson (2009). 

To create this model:

.. code-block:: python

    from cratermaker.components.scaling import Scaling
    scaling_model = Scaling.make("richardson2009")

Available Methods
------------------

.. autosummary::
    :toctree: generated/

    Richardson2009.final_to_transient
    Richardson2009.transient_to_final
    Richardson2009.projectile_to_transient
    Richardson2009.transient_to_projectile
    Richardson2009.get_morphology_type

Attributes
----------

.. autosummary::
    :toctree: generated/

    Richardson2009.target
    Richardson2009.impactor
    Richardson2009.target_density
    Richardson2009.material_catalogue
    Richardson2009.material_name
    Richardson2009.K1
    Richardson2009.mu
    Richardson2009.Ybar
    Richardson2009.catalogue_key
    Richardson2009.transition_diameter
    Richardson2009.transition_nominal
    Richardson2009.complex_enlargement_factor
    Richardson2009.simple_enlargement_factor
    Richardson2009.final_exp

Usage example
-------------

.. code-block:: python

    from cratermaker.components.scaling.richardson2009 import Richardson2009
    scaling = Richardson2009()


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

Usage example
-------------

.. code-block:: python

    from cratermaker.components.morphology import Morphology
    morphology = Morphology.make()


.. currentmodule:: cratermaker.components.morphology.simplemoon

SimpleMoon
==========

.. autosummary::
    :toctree: generated/

    SimpleMoon.make
    SimpleMoon.crater_profile
    SimpleMoon.ejecta_profile
    SimpleMoon.ejecta_distribution
    SimpleMoon.form_crater
    SimpleMoon.form_ejecta
    SimpleMoon.ray_intensity
    SimpleMoon.make

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

    from cratermaker.components.morphology.simplemoon import SimpleMoon
    simplemoon = SimpleMoon.make()


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

Usage example
-------------

.. code-block:: python

    from cratermaker.components.impactor import Impactor
    impactor = Impactor.make()


.. currentmodule:: cratermaker.components.impactor.asteroids

Asteroid impactors
================== 

.. autosummary::
    :toctree: generated/

    AsteroidImpactors.make
    AsteroidImpactors.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    AsteroidImpactors.target_name
    AsteroidImpactors.target
    AsteroidImpactors.sample_angles
    AsteroidImpactors.sample_velocities
    AsteroidImpactors.sample_directions
    AsteroidImpactors.mean_velocity
    AsteroidImpactors.angle
    AsteroidImpactors.direction
    AsteroidImpactors.velocity
    AsteroidImpactors.density
    AsteroidImpactors.vertical_velocity

Usage example
-------------

.. code-block:: python

    from cratermaker.components.impactor.asteroids import AsteroidImpactors
    asteroids = AsteroidImpactors.make()


.. currentmodule:: cratermaker.components.impactor.comets

Comet impactors
===============

.. autosummary::
    :toctree: generated/

    CometImpactors.make
    CometImpactors.new_projectile

Attributes
----------

.. autosummary::
    :toctree: generated/

    CometImpactors.target_name
    CometImpactors.target
    CometImpactors.sample_angles
    CometImpactors.sample_velocities
    CometImpactors.sample_directions
    CometImpactors.mean_velocity
    CometImpactors.angle
    CometImpactors.direction
    CometImpactors.velocity
    CometImpactors.density
    CometImpactors.vertical_velocity

Usage example
-------------

.. code-block:: python

    from cratermaker.components.impactor.comets import CometImpactors
    comets = CometImpactors.make()


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
