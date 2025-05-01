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

    Simulation.surface
    Simulation.production
    Simulation.scaling
    Simulation.crater
    Simulation.projectile
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

.. autosummary::
    :toctree: generated/

    Surface.maker
    Surface.available
    Surface.generate_grid
    Surface.check_if_regrid




.. _api-Surface:

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

Attributes
----------

.. autosummary::
    :toctree: generated/

    IcosphereSurface.name
    IcosphereSurface.radius
    IcosphereSurface.uxgrid
    IcosphereSurface.file
    IcosphereSurface.regrid
    IcosphereSurface.gridlevel

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import IcosphereSurface
    grid = IcosphereSurface.maker(3, radius=1737.4)


.. currentmodule:: cratermaker.components.grid.arbitrary_resolution

Arbitrary resolution grid
=========================

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionSurface.maker
    ArbitraryResolutionSurface.generate_face_distribution
    IcosphereSurface.generate_face_distribution
    IcosphereSurface.generate_grid
    IcosphereSurface.check_if_regrid
    IcosphereSurface.to_config

Attributes
----------

.. autosummary::
    :toctree: generated/

    ArbitraryResolutionSurface.pix
    ArbitraryResolutionSurface.radius

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import ArbitraryResolutionSurface
    grid = ArbitraryResolutionSurface.maker(pix=100, radius=1737.4)


.. currentmodule:: cratermaker.components.grid.hireslocal

Hi-res local grid
=================

.. autosummary::
    :toctree: generated/

    HiResLocalSurface.maker
    HiResLocalSurface.generate_face_distribution

Attributes
----------

.. autosummary::
    :toctree: generated/

    HiResLocalSurface.pix
    HiResLocalSurface.radius
    HiResLocalSurface.local_radius
    HiResLocalSurface.local_location
    HiResLocalSurface.superdomain_scale_factor

Usage example
-------------

.. code-block:: python

    from cratermaker.components.grid import HiResLocalSurface
    grid = HiResLocalSurface.maker(pix=50, radius=1737.4, local_radius=100)


.. _api-Production:

.. currentmodule:: cratermaker.components.production

Production
==========

.. autosummary::
    :toctree: generated/

    Production.maker

Usage example
-------------

.. code-block:: python

    from cratermaker import Production
    production = Production.maker("neukum")


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

    from cratermaker.components.production.neukum import NeukumProduction
    neukum = NeukumProduction.maker()


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

    from cratermaker.components.production.powerlaw import PowerLawProduction
    powerlaw = PowerLawProduction.maker()


.. _api-Scaling:

.. currentmodule:: cratermaker.components.scaling

Scaling
=======

The `Scaling` class is an abstract base class for crater scaling relationships. 
Use the `Scaling.maker` method to create a specific scaling model. 
Available models: see `Scaling.available()`.

Example:

.. code-block:: python

    from cratermaker import Scaling
    scaling_model = Scaling.maker("richardson2009")

Main Method
-----------

.. autosummary::
    :toctree: generated/

    Scaling.maker

Attributes
----------

.. autosummary::
    :toctree: generated/

    Scaling.target
    Scaling.projectile
    Scaling.target_density

Usage example
-------------

.. code-block:: python

    from cratermaker.components.scaling import Scaling
    scaling_model = Scaling.maker("richardson2009")


.. _api-Richardson2009:

.. currentmodule:: cratermaker.components.scaling.richardson2009

Richardson 2009 Scaling Model
=============================

The `Richardson2009` model implements crater scaling based on Richardson (2009). 

To create this model:

.. code-block:: python

    from cratermaker.components.scaling import Scaling
    scaling_model = Scaling.maker("richardson2009")

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
    Richardson2009.projectile
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

    Morphology.maker

Attributes
----------

.. autosummary::
    :toctree: generated/

    Morphology.crater

Usage example
-------------

.. code-block:: python

    from cratermaker.components.morphology import Morphology
    morphology = Morphology.maker()


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

    from cratermaker.components.morphology.simplemoon import SimpleMoon
    simplemoon = SimpleMoon.maker()


.. currentmodule:: cratermaker.components.projectile

.. _api-Projectile:

Projectile
==========

The Projectile class is an operations class defining the interface for generating projectile velocities, angles, and densities for a given target body.

.. autosummary::
    :toctree: generated/

    Projectile.maker

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

Usage example
-------------

.. code-block:: python

    from cratermaker.components.projectile import Projectile
    projectile = Projectile.maker()


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

    from cratermaker.components.projectile.asteroids import AsteroidProjectiles
    asteroids = AsteroidProjectiles.maker()


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

    from cratermaker.components.projectile.comets import CometProjectiles
    comets = CometProjectiles.maker()


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
    component_utils.ComponentBase.maker
    component_utils.ComponentBase.name
    component_utils.ComponentBase.register
    component_utils.ComponentBase.available
    component_utils.import_components
