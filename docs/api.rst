.. currentmodule:: cratermaker

.. _api:

#########################
Cratermaker API reference
#########################

This section of the documentation provides a detailed reference for the Production classes in the Cratermaker project.


Simulation
==========

The Simulation class orchestrates the processes involved in running a crater simulation.


.. autoclass:: cratermaker.core.simulation.Simulation
   :members:
   :undoc-members:
   :show-inheritance:

Surface
=======

This Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

.. autoclass:: cratermaker.core.surface.Surface
   :members:
   :undoc-members:
   :show-inheritance:

Production
==========

The Production class serves as a base class for computing the production function for craters and impactors.


.. autoclass:: cratermaker.core.production.Production
   :members:
   :undoc-members:
   :show-inheritance:

NeukumProduction
================

The NeukumProduction class extends the Production class to implement the Neukum production function for the Moon and Mars.

.. autoclass:: cratermaker.core.production.NeukumProduction
   :members:
   :undoc-members:
   :show-inheritance:

Target
======

The Target class represents the target body for the simulation and its physical properties, such as radius, gravity, and surface composition.


.. autoclass:: cratermaker.core.target.Target
   :members:
   :undoc-members:
   :show-inheritance:


Material
========

The Material class represents a material that can be used in the simulation. The properties defined in a Material object are used in crater scaling calculations.


.. autoclass:: cratermaker.core.target.Material
   :members:
   :undoc-members:
   :show-inheritance:


Crater
======

The Crater class represents a single crater in the simulation. 

.. autoclass:: cratermaker.core.crater.Crater
   :members:
   :undoc-members:
   :show-inheritance:

Projectile
==========

The Projectile class represents a single projectile in the simulation. It is used to generate a crater of a given size.


.. autoclass:: cratermaker.core.crater.Projectile
   :members:
   :undoc-members:
   :show-inheritance:


Scale
=====

Scale is an operations class for computing the scaling relationships between impactors and craters.  This class encapsulates the logic for converting between projectile properties and crater properties, as well as determining crater morphology based on size and target properties.

.. autoclass:: cratermaker.core.scale.Scale
   :members:
   :undoc-members:
   :show-inheritance:


Morphology
==========

The Morphology class represents the dimensions of different morphological characteristics of a crater.

.. autoclass:: cratermaker.core.morphology.Morphology
   :members:
   :undoc-members:
   :show-inheritance:



Utility function
================

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

