.. _api-projectile:

##########
Projectile
##########

The Projectile class is an operations class defining the interface for generating projectile velocities, angles, and densities for a given target body.

*************************
Available Implementations
*************************

+------------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------------------------------------------+
| Class                                                                        | Argument name     | Example Usage                                                                                             |
+==============================================================================+===================+===========================================================================================================+
| :py:class:`~cratermaker.components.projectile.asteroids.AsteroidProjectiles` | "asteroids"       | asteroids = Projectile.maker("asteroids")                                                                 |
+------------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------------------------------------------+
| :py:class:`~cratermaker.components.projectile.comets.CometProjectiles`       | "comets"          | comets = Projectile.maker("comets")                                                                       |
+------------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------------------------------------------+
| :py:class:`~cratermaker.components.projectile.generic.GenericProjectiles`    | "generic"         | generic = Projectile.maker("generic",sample=True,projectile_mean_velocity=20e3,projectile_density=1500.0) |
+------------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------------------------------------------+

.. autoclass:: cratermaker.components.projectile.Projectile
   :members:
   :undoc-members:
   :no-index-entry:

.. currentmodule:: cratermaker.components.projectile.asteroids

*******************
AsteroidProjectiles
*******************

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.asteroids.AsteroidProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(projectile="asteroids")
    print(sim.projectile)


.. currentmodule:: cratermaker.components.projectile.comets

****************
CometProjectiles
****************

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.comets.CometProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(projectile="comets")
    print(sim.projectile)



.. currentmodule:: cratermaker.components.projectile.generic

******************
GenericProjectiles
******************

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.generic.GenericProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(projectile="generic", sample=True, projectile_mean_velocity=20e3, projectile_density=1500.0)
    print(sim.projectile)

