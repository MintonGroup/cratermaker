.. _api-projectile:

##########
Projectile
##########

The Projectile class is an operations class defining the interface for generating projectile velocities, angles, and densities for a given target body.

Available Projectile Implementations
------------------------------------

+------------------------------------------------------------------------------+-------------------+-------------------------------------------------------------------------------------+
| Class                                                                        | Instantiation     | Example Usage                                                                       |
+==============================================================================+===================+=====================================================================================+
| :py:class:`~cratermaker.components.projectile.asteroids.AsteroidProjectiles` | "asteroids"       | asteroids = Projectile.maker("asteroids")                                           |
+------------------------------------------------------------------------------+-------------------+-------------------------------------------------------------------------------------+
| :py:class:`~cratermaker.components.projectile.comets.CometProjectiles`       | "comets"          | comets = Projectile.maker("comets")                                                 |
+------------------------------------------------------------------------------+-------------------+-------------------------------------------------------------------------------------+
| :py:class:`~cratermaker.components.projectile.generic.GenericProjectiles`    | "generic"         | generic = Projectile.maker("generic",sample=True,mean_velocity=20e3,density=1500.0) |
+------------------------------------------------------------------------------+-------------------+-------------------------------------------------------------------------------------+

.. autoclass:: cratermaker.components.projectile.Projectile
   :members:
   :undoc-members:
   :no-index-entry:

.. currentmodule:: cratermaker.components.projectile.asteroids

Asteroid projectiles
==================== 

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.asteroids.AsteroidProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Projectile
    asteroids = Projectile.maker("asteroids")


.. currentmodule:: cratermaker.components.projectile.comets

Comet projectiles
=================

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.comets.CometProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Projectile
    comets = Projectile.maker("comets")



.. currentmodule:: cratermaker.components.projectile.generic

Generic projectiles
===================

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.generic.GenericProjectiles
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Projectile
    generic = Projectile.maker("generic", sample=True, mean_velocity=20e3, density=1500.0)

