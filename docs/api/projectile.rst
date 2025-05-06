.. _api-projectile:

Projectile
==========

The Projectile class is an operations class defining the interface for generating projectile velocities, angles, and densities for a given target body.

Available Projectile Implementations
------------------------------------

+------------------------+-------------------+----------------------------------------------+
| Class                  | Instantiation     | Example Usage                                |
+========================+===================+==============================================+
| AsteroidProjectiles    | "asteroids"       | asteroids = Projectile.maker("asteroids")    |
+------------------------+-------------------+----------------------------------------------+
| CometProjectiles       | "comets"          | comets = Projectile.maker("comets")          |
+------------------------+-------------------+----------------------------------------------+

.. autoclass:: cratermaker.components.projectile.Projectile
   :members:
   :undoc-members:
   :no-index:

.. currentmodule:: cratermaker.components.projectile.asteroids

Asteroid projectiles
==================== 

See `Projectile`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.projectile.asteroids.AsteroidProjectiles
   :members:
   :undoc-members:
   :no-index:

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
   :no-index:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Projectile
    comets = Projectile.maker("comets")

