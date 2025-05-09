.. _api-scaling:

Scaling
=======

The `Scaling` class is an abstract base class for crater scaling relationships. 
Use the `Scaling.maker` method to create a specific scaling model.

Available Scaling Implementations
---------------------------------

+-------------------+----------------+--------------------------------------------------------+
| Class             | Instantiation  | Example Usage                                          |
+===================+================+========================================================+
| MonteCarloScaling | "montecarlo"   | scaling_model = Scaling.maker("montecarlo",            |
|                   |                |                                target="Mars",          |
|                   |                |                                projectile="asteroids") |
+-------------------+----------------+--------------------------------------------------------+
| CTEMScaling       | "ctem"         | scaling_model = Scaling.maker("ctem",                  |
|                   |                |                               target="Mars",           |
|                   |                |                               projectile="asteroids")  |
+-------------------+----------------+--------------------------------------------------------+

.. autoclass:: cratermaker.components.scaling.Scaling
   :members:
   :undoc-members:
   :no-index:


.. _api-MonteCarloScaling:

.. currentmodule:: cratermaker.components.scaling.montecarlo

Monte Carlo scaling model
=========================

See `Scaling`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.scaling.montecarlo.MonteCarloScaling
   :members:
   :undoc-members:
   :no-index:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Scaling
    scaling_model = Scaling.maker("montecarlo", target="Mars", projectile="asteroids")

.. _api-CTEMScaling:

.. currentmodule:: cratermaker.components.scaling.ctem


CTEM Scaling model
==================

See `Scaling`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.scaling.ctem.CTEMScaling
   :members:
   :undoc-members:
   :no-index:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Scaling
    scaling_model = Scaling.maker("ctem", target="Mars", projectile="asteroids")

