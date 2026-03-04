.. _api-scaling:

#######
Scaling
#######

The `Scaling` class is an abstract base class for crater scaling relationships. 
Use the `Scaling.maker` method to create a specific scaling model.

*************************
Available Implementations
*************************

+--------------------------------------------------------------------------+----------------+------------------------------------------------------------------------------------+
| Class                                                                    | Argument name  | Example Usage                                                                      |
+==========================================================================+================+====================================================================================+
| :py:class:`~cratermaker.components.scaling.montecarlo.MonteCarloScaling` | "montecarlo"   | scaling_model = Scaling.maker("montecarlo", target="Mars", projectile="asteroids") |
+--------------------------------------------------------------------------+----------------+------------------------------------------------------------------------------------+
| :py:class:`~cratermaker.components.scaling.ctem.CTEMScaling`             | "ctem"         | scaling_model = Scaling.maker("ctem", target="Mars", projectile="asteroids")       |
+--------------------------------------------------------------------------+----------------+------------------------------------------------------------------------------------+

.. autoclass:: cratermaker.components.scaling.Scaling
   :members:
   :undoc-members:
   :no-index-entry:


.. _api-MonteCarloScaling:

.. currentmodule:: cratermaker.components.scaling.montecarlo

*****************
MonteCarloScaling
*****************

See `Scaling`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.scaling.montecarlo.MonteCarloScaling
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python
   :linenos:

    from cratermaker import Scaling
    scaling_model = Scaling.maker("montecarlo", target="Mars", projectile="asteroids")

.. _api-CTEMScaling:

.. currentmodule:: cratermaker.components.scaling.ctem

***********
CTEMScaling
***********

See `MonteCarloScaling`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.scaling.ctem.CTEMScaling
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python
   :linenos:

    from cratermaker import Scaling
    scaling_model = Scaling.maker("ctem", target="Mars", projectile="asteroids")

