.. _api-production:

##########
Production
##########

Available Production Implementations
------------------------------------

+----------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------+
| Class                                                                      | Instantiation     | Example Usage                                                         |
+============================================================================+===================+=======================================================================+
| :py:class:`~cratermaker.components.production.neukum.NeukumProduction`     | "neukum"          | production = Production.maker("neukum",version="projectile")          |
+----------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------+
| :py:class:`~cratermaker.components.production.powerlaw.PowerLawProduction` | "powerlaw"        | production = Production.maker("powerlaw", slope=-4.0, N1_coef=1.0e-6) |
+----------------------------------------------------------------------------+-------------------+-----------------------------------------------------------------------+

.. autoclass:: cratermaker.components.production.Production
   :members:
   :undoc-members:
   :no-index-entry:


.. _api-NeukumProduction:

.. currentmodule:: cratermaker.components.production.neukum

NeukumProduction
================

See `Production`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.production.neukum.NeukumProduction
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Production
    production = Production.maker("neukum", version="projectile")


.. _api-PowerLawProduction:

.. currentmodule:: cratermaker.components.production.powerlaw

Power law production function
=============================

See `Production`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.production.powerlaw.PowerLawProduction
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Production
    production = Production.maker("powerlaw", slope=-4.0, N1_coef=1.0e-6)

