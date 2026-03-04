.. _api-morphology:

##########
Morphology
##########

The Morphology class is an operations class for computing the morphology of a crater based on its size and target properties. It encapsulates the logic for altering the topography of the surface based on the crater properties.

*************************
Available Implementations
*************************

+----------------------------------------------------------------------+----------------+---------------------------------------------+
| Class                                                                | Argument name  | Example Usage                               |
+======================================================================+================+=============================================+
| :py:class:`~cratermaker.components.morphology.simplemoon.SimpleMoon` | "simplemoon"   | morphology = Morphology.maker("simplemoon") |
+----------------------------------------------------------------------+----------------+---------------------------------------------+

.. autoclass:: cratermaker.components.morphology.Morphology
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.MorphologyCrater
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.MorphologyCraterVariable
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.CraterQueueManager
   :members:
   :undoc-members:
   :no-index-entry:


.. currentmodule:: cratermaker.components.morphology.simplemoon

**********
SimpleMoon
**********

See `Morphology`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.morphology.simplemoon.SimpleMoon
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.simplemoon.SimpleMoonCrater
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.simplemoon.SimpleMoonCraterFixed
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python
   :linenos:

    from cratermaker import Morphology
    morphology = Morphology.maker("simplemoon")