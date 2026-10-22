.. _api-morphology:

##########
Morphology
##########

The Morphology class is an operations class for computing the morphology of a crater based on its size and target properties. It encapsulates the logic for altering the topography of the surface based on the crater properties.

*************************
Available Implementations
*************************

+-------------+----------------+---------------------------------------------+
| Class       | Argument name  | Example Usage                               |
+=============+================+=============================================+
| |BasicMoon| | "basicmoon"    | morphology = Morphology.maker("basicmoon")  |
| |RealMoon|  | "realmoon"     | morphology = Morphology.maker("realmoon")   |
+-------------+----------------+---------------------------------------------+

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


.. currentmodule:: cratermaker.components.morphology.basicmoon

*******************
BasicMoonMorphology
*******************

See `Morphology`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.morphology.basicmoon.BasicMoonMorphology
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.basicmoon.BasicMoonCrater
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.basicmoon.BasicMoonCraterFixed
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(morphology="basicmoon")
    print(sim.morphology)



******************
RealMoonMorphology
******************

See `Morphology`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.morphology.realmoon.RealMoonMorphology
   :members:
   :undoc-members:
   :no-index-entry:


.. autoclass:: cratermaker.components.morphology.realmoon.RealMoonCrater
   :members:
   :undoc-members:
   :no-index-entry:

.. autoclass:: cratermaker.components.morphology.realmoon.PSD1D
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(morphology="realmoon")
    print(sim.morphology)