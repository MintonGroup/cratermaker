.. _api-target:

######
Target
######

The Target class contains the properties of the celestial body being impacted. 

************************
Available Implementation
************************

Currently, there is only one implementation of Target. Passing the name of the target body, such as "Moon" or "Mars" to the :py:func:`Target.maker() <cratermaker.components.target.Target.maker>` factory method will return a Target object with the properties of the specified body. The 

.. autoclass:: cratermaker.components.target.Target
   :members:
   :undoc-members:
   :no-index-entry:

Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(target="Titan")
    print(sim.target)
