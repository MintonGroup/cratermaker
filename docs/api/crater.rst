.. _api-crater:

######
Crater
######

The ``Crater`` class represents a single crater in the simulation. It is used to model the crater resulting from an impact, including its size, shape, depth, and other morphological features. It also defines the properties of the projectile, such as its size, velocity, material, and angle of impact.


Base Crater Class
=================

The ``Crater`` class is the base class for all crater objects in the simulation. It contains the core properties and methods that are common to all crater types. It is a composite object that contains both fixed and variable properties of the crater, as well as a reference to the target body. Fixed and variable properties are stored internally as subclasses, but all properties are accessible directly from the crater object for ease of use.

.. autoclass:: cratermaker.components.crater.Crater
    :members:
    :undoc-members:
    :no-index-entry:

Fixed Properties
================

The ``CraterFixed`` class represents the properties of the crater that are fixed at when the crater object is instanteated. These properties include things such as the original values of morphological parameters and location. 

.. autoclass:: cratermaker.components.crater.CraterFixed
    :members:
    :undoc-members:
    :no-index-entry:

Variable Properties
===================

The ``CraterVariable`` class represents the properties of the crater that can be altered after the crater object is instantiated. These properties include things such as measured values of morphological parameters, which can be updated as the crater evovles on the surface.

.. autoclass:: cratermaker.components.crater.CraterVariable
    :members:
    :undoc-members:
    :no-index-entry:

