.. _api-crater:

######
Crater
######

The ``Crater`` class represents a single crater in the simulation. It is used to model the crater resulting from an impact, including its size, shape, depth, and other morphological features. It also defines the properties of the projectile, such as its size, velocity, material, and angle of impact.

************************
Available Implementation
************************


The ``Crater`` class is the base class for all crater objects in the simulation. It contains the core properties and methods that are common to all crater types. It is a composite object that contains both fixed and variable properties of the crater, as well as a reference to the target body. Fixed and variable properties are stored internally as subclasses, but all properties are accessible directly from the crater object for ease of use.

.. autoclass:: cratermaker.components.crater.Crater
    :members:
    :undoc-members:
    :no-index-entry:

***********
CraterFixed
***********

The ``CraterFixed`` class represents the properties of the crater that are fixed at when the crater object is instantiated. These properties include things such as the original values of morphological parameters and location. This class is not meant to be instantiated directly, but rather it is used internally by the `Crater` class to store fixed properties of the crater. 

.. autoclass:: cratermaker.components.crater.CraterFixed
    :members:
    :undoc-members:
    :no-index-entry:

**************
CraterVariable
**************

The ``CraterVariable`` class represents the properties of the crater that can be altered after the crater object is instantiated. These properties include things such as measured values of morphological parameters, which can be updated as the crater evovles on the surface. This class is not meant to be instantiated directly, but rather it is used internally by the `Crater` class to store variable properties of the crater.

.. autoclass:: cratermaker.components.crater.CraterVariable
    :members:
    :undoc-members:
    :no-index-entry:

