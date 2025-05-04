.. _api-simulation:

Simulation
==========

The Simulation class is the main class for the cratermaker project. It is used to create a simulation of a crater on a given target 
body. The Simulation class is used to generate craters of a given size and morphology based on the production function, morphology
function, and crater scaling relationship model. The surface of the target body is represented by a Surface attribute called
`surface`, which contains a UxDataset object called `surface.uxds`. This is an unstructured grid dataset that contains data for the target body surface.

Creating a Simulation
---------------------

.. autoclass:: cratermaker.core.simulation.Simulation
   :members:
   :undoc-members:
   :no-index:
