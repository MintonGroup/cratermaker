.. _api-simulation:

##########
Simulation
##########

The Simulation class is the main class for the cratermaker project. It is used to create a simulation of a crater on a given target 
body. The Simulation class is used to generate craters of a given size and morphology based on the production function, morphology
function, and crater scaling relationship model. The surface of the target body is represented by a Surface attribute called
`surface`, which contains a UxDataset object called `surface.uxds`. This is an unstructured grid dataset that contains data for the target body surface.
Arguments to any of the component models can be passed in as arguments to a Simulation instantiation.

.. autoclass:: cratermaker.core.simulation.Simulation
   :members:
   :undoc-members:
   :no-index-entry:


Usage example
=============

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(simdir="ceres_local_sim",
                     target="Ceres", 
                     surface="hireslocal", 
                     pix=50, 
                     local_radius=50e3, 
                     local_location=(58.2,25.4), 
                     projectile="generic", 
                     sample=True, 
                     projectile_mean_velocity=5e3, 
                     projectile_density=1500.0,
                     ask_overwrite=False,
                     reset=True,
                     )
    print(sim)

