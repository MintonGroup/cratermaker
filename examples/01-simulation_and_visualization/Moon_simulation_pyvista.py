"""
Run a simulation of the Moon and visualize with PyVista
=======================================================

.. rubric:: By Dennise Valadez and David Minton

This example demonstrates how to run a crater population simulation on the Moon and visualize the final surface using PyVista.

We will use the ``run`` method in a Cratermaker Simulation object for an age of 4.31 billion years using the default Neukum production function [#]_. This will produces a cratering history that is similar to that seen on the most ancient terrains on the Moon. However, keep in mind that this is not necessarily the "true" age of the surface, as the Neukum production function is very poorly constrained prior to 3.9 billion years ago. Also, because we are not setting a random seed, the exact crater population will be different each time the simulation is run.

We also reduce the gridlevel to 6 to speed up the simulation for this example. Even at low resolution, it will take several minutes to run the simulation. We will plot them using a reversed grayscale colormap to somewhat mimic the appearence of the lunar surface, or at least what the Moon would look like if every deep crater was filled with mare basalt!



References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., (2001) Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi: 10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`
"""

import pyvista as pv

import cratermaker as cm

# Initialize a quick Moon simulation. We will reduce the resolution to gridlevel 6 to speed up the simulation for this example. We set `ask_overwrite` to False to avoid being prompted during the example.
sim = cm.Simulation(target="Moon", gridlevel=6, ask_overwrite=False)

sim.run(age=4310)
sim.export(driver="VTK")
# Load and visualize the resulting surface using Pyvista
mesh = pv.read(sim.surface.output_dir / "surface000001.vtp")

plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars="face_elevation", cmap="Greys_r", show_edges=False)
plotter.show()
