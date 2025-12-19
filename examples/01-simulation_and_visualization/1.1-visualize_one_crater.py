"""
Manually emplace a single crater on the Moon and visualize it with PyVista
==========================================================================

This example shows how to emplace a single large crater on the Moon using Cratermaker's `Simulation.emplace()` method. The crater is defined with a specific diameter and location. The final surface is exported and visualized using PyVista.

"""

import pyvista as pv

import cratermaker as cm

sim = cm.Simulation(ask_overwrite=False)
sim.emplace(final_diameter=500e3, location=(45, 60))
sim.show(cmap="cividis")
