"""
Manually emplace a single crater on the Moon and visualize it with PyVista
==========================================================================

This example shows how to emplace a single large crater on the Moon using Cratermaker's `Simulation.emplace()` method. The crater is defined with a specific diameter and location. The final surface is exported and visualized using PyVista.

"""

import cratermaker as cm

simdir = "simdata-1_1"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages. Alternatively,
# passing ask_overwrite=False and reset=True to Simulation will also allow the example to run without requiring any prompts.
cm.cleanup(simdir)

sim = cm.Simulation(gridlevel=6, simdir=simdir)
sim.emplace(final_diameter=500e3, location=(45, 60))
sim.show(variable_name="face_elevation", cmap="cividis")
