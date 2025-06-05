
"""
Run a simulation of the Moon and visualize with PyVista
===============================================

This example demonstrates how to run a short crater population simulation on the Moon
and visualize the final surface using PyVista.

The simulation uses generates just a few large craters on the Lunar surface to keep the run time short. Cratermaker automatically exports the surface as a VTK PolyData
file (`.vtp`), which is then visualized using PyVista.

"""

import cratermaker as cm
import pyvista as pv

# Initialize a quick Moon simulation
sim = cm.Simulation(target="Moon")

# Simulate 10 craters > 50 km
sim.run(diameter_number=(50e3, 10), ninterval=1)

# Load and visualize the resulting surface using Pyvista
mesh = pv.read("export/surface000001.vtp")
scalars = "node_elevation" if "node_elevation" in mesh.array_names else mesh.array_names[0]

plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars=scalars, cmap="gray", show_edges=False)
plotter.show()
