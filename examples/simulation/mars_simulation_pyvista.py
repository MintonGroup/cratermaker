"""
Visualize a Mars crater simulation with PyVista
===============================================

This example demonstrates how to run a short crater population simulation on Mars
and visualize the final surface state using PyVista.

The simulation uses the default Mars configuration and runs for 1 billion years.
Cratermaker automatically exports the surface as a VTK PolyData file (`.vtp`),
which is then visualized using pyvista.

"""

import cratermaker as cm
import pyvista as pv

# Run a simple Mars sim
sim = cm.Simulation(target="Mars")
sim.run(age=1000)  # Simulate 1 billion years

# Load surface
mesh = pv.read("export/surface000001.vtp")

# Visualize using PyVista
scalars = "node_elevation" if "node_elevation" in mesh.array_names else mesh.array_names[0]
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars=scalars, cmap="terrain", show_edges=False)
plotter.show()
