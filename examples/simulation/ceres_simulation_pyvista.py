"""
Manually emplace craters on Ceres and visualize
===============================================

This example shows how to emplace a few craters on Ceres using Cratermaker’s
`Simulation.emplace()` method. Craters are defined with specific diameters and optional
locations. The final surface is exported and visualized using PyVista.

"""

import cratermaker as cm
import pyvista as pv

# Create a simulation on Ceres with default surface and morphology
sim = cm.Simulation(target="Ceres")

# Emplace craters manually with fixed diameters (in meters)
sim.emplace(final_diameter=40e3, location=(0, 0))       # Equatorial crater
sim.emplace(final_diameter=60e3, location=(45, 120))    # Mid-lat crater
sim.emplace(final_diameter=20e3, location=(-30, -60))   # Southern hemisphere

# Save and export the result
sim.export(format="vtp")

# Load and display in greyscale
mesh = pv.read("export/surface000000.vtp")
scalars = "node_elevation" if "node_elevation" in mesh.array_names else mesh.array_names[0]

plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars=scalars, cmap="gray", show_edges=False)
plotter.show()
