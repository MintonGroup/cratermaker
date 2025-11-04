"""
Visualize a Mars crater simulation with PyVista
===============================================

This example demonstrates how to run a short crater population simulation on Mars and visualize the final surface state using PyVista.

The simulation uses the default Mars configuration and runs for 1 billion years.  Cratermaker automatically exports the surface as a VTK PolyData file (`.vtp`), which is then visualized using pyvista.

"""

import pyvista as pv

import cratermaker as cm

# Run a simple Mars sim
sim = cm.Simulation(target="Mars", ask_overwrite=False)
sim.run(age=1000)  # Simulate 1 billion years
sim.export(format="vtp")
# Load surface
mesh = pv.read(sim.surface.output_dir / "surface000001.vtp")

# Visualize using PyVista
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars="face_elevation", cmap="terrain", show_edges=False)
plotter.show()
