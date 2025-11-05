"""
Run a simulation of a local region on Mars.
===========================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates how to use the HiResLocal Surface to run a crater population over a small region of a planet, in this case Mars.

For this example, we will The simulation uses the default Mars configuration and runs for 1 billion years.  Cratermaker automatically exports the surface as a VTK PolyData file (`.vtp`), which is then visualized using pyvista.

"""

import pyvista as pv

import cratermaker as cm

# Run a simple Mars sim
sim = cm.Simulation(
    target="Mars", surface="hireslocal", local_location=(0, 0), pix=100.0, local_radius=50.0e3, ask_overwrite=False, reset=True
)
sim.run(age=1000)  # Simulate 1 billion years
sim.export(format="vtp")
# Load surface
mesh = pv.read(sim.surface.output_dir / "local_surface000001.vtp")

# Visualize using PyVista
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars="face_elevation", cmap="Oranges", show_edges=False)
plotter.show()
