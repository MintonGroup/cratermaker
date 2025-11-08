"""
Run a simulation of a local region on Mars.
===========================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates how to use the HiResLocal Surface to run a crater population over a small region of a planet, in this case Mars.

For this example, we will The simulation uses the default Mars configuration and runs for 2 billion years on a local region with a radius of 20 km and resolution of 100 m/pix.  We then visualize the surface using PyVista and a Mars-like colormap.


"""

import pyvista as pv

import cratermaker as cm

# Run a simple Mars sim
sim = cm.Simulation(
    target="Mars",
    surface="hireslocal",
    local_location=(0, 0),
    pix=100.0,
    local_radius=20.0e3,
    ask_overwrite=False,
    rng_seed=86186233406,  # This will ensure we get the same crater population each time we run the example
)
sim.run(age=2000)

sim.export(driver="VTK")
mesh = pv.read(sim.surface.output_dir / "local_surface000001.vtp")
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars="face_elevation", cmap="pink", show_edges=False)
# Orient the camera to look down on the surface
plotter.camera_position = "yz"
plotter.show()
