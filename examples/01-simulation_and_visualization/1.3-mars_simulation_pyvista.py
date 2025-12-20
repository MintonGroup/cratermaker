"""
Run a simulation of a local region on Mars.
===========================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates how to use the HiResLocal Surface to run a crater population over a small region of a planet, in this case Mars.

For this example, we will The simulation uses the default Mars configuration and runs for 2 billion years on a local region with a radius of 20 km and resolution of 100 m/pix.  We then visualize the surface using PyVista and a Mars-like colormap.
We also pass an option that will automatically generate hillshade plots each time the surface is saved. You can also pass in "elevation" to generate elevation plots instead of hillshade plots.


"""

import pyvista as pv
from IPython.display import Image

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
sim.run(age=2000, plot_style="elevation", cmap="pink", scalebar=True, label="Mars region simulation")
sim.show(cmap="pink")

# Alternatively, this will generate hillshade images with a default time stamp
# sim.run(age=2000, plot_style="hillshade")

# We can also display the saved image directly. The name will depend on the plot_style option used above.
Image(filename=sim.surface.plot_dir / "elevation000001.png")
