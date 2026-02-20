"""
Run a simulation of a local region on Mars.
===========================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates how to use the HiResLocal Surface to run a crater population over a small region of a planet, in this case Mars.

For this example, we will The simulation uses the default Mars configuration and runs for 1 billion years on a local region with a radius of 20 km and resolution of 100 m/pix.  We then visualize the surface using PyVista and a Mars-like colormap.
We also pass an option that will automatically generate hillshade plots each time the surface is saved. You can also pass in "elevation" to generate elevation plots instead of hillshade plots.


"""

from IPython.display import Image

import cratermaker as cm

simdir = "simdata-1_3"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages. Alternatively,
# passing ask_overwrite=False and reset=True to Simulation will also allow the example to run without requiring any prompts.
cm.cleanup(simdir)


# Run a simple Mars sim
sim = cm.Simulation(
    target="Mars",
    surface="hireslocal",
    local_location=(0, 0),
    pix=200.0,
    local_radius=20.0e3,
    rng_seed=86186233405,  # This will ensure we get the same crater population each time we run the example
    simdir=simdir,
)
# This will automatically generate a hillshade plot every time the simulation is saved. All components of Cratermaker have a save_actions property that can be used to specify actions to be performed when the save function is called. This is useful for automatically generating plots or other outputs at specified intervals during the simulation.
sim.save_actions = [
    {
        "plot": {
            "plot_style": "hillshade",
            "cmap": "pink",
            "scalebar": True,
            "label": "Mars region simulation",
            "show": True,
            "save": True,
        }
    }
]
sim.run(age=1000)
sim.show3d(variable_name="face_elevation", cmap="pink")

# Alternatively, this will generate hillshade images with a default time stamp
# sim.run(age=1000)

# We can also display the saved image directly. The name will depend on the plot_style option used above.
Image(filename=sim.surface.plot_dir / "local_surface_hillshade000001.png")
