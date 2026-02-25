"""
Emplace a crater a specific distance and bearing from the center of a HiResLocalSurface
=======================================================================================

This example shows how to emplace a crater at a specific distance and bearing from the center of a HiResLocalSurface using the new `relative_location` arguments that can be passed to the `emplace` method. This is useful to place a crater on a HiResLocalSurface at a specific location relative to the center without trying to figure out what lat,lon coordinates to use. In this example, we will emplace a sequence of craters in a spiral pattern. This not meant to be a realistic planetary surface, but it looks pretty cool!

"""

import cratermaker as cm

simdir = "simdata-1_6"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages. Alternatively,
# passing ask_overwrite=False and reset=True to Simulation will also allow the example to run without requiring any prompts.
cm.cleanup(simdir)

sim = cm.Simulation(
    surface="hireslocal",
    simdir=simdir,
    local_location=(0, 0),
    pix=10.0,
    local_radius=2000.0,
    do_counting=False,
)

# Emplace craters in a spiral pattern with craters getting slightly largers as they get farther from the center
total_rotation = 3 * 360.0
ncraters = 36
bearing_spacing = int(total_rotation / ncraters)
distance_spacing = sim.surface.local_radius / ncraters
for i in range(ncraters):
    bearing = i * bearing_spacing
    distance = i * distance_spacing
    diameter = 50.0 + i * 2.5
    sim.emplace(diameter=diameter, relative_location={"distance": distance, "bearing": bearing})

# Use the "c" key to show the circle representing the crater
sim.show3d(variable_name="face_elevation", cmap="cividis")
