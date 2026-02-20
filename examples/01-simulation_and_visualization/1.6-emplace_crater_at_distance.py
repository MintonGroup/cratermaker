"""
Emplace a crater a specific distance and bearing from the center of a HiResLocalSurface
=======================================================================================

This example shows how to emplace a crater at a specific distance and bearing from the center of a HiResLocalSurface using the new `compute_location_from_distance_bearing` function of the Surface class. This is useful to place a crater on a HiResLocalSurface at a specific location relative to the center without trying to figure out what lat,lon coordinates to use.

"""

import cratermaker as cm

simdir = "simdata-1_6"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages. Alternatively,
# passing ask_overwrite=False and reset=True to Simulation will also allow the example to run without requiring any prompts.
cm.cleanup(simdir)

sim = cm.Simulation(
    surface="hireslocal",
    simdir="craterstats_compare_lores",
    local_location=(0, 0),
    pix=10.0,
    local_radius=2000.0,
    ask_overwrite=False,
    reset=True,
    rng_seed=349572341965,
)

# Put a 500 m crater 1 km due north of the center
location = sim.surface.local.compute_location_from_distance_bearing(distances=1000.0, bearings=0)
sim.emplace(diameter=500.0, location=location)

# Use the "c" key to show the circle representing the crater
sim.show(variable_name="face_elevation", cmap="cividis")
