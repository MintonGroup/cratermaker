""""
Run a simulation of the Moon with basins emplaced using QuasiMC mode
====================================================================

.. rubric:: By Austin Blevins and David Minton

This example shows how to run a lunar simulation in 'Quasi-Monte Carlo' mode, where the basins are read from a csv file that gives their diameters and locations, along with either an age range or N(20) range.

We will use the ``run`` method in a Cratermaker Simulation object for an age of 4.31 billion years using the default Neukum production function [#]_. The largest craters (aka basins) will not follow this production function, but instead will be emplaced according to the csv file. The age of South Pole Aitken basin is set to simulation start, and Imbrium is set to 3.9 billion years ago. The other basins can vary based on N(20) values measured in Orgel et al. (2018) [#]_. The 74 basins catalogued by Neumann et al. (2015) [#]_ are included in this file.

Like example 1.2, we also reduce the gridlevel to 6 to speed up the simulation for this example. Even at low resolution, it will take several minutes to run the simulation. We will plot using a reversed grayscale colormap to somewhat mimic the appearence of the lunar surface.
"""

import cratermaker as cm

simdir = "simdata-1_7"

# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.


# Initialize a quick Moon simulation. We will reduce the resolution to gridlevel 6 and turn off counting to speed up the simulation for this example.
sim = cm.Simulation(target="Moon", gridlevel=6, do_counting=False, simdir=simdir, ask_overwrite=False, reset=True, quasimc_file="qmc_input.csv")

sim.run(age=4310)
sim.show3d(variable_name="face_elevation", cmap="Greys_r")