"""
Load a DEM of phobos
====================

This example shows how to load arbitrary DEM data using the data composer.

"""

from cratermaker import Simulation

simdir = "simdata-1_7"

# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cratermaker.cleanup(simdir) will remove all pre-existing output files.
sim = Simulation(
    target="Phobos",
    surface="icosphere",
    simdir=simdir,
    gridlevel=7,
    ask_overwrite=False,
    reset=True,
)

with sim.surface.data_composer() as composer:
    composer.add_data("https://planetarymaps.usgs.gov/mosaic/Phobos_ME_HRSC_DEM_Global_2ppd.tif")

sim.show3d()
