"""
Create a DataSurface centered on Kepler crater
==============================================

.. rubric:: By David Minton

This example demonstrates how to use the DataSurface to fetch real DEM data for a local region on the Moon. In this case, we center the surface on Kepler crater (321.9913E, 8.121N) with a radius of 50 km and a resolution of 200 m/pix. We then visualize the surface using PyVista both with and without the superdomain.

"""

from cratermaker import Surface

simdir = "simdata-1_5"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages.Alternatively,
# passing ask_overwrite=False to Surface.maker() also allow the example to run without requiring any prompts.
from cratermaker import cleanup

cleanup(simdir)

surface = Surface.maker("datasurface", local_location=(321.9913, 8.121), local_radius=50.0e3, pix=200.0, simdir=simdir)
surface.show(superdomain=False)
surface.show(superdomain=True)
