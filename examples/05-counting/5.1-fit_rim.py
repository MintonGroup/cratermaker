"""
Fit a crater rim given a DEM and approximate crater size and location
=====================================================================

.. rubric:: By David Minton

In this example, we will create a DataSurface centered on a region of the Moon that contains the Lansberg B crater. We will then supply a slightly wrong crater size and location, and use the Counting class to fit the crater rim to the DEM data.

"""

import os

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Crater, Simulation

simdir = "simdata-5_1"
# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.


# Lansberg B is a 9 km crater relatively fresh simple crater located at (28.14°W, 2.493°S).
# Start by creating a (slightly) incorrect Crater object representing our initial guess for Lansberg B.
lansberg_b = Crater.maker(name="Lansberg B", diameter=9.5e3, location=(-28.1, -2.45))

# Next, we will create a DataSurface that should be large enough to encompass the correct crater rim.
sim = Simulation(
    surface="datasurface",
    local_location=lansberg_b.location,
    local_radius=lansberg_b.radius * 3.0,
    simdir=simdir,
    ask_overwrite=False,
    reset=True,
)

lansberg_b = sim.counting.add(lansberg_b)
# Now refine the fit of the crater rim using the Counting class.
lansberg_b = sim.counting.fit_rim(crater=lansberg_b, fit_ellipse=False, fit_center=True)

# If we print the crater object, we will see that the original parameters are retained, but the values from the fit are prepended by `measured_`
print(lansberg_b)

# We can plot the surface with the initial (cyan dashed line) and fitted (white solid line) crater rims overlaid.
sim.plot(
    plot_style="hillshade",
    variable_name="face_elevation",
    save=False,
    show=True,
    include_counting=True,
    observed_original_color="cyan",
    observed_color="white",
    close_when_done=False,  # This is normally not necessary, but needed to render in the documentation build process
)
plt.show()

# If you want to see the score that the rim finder used, just pass `plot_score=True` to the plotting function above
sim.plot(
    plot_style="hillshade",
    variable_name="rimscore",
    cmap="magma",
    save=False,
    show=True,
    include_counting=True,
    observed_color="white",
    close_when_done=False,
)
plt.show()
