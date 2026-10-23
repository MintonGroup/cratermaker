"""
Create a shaded topographic representation of a crater
======================================================

.. rubric:: By David Minton

This example showcases how to create a crater and ejecta profile using the "basicmoon" morphology model from the Cratermaker and visual its topography. This will mimic how CTEM generates a test crater, though it is much simpler to run than that venerable old Fortran-based beast of a code!. The crater is created with a radius of 1 km. The hill shade uses the same settings that CTEM uses.

Note, that this example is rather more complex than it really needs to be. In practice, the Morphology object is always associated with a Surface object, which contains its own set of methods for visualizing the surface morphology. For this example, we are bypassing the Surface object and its unstructured mesh functionality entirely and generating our own flat grid.  Example 1.1-visualize_one_crater.py demonstrates a more "typical" approach.

"""

import numpy as np

from cratermaker import Simulation

# Generate 1000x1000 grid centered at (0, 0)
gridsize = 1000
extent = 5e3  # +/- extent in meters
pix = 2 * extent / gridsize  # pixel size in meters
x = np.linspace(-extent, extent, gridsize)
y = np.linspace(-extent, extent, gridsize)
xx, yy = np.meshgrid(x, y)

radial_distances = np.sqrt(xx**2 + yy**2)
bearings = (
    np.degrees(np.arctan2(xx, yy)) + 180
)  # Bearing angles aren't strictly necessary for the basicmoon morphology craters, but we include them here anyway to show you how to define them for a custom grid.

# Because we are not going to be using the Surface object, we will initialize it with a minimum icosphere surface so that we aren't generating a lot of output data we don't need.
sim = Simulation(simdir="simdata-4_1", gridlevel=0, morphology="basicmoon", ask_overwrite=False, reset=True)
crater = sim.Crater.maker(radius=1.0e3)

# Compute crater and ejecta profiles
dem = sim.morphology.crater_profile(
    crater=crater, radial_distances=radial_distances, bearings=bearings, include_crater=True, include_ejecta=True
)


def plot(dem, pix):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource

    dpi = 300
    gridsize = dem.shape[0]
    # Create hillshade
    azimuth = 300.0  # user['azimuth']
    solar_angle = 20.0  # user['solar_angle']
    height = gridsize / dpi
    width = gridsize / dpi
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = plt.axes([0, 0, 1, 1])

    ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
    hillshade = ls.hillshade(dem, dx=pix, dy=pix, fraction=1)

    # Plot the shaded relief
    ax.imshow(
        hillshade,
        interpolation="nearest",
        cmap="gray",
        vmin=0.0,
        vmax=1.0,
        extent=(-extent, extent, -extent, extent),
    )
    plt.axis("off")
    plt.show()


plot(dem, pix)
