"""
Fit a crater rim given a DEM and approximate crater size and location
=====================================================================

.. rubric:: By David Minton

In this example, we will create a DataSurface centered on a region of the Moon that contains the Lansberg B crater. We will then supply a slightly wrong crater size and location, and use the Counting class to fit the crater rim to the DEM data.

"""

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Counting, Crater, Surface
from cratermaker.utils.general_utils import format_large_units


def plot_fits(surface, crater=None, plot_score=False, imagefile=None):
    """
    Plot the surface with crater fits overlaid.

    Parameters
    ----------
    surface : Surface
        The DataSurface object to plot.
    crater : Crater, optional
        A Crater object containing the initial and/or fit crater rims to plot.
    plot_score : bool, optional
        Whether to plot the rimscore variable on the surface.
    imagefile : str, optional
        If provided, save the plot to this file instead of showing it.
    """
    W = int(max(1, np.ceil((2.0 * surface.local_radius) / surface.pix)))
    fig, ax = plt.subplots(figsize=(1, 1), dpi=W, frameon=False)
    surface.plot(show=False, style="hillshade", ax=ax)
    if crater:
        # Plot the initial guess
        crater.to_geoseries(surface=surface, split_antimeridian=False, measured=False).to_crs(surface.local.crs).plot(
            ax=ax, facecolor="none", edgecolor="cyan", linewidth=0.1, linestyle=":"
        )

        # Plot the fit
        crater.to_geoseries(surface=surface, split_antimeridian=False, measured=True).to_crs(surface.local.crs).plot(
            ax=ax, facecolor="none", edgecolor="white", linewidth=0.2
        )

    # Plot the rimscore
    if plot_score and "rimscore" in surface.uxds.data_vars:
        surface.plot(
            variable="rimscore",
            show=False,
            style="elevation",
            cmap="magma",
            ax=ax,
            imagefile=imagefile,
        )
    if imagefile is not None:
        plt.savefig(imagefile, bbox_inches="tight", pad_inches=0)
    else:
        plt.show()
    plt.close(fig)
    return


# Lansberg B is a 9 km crater relatively fresh simple crater located at (28.14°W, 2.493°S).
# Start by creating a (slightly) incorrect Crater object representing our initial guess for Lansberg B.

lansberg_b = Crater.maker(diameter=9.5e3, location=(-28.1, -2.45))

# Next, we will create a DataSurface that should be large enough to encompass the correct crater rim.
surface = Surface.maker(
    surface="datasurface",
    local_location=lansberg_b.location,
    local_radius=lansberg_b.radius * 3.0,
    ask_overwrite=False,
)

# Now refine the fit of the crater rim using the Counting class.
lansberg_b = Counting.maker(surface=surface, ask_overwrite=False).fit_rim(crater=lansberg_b, fit_ellipse=False, fit_center=True)

# If we print the crater object, we will see that the original parameters are retained, but the values from the fit are prepended by `measured_`
print(lansberg_b)

# We can plot the surface with the initial (cyan dashed line) and fitted (white solid line) crater rims overlaid.
plot_fits(surface=surface, crater=lansberg_b)

# If you want to see the score that the rim finder used, just pass `plot_score=True` to the plotting function above
plot_fits(surface=surface, crater=lansberg_b, plot_score=True)
