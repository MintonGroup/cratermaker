"""
Create a DataSurface centered on Kepler crater
==============================================

.. rubric:: By David Minton

This example demonstrates how to use the DataSurface to fetch real DEM data for a local region on the Moon. In this case, we center the surface on Kepler crater (321.9913E, 8.121N) with a radius of 50 km and a resolution of 200 m/pix. We then visualize the surface using PyVista both with and without the superdomain.

"""

import os

# This example can take some time due to the use of DataSurface, which must download data from the NASA PDS. For document building, we will simple display a pre-rendered image rather than running the script to save resources.
simdir = "simdata-1_5"
if not os.environ.get("READTHEDOCS", False):
    from cratermaker import Surface

    # Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
    # prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
    # own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.
    surface = Surface.maker(
        "datasurface",
        local_location=(321.9913, 8.121),
        local_radius=50.0e3,
        pix=200.0,
        simdir=simdir,
        ask_overwrite=False,
        reset=True,
    )
    if os.environ.get("RENDERFIGS", False):
        surface.plot(plot_style="hillshade", show=False, save=True)
    surface.show3d(superdomain=False)
    surface.show3d(superdomain=True)

else:
    from pathlib import Path

    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    plt.imshow(mpimg.imread(Path(simdir) / "surface_images" / "local_surface_hillshade.png"))
    plt.axis("off")
    plt.show()
