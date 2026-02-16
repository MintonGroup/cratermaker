"""
Plot the ray intensity map for a crater
=======================================

.. rubric:: By David Minton

This example demonstrates how to create a ray intensity map for a crater using the "simplemoon" morphology model.

The ray intensity model simulates the spatial modulation of ejecta distribution in the form of rays, a feature observed around many planetary impact craters. These rays are modeled using a layered pattern of radially aligned Gaussian contributions, whose intensity decays with radial distance from the crater rim. The number and angular spread of rays vary stochastically to produce a realistic, azimuthally varying intensity field.

The model uses the ray intensity function compiled from Rust code for performance. This function computes the contribution of each ray at a given point based on its radial distance and angular bearing from the crater center. The output intensity is plotted in a 2D space normalized to the crater radius and visualized on a logarithmic color scale to emphasize the structure of the rays.

Note, that this example is rather more complex than it really needs to be. In practice, the Morphology object is always associated with a Surface object, which contains its own set of methods for visualizing the surface morphology. For this example, we are bypassing the Surface object functionality entirely and generating our own grid for illustration purposes. Example 1.1-visualize_one_crater.py demonstrates a more "typical" approach.

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from cratermaker import Crater, Morphology

simdir = "simdata-4_4"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages.Alternatively,
# passing ask_overwrite=False to Morphology.maker() also allow the example to run without requiring any prompts.
from cratermaker import cleanup

cleanup(simdir)


crater = Crater.maker(radius=10.0e3)
# Because we are not explicitly passing a Surface object, the Morphology constructor will generate a default surface. We pass the "simdir" and "gridlevel" arguments to control the Surface generation, even though we don't make use of it directly here.
morphology = Morphology.maker(simdir=simdir, gridlevel=4)
rc = crater.radius

grid_size = 1000
rmax = 20
morphology.ejecta_truncation = rmax
extent = rmax * rc
x = np.linspace(-extent, extent, grid_size)
y = np.linspace(-extent, extent, grid_size)
xx, yy = np.meshgrid(x, y)

# Convert Cartesian to polar coordinates
r = np.sqrt(xx**2 + yy**2)
theta = np.degrees(np.arctan2(yy, xx)) + 180.0

# Compute ray intensity
r_flat = r.ravel()
theta_flat = theta.ravel()
intensity_flat = morphology.ray_intensity(crater, r_flat, theta_flat)
intensity = intensity_flat.reshape(r.shape)
# Plot the intensity colormap
fig, ax = plt.subplots(figsize=(6, 5))
c = ax.imshow(
    intensity,
    norm=LogNorm(vmin=1e-4, vmax=intensity.max()),
    origin="lower",
    extent=[-rmax, rmax, -rmax, rmax],
    cmap="inferno",
)
ax.set_xlabel("x (r_c)")
ax.set_ylabel("y (r_c)")
ax.set_title("Ray Intensity Map")
fig.colorbar(c, ax=ax, label="Intensity")

plt.tight_layout()
plt.show()
print("Done")
