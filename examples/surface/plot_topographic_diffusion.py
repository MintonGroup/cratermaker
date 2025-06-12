"""
Topographic Diffusion
=====================

.. rubric:: By Dennise Valadez and David Minton

This example demonstrates how to use the :meth:`apply_diffusion` method in the :ref:`Surface <ug-surface>` class to model topographic diffusion. In this example, we will simulate the change in elevation over time of a hill with a Gaussian profile. This example has an analytical solution, which we will compare against the numerical solution provided by the diffusion method.

For this example we will use the "hireslocal" surface type, which will give us a high-resolution local region of the surface that will be approximately a flat surface.

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from cratermaker import Surface


def analytical_elevation(r, h0, sigma, kappa, t):
    """Analytical solution for elevation profile after diffusion."""
    variance_t = sigma**2 + 2 * kappa * t
    return h0 * (sigma**2 / variance_t) * np.exp(-(r**2) / (2 * variance_t))


# We will set up a surface with a 2 km radius high resolution region at 10 m/pix. The location is arbitrary, but we will set it to (0, 0) for simplicity (the equator/prime meridian).
local_radius = 2000.0
pix = 10.0
local_location = (0, 0)
superdomain_scale_factor = 10000

# The initial elevation will be a Gaussian hill with a height of 100 m and a standard deviation of 200 m.
h0 = 100.0
sigma = 200.0

# The diffusivity κ will be set to 5000 m² per time step (in arbitrary time units), and we will run the simulation for a total time of 10 time steps.
kappa = 5000.0
nsteps = 10
total_time = 10.0  # total simulated time in arbitrary units
dt = total_time / nsteps

surface = Surface.maker(
    "hireslocal",
    local_location=local_location,
    pix=pix,
    local_radius=local_radius,
    superdomain_scale_factor=superdomain_scale_factor,
)

# Generate the initial elevation profile and apply it to the surface
r = surface.local.face_distance
h = analytical_elevation(r, h0, sigma, kappa, 0.0)
isort = np.argsort(r)
surface.local.update_elevation(h)

# We will loop over the number of steps, applying diffusion at each step and plotting the results, starting from the initial condition
for i in range(nsteps + 1):
    t = i * dt
    label_num = "Initial" if i == 0 else f"Step {i}"
    label_ana = "Analytical" if i == nsteps else None
    h_numerical = surface.local.face_elevation
    h_analytical = analytical_elevation(r, h0, sigma, kappa, t)
    plt.plot(r[isort], h_numerical[isort], label=label_num)
    plt.plot(r[isort], h_analytical[isort], linestyle="--", color="black", alpha=0.5, label=label_ana)
    # Compare difference between analytical and numerical statistically
    rmse = np.sqrt(np.mean((h_numerical - h_analytical) ** 2))
    nrmse = rmse / np.ptp(h_analytical)
    print(f"{label_num:<7}: Normalized root mean square error = {nrmse:.4%}")
    surface.local.apply_diffusion(kappa)

plt.xlabel("Distance (m)")
plt.ylabel("Elevation (m)")
plt.xlim(0, local_radius)

plt.title("Elevation Profile with Topographic Diffusion")
plt.legend()
plt.tight_layout()
plt.show()
