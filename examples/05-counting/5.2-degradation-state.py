"""
Compute the degradation state of a crater
=========================================

.. rubric:: By David Minton

In this example, we will emplace a small simple crater on a HiResLocal surface, and then apply topographic diffusion to the surface in order to simulate diffusive degradation of the crater. We will then compare the applied diffusive degradation amount with the estimated degradation state of the crater using the DepthCount model, which uses depth-to-diameter as a proxy for degradation state. If the proxy is good, the values of degradation state computed by Cratermaker vs applied degradation should plot along a 1:1 line.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Simulation

simdir = "simdata-5_2"

# Note, that for these examples we pass ask_overwrite=False and reset=True to the Simulation constructor. This will suppress
# prompts that ask the user if they want to overwrite existing files, which would would prevent these examples from running on their
# own when building the documentation pages. Alternatively, calling cm.cleanup(simdir) will remove all pre-existing output files.

sim = Simulation(
    surface="hireslocal",
    simdir=simdir,
    local_location=(0, 0),
    pix=5.0,
    local_radius=500.0,
    ask_overwrite=False,
    reset=True,
)

# Make a 100 m crater
diameter = 100.0
crater = sim.emplace(diameter=diameter, location=(0, 0))[0]
K_func = sim.counting.measure_degradation_state

# We'll use the crater's visibility function to help set the diffusion amount. We'll go to 1.5x the value of Kv
niter = 10
Kv = sim.counting.visibility_function(crater)
K_max = 1.5 * Kv
dK = K_max / niter

sim.save()
K_measured_vals = [K_func(crater)]
K_applied_vals = [0.0]
for K in np.linspace(dK, K_max, niter, endpoint=True):
    # Apply dK topographic diffusion to the surface
    sim.surface.local.apply_diffusion(dK)

    # Record what Cratermaker computes the degradation state to be
    K_measured_vals.append(K_func(crater))

    # Record what the degradation state ought to be
    K_applied_vals.append(K)

# Normalize by radius^2 to make the plot unitless
K_r2_measured = np.array(K_measured_vals) / crater.radius**2
K_r2_applied = np.array(K_applied_vals) / crater.radius**2
fig, ax = plt.subplots(layout="constrained")
ax.scatter(K_r2_applied, K_r2_measured, label="Cratermaker simulation")
ax.plot(K_r2_applied, K_r2_measured, color="k", linestyle="--", label="1:1")
ax.set_xlabel("Applied diffusive degradation ($K/r^2)$")
ax.set_ylabel("Estimated degradation state ($K/r^2$)")
ax.legend()
plt.show()
