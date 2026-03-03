"""
Compute the degradation state of a crater
=========================================

.. rubric:: By David Minton

In this example, we will emplace a small simple crater on a HiResLocal surface, and then apply topographic diffusion to the surface in order to simulate diffusive degradation of the crater. We will then compare the applied diffusive degradation amount with the estimated degradation state of the crater using the SimpleCount model, which uses depth-to-diameter as a proxy for degradation state. If the proxy is good, this should plot as a 1:1 line.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Simulation

simdir = "simdata-5_2"

# Note, that for these examples we use cratermaker's cleanup function to prepare a fresh directory for the example to run. This will
# prevent prompts that will prevent these examples from running on their own when building the documentation pages. Alternatively,
# passing ask_overwrite=False to Surface.maker() also allow the example to run without requiring any prompts.
from cratermaker import cleanup

cleanup(simdir)


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
Kd_func = sim.counting.measure_degradation_state

# We'll use the crater's visibility function to help set the diffusion amount. We'll go to 1.5x the value of Kv
niter = 10
Kv = sim.counting.visibility_function(crater)
Kd_max = 1.5 * Kv
dK = Kd_max / niter

sim.save()
Kd_vals = [Kd_func(crater)]
K_vals = [0.0]
for K in np.linspace(dK, Kd_max, niter, endpoint=True):
    # Apply dK topographic diffusion tot he surface
    # sim.interval += 1
    sim.surface.local.apply_diffusion(dK)
    # sim.save()
    # Estimate the degradation state of the crater.

    # These should ideally be the same
    Kd_vals.append(Kd_func(crater))
    K_vals.append(K)

# Normalize by radius^2 to make the plot unitless
Kd_r2 = np.array(Kd_vals) / crater.radius**2
K_r2 = np.array(K_vals) / crater.radius**2
fig, ax = plt.subplots(layout="constrained")
ax.scatter(K_r2, Kd_r2)
ax.plot(K_r2, K_r2, color="k", linestyle="--")
ax.set_ylabel("Estimated degradation state ($K_d/r^2$)")
ax.set_xlabel("Applied diffusive degradation ($K/r^2)$")
plt.show()
