"""
Plot a distribution of random impact angles
===========================================

.. rubric:: By David Minton

This example demonstrates the generation of random impact angles.

"""

import matplotlib.pyplot as plt
import numpy as np

from cratermaker.utils.montecarlo_utils import get_random_impact_angle

# Sample data generation
size = 10000
angles = get_random_impact_angle(size=size)

# Number of bins
bins = 50
observed_counts, bins_ang = np.histogram(angles, bins=bins, range=(0.0, 90.0))

# Calculate expected distribution
uniform_dist = np.linspace(0, 1, size)
transformed_angles = np.rad2deg(np.arcsin(np.sqrt(uniform_dist)))
expected_counts, _ = np.histogram(transformed_angles, bins=bins, range=(0.0, 90.0))

# Plotting
fig, ax = plt.subplots(figsize=(8, 4))

# Observed counts
ax.bar(
    bins_ang[:-1],
    observed_counts,
    width=np.diff(bins_ang),
    align="edge",
    label="Observed",
    alpha=0.5,
)

# Expected counts
ax.plot(bins_ang[:-1], expected_counts, label="Expected", color="red")

ax.set_xlabel("Impact Angle (deg)")
ax.set_ylabel("Count")
ax.legend()
ax.set_title("Impact Angle Distribution")

plt.tight_layout()
plt.show()
