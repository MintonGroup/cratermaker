"""
Plot a distribution of random velocities given a mean
=====================================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates the generation of random impact velocities using the :ref:`Projectile <ug-projectile>` component. We use the "generic" projectile model with a fixed mean velocity, and draw new projectiles with :py:meth:`~cratermaker.components.projectile.Projectile.new_projectile`, which samples each impact velocity from a modified Maxwell-Boltzmann distribution that accounts for the escape velocity of the target. If vmean > vescape, then the distribution will be adjusted so that the escape velocity is the minimum velocity by computing an encounter velocity then summing the encounter and escape velocities in quadrature. If vmean < vescape, then the distribution will be a truncated Maxwellian.

Here we sample 10000 projectile velocities on the Moon with a mean of 0.75 and 2 times the lunar escape velocity. Then we compare a histogram of the generated velocities with the unmodified Maxwell-Boltzmann distribution that has the same mean.

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import maxwell

from cratermaker import Projectile, Target

# Use Moon as target
target = Target.maker("Moon")
vescape = target.escape_velocity
size = 10000

factors = np.array([0.75, 2.0])
vmean_values = factors * vescape
titles = [f"vmean = {float(factor):0.2f} $\\times$ vescape" for factor in factors]

fig, axs = plt.subplots(2, 1, figsize=(8, 8), layout="constrained")

for row, vmean in enumerate(vmean_values):
    # Initialize a generic projectile population with a fixed mean velocity. New projectiles always account for the target's escape velocity.
    proj = Projectile.maker("generic", projectile_mean_velocity=vmean, projectile_density=3000, target=target)

    velocities = np.array([proj.new_projectile().velocity for _ in range(size)])

    bins = np.linspace(0, vmean * 3 / vescape, 50)
    histogram, bins = np.histogram(velocities / vescape, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    theoretical_dist = maxwell.pdf(bin_centers * vescape, scale=vmean / np.sqrt(8 / np.pi)) * vescape

    ax = axs[row]
    ax.bar(
        bin_centers,
        histogram,
        width=bins[1] - bins[0],
        alpha=0.5,
        label="Observed",
    )
    ax.plot(
        bin_centers,
        theoretical_dist,
        label="Maxwell-Boltzmann",
        color="red",
    )
    rms = np.sqrt(np.mean(velocities**2)) / vescape
    ax.axvline(rms, color="orange", linestyle="--", label="RMS Velocity")
    ax.annotate(
        f"RMS = {rms:.2f}",
        xy=(rms, ax.get_ylim()[1] * 0.65),
        xytext=(5, 0),
        textcoords="offset points",
        rotation=90,
        va="top",
        ha="left",
        fontsize=12,
    )
    ax.set_ylabel("Probability Density")
    ax.set_xlabel("$V / V_{{escape}}$")
    ax.legend()
    ax.set_title(titles[row])

plt.show()
