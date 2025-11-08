"""
Plot a distribution of random velocities given a mean
=====================================================

.. rubric:: By David Minton and Dennise Valadez

This example demonstrates the generation of random velocities using the :func:`get_random_velocity` function from the :ref:`montecarlo_utils <api-utils-montecarlo>` module. The expected distribution of imapct velocities is a modified Maxwell-Boltzmann distribution. The modification is related to the relationship between the mean velocity and the escape velocity. If vmean > vescape, then the distribution will adjusted so that the escape velocity is the minimum velocity by computing an encounter velocity then summing the encounter and escape velocities in quadrature. If vmean < vescape, then the distirbution will be a truncated maxwellian.

Here we sample 10000 velocities with a mean of 0.5 and 2x the escape velocity of the Moon, both with and without the escape velocity argument. Then we will compare a histogram of generated velocities with .

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

fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex="col", sharey="row")

for col, vmean in enumerate(vmean_values):
    # Initialize projectiles with and without escape velocity
    proj = Projectile(mean_velocity=vmean, density=3000)

    velocities = np.array([proj.new_projectile().velocity for _ in range(size)])

    bins = np.linspace(0, vmean * 3 / vescape, 50)
    histogram, bins = np.histogram(velocities / vescape, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    theoretical_dist = (
        maxwell.pdf(bin_centers * vescape, scale=vmean / np.sqrt(8 / np.pi)) * vescape
    )

    ax = axs[col]
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
    if col == 0:
        ax.set_ylabel("Probability Density")
        ax.set_xlabel("$V / V_{{escape}}$")
        ax.legend()
        ax.set_title(titles[col])

plt.tight_layout()
plt.show()
