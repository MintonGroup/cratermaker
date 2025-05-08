"""
Plot a distribution of random sizes
===================================

.. rubric:: By David Minton

This example demonstrates the generation of random impact sizes, which is used to sample production function SFDs. In this example, we just draw from a simple power law.

"""

import matplotlib.pyplot as plt
import numpy as np

from cratermaker.utils.montecarlo_utils import get_random_size

num_realizations = 100
nbins = 10
size = 100000
Dhi = 100.0
p = 2.0
C = Dhi**p
Dlo = (size / C) ** (-1.0 / p)
diameters = np.exp(np.linspace(np.log(Dlo), np.log(100 * Dhi), nbins))
# Test a simple power law SFD
cdf = C * diameters**-p

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)

ax.set_xlabel("$D$")
ax.set_ylabel("$N_{>D}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim([1.0, 1.1 * size])

for i in range(num_realizations):
    new_diameters = get_random_size(diameters, cdf, mu=size)
    new_diameters = np.sort(new_diameters)[::-1]
    nval = np.arange(1, new_diameters.size + 1)
    if i == 0:
        label = "Sampled SFD"
    else:
        label = None

    ax.plot(new_diameters, nval, label=label, alpha=0.25, color="cornflowerblue")

ax.plot(diameters, cdf, label="Model SFD", color="black")
ax.legend(loc="upper right")

plt.tight_layout()
plt.show()
