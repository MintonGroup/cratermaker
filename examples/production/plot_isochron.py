"""
Plot isochrons for the Moon and Mars for 1 Ma, 1 Ga, and 4 Ga using the Neukum Production Function
==================================================================================================

.. rubric:: By Austin Blevins and David Minton

In this example, we will be using the "neukum" prodcution model in Cratermaker to plot isochrons for three different age surfaces. We will also format the plots so with similar axes as Figure 2 of Neukum, Ivanov, and Hartmann (2001) [#]_, but with a handy logarithmic grid.

The Neukum production function (NPF) is only defined over a limited range, from crater diameters of 0.01 km to 300 km. Because Cratermaker is intended to primarily be used for forward modeling of landscape evolution where regolith mixing by small impactors plays a significant role, we require a production function that can extend to much smaller diameters. Therefore we extrapolate the lower end of the NPF with a simple power law with a slope equal to that of the lower end.

The upper range of the NPF must also be extrapolated, but for a different reason. If we truncated the NPF at 300 km, this would artificially limit the size of the largest craters that could be produced in a simulation and simulations of the ancient lunar or martian surface would contain dozens of 300 km craters. To avoid this, we extrapolate the upper end of the NPF, however we cannot simply extrapolate the upper end with its own slope, as this would result in far too many very large craters. This is similar to the issue that was discussed in Minton et al. (2015) [#]_. Therefore we steepen the slope of the upper end by shifting the exponent by -2. Futher work would be needed to determine the best fit to the upper end of the NPF.

For this example, we plot the valid range of the NPF in solid black, and the extrapolations in orange with a dash-dot line style.

References
----------

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., (2001) Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55-86. `doi: 10.1023/A:1011989004263 <https://doi.org/10.1023/A:1011989004263>`
.. [#] Minton, D.A., Richardson, J.E., Fassett, C.I., (2015) Re-examining the main asteroid belt as the primary source of ancient lunar craters. Icarus 247, 172-190. `doi: 10.1016/j.icarus.2014.10.018 <https://doi.org/10.1016/j.icarus.2014.10.018>`

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

from cratermaker import Production

fig, axs = plt.subplots(1, 2, figsize=(8, 7))
axes = dict(zip(["Moon", "Mars"], axs))
tvals = [0.01, 1.0, 4.0]
x_min = 1e-3
x_max = 1e4
y_min = 1e-12
y_max = 1e6
nD = 1000
Dvals = np.logspace(np.log10(x_min), np.log10(x_max), num=nD)
for version, ax in axes.items():
    production = Production.maker("neukum", version=version)
    ax.title.set_text(version)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("$\\mathregular{N_{>D} (km^{-2})}$")
    ax.set_xlabel("Diameter (km)")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # This section will create a proper log-log grid
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100)
    )
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(
        ticker.LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100)
    )
    ax.grid(True, which="minor", ls="-", lw=0.5, zorder=5)
    ax.grid(True, which="major", ls="-", lw=1, zorder=10)

    # Split the valid and extrapolated range to plot them in different line styles and colors
    # Cratermaker uses meters for length and My for time, so we need to convert to km and Ga
    # for the plots to match the figure in Neukum et al. (2001)
    inrange = (Dvals > production.sfd_range[0] * 1e-3) & (
        Dvals < production.sfd_range[1] * 1e-3
    )
    lo = Dvals < production.sfd_range[0] * 1e-3
    hi = Dvals > production.sfd_range[1] * 1e-3
    for t in tvals:
        Nvals = production.function(diameter=Dvals * 1e3, age=t * 1e3)
        Nvals = Nvals * 1e6
        ax.plot(
            Dvals[inrange], Nvals[inrange], "-", color="black", linewidth=1.0, zorder=50
        )
        ax.plot(Dvals[lo], Nvals[lo], "-.", color="orange", linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], "-.", color="orange", linewidth=2.0, zorder=50)
        labeli = int(0.25 * nD)
        ax.text(
            Dvals[labeli],
            3 * Nvals[labeli],
            f"{t:.2f} ",
            ha="left",
            va="top",
            rotation=-72,
        )

plt.tight_layout()
plt.show()
