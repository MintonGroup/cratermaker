"""
Plot isochrons for the Moon and Mars for 1 Ma, 1 Ga, and 4 Ga using the Neukum Production Function
==================================================================================================

.. rubric:: By Austin Blevins and David Minton

In this example, we will be using the "neukum" prodcution model in Cratermaker to plot isochrons for three different age surfaces. We will also format the plots so with similar axes as Figure 2 of Neukum, Ivanov, and Hartmann (2001) [#]_. 

.. [#] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55–86. https://doi.org/10.1023/A:1011989004263
"""

from cratermaker import Production
import matplotlib.pyplot as plt
from matplotlib import ticker 
import numpy as np

fig, axs = plt.subplots(1, 2, figsize=(8, 7))
axes = dict(zip(['Moon','Mars'],axs))
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
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel('$\\mathregular{N_{>D} (km^{-2})}$')
    ax.set_xlabel('Diameter (km)')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2,10), numticks=100))
    ax.grid(True,which="minor",ls="-",lw=0.5,zorder=5)
    ax.grid(True,which="major",ls="-",lw=1,zorder=10)
    inrange = (Dvals > production.sfd_range[0]*1e-3) & (Dvals < production.sfd_range[1]*1e-3)
    lo = Dvals < production.sfd_range[0]*1e-3
    hi = Dvals > production.sfd_range[1]*1e-3
    for t in tvals:
        Nvals = production.function(diameter=Dvals*1e3,age=t*1e3)
        Nvals = Nvals * 1e6 
        ax.plot(Dvals[inrange], Nvals[inrange], '-', color='black', linewidth=1.0, zorder=50)
        ax.plot(Dvals[lo], Nvals[lo], '-.', color='orange', linewidth=2.0, zorder=50)
        ax.plot(Dvals[hi], Nvals[hi], '-.', color='orange', linewidth=2.0, zorder=50)
        labeli = int(0.25*nD)
        ax.text(Dvals[labeli],3*Nvals[labeli],f"{t:.2f} ", ha="left", va="top",rotation=-72)

plt.tight_layout()
plt.show()