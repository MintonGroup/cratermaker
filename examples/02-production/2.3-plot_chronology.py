"""
Plot the NPF Chronology function for the Moon and Mars
=======================================================

.. rubric:: By Austin Blevins and David Minton

Here we plot the chronology component of the Neukum production function for the Moon and Mars. The chronology function gives the cumulative number of craters larger than 1 km in diameter per unit area, N(1), that accumulate on a surface of a given age. The plot shows N(1) per square kilometer as a function of surface age.
"""

import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Production

moon = Production.maker("neukum", version="Moon")
mars = Production.maker("neukum", version="Mars")
age_vals = np.linspace(4.5e3, 0.0, num=1000)

N1_moon = moon.chronology(age=age_vals)
N1_mars = mars.chronology(age=age_vals)

fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
ax.set_yscale("log")
ax.set_ylabel("$\\mathregular{N(1) (km^{-2})}$")
ax.set_xlabel("Time (Gy ago)")
ax.set_xlim(4.5, 0)

# Convert age from My to Gy and N1 from m² to km² for plotting
ax.plot(
    age_vals * 1e-3,
    N1_moon * 1e6,
    "-",
    color="dimgrey",
    linewidth=2.0,
    zorder=50,
    label="Moon",
)
ax.plot(
    age_vals * 1e-3,
    N1_mars * 1e6,
    "-",
    color="orange",
    linewidth=2.0,
    zorder=50,
    label="Mars",
)

ax.legend()

plt.show()
