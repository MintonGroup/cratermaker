"""
Depth vs diameter of fresh lunar craters
========================================

.. rubric:: By David Minton

This example shows the initial depth-to-diameter ratio of fresh lunar craters in Cratermaker, and compares them with measurements from the literature.

Sources
-------

- The model for the rim floor elevations as a function of crater diameter for large craters (diameters larger than ~5km), including their scatter, is taken from Pike (1977) [#]_
- For small craters (diameters smaller than ~5km), the depth-to-diameter model is taken as a weighted average between two models: 1) A modified version of the profile function given in Fassett and Thomson (2014) [#]_, which takes depth as a variable quantity. 2) The suite of profile functions given in Yang et al. (2021) [#]_. The weighting factor is set such that for craters larger than 5 km, all craters use the Fassett and Thomson profile, and all craters smaller than 50 m use one of the Yang et al. profiles (chosen at random for diameters less than 180 m, otherwise the "normal" version is used). In between, the weighting of the FT14 to Y21 profile models (and their depth) is chosen from a weighted gaussian function that biases the profiles to the shallower Y21 model for small craters and to the deeper FT14 model for large craters within the range.
- Measurements of depth-to-diameter values for fresh craters are drawn from Pike (1974) [#]_, Stopar et al. (2017) [#]_ (only the Class A craters), and Hoover et al. (2024) [#]_.


References
----------

.. [#] Pike, R.J., 1977. Size-dependence in the shape of fresh impact craters on the moon. Presented at the In: Impact and explosion cratering: Planetary and terrestrial implications; Proceedings of the Symposium on Planetary Cratering Mechanics, pp. 489-509.
.. [#] Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. `doi:10.1002/2014JE004698 <https://doi.org/10.1002/2014JE004698>`_
.. [#] Yang, X., Fa, W., Du, J., Xie, M., Liu, T., 2021. Effect of Topographic Degradation on Small Lunar Craters: Implications for Regolith Thickness Estimation. Geophysical Research Letters 48, e2021GL095537. `doi:10.1029/2021GL095537 <https://doi.org/10.1029/2021GL095537>`_
.. [#] Pike, R.J., 1974. Depth/diameter relations of fresh lunar craters: Revision from spacecraft data. Geophysical Research Letters 1, 291-294. `doi:10.1029/GL001i007p00291 <https://doi.org/10.1029/GL001i007p00291>`_
.. [#] Stopar, J.D., Robinson, M.S., Barnouin, O.S., McEwen, A.S., Speyerer, E.J., Henriksen, M.R., Sutton, S.S., 2017. Relative depths of simple craters and the nature of the lunar regolith. Icarus. `doi:10.1016/j.icarus.2017.05.022 <https://doi.org/10.1016/j.icarus.2017.05.022>`_
.. [#] Hoover, R.H., Robbins, S.J., Hynek, B.M., Hayne, P.O., 2024. Depth-to-diameter Ratios of Fresh Craters on the Moon and Implications for Surface Age Estimates. Planet. Sci. J. 5, 26. `doi:10.3847/PSJ/ad18d4 <https://doi.org/10.3847/PSJ/ad18d4>`_



Depth vs diameter measurement files
-----------------------------------

  :download:`Hoover2024Table2.csv </_static/Hoover2024Table2.csv>`
  :download:`Pike1974-Mare.csv </_static/Pike1974-Mare.csv>`
  :download:`Pike1974-Highlands.csv </_static/Pike1974-Highlands.csv>`
  :download:`S2-Stopar et al. Simple Craters.csv </_static/S2-Stopar_et_al._Simple_Craters.csv>`


"""

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from cratermaker import Simulation

simdir = "simdata-4.5"

xmin, xmax = 0.01, 1000.0
ymin, ymax = 0.0, 0.30


# Compute depth-to-diameter for 1000 craters in Cratermaker across the size range of interest
sim = Simulation(simdir=simdir, gridlevel=0, ask_overwrite=False)
cm_diameter = np.exp(sim.rng.uniform(low=np.log(xmin), high=np.log(xmax), size=1000))
cm_dD = []
for diameter in cm_diameter * 1e3:
    crater = sim.Crater.maker(diameter=diameter)
    cm_dD.append(crater.depth_to_diameter)


# Read in observational data
df_pike_mare = pl.read_csv("Pike1974-Mare.csv")
df_pike_highlands = pl.read_csv("Pike1974-Highlands.csv")
df_hoover = pl.read_csv("Hoover2024Table2.csv")
df_stopar = pl.read_csv("S2-Stopar_et_al._Simple_Craters.csv")
diam_stopar = df_stopar["Diameter, crater (m)"].to_numpy() * 1e-3
dD_stopar = df_stopar["depth/ Diameter"].to_numpy()
class_stopar = df_stopar["Class (FINAL)"].to_numpy()

# Plot depth/diameter vs diameter
fig, ax = plt.subplots(figsize=(15, 8), layout="constrained")

ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)
ax.set_xscale("log")

ax.set_xlabel("Diameter (km)")
ax.set_ylabel("Depth/Diameter")

ax.scatter(
    x=df_pike_mare["Diameter (km)"],
    y=df_pike_mare["Depth (km)"] / df_pike_mare["Diameter (km)"],
    label="Pike (1974) Mare",
    zorder=10,
)
ax.scatter(
    x=df_pike_highlands["Diameter (km)"],
    y=df_pike_highlands["Depth (km)"] / df_pike_highlands["Diameter (km)"],
    label="Pike (1974) Highlands",
    zorder=10,
)

ax.errorbar(
    x=df_hoover["D(m)"] * 1e-3,
    y=df_hoover["d/D"],
    yerr=df_hoover["err"],
    fmt="s",
    label="Hoover et al. (2024)",
    zorder=10,
)

ax.scatter(
    x=diam_stopar[class_stopar == "A"],
    y=dD_stopar[class_stopar == "A"],
    label="Stopar et al. (2017) Class A",
    zorder=10,
)

ax.scatter(
    x=cm_diameter,
    y=cm_dD,
    label="Cratermaker",
    marker="o",
    color="black",
    s=5,
    zorder=1,
)
ax.legend(loc="upper right")
plt.show()


# Plot depth vs diameter
fig, ax = plt.subplots(figsize=(10, 8), layout="constrained")
ymin, ymax = 0.001, 10.0
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)
ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Diameter (km)")
ax.set_ylabel("Depth (km)")

ax.scatter(
    x=df_pike_mare["Diameter (km)"],
    y=df_pike_mare["Depth (km)"],
    label="Pike 1974 Mare",
    zorder=10,
)
ax.scatter(
    x=df_pike_highlands["Diameter (km)"],
    y=df_pike_highlands["Depth (km)"],
    label="Pike 1974 Highlands",
    zorder=10,
)
ax.errorbar(
    x=df_hoover["D(m)"] * 1e-3,
    y=df_hoover["d/D"] * df_hoover["D(m)"] * 1e-3,
    yerr=df_hoover["D(m)"] * df_hoover["err"] * 1e-3,
    fmt="s",
    label="Hoover et al. (2024)",
    zorder=10,
)
ax.scatter(
    x=diam_stopar[class_stopar == "A"],
    y=dD_stopar[class_stopar == "A"] * diam_stopar[class_stopar == "A"],
    label="Stopar et al. (2017) Class A",
    zorder=10,
)

ax.scatter(
    x=cm_diameter,
    y=cm_dD * cm_diameter,
    label="Cratermaker",
    marker="o",
    color="black",
    s=5,
    zorder=1000,
)

ax.legend(loc="lower right")
plt.show()
