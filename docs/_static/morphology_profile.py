import matplotlib.pyplot as plt
import numpy as np

from cratermaker import Crater, Morphology

crater = Crater.maker(final_radius=1.0e3)
morphology = Morphology.maker("simplemoon")

rc = crater.final_radius
rvals = np.linspace(0, 3.0 * rc, 1000)
hvals = morphology.crater_profile(crater, rvals)
ejvals = morphology.ejecta_profile(crater, rvals)

# Filter out the interior of the crater and add the ejecta to the crater profile to draw the combined profile
rej = rvals[ejvals > 0]
ejvals = hvals[ejvals > 0] + ejvals[ejvals > 0]

rnorm = rvals / rc
size_ratio = (hvals.max() - hvals.min()) / (2 * rvals.max())

w = 200
h = w * size_ratio

fig, ax = plt.subplots(figsize=(w, h), layout="constrained")


ax.set_aspect("equal")
linewidth = 16
ax.plot(-rnorm[rnorm <= 1], hvals[rnorm <= 1] / rc, color="black", linewidth=linewidth)
ax.plot(-rnorm[rnorm >= 1], ejvals / rc, color="black", linewidth=linewidth)
ax.plot(rnorm[rnorm <= 1], hvals[rnorm <= 1] / rc, color="black", linewidth=linewidth)
ax.plot(rnorm[rnorm >= 1], ejvals / rc, color="black", linewidth=linewidth)
ax.axis("off")
plt.savefig("morphology_profile.svg", transparent=True)
