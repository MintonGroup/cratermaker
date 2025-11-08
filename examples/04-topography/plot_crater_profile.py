"""
Create a crater and ejecta profile with the "simplemoon" morphology model
=========================================================================

.. rubric:: By David Minton

This example showcases how to create a crater and ejecta profile using the "simplemoon" morphology model from the CraterMaker package. The crater is created with a radius of 1 km, and the profiles are plotted in a 2D space normalized to the crater radius.

"""

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

fig, ax = plt.subplots(figsize=(8, 2))

ax.set_xlabel("Distance ($r/r_c$)")
ax.set_ylabel("Height ($h/r_c$)")

ax.set_aspect("equal")
ax.plot(rvals / rc, hvals / rc)
ax.plot(rej / rc, ejvals / rc)

plt.show()
