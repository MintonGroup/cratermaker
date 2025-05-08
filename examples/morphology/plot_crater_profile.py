"""
Create a crater and ejecta profile with the "simplemoon" morphology model
=========================================================================
This example showcases how to create a crater and ejecta profile using the "simplemoon" morphology model from the CraterMaker package. The crater is created with a radius of 1 km, and the profiles are plotted in a 2D space normalized to the crater radius.

.. rubric:: By David Minton

"""

import numpy as np
import matplotlib.pyplot as plt
from cratermaker import Crater, Morphology

morphology = Morphology.maker(crater=Crater.maker(final_radius=1.0e3))

rc = morphology.crater.final_radius
rvals = np.linspace(0, 3.0 * rc, 1000)
hvals = morphology.crater_profile(rvals)
ejvals = morphology.ejecta_profile(rvals)

# Filter out the interior of the crater and add the ejecta to the crater profile to draw the combined profile
rej = rvals[ejvals > 0]
ejvals = hvals[ejvals > 0] + ejvals[ejvals > 0]

fig, ax = plt.subplots(figsize=(8, 2))

ax.set_xlabel('Distance ($r/r_c$)')
ax.set_ylabel('Height ($h/r_c$)')

ax.set_aspect('equal')
ax.plot(rvals/rc, hvals/rc)
ax.plot(rej/rc, ejvals/rc)

plt.show()