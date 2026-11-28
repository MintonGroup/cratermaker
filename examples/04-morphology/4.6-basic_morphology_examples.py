"""
Crater morphology across sizes
==============================

.. rubric:: By David Minton

This example shows fresh craters made with the "basicmoon" morphology model across the range of lunar crater sizes: a simple crater, a crater at the simple-to-complex transition, a complex crater, and a multiring basin.

Each crater is emplaced at the center of its own simulation using the "hireslocal" surface, which has a high resolution local region surrounded by a coarser superdomain covering the rest of the Moon. The local region is scaled to the size of each crater, with a radius equal to the crater diameter, so each crater and its continuous ejecta blanket fill the view at the same relative resolution. The images are generated with the simulation's own ``plot`` method, which draws a hillshade of the local region.
"""

import matplotlib.pyplot as plt

from cratermaker import Scaling, Simulation

# Number of pixels across the radius of the local region. This sets the resolution of each crater relative to its size.
pixels_per_radius = 200

# Each crater gets its own simulation directory. The transitional crater diameter is the nominal simple-to-complex transition
# diameter of the lunar scaling model.
examples = {
    "Simple": 4e3,
    "Transitional": Scaling.maker(target="Moon").transition_nominal,
    "Complex": 100e3,
    "Multiring basin": 700e3,
}

fig, axes = plt.subplots(2, 2, figsize=(10, 10), layout="constrained")
for ax, (name, diameter) in zip(axes.flat, examples.items(), strict=True):
    sim = Simulation(
        surface="hireslocal",
        morphology="basicmoon",
        simdir=f"simdata-4.6/{name.lower().replace(' ', '_')}",
        local_radius=diameter,
        pix=diameter / pixels_per_radius,
        rng_seed=4600,
        reset=True,
        ask_overwrite=False,
    )
    if name == "Transitional":
        sim.emplace(diameter=diameter, location=(0, 0), morphology_type="transitional")
    else:
        sim.emplace(diameter=diameter, location=(0, 0))
    sim.plot(label=f"{name} ({diameter / 1e3:.0f} km)", ax=ax, show=False, save=False, close_when_done=False)

plt.show()
