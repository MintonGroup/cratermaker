.. currentmodule:: cratermaker

.. _scaling:

#######
Scaling
#######

.. rubric:: by Dennise Valadez and David A. Minton

Cratermaker's ``scaling`` module provides tools to convert projectile parameters (e.g., projectile diameter, velocity, and target properties) into a final crater diameter. It also categorizes craters by their morphology: **simple**, **transitional**, or **complex**.

This example simulates how crater diameters vary with projectile size across different planetary bodies — specifically **Mars**, **Venus**, **Mercury**, and the **Moon** — and visualizes this relationship.

Example: Crater scaling for various planetary surfaces
=======================================================

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   from cratermaker import (
       Projectile, Target, scaling
   )

   def plot_projectile_to_crater_scaling():
       bodies = ["Mars", "Venus", "Mercury", "Moon"]
       proj_diameters = np.logspace(1, 4, 100)
       markersize = 5
       markerstyle = 'D'

       fig, ax = plt.subplots(figsize=(9, 7))
       ax.set_xscale("log")
       ax.set_yscale("log")
       ax.set_xlabel("Projectile Diameter [m]")
       ax.set_ylabel("Final Crater Diameter [m]")
       ax.set_title("Projectile to Crater Scaling for Mars, Venus, Mercury, and the Moon")

       colors = [
           "blue", "green", "red", "orange", "purple", "brown",
           "cyan", "magenta", "lime", "olive", "pink", "gray"
       ]

       color_index = 0
       label_color_map = {}

       for body in bodies:
           target = Target(body)
           scaling = scaling(target=target)
           crater_diams = []
           morphs = []

           for d in proj_diameters:
               proj = Projectile(
                   diameter=d,
                   mean_velocity=target.escape_velocity + 1000,
                   target=target
               )
               crater = scaling.projectile_to_crater(proj)
               crater_diams.append(crater.diameter)
               morphs.append(scaling.get_morphology_type(crater.diameter))

           for morph in np.unique(morphs):
               mask = np.array(morphs) == morph
               label = f"{body} – {morph.capitalize()}"

               if label not in label_color_map:
                   color = colors[color_index % len(colors)]
                   label_color_map[label] = color
                   color_index += 1

               ax.plot(
                   proj_diameters[mask],
                   np.array(crater_diams)[mask],
                   markerstyle,
                   markersize=markersize,
                   linestyle='None',
                   label=label,
                   color=label_color_map[label],
                   markeredgewidth=0.5
               )

       handles, labels = ax.get_legend_handles_labels()
       unique = dict(zip(labels, handles))
       ax.legend(unique.values(), unique.keys(), fontsize=9)

       ax.grid(True, which='both', linestyle='--', alpha=0.3)
       plt.tight_layout()
       plt.show()

   if __name__ == "__main__":
       plot_projectile_to_crater_scaling()

Explanation
-----------

- The function creates a range of projectile diameters from 10 to 10,000 meters.
- For each body, a ``Target`` and corresponding ``scaling`` object are created.
- Projectiles are launched at a velocity just above escape velocity.
- The resulting crater diameter and morphology are computed using ``scaling.projectile_to_crater`` and ``scaling.get_morphology_type``.
- Crater types are plotted using basic colors, categorized by body and morphology.

This visualization helps compare how gravity and crustal properties of each planetary body influence the final crater size and morphology classification.

.. image:: ../_static/scaling_example.png
   :alt: Scaling relationships for Mars, Venus, Mercury, and the Moon.
   :scale: 50 %
   :align: center
