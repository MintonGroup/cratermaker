.. currentmodule:: cratermaker

Scaling
=======

Cratermaker's ``scaling`` module provides tools to convert projectile parameters (e.g., projectile diameter, velocity, and target properties) into a final crater diameter. It also categorizes craters by their morphology: **simple**, **transitional**, or **complex**.

This example simulates how crater diameters vary with projectile size across different planetary bodies — specifically **Mars**, **Venus**, **Mercury**, and the **Moon** — and visualizes this relationship.


Explanation
-----------

- The function creates a range of projectile diameters from 10 to 10,000 meters.
- For each body, a ``Target`` and corresponding ``scaling`` object are created.
- Projectiles are launched at a velocity just above escape velocity.
- The resulting crater diameter and morphology are computed using ``scaling.projectile_to_crater`` and ``scaling.get_morphology_type``.
- Crater types are plotted using basic colors, categorized by body and morphology.

More Scaling examples
---------------------

See more examples at  :ref:`gal-scaling`

