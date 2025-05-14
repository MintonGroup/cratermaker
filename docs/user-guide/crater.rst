.. currentmodule:: cratermaker

.. image:: ../_static/crater_icon.svg
    :alt: Production
    :align: center
    :width: 300px

.. _ug-crater:

Crater
======

The ``Crater`` method is used to represent craters and the relationship properties of the projectile.

The following is a simple example that creates a crater on the Moon from a 100 meter diameter projectile:

.. code-block:: python

    from cratermaker.components.target import Target
    from cratermaker.core.crater import Crater

    target = Target.maker("Moon")

    crater = Crater.maker(
        target=target,
        projectile_diameter=100,
        projectile_density=3000,
        projectile_velocity=20000,
        projectile_angle=45,
        location=(0, 0),
    )

    print(f"Crater final diameter: {crater.final_diameter:.2f} m")
