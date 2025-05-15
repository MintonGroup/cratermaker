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

.. ipython:: python
    from cratermaker import Crater
    #By default the crater method will default to the Moon 
    crater = Crater.maker(
        projectile_diameter=100,
        projectile_density=3000,
        projectile_velocity=20000,
        projectile_angle=45,
        location=(0, 0),
    )

    print(f"Crater final diameter: {crater.final_diameter:.2f} m")


Examples
-------------------

.. ipython:: python

    # Example 1 — Asteroid impact on the Moon
    from cratermaker import Crater

    crater = Crater.maker(
        target="Moon",                 # Built-in target
        projectile="asteroids",        # Built-in projectile model
        projectile_diameter=1500,      # Diameter in meters
        projectile_velocity=20000,     # Velocity in m/s
        projectile_angle=45,           # Angle in degrees
        location=(0, 0)                
    )

    print(crater)


.. ipython:: python

    # Example 2 — Icy comet impact on Europa
    from cratermaker import Crater

    crater = Crater.maker(
        target="Europa",
        projectile="comets",           
        projectile_density=900,        # for icy projectile (kg/m^3)
        projectile_diameter=2000,
        projectile_velocity=18000,
        projectile_angle=25,
        location=(150, -45),
    )

    print(crater)
