.. currentmodule:: cratermaker

.. image:: ../_static/crater_icon.svg
    :alt: Production
    :align: center
    :width: 300px

.. _ug-crater:

Crater
======

The ``Crater`` method is used to represent craters and the relationship properties of the projectile.

To start, the following is a simple example that creates a crater on the Moon from a 100 meter diameter projectile:

.. ipython:: python

    from cratermaker import Crater

    # By default the crater method will default to the Moon 
   
    crater = Crater.maker(
        projectile_diameter=100,
        projectile_density=3000,
        projectile_velocity=20000,
        projectile_angle=45,
        location=(0, 0),
    )

    print(f"Crater final diameter: {crater.final_diameter:.2f} m")


Additional Examples
-------------------

.. ipython:: python

    # Asteroid impact on the Moon
    
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

    # Icy comet impact on Europa
    
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


.. ipython:: python

    # High-speed asteroid impact on Mercury
    
    from cratermaker import Crater

    crater = Crater.maker(
        target="Mercury",
        projectile="asteroids",
        projectile_diameter=3000,      # 3 km diameter
        projectile_velocity=35000,     # high-speed
        projectile_angle=55,
        location=(30, 10),
    )

    print(crater)


.. note::

    The following example uses a custom-defined ``Target`` object. 
    To learn how to define targets manually, see the :ref:`ug-target` documentation.

.. ipython:: python

    # Custom icy target (Kuiper Belt Object)
    
    from cratermaker import Crater, Target

    kuip_belt_obj = Target(
        name="KBO-1",
        diameter=300e3,                # 300 km
        mass=1e19,                     # Approximate small icy body
        material="Ice",
        transition_scale_type="ice"
    )

    crater = Crater.maker(
        target=kuip_belt_obj,
        projectile_density=2800,       # rocky projectile
        projectile_diameter=500,       # meters
        projectile_velocity=15000,     # meters/second
        projectile_angle=20,           # small angle
        location=(75, -10),
    )

    print(crater)
