.. currentmodule:: cratermaker

.. image:: ../_static/crater_icon.svg
    :alt: Production
    :align: center
    :width: 300px

.. _ug-crater:

Crater
======

The :ref:`Crater <api-crater>` class is one of the core components of Cratermaker. It is a dataclass that represents the properties of a crater and the projectile that formed it. Like all other components, it is instantiated with its `maker` factory method. You can create a crater by either specifying either its size or the size of a projectile. In either case, the crater and projectile properties will be computed using :ref:`Target <ug-scaling>`, :ref:`Scaling <ug-scaling>`, and :ref:`Projectile <ug-projectile>` objects. These can be provided to :ref:`Crater.maker() <api-crater>`, but if not, then the :ref:`default <ug-defaults>` values will be used. 

At a minimum, you need to specify exactly one size input, which can be one of the following:

- ``final_diameter``: The final rim-to-rim diameter (post collapse phase) of the crater in meters.
- ``final_radius``: The final rim-to-rim radius (post collapse phase) of the crater in meters.
- ``transient_diameter``: The transient diameter (pre-collapse phase) of the crater in meters.  
- ``transient_radius``: The transient radius (pre-collapse phase) of the crater in meters.  
- ``projectile_diameter``: The diameter of the projectile that formed the crater in meters.
- ``projectile_radius``: The radius of the projectile that formed the crater in meters.
- ``projectile_mass``: The mass of the projectile that formed the crater in kilograms.

All other parameters are optional, and if not specified, will be determined by the provided (or default) Target, Scaling, and Projectile objects. We can demonstrate this behavior by specifiying the minimum set of arguments for a Crater, in which we define a crater with a final rim-to-rim diameter of 100 meters and then printing the resulting crater object to show its full set of properties:

.. ipython:: python

    from cratermaker import Crater

    crater = Crater.maker(final_diameter=100)
    print(crater)

As we can see, the Crater object contains values for the transient crater diameter, projectile properties, including size, mass, and impact velocity, angle and direction. It also contains a set of location coordinates on the target body, and a string value of the morphology type indicating that this is a simple crater. 

Additional Examples
-------------------

.. ipython:: python

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
