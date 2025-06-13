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

- `final_diameter`: The final rim-to-rim diameter (post collapse phase) of the crater in meters.
- `final_radius`: The final rim-to-rim radius (post collapse phase) of the crater in meters.
- `transient_diameter`: The transient diameter (pre-collapse phase) of the crater in meters.  
- `transient_radius`: The transient radius (pre-collapse phase) of the crater in meters.  
- `projectile_diameter`: The diameter of the projectile that formed the crater in meters.
- `projectile_radius`: The radius of the projectile that formed the crater in meters.
- `projectile_mass`: The mass of the projectile that formed the crater in kilograms.

All other parameters are optional, and if not specified, will be determined by the provided (or default) Target, Scaling, and Projectile objects. We can demonstrate this behavior by specifiying the minimum set of arguments for a Crater, in which we define a crater with a final rim-to-rim diameter of 100 meters and then printing the resulting crater object to show its full set of properties:

.. ipython:: python

    from cratermaker import Crater

    crater = Crater.maker(final_diameter=100.0)
    print(crater)

As we can see, the Crater object contains values for the transient crater diameter, projectile properties, including size, mass, and impact velocity, angle and direction. It also contains a set of location coordinates on the target body, and a string value of the morphology type indicating that this is a simple crater. 

Because we didn't specify them, the :ref:`default values <ug-defaults>` of :ref:`Target <ug-target>`, :ref:`Scaling <ug-scaling>`, and :ref:`Projectile <ug-projectile>` were used. Therefore the above is quivalent to:

.. ipython:: python

    crater = Crater.maker(final_diameter=100.0, 
                          target="Moon", 
                          projectile="asteroids", 
                          scaling="montecarlo")
    print(crater)

Crater from Projectile
----------------------

You can also create a crater by specifying the projectile size, which will then compute the crater size based on the projectile properties and the target body. For example, if we want to create a crater on Europa with a projectile diameter of 1.5 km , we can do so as follows:

.. ipython:: python

    crater = Crater.maker(target="Europa", projectile_diameter=1500)
    print(crater)

Because we chose Europa as a target, the default projectile population is "comets" instead of "asteroids."


Specifying Crater Properties 
--------------------------------

Besides the crater or projectile sizes values, which are required, there are a number of optional arguments you can pass to :meth:`Crater.maker()`, including

- `projectile_density`: The density of the projectile in kg/m\ :sup:`3`. If not provided, it will be defined through the Projectile component provided. (e.g., "asteroids" or "comets").
- `projectile_mean_velocity`: The mean velocity in m/s from which to sample a projectile velocity. This will override the value from the Projectile model.
- `projectile_velocity`: The total impact velocity of the projectile in m/s. If not provided, it is drawn from a distribution based on the target body and projectile populations.
- `projectile_vertical_velocity`: The vertical component of the velocity in m/s. The total impact velocity will be determined by the projectile angle (either provided or drawn from a distribution).
- `projectile_angle`: The impact angle in degrees (0-90). If not provided, it is drawn from a distribution that peaks at 45 degrees.
- `projectile_direction`: The direction of the impact in degrees (0-360) relative to true north. If not provided, it is drawn from a uniform distribution.
- `location`: The (longitude, latitude) location of the crater in degrees.
- `projectile_location`: The (longitude, latitude) location of the projectile impact. This is equivalent to `location`, which takes precedence
- `age`: The age of the crater in My. If not provided, it is not set. This is used by the :ref:`Simulation <ug-simulation>` for tracking the age at which a crater formed during the simulation. Setting it has no effect on the crater properties.


Suppose we wish to create a crater with a 10km transient diameter on Europa, but want to specify that it formed from an impactor with a density of 900 kg/m\ :sup:`3`, an impact velocity of 30 km/s, an impact angle of 25°, and a location at 45° S, 150° E. We can do so as follows:

.. ipython:: python

    crater = Crater.maker(
        target="Europa",
        transient_diameter=10e3,
        projectile_density=900,        
        projectile_velocity=30e3,
        projectile_angle=25.0,
        location=(150, -45),
    )
    print(crater)

We can specify the projectile angle from a combination of the vertical velocity and total impact velocity. In this example we can create an impact on Mercury using a 3 km projectile that indirectly specifies the projectile angle from a 35 km/s velocity and a vertical velocity of 20 km/s. 

.. ipython:: python

    crater = Crater.maker(
        target="Mercury",
        projectile_diameter=3000.0,
        projectile_velocity=35.0e3,
        projectile_vertical_velocity=20e3,
    )
    print(crater)


Crater on a Custom Target
-------------------------
The following example uses a custom-defined :ref:`Target <ug-target>` object. Here we will create a 300 km diameter body in the Kuiper belt with a mass of 10 :sup:`19` kg (about 700 kg/m\ :sup:`3` bulk density). It will use both the ice model for its surface material and simple-to-complex scaling model. 

.. ipython:: python

    from cratermaker import Crater, Target

    kuip_belt_obj = Target(
        name="KBO-1",
        diameter=300e3,                
        mass=1e19,                     
        material="Ice",
        transition_scale_type="ice"
    )
    print(kuip_belt_obj)

Next we will feed our custom target body to the :meth:`Crater.maker()` method to create a crater on it. We will use a rocky projectile with a diameter of 500 m, an impact velocity of 5 km/s. All other properties will be computed.

.. ipython:: python

    crater = Crater.maker(
        target=kuip_belt_obj,
        projectile_density=2800,       
        projectile_diameter=500,       
        projectile_velocity=5000,     
    )
    print(crater)
