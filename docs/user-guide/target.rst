.. currentmodule:: cratermaker

.. image:: ../_static/target_icon.svg
    :alt: Target
    :align: center
    :width: 300px


.. _ug-target:

Target
======

The Target class contains basic information about a celestial body on which craters are emplaced. It contains information about the body's size, material properties, and surface gravity. To create a standalone target body, you use the :func:`cratermaker.Target.maker` method. 

.. ipython:: python
    :okwarning:

    from cratermaker import Target
    target = Target.maker("Mars")
    print(target)

.. _ug-target-examples:

Examples
--------

Below are examples of how to use the Target class to explore planetary properties.

**Example 1: Defining the Moon as a target using the built-in catalogue**

.. ipython:: python
    :okwarning:

    from cratermaker import Target
    moon = Target.maker("Moon")

    # gravity and escape velocity
    gravity = moon.gravity              # in m/s²
    escape_velocity = moon.escape_velocity  # in m/s

    print(f"{moon.name}:")
    print(f"  Gravity = {gravity:.2f} m/s²")
    print(f"  Escape velocity = {escape_velocity:.2f} m/s")

**Example 2: Viewing the Target catalogue**

.. ipython:: python
    :okwarning:

    from cratermaker import Target
    # the classmethod .maker() is needed to get a temporary Target
    catalogue = Target.maker("Moon").catalogue

    # Print all available targets in catalogue
    print("Available targets:")
    for name in catalogue:
        print("-", name)

**Example 3: Defining a body not in the catalogue (Eris)**

.. ipython:: python
    :okwarning:

    from cratermaker import Target

    # Define Eris 
    eris = Target(
        name="Eris",
        diameter=2326e3,   # in meters
        mass=1.66e22,      # in kg
        material="Ice",
        transition_scale_type="ice"
    )

    print(f"{eris.name}: Gravity = {eris.gravity:.2f} m/s², Escape velocity = {eris.escape_velocity:.2f} m/s")
