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
    print(moon)

**Example 2: Viewing the Target catalogue**

.. ipython:: python
    :okwarning:

    from cratermaker import Target
    # Access the target catalogue options 
    catalogue = Target.maker().catalogue

    # Print available targets
    print(f"Available targets:\n{catalogue}")


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
    print(eris)
