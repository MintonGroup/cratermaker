.. currentmodule:: cratermaker

.. image:: ../_static/simulation_icon.svg
    :alt: Simulation
    :align: center
    :width: 300px

.. _ug-simulation:

Simulation
==========

The :class:`Simulation` class manages the process of evolving cratered planetary surfaces through time. It supports configuring a target body, choosing a production function, and generating craters through emplacement or time evolution.

To run a simulation, use the :class:`cratermaker.Simulation` constructor and call :meth:`Simulation.run` or :meth:`Simulation.emplace`.

.. ipython:: python
    :okwarning:

    from cratermaker import Simulation
    sim = Simulation(target="Moon")
    sim.run(age=100)

    print(f"Craters emplaced: {len(sim.true_crater_list)}")

.. _ug-simulation-examples:

Examples
--------

Below are examples of how to use the Simulation class.

**Example 1: Run a short simulation on Mars**

.. ipython:: python
    :okwarning:

    from cratermaker import Simulation
    sim = Simulation(target="Mars")
    sim.run(age=500)
    print(f"Crater count: {len(sim.true_crater_list)}")

**Example 2: Simulate default Moon configuration**

.. ipython:: python
    :okwarning:

    from cratermaker import Simulation
    sim = Simulation()
    sim.run(age=1000)
    print(f"Crater count: {len(sim.true_crater_list)}")

**Example 3: Emplacing a single crater on the Moon**

You can manually emplace a crater with a specified size and location using the :meth:`Simulation.emplace` method.

.. ipython:: python
    :okwarning:

    from cratermaker import Simulation

    # Create a simulation on the Moon
    sim = Simulation(target="Moon")

    # Emplace a single crater with 20 km diameter 
    sim.emplace(final_diameter=20000, location=(0, 0))

    # Inspect crater
    crater = sim.true_crater_list[0]
    print(f"Crater diameter: {crater.final_diameter:.1f} m")
    print(f"Location (lon, lat): {crater.location}")
