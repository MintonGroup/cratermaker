.. currentmodule:: cratermaker

.. image:: ../_static/full_simulation.png
    :alt: Simulation
    :align: center
    :width: 300px

.. _ug-simulation:

Simulation
==========

The :class:`Simulation` class manages the evolution of cratered planetary surfaces. It contains components such as target bodies, impactor distributions, and emplacement to simulate crater formation over time.

A `Simulation` instance can be initiated with keyword arguments including:

- :class:`Target`
- Production function
- Scaling laws
- Ejecta and saturation models

Any valid arguments for these components can be passed directly to the :class:`Simulation` constructor.

**Basic usage:**

The following example configures a simulation targeting the Moon and runs it for 100 Myr. Then it accesses the list of emplaced craters.

.. code-block:: python

    from cratermaker import Simulation
    sim = Simulation(target="Moon")
    sim.run(age=100)
    craters = sim.true_crater_list

**Available parameters include** :

- ``target`` (str or Target): Name of target body (e.g., "Moon", "Mars"), or a `Target` instance.
- ``radius``, ``diameter``, ``mass``, ``density``, ``material``, ``transition_scale_type``: Parameters forwarded to `Target`.
- ``production`` (str): Choose the crater production function.
- ``scaling`` (str): Select the impact scaling law.
- ``morphology``, ``surface``, ``projectile`` (str): Other model component selectors.
- ``rng``, ``rng_seed``, ``rng_state``: Random number generator configuration.
- ``resume_old`` (bool): Resume a previously saved simulation.

**Emplacement-specific arguments** (for use with :meth:`Simulation.emplace`):

- ``final_diameter`` (float): Diameter of the crater in meters.
- ``location`` (tuple): Longitude and latitude in degrees where the crater should be emplaced.
- Additional projectile/crater attributes may also be passed (e.g., ``projectile_mass``, ``velocity``).

**Common methods:**

- :meth:`Simulation.run`: Evolve a surface over a specified time or crater production interval.
- :meth:`Simulation.emplace`: Emplace one or more specific craters manually.

More detailed component examples are provided in the Gallery section.

.. seealso::

    - :ref:`api-simulation` for the API reference
    - :ref:`gal-simulation` for example simulations, including visualizations
