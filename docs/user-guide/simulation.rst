.. currentmodule:: cratermaker

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()

.. image:: ../_static/full_simulation.png
    :alt: Simulation
    :align: center
    :width: 300px
    :class: dark-light

.. _ug-simulation:

Simulation
==========

The :class:`~cratermaker.core.simulation.Simulation` class manages the evolution of cratered planetary surfaces. It contains components such as target bodies, impactor distributions, and emplacement to simulate crater formation over time.

A `Simulation` instance can be initiated with keyword arguments specifying the various component models to be used in the simulation. These components include:

- :ref:`Target <ug-target>`
- :ref:`Scaling <ug-scaling>`
- :ref:`Production <ug-production>`
- :ref:`Morphology <ug-morphology>`
- :ref:`Surface <ug-surface>`
- :ref:`Projectile <ug-projectile>`
- :ref:`Counting <ug-counting>`

Any valid arguments for these components can be passed directly to the :class:`~cratermaker.core.simulation.Simulation` constructor, otherwise default values will be used (see the guide for :ref:`Default Behavior <ug-defaults>`).


Creating craters and populations of craters
-------------------------------------------

The :class:`~cratermaker.core.simulation.Simulation` class has set of core methods that can be used to generate craters, using the requested (or default) Production and Morphology components onto the requested Surface component. These are :meth:`~cratermaker.core.simulation.Simulation.emplace`, :meth:`~cratermaker.core.simulation.Simulation.populate`, and :meth:`~cratermaker.core.simulation.Simulation.run`. Each successive method makes use of the others (that is, `populate` calls `emplace` on the craters in generates, and `run` calls `populate` for each interval of time in the simulation). We will discuss each of these in turn:

Emplace one or more specific craters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~cratermaker.core.simulation.Simulation.emplace` is used to emplace one or more craters onto a surface. The arguments are either the parameters that you would pass to :meth:`~cratermaker.components.crater.Crater.maker` to generate a single crater, a single :class:`~cratermaker.components.crater.Crater` object, or a list of :class:`~cratermaker.components.crater.Crater` objects. So for instance, you can emplace a single 300 km crater onto a random location on a surface, you would do:

.. code-block:: python

    from cratermaker import Simulation

    sim = Simulation()
    sim.emplace(diameter=300e3)

To place 3 craters onto the surface, you could do the following:

.. code-block:: python

    from cratermaker import Crater

    craters = [Crater.maker(diameter=d) for d in [100e3, 200e3, 300e3]]
    sim.emplace(craters)

Emplace a population of craters drawn from a production function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next level of methods is :meth:`~cratermaker.core.simulation.Simulation.populate`, which is used to generate a population of craters drawn from whichever :ref:`Production <ug-production>` component module is loaded into the simulation. The arguments for this method are the same as those for :meth:`~cratermaker.components.production.Production.populate`, which can either be the time interval over which to draw from the production function, or a pair of diameter/number values that specify the position in the production function to draw from. For instance, to emplace a population of craters corresponding to the production function between 200 My to 100 My before present, you would do:

.. code-block:: python

    sim.populate(time_start=200, time_end=100)


The :meth:`~cratermaker.core.simulation.Simulation.populate` method is really meant to be a helper method for the :meth:`~cratermaker.core.simulation.Simulation.run` method, which is used to run a full simulation over a specified time interval. For instance, the 
The :meth:`~cratermaker.core.simulation.Simulation.populate` method does not track the time of emplacement or advance the time interval of the simulation.



Running a Simulation
--------------------

Once a :class:`~cratermaker.core.simulation.Simulation` instance has been created, the simulation can be run using the :meth:`~cratermaker.core.simulation.Simulation.run` method. This method requires an `age` argument specifying the duration of the simulation in million years (My).  

The following example configures a simulation targeting the Moon and runs it for 4 Gy, then we plot the number of true emplaced craters as well as the number of observed craters that the crater counting algorithm has determined are observable. To speed up the example, we have reduced the resolution of the surface grid from the default value of 8 to 6 using the ``gridlevel`` argument.


.. ipython:: python
    :okwarning:
    :suppress:

    cleanup()

.. ipython:: python
    :okwarning:

    from cratermaker import Simulation
    sim = Simulation(target="Moon", gridlevel=6)
    sim.run(age=4000)
    print(f"Number of true emplaced craters: {len(sim.emplaced)}")
    print(f"Number of observed craters: {len(sim.observed)}")

More detailed component examples are provided in the Gallery section.

.. seealso::

    - :ref:`api-simulation` for the API reference
    - :ref:`gal-simulation` for example simulations, including visualizations

.. ipython:: python
    :okwarning:
    :suppress:

    from cratermaker import cleanup
    cleanup()