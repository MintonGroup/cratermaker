.. currentmodule:: cratermaker

################################
Getting Started with Cratermaker
################################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs. It is built to be highly modular, and so many of the core components of Cratermaker can be run independently, outside of its core landscape evolution simulation. Each component of Cratermaker is built around a set of classes, which include:

- **Simulation**: This is main class that can be used to run a comprehensive cratered landscape evolution model.
- **Scaling**: Contains the model for relating projectile size and crater size.
- **Projectile**: Models the properties of projectiles, including their size, velocity, and density. 
- **Production**: Manages the production function, which determines the number of craters/projectiles produced over time and surface area.
- **Target**: Models the properties of the target body, including its size, material properties, and surface gravity.
- **Surface**: Represents the surface of the target body, including its topography and material properties.
- **Grid**: Handles the construction of the unstructured grid used for the modeling the surface, including its resolution and geometry.
- **Morphology**: Used to construct the three-dimensional representation of a crater onto the surface.
- **Crater**: Represents a single crater and the properties of the projectile that formed it.

While some components require others to function (for instance, a Scaling model needs both a Target and Projectile component in order to perform the scaling computations), such components can be invoked by themselves, and any missing components are built internally using a resonable default. The class definition of each component contains a special construction method called ``.maker()``, which is used to create the component. This method can be used to create a component without having to create the entire simulation.

In thhis guide, we will walk through the basic process of installing Cratermaker and running through some basic examples to demonstrate its functionality.

Installation
============
To begin, simply install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: bash 

      pip install cratermaker


Simulating a single crater
==========================

Here’s a simple example demonstrating how to initialize a simulation for the Moon, then emplace a single 100 km diameter crater at a random location on the surface:

.. code-block:: python

    import cratermaker
    sim = cratermaker.Simulation()
    sim.emplace_crater(final_diameter=100e3)


This example uses the default Moon target and emplaces a crater of 100 km in diameter. The ``emplace_crater`` method handles the complexities of crater formation, including the selection of a random location and the necessary calculations to model the impact.

.. note::

      The first time you run a simulation in a given directory, Cratermaker generates a configuration file called ``cratermaker.yaml`` and a surface mesh in ``surface_data\grid.nc``. By default, these are stored in the current working directory from where it was invoked, but this can be changed by passing the argument ``simdir="your\custom\path"`` to ``Simulation``.  


Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
=====================================================================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. The following will export the 
surface mesh to two files inside a folder called ``vtk_files``::

      sim.export_vtk()

This will generate a series of files with the patern ``surfXXXXXX.vtp``, which represent the grid and its associated data at each output time interval of the simulation. If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView. 

Then, open ParaView and select ``File -> Open``, navigate to the ``vtk_files`` directory, and select the file group. Clicking ``Apply`` will load the file. You can then select ``File -> Open`` again and load the other file. Under ``Coloring`` select ``face_elevation`` or  ``node_elevation`` selected to visualize the surface elevation. 

 The following image shows the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_100km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location.
      :align: center

Simulating a population of craters
==================================

Simulating a single crater is useful for testing, but Cratermaker is designed to simulate populations of craters over time. The following example demonstrates how to initialize a simulation of the Moon and emplace a population of craters using the default production function, which is the Neukum production function. The simulation will run for 1 billion years.

.. code-block:: python

    import cratermaker
    sim = cratermaker.Simulation()
    sim.run(age=1000)

.. note::

      The default units for Cratermaker are meters for length and million years for time.


Default Behavior
================
As Cratermaker is designed to be easy to use, all of its component classes are built to be invoked without any arguments. It is important to understand what the default behavior of these classes is, as this will affect the results of your simulation.

- **Scaling**: The default scaling model is called ``richardson2009``, as it is similar to the one used by the Cratered Terrain Evolution Model (CTEM) that is a progenitor to Cratermaker (see [1]_). The projectile to transient scaling model is mostly based on Holsapple (1993) ([2]_) with some additional scaling parameters for ice given by Kraus et al. (2011) ([3]_).
- **Projectile**: Reference Minton et al. (2010) and Zahnle et al. (2002).
- **Production**: Neukum et al. (2001)
- **Target**: The default target is the Moon.
- **Grid**: Level 8 icosphere
- **Morphology**: Simple moon (CTEM-based)

.. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
.. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
.. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016    