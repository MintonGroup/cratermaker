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


Creating a Simulation
===================================

Here’s a simple example demonstrating how to initialize a simulation for the Moon, then emplace a single 100 km diameter crater at a random location on the surface:

.. code-block:: python

    import cratermaker
    sim = cratermaker.Simulation()
    sim.emplace_crater(diameter=100e3)


This example uses the default Moon target and emplaces a crater of 100 km in diameter. The ``emplace_crater`` method handles the complexities of crater formation, including the selection of a random location and the necessary calculations to model the impact.

.. note::

      The first time you run a simulation in a given directory, Cratermaker generates a configuration file called ``cratermaker.yaml`` and a surface mesh in ``surface_data\grid.nc``. By default, these are stored in the current working directory from where it was invoked, but this can be changed by passing the argument ``simdir="your\custom\path"`` to ``Simulation``.  


#####################################################################
Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
#####################################################################

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. The following will export the 
surface mesh to two files inside a folder called ``vtk_files``::

      sim.export_vtk()

This will generate a series of files with the patern ``surfXXXXXX.vtp``, which represent the grid and its associated data at each output time interval of the simulation. If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView. 

Then, open ParaView and select ``File -> Open``, navigate to the ``vtk_files`` directory, and select the file group. Clicking ``Apply`` will load the file. You can then select ``File -> Open`` again and load the other file. Under ``Coloring`` select ``face_elevation`` or  ``node_elevation`` selected to visualize the surface elevation. 

 The following image shows the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_100km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location.
      :align: center



Default Behavior
================

By default, Cratermaker simulates impacts on the Moon's surface. The ``Simulation`` class initializes with the following default parameters:

- **Target**: Moon, represented by the ``target`` parameter. Other bodies can be specified. See :ref:`api-Target`.
- **Production**: The default production function depends on the target body. For inner solar system bodies (Mercury, Venus, Earth, Moon, and Mars) the ``NeukumProduction``, and a simple power law production function is used. For details, see :ref:`api-Production`. 
- **Scaling**: The default projectile to crater scaling relationship model is based on Richardson (2009), with some modifications. See :ref:`api-Scaling` for details.

These defaults provide a balanced starting point for typical crater simulations, making it straightforward for new users to begin without extensive configuration. See :ref:`api-Simulation` for details on the ``Simulation`` class and its parameters.
