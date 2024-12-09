.. currentmodule:: cratermaker

################
Getting Started
################

First, set up a Python environment and set the system environment variables using ```buildscripts/set_environment.sh```. This script sets the necessary environment variables for the Fortran compiler, CMake, and other dependencies.  You can install cratermaker with pip::
   
      pip install -e .

You can verify the installation was successful by running the tests::

      pytest tests/


##################
Emplacing a Crater
##################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs.

Default Behavior
================

By default, Cratermaker emulates impacts on the Moon's surface. The ``Simulation`` class initializes with the following default parameters:

- **Target Body**: Moon, represented by the ``target`` parameter. Other bodies can be specified. See :ref:`api-Target`.
- **Material**: The default material for the Moon is ``soft rock``. Custom materials can be specified. See :ref:`api-Material`.
- **Production**: The default production function depends on the target body. For inner solar system bodies (Mercury, Venus, Earth, Moon, and Mars) the ``NeukumProduction``, and a simple power law production function is used. For details, see :ref:`api-Production` and :ref:`api-NeukumProduction`.
- **Scaling**: The default projectile to crater scaling relationship model is based on Holsapple (1993), with some modifications. See :ref:`api-Scale` for details.
- **Pixel Resolution**: Cratermaker uses an unstructured mesh to represent the surface, and so does not have a fixed resolution. When generating the mesh, you can set an approximate size scale for the mesh faces using the ``pix`` parameter. The default value is 10\ :sup:`-3`` times the surface area of the target body.  

These defaults provide a balanced starting point for typical crater simulations, making it straightforward for new users to begin without extensive configuration. See :ref:`api-Simulation` for details on the ``Simulation`` class and its parameters.

Customization
=============

All default parameters can be customized. Users can specify different target bodies, materials, resolution levels, and models for scaling, morphology, and production functions. This flexibility allows advanced users to tailor the simulation to specific research needs or to simulate conditions on other celestial bodies.

Example: Emplacing a Single Crater
===================================

Here’s a simple example demonstrating how to initialize a simulation for the Moon, then emplace a single 100 km diameter crater at a random location on the surface:

.. code-block:: python

    import cratermaker
    sim = cratermaker.Simulation()
    sim.emplace_crater(diameter=100e3)


This example uses the default Moon target and emplaces a crater of 100 km in diameter. The ``emplace_crater`` method handles the complexities of crater formation, including the selection of a random location and the necessary calculations to model the impact.

.. note::

      The first time you run a simulation in a given directory, Cratermaker generates a surface mesh using Jigsaw, which can take several minutes at the default resolution. The generated mesh is stored in the ``surface_data`` directory for quick loading in subsequent runs, provided that you don't change the grid parameters. 


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
