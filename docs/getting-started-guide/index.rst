.. currentmodule:: cratermaker

.. _getting-started-guide:

################################
Getting Started with Cratermaker
################################

Cratermaker is designed to be both user-friendly and highly configurable, ensuring that users can start simulations quickly with reasonable defaults or tailor their simulations to specific needs. It is built to be highly modular, and so many of the core components of Cratermaker can be run independently, outside of its core landscape evolution simulation.

.. _getting-started-guide-normal-installation:

Normal installation
===================
To begin, simply install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: bash 

      pip install cratermaker

.. _getting-started-guide-source-installation:

Installing Cratermaker from Source
==================================

Cratermaker is mostly Python, but it does have some Rust components. To install Cratermaker from source, you will need to have Rust installed. You can install Rust using the following command:

.. code-block:: bash

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

The build and install an editable version  you can install Cratermaker from source using the following command:

.. code-block:: bash

    pip install -e .      

.. _getting-started-guide-single-crater:

Simulating a single crater
==========================

Here’s a simple example demonstrating how to initialize a simulation for the Moon, then emplace a single 100 km diameter crater at a random location on the surface:

.. ipython:: python
    :okwarning:
    
    import cratermaker as cm
    sim = cm.Simulation()
    sim.emplace_crater(final_diameter=100e3)


This example uses the default Moon target and emplaces a crater of 100 km in diameter. The ``emplace_crater`` method handles the complexities of crater formation, including the selection of a random location and the necessary calculations to model the impact.

.. note::

      The first time you run a simulation in a given directory, Cratermaker generates a configuration file called ``cratermaker.yaml`` and a surface mesh in ``surface_data\grid.nc``. By default, these are stored in the current working directory from where it was invoked, but this can be changed by passing the argument ``simdir="your\custom\path"`` to ``Simulation``.  

.. _getting-started-guide-visualizing-paraview:

Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
=====================================================================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. The following will export the 
surface mesh to two files inside a folder called ``vtk_files``:

.. ipython:: python
    :okwarning:
    
      sim.export("vtk")

This will generate a series of files with the patern ``surfXXXXXX.vtp``, which represent the grid and its associated data at each output time interval of the simulation. If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView. 

Then, open ParaView and select ``File -> Open``, navigate to the ``export`` directory, and select the file group. Clicking ``Apply`` will load the file. You can then select ``File -> Open`` again and load the other file. Under ``Coloring`` select ``face_elevation`` or  ``node_elevation`` selected to visualize the surface elevation. 

 The following image shows the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_100km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location.
      :align: center


Simulating a population of craters
==================================

Simulating a single crater is useful for testing, but Cratermaker is designed to simulate populations of craters over time. The following example demonstrates how to initialize a simulation of the Moon and emplace a population of craters using the default Neukum production function. The simulation will run for 4.31 Gy, and save the state of the surface in intervals of 50 My.

.. code-block:: python

    import cratermaker as cm
    sim = cm.Simulation()
    sim.run(age=4300, age_interval=50)

.. note::

      The default units for Cratermaker are meters for length and million years for time.


