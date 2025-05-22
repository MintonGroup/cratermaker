.. currentmodule:: cratermaker


.. image:: ../_static/full_simulation.png
    :alt: Simulation
    :align: center
    :width: 300px

.. _ug-simulation:

##########
Simulation
##########



The Simulation class is the main orchestration tool used to manage and run landscape evolution simulations.


Simulating a population of craters
==================================

Simulating a single crater is useful for testing, but Cratermaker is designed to simulate populations of craters over time. The following example demonstrates how to initialize a simulation of the Moon and emplace a population of craters using the default Neukum production function. The simulation will run for 4.31 Gy, and save the state of the surface in intervals of 50 My.

.. ipython:: python

    import cratermaker as cm
    sim = cm.Simulation()
    sim.run(age=4300, age_interval=50)

.. note::

      The default units for Cratermaker are meters for length and million years for time.



Optional:  Simple Mars simulation visualization with PyVista
============================================================

Once you've simulated a population of craters, you can export the surface mesh and view it in 3D using the `pyvista` Python package.

The example below runs a short simulation on Mars with a population of impactors over a short time scale.

.. ipython:: python

    import cratermaker as cm
    import pyvista as pv
    import glob
    from cratermaker.utils.export import to_vtk

    # simulation of Mars
    sim = cm.Simulation(target="Mars")
    sim.run(age=1000, age_interval=100)  #1 billion years

    # Export the cratered surface to VTK format
    to_vtk(sim.surface)
    print("Exported VTK file")

    # Visualize with PyVista
    vtp_file = sorted(glob.glob("export/*.vtp"))[-1]
    mesh = pv.read(vtp_file)

    scalars = "node_elevation" if "node_elevation" in mesh.array_names else mesh.array_names[0]
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, scalars=scalars, cmap="terrain", show_edges=False)
    plotter.show()

.. note::

    This example requires the ``pyvista`` package. Install it with ``pip install pyvista``.
    You can also open the exported ``.vtp`` file in Paraview for advanced visualization.

``node_elevation`` is a scalar field representing surface height. If it's unavailable, the script will default to the first available scalar for coloring.
