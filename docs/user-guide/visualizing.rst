.. currentmodule:: cratermaker

.. _ug-visualizing:

Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
=====================================================================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. In this example, we will emplace a 100 km crater on the Moon at latitude and longitude of 0°, and then visualize the surface mesh using ParaView.

.. code-block:: python

    import cratermaker as cm
    sim = cm.Simulation()
    sim.emplace_crater(final_diameter=100e3, location=(0,0))
    sim.export("vtk")


The simulation will generate several files in a folder called ``surface_data``, including ``grid.nc`` and ``surf000000.nc``. When exported to vtk format, a file called ``surf000000.vtp`` will be placed in the ``export``. In this example, the simulation only contains one interval, so only one file is created (see :ref:`ug-Simulation` for how to run multi-interval simulations). 

If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView.  Then, open ParaView and select ``File -> Open``, navigate to the ``export`` directory, and select the file group. Clicking ``Apply`` will load the file. Under ``Coloring`` select either``face_elevation`` or  ``node_elevation`` selected to visualize the surface elevation. 

 The following image shows what you should see. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_100km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location.
      :align: center

