.. currentmodule:: cratermaker

.. _ug-visualizing:

Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
=====================================================================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. In this example, we will emplace a 500 km crater on the Moon at a location of 60° N, 45°E, and then visualize the surface mesh using ParaView.

.. code-block:: python

    import cratermaker as cm
    sim = cm.Simulation()
    sim.emplace(final_diameter=500e3, location=(45,60))
    sim.export("vtp")


The simulation will generate several files in a folder called ``surface_data``, including ``grid.nc`` and ``surf000000.nc``. When exported to vtk format, a file called ``surf000000.vtp`` will be placed in the ``export``. In this example, the simulation only contains one interval, so only one file is created (see :ref:`ug-Simulation` for how to run multi-interval simulations). 

If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView.  Then, open ParaView and select ``File -> Open``, navigate to the ``export`` directory where you ran the simulation, and select ``surf000000.vtp``. Click ``Apply`` in the Properties windows (or the closed eyeball next to the filename in the Pipeline Browser) and it will render your surface mesh with the default settings.

 The following image shows what you should see. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_500km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 500 km diameter crater.
      :align: center

