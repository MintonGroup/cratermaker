.. currentmodule:: cratermaker

.. _ug-visualizing:

Visualizing the surface 
=======================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with tools like `PyVista <https://docs.pyvista.org/>`__  and `ParaView <https://www.paraview.org/>`__. In this example, we will emplace a 500 km crater on the Moon at a location of 60° N, 45°E, and then visualize the surface mesh using Pyvista.

The simulation will generate several files in a folder called ``surface``, including ``grid.nc`` and ``surf000000.nc``. When exported to vtk format, a file called ``surf000000.vtp`` will also be placed in the ``surface`` folder. In this example, the simulation only contains one interval, so only one file is created (see :ref:`ug-Simulation` for how to run multi-interval simulations). 

We can then open up the mesh in PyVista for visualization

.. ipython:: python
    :okwarning:
    :suppress:

    # Remove any existing output directory for a clean test
    from pathlib import Path
    out_dirs = ["surface","craters","export"]
    for d in out_dirs:
        if Path(d).exists():
            import shutil
            shutil.rmtree(d)

.. pyvista-plot::

    import cratermaker as cm
    sim = cm.Simulation()
    sim.emplace(final_diameter=500e3, location=(45,60))
    sim.export(driver="VTK")

    import pyvista as pv
    mesh = pv.read("surface/surface000000.vtp")
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, scalars="face_elevation", cmap="cividis", show_edges=False)
    plotter.show()
