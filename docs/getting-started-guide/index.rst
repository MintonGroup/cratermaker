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

Normal installation
===================
To begin, simply install Cratermaker using pip into the Python environment of your choice, provided it is at least Python 3.10 or above.

.. code-block:: bash 

      pip install cratermaker


Installing Cratermaker from Source
==================================

Cratermaker is mostly Python, but it does have some Rust components. To install Cratermaker from source, you will need to have Rust installed. You can install Rust using the following command:

.. code-block:: bash

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

The build and install an editable version  you can install Cratermaker from source using the following command:

.. code-block:: bash

    pip install -e .      


Simulating a single crater
==========================

Here’s a simple example demonstrating how to initialize a simulation for the Moon, then emplace a single 100 km diameter crater at a random location on the surface:

.. ipython:: python
    :okwarning:
    
    import cratermaker
    sim = cratermaker.Simulation()
    sim.emplace_crater(final_diameter=100e3)


This example uses the default Moon target and emplaces a crater of 100 km in diameter. The ``emplace_crater`` method handles the complexities of crater formation, including the selection of a random location and the necessary calculations to model the impact.

.. note::

      The first time you run a simulation in a given directory, Cratermaker generates a configuration file called ``cratermaker.yaml`` and a surface mesh in ``surface_data\grid.nc``. By default, these are stored in the current working directory from where it was invoked, but this can be changed by passing the argument ``simdir="your\custom\path"`` to ``Simulation``.  


Visualizing the surface with `ParaView <https://www.paraview.org/>`__ 
=====================================================================

Cratermaker can export the surface mesh to a VTK file, which can be visualized with ParaView. The following will export the 
surface mesh to two files inside a folder called ``vtk_files``:

.. ipython:: python
    :okwarning:
    
      sim.export_vtk()

This will generate a series of files with the patern ``surfXXXXXX.vtp``, which represent the grid and its associated data at each output time interval of the simulation. If you haven't already, be sure to `download <https://www.paraview.org/download/>`__ and install ParaView. 

Then, open ParaView and select ``File -> Open``, navigate to the ``vtk_files`` directory, and select the file group. Clicking ``Apply`` will load the file. You can then select ``File -> Open`` again and load the other file. Under ``Coloring`` select ``face_elevation`` or  ``node_elevation`` selected to visualize the surface elevation. 

 The following image shows the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location. The crater is visible as a depression in the surface mesh.

 .. image:: ../_static/paraview_100km_crater_moon.png
      :alt: ParaView visualization of the surface mesh of the Moon with a single 100 km diameter crater emplaced at a random location.
      :align: center

Simulating a population of craters
==================================

Simulating a single crater is useful for testing, but Cratermaker is designed to simulate populations of craters over time. The following example demonstrates how to initialize a simulation of the Moon and emplace a population of craters using the default production function, which is the Neukum production function. The simulation will run for 4.31 Gy, and save the state of the surface in intervals of 50 My.

.. code-block:: python

    sim = cratermaker.Simulation()
    sim.run(age=4300, age_interval=50)

.. note::

      The default units for Cratermaker are meters for length and million years for time.


Default Behavior
================
As Cratermaker is designed to be easy to use, all of its component classes are built to be invoked without any arguments. When a ``Simulation`` object is created as in the examples above, a set of component objects are created with their default parameters. It is important to understand what the default behavior of these classes is, as this will affect the results of your simulation.

- **Target**: The default target is the Moon. There are a number of known bodies that can be selected as targets, including some jovian satellites, and small bodies, but currently few of them have been tested. 
- **Scaling**: The default scaling model is called ``richardson2009``, as it is similar to the one used by the Cratered Terrain Evolution Model (CTEM) that is a progenitor to Cratermaker [1]_. The projectile to transient scaling model is mostly based on Holsapple (1993) [2]_ with some additional scaling parameters for ice given by Kraus et al. (2011) [3]_.
- **Production**: There are two production function models available: ``neukum`` and ``powerlaw``. The default is ``neukum``, and there are three versions available of this model: ``Moon`` [4]_, ``Mars`` [5]_, and ``Projectile`` [6]_.  The version will be selected based on the target body. When the target body is either ``Moon`` or ``Mars``, then one of the two crater-based production functions are selected, and projectile sizes are determined by the Scaling model. For other bodies, ``Projectile`` is chosen, and the crater size is determined by the Scaling model.  
- **Projectile**: There are two projectile velocity and density models available: ``asteroids`` and ``comets``. The default projectile model is determined by the target body. If the body is an inner solar system body, then ``asteroids`` is chosen, otherwise ``comets`` is chosen. The asteroid velocity distributions for the Moon are from Yue et al. (2013) [7]_, and for all other inner planets from Minton et al. (2010) [8]_. Comet velocities are from Zahnle et al. (2003) [9]_.
- **Grid**: There are three grid models available: ``icosphere``, ``arbitrary_resolution``, and ``hireslocal``. The default is ``icosphere``, which builds fast an efficient representation of a sphere. The *resolution* of the grid (the number of faces of the mesh) is determined by the formula :math:`20 \times 4^n`, where :math:`n` is given by the argument ``gridlevel`` with a default value of 8.
- **Morphology**: Currently one morphology model is available: ``simplemoon``. This is a model that similar to that used by CTEM. Most of the parameters are taken from Pike (1977) [10]_, except for simple crater profiles, which use a model from Fassett and Thomson (2014) [11]_. Ejecta blanket scaling is from McGetchin et al. (1973) [12]_.  

Each of the components can be accessed through the ``sim`` object. For example, to access the scaling model:

.. ipython:: python
    :okwarning:

    print(sim.scaling)


References
==========

.. [1] Richardson, J.E., 2009. Cratering saturation and equilibrium: A new model looks at an old problem. Icarus 204, 697-715. https://doi.org/10.1016/j.icarus.2009.07.029
.. [2] Holsapple, K.A., 1993. The scaling of impact processes in planetary sciences 21, 333-373. https://doi.org/10.1146/annurev.ea.21.050193.002001
.. [3] Kraus, R.G., Senft, L.E., Stewart, S.T., 2011. Impacts onto H2O ice: Scaling laws for melting, vaporization, excavation, and final crater size. Icarus 214, 724-738. https://doi.org/10.1016/j.icarus.2011.05.016
.. [4] Neukum, G., Ivanov, B.A., Hartmann, W.K., 2001. Cratering Records in the Inner Solar System in Relation to the Lunar Reference System. Space Science Reviews 96, 55–86. https://doi.org/10.1023/A:1011989004263
.. [5] Hartmann, W.K., Neukum, G., 2001. Cratering Chronology and the Evolution of Mars. Space Science Reviews 96, 165–194.
.. [6] Ivanov, B.A., Neukum, G., Wagner, R., 2001. Size-Frequency Distributions of Planetary Impact Craters and Asteroids, in: Collisional Processes in the Solar System. Springer Netherlands, Dordrecht, pp. 1–34. https://doi.org/10.1007/978-94-010-0712-2_1
.. [8] Yue, Z., Johnson, B.C., Minton, D.A., Melosh, H.J., Di, K., Hu, W., Liu, Y., 2013. Projectile remnants in central peaks of lunar impact craters. Nature Geosci 6, 435 EP-. https://doi.org/10.1038/ngeo1828
.. [7] Minton, D.A., Malhotra, R., 2010. Dynamical erosion of the asteroid belt and implications for large impacts in the inner Solar System. Icarus 207, 744–757. https://doi.org/10.1016/j.icarus.2009.12.008
.. [9] Zahnle, K., Schenk, P., Levison, H., Dones, L., 2003. Cratering rates in the outer Solar System. Icarus 163, 263–289. https://doi.org/10.1016/S0019-1035(03)00048-4
.. [10] Pike, R.J., 1977. Size-dependence in the shape of fresh impact craters on the moon. Presented at the In: Impact and explosion cratering: Planetary and terrestrial implications; Proceedings of the Symposium on Planetary Cratering Mechanics, pp. 489–509.
.. [11] Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. https://doi.org/10.1002/2014JE004698
.. [12] McGetchin, T.R., Settle, M., Head, J.W., 1973. Radial thickness variation in impact crater ejecta: Implications for lunar basin deposits. Earth Planet. Sci. Lett. 20, 226–236. https://doi.org/10.1016/0012-821X(73)90162-3
 

