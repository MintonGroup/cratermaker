.. currentmodule:: cratermaker

.. _ug-components-overview:

#########################
Components of Cratermaker
#########################

Each component of Cratermaker is built around a set of classes, which include:

- **Target**: Models the properties of the target body, including its size, material properties, and surface gravity. 
- **Projectile**: Models the properties of projectiles, including their size, velocity, and density. 
- **Scaling**: Contains the model for relating projectile size and crater size. 
- **Production**: Manages the production function, which determines the number of craters/projectiles produced over time and surface area. 
- **Surface**: Represents the surface of the target body, including its resolution, topography, material properties. 
- **Morphology**: Used to construct the three-dimensional representation of a crater onto the surface. 
- **Crater**: Represents a single crater and the properties of the projectile that formed it. 
- **Simulation**: This is main class that can be used to run a comprehensive cratered landscape evolution model. 

The class definition of each component contains a special construction method called ``.maker()``, which is used to create the component. This method can be used to create a component without having to create the entire simulation. While some components require others to function (for instance, a Scaling model needs both a Target and Projectile component in order to perform the scaling computations), when creating such an object, any missing components are built internally using resonable defaults (see :ref:`ug-defaults` for details).

