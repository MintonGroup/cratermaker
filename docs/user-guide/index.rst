.. currentmodule:: cratermaker

.. _user-guide:

.. toctree::
   :maxdepth: 3

###########
User Guide
###########

This section provides an overview of the components of the Cratermaker package, with several examples to demonstrate how each component works.


Components of Cratermaker
=========================


.. grid:: 1 
    :gutter: 2

    .. grid-item-card:: Simulation
        :columns: 12
        :img-top: ../_images/full_simulation.png
        :class-img-top: dark-light
        :link: simulation
        :link-type: doc

        This is main class that can be used to run a comprehensive cratered landscape evolution model. 


    .. grid-item-card::  Surface
        :columns: 12 6 4 3
        :img-top: ../_images/surface_icon_light.svg
        :class-img-top: only-light
        :class-item: only-light
        :link: surface
        :link-type: doc

        Represents the surface of the target body, including its resolution, topography, material properties. 

    .. grid-item-card::  Surface
        :columns: 12 6 4 3
        :img-top: ../_images/surface_icon_dark.svg
        :class-img-top: only-dark
        :class-item: only-dark
        :link: surface
        :link-type: doc

        Represents the surface of the target body, including its resolution, topography, material properties.


    .. grid-item-card::  Morphology
        :columns: 12 6 4 3
        :img-top: ../_images/morphology_icon_light.svg
        :class-img-top: only-light
        :class-item: only-light
        :link: morphology
        :link-type: doc

        Used to construct the three-dimensional representation of a crater onto the surface. 

    .. grid-item-card::  Morphology
        :columns: 12 6 4 3
        :img-top: ../_images/morphology_icon_dark.svg
        :class-img-top: only-dark
        :class-item: only-dark
        :link: morphology
        :link-type: doc

        Used to construct the three-dimensional representation of a crater onto the surface. 


    .. grid-item-card::  Production
        :columns: 12 6 4 3
        :img-top: ../_images/production_icon.svg
        :class-img-top: dark-light
        :link: production
        :link-type: doc

        Manages the production function, which determines the number of craters/projectiles produced over time and surface area.

    .. grid-item-card:: Counting
        :columns: 12 6 4 3
        :img-top: ../_images/counting_icon_light.svg
        :class-img-top: only-light
        :class-item: only-light
        :link: counting
        :link-type: doc

        Contains methods for counting craters on the surface.

    .. grid-item-card:: Counting
        :columns: 12 6 4 3
        :img-top: ../_images/counting_icon_dark.svg
        :class-img-top: only-dark
        :class-item: only-dark
        :link: counting
        :link-type: doc

        Contains methods for counting craters on the surface.

    .. grid-item-card::  Target
        :columns: 12 6 4 3
        :img-top: ../_images/target_icon.svg
        :class-img-top: dark-light
        :link: target
        :link-type: doc

        Models the properties of the target body, including its size, material properties, and surface gravity.  


    .. grid-item-card::  Projectile
        :columns: 12 6 4 3
        :img-top: ../_images/projectile_icon.svg
        :class-img-top: dark-light
        :link: projectile
        :link-type: doc

        Models the properties of projectiles, including their size, velocity, and density. 


    .. grid-item-card::  Scaling
        :columns: 12 6 4 3
        :img-top: ../_images/scaling_icon_light.svg
        :class-img-top: only-light
        :class-item: only-light
        :link: scaling
        :link-type: doc

        Contains the model for relating projectile size and crater size. 

    .. grid-item-card::  Scaling
        :columns: 12 6 4 3
        :img-top: ../_images/scaling_icon_dark.svg
        :class-img-top: only-dark
        :class-item: only-dark
        :link: scaling
        :link-type: doc

        Contains the model for relating projectile size and crater size. 


    .. grid-item-card::  Crater
        :columns: 12 6 4 3
        :img-top: ../_images/crater_icon.svg
        :class-img-top: dark-light
        :link: crater
        :link-type: doc

        Contains the model for relating projectile size and crater size. 


    .. grid-item-card::  Exporting and Visualization
        :columns: 12 6 4 4
        :img-top: ../_images/visualizing_icon.svg
        :class-img-top: dark-light
        :link: visualizing
        :link-type: doc

        Visualizing data and exporting it to different formats for use with other tools.

    .. grid-item-card::  Example gallery
        :columns: 12 6 4 4
        :img-top: ../_images/gallery_icon.png
        :class-img-top: dark-light
        :link: ../auto_examples/index
        :link-type: doc

        See full working examples of how to use Cratermaker.

    .. grid-item-card::  API reference
        :columns: 12 6 4 4
        :img-top: ../_images/api_icon_light.svg
        :class-img-top: only-light
        :class-item: only-light
        :link: ../api/index
        :link-type: doc

        The API guide contains a detailed description of the Cratermaker API.  The reference describes how the methods work and which parameters can be used. It assumes that you have an understanding of the key concepts.

    .. grid-item-card::  API reference
        :columns: 12 6 4 4
        :img-top: ../_images/api_icon_dark.svg
        :class-img-top: only-dark
        :class-item: only-dark
        :link: ../api/index
        :link-type: doc


        The API guide contains a detailed description of the Cratermaker API.  The reference describes how the methods work and which parameters can be used. It assumes that you have an understanding of the key concepts.

.. toctree::
   :maxdepth: 1

   simulation
   morphology
   production
   counting
   projectile
   scaling
   surface
   target
   crater
   defaults
   visualizing


