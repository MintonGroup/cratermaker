.. currentmodule:: cratermaker

.. _user-guide:

###########
User Guide
###########

This section provides an overview of the components of the Cratermaker package, with several examples to demonstrate how each component works.


.. toctree::
   :maxdepth: 1

   installation


Components of Cratermaker
=========================


.. grid:: 1 
    :gutter: 2

    .. grid-item-card:: Simulation
        :columns: 12
        :img-top: ../_static/full_simulation.png
        :link: simulation
        :link-type: doc

        This is main class that can be used to run a comprehensive cratered landscape evolution model. 


    .. grid-item-card::  Surface
        :columns: 3
        :img-top: ../_static/surface_icon.svg
        :link: surface
        :link-type: doc

        Represents the surface of the target body, including its resolution, topography, material properties. 


    .. grid-item-card::  Morphology
        :columns: 3
        :img-top: ../_static/morphology_icon.svg
        :link: morphology
        :link-type: doc

        Used to construct the three-dimensional representation of a crater onto the surface. 

    .. grid-item-card::  Production
        :columns: 3
        :img-top: ../_static/production_icon.svg
        :link: production
        :link-type: doc

        Manages the production function, which determines the number of craters/projectiles produced over time and surface area.

    .. grid-item-card:: Counting
        :columns: 3
        :img-top: ../_static/counting_icon.svg
        :link: counting
        :link-type: doc

        Contains methods for counting craters on the surface.

    .. grid-item-card::  Target
        :columns: 3
        :img-top: ../_static/target_icon.svg
        :link: target
        :link-type: doc

        Models the properties of the target body, including its size, material properties, and surface gravity.  


    .. grid-item-card::  Projectile
        :columns: 3
        :img-top: ../_static/projectile_icon.svg
        :link: projectile
        :link-type: doc

        Models the properties of projectiles, including their size, velocity, and density. 


    .. grid-item-card::  Scaling
        :columns: 3
        :img-top: ../_static/scaling_icon.svg
        :link: scaling
        :link-type: doc

        Contains the model for relating projectile size and crater size. 


    .. grid-item-card::  Crater
        :columns: 3
        :img-top: ../_static/crater_icon.svg
        :link: crater
        :link-type: doc

        Contains the model for relating projectile size and crater size. 



.. grid:: 1 
    :gutter: 2

    .. grid-item-card::  Example gallery
        :columns: 6
        :img-top: ../_static/index_gallery.svg
        :link: ../auto_examples/index
        :link-type: doc

        See full working examples of how to use Cratermaker.

    .. grid-item-card::  API reference
        :columns: 6
        :img-top: ../_static/index_api_reference.svg
        :link: ../api/index
        :link-type: doc

        The API guide contains a detailed description of the Cratermaker API.  The reference describes how the methods work and which parameters can be used. It assumes that you have an understanding of the key concepts.

.. toctree::
   :hidden:

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


