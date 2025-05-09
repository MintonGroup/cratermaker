.. _api-surface:

Surface
=======

The Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

Available Surface Implementations
---------------------------------

+-----------------------------+------------------------+-------------------------------------------------+
| Class                       | Instantiation          | Example Usage                                   |
+=============================+========================+=================================================+
| IcosphereSurface            | "icosphere"            | surface = Surface.maker("icosphere",            |
|                             |                        |                       gridlevel=7)              |
+-----------------------------+------------------------+-------------------------------------------------+
| ArbitraryResolutionSurface  | "arbitrary_resolution" | surface = Surface.maker("arbitrary_resolution", |
|                             |                        |                       pix=100)                  |
+-----------------------------+------------------------+-------------------------------------------------+
| HiResLocalSurface           | "hireslocal"           | surface = Surface.maker("hireslocal",           |
|                             |                        |                       pix=50, local_radius=1e3, |
|                             |                        |                       local_location=(0,9))     |
+-----------------------------+------------------------+-------------------------------------------------+

.. autoclass:: cratermaker.components.surface.Surface
   :members:
   :undoc-members:
   :no-index:

.. _api-IcosphereSurface:

.. currentmodule:: cratermaker.components.surface.icosphere

Icosphere grid
==============

See `Surface`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.surface.icosphere.IcosphereSurface
   :members:
   :undoc-members:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Surface
    surface = Surface.maker("icosphere",gridlevel=7)

.. _api-ArbitraryResolutionSurface:

.. currentmodule:: cratermaker.components.surface.arbitrary_resolution

Arbitrary resolution grid
=========================

See `Surface`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.surface.arbitrary_resolution.ArbitraryResolutionSurface
   :members:
   :undoc-members:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Surface
    surface = Surface.maker("arbitrary_resolution",pix=100)

.. _api-HiResLocalSurface:

.. currentmodule:: cratermaker.components.surface.hireslocal

Hi-res local grid
=================

See `Surface`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.surface.hireslocal.HiResLocalSurface
   :members:
   :undoc-members:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Surface
    surface = Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9))
