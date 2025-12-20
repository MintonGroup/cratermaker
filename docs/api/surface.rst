.. _api-surface:

Surface
=======

The Surface class is used for handling surface-related data and operations in the cratermaker project. It provides methods for setting elevation data, calculating distances and bearings, and other surface-related computations. The Surface class extends UxDataset for the cratermaker project.

Available Surface Implementations
---------------------------------

+-----------------------------+------------------------+-----------------------------------------------------------+
| Class                       | Instantiation          | Example Usage                                             |
+=============================+========================+===========================================================+
| IcosphereSurface            | "icosphere"            | surface = Surface.maker("icosphere",                      |
|                             |                        |                          gridlevel=7)                     |
+-----------------------------+------------------------+-----------------------------------------------------------+
| ArbitraryResolutionSurface  | "arbitrary_resolution" | surface = Surface.maker("arbitrary_resolution",           |
|                             |                        |                         pix=100)                          |
+-----------------------------+------------------------+-----------------------------------------------------------+
| HiResLocalSurface           | "hireslocal"           | surface = Surface.maker("hireslocal",                     |
|                             |                        |                         pix=50,                           |  
|                             |                        |                         local_radius=1e3,                 |
|                             |                        |                         local_location=(0,9))             |
+-----------------------------+------------------------+-----------------------------------------------------------+
| DataSurface                 | "datasurface"          | surface = Surface.maker("datasurface",                    |
|                             |                        |                         local_location=(321.9913, 8.121), |
|                             |                        |                         local_radius=50.0e3,              | 
|                             |                        |                         pix=200.0)                        |
+-----------------------------+------------------------+-----------------------------------------------------------+

.. autoclass:: cratermaker.components.surface.Surface
   :members:
   :undoc-members:
   :no-index:


.. autoclass:: cratermaker.components.surface.LocalSurface
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

.. autoclass:: cratermaker.components.surface.hireslocal.LocalHiResLocalSurface
   :members:
   :undoc-members:

Usage example
-------------

.. code-block:: python
   :linenos:

    from cratermaker import Surface
    surface = Surface.maker("hireslocal", pix=50, local_radius=1e3, local_location=(0,9))

.. _api-DataSurface:

.. currentmodule:: cratermaker.components.surface.datasurface

DataSurface grid
=================

See `api-HiResLocalSurface`_ for inherited methods and attributes.

.. autoclass:: cratermaker.components.surface.datasurface.DataSurface
   :members:
   :undoc-members:


Usage example
-------------

This creates a DataSurface with a local region containing Kepler crater. The resolution of the local surface will be approxiamtely half of what it would be if ``pix`` was not set.

.. code-block:: python
   :linenos:

    from cratermaker import Surface
    surface = Surface.maker("datasurface", pix=200.0, local_location=(321.9913, 8.121), local_radius=50.0e3)
