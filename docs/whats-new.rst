.. currentmodule:: cratermaker

What's New
==========

.. _whats-new.2025.11.2-alpha:

:release:`v2025.11.2-alpha`
---------------------------

- Added export methods to Surface, LocalSurface, and Simulation. Surface data can now be exported to VTK, GeoPackage, ESRI Shapefile, and more. :pull:67 `David Minton`_


.. _whats-new.2025.11.1-alpha:

:release:`v2025.11.1-alpha`
---------------------------

- Minor fix that allows example scripts to show interactive output in the Sphinx Gallery.



.. _whats-new.2025.11.0-alpha:

:release:`v2025.11.0-alpha`
--------------------------


- Major updates to documentation and GIS exporting tools. :pull:66 `David Minton`_



.. _whats-new.2025.10.1-alpha:

:release:`v2025.10.1-alpha`
--------------------------


- Minor updates to exporting tools.


.. _whats-new.2025.10.0-alpha:

:release:`v2025.10.0-alpha`
--------------------------


- Added crater user guide by :pull:55 `Dennise Valadez`_



.. _whats-new.2025.11.0-alpha:

:release:`v2025.11.0-alpha`
---------------------------


- Major updates to documentation and GIS exporting tools. :pull:66 `David Minton`_


.. _whats-new.2025.10.1-alpha:

:release:`v2025.10.1-alpha`
---------------------------


- Minor updates to exporting tools.


.. _whats-new.2025.10.0-alpha:

:release:`v2025.10.0-alpha`
---------------------------


- Added crater user guide by :pull:55 `Dennise Valadez`_


.. _whats-new.2025.4.1-alpha:

:release:`v2025.4.1-alpha`
--------------------------

Bug Fixes
~~~~~~~~~

- Fixed broken dependencies in environment definitions.

Contributors
~~~~~~~~~~~~
- `David Minton`_


.. _whats-new.2025.4.0-alpha:

:release:`v2025.4.0-alpha`
--------------------------

What’s Changed
~~~~~~~~~~~~~~

-  Wheels built for Linux and Mac are now deployed to PyPI.
-  Restructure code to remove the strict dependency on the conda only
   packages by `David Minton`_ :pull:16

Contributors
~~~~~~~~~~~~
- `David Minton`_


.. _whats-new.2024.11.0-alpha:

:release:`v2024.11.0-alpha`
---------------------------

- Updated the requirements to use newer versions of xarray and uxarray. It should be compatible with Numpy 2 now. 
- UxArray changed the `face_area` calculation to use the true dimensions rather than assuming a unit sphere. Changed it so that the area of faces is correct and added a new test.`

Contributors
~~~~~~~~~~~~
- `David Minton`_

.. _whats-new.2024.8.0-alpha:

:release:`v2024.8.0-alpha`
--------------------------

- Prerelease for development and testing.

Contributors
~~~~~~~~~~~~
- `David Minton`_
- `Austin Blevins`_
- `Jun Du`_
- `Dennise Valadez`_
- `Elizabeth Norman`_

.. _David Minton: https://github.com/profminton
.. _Austin Blevins: https://github.com/austinblevins
.. _Jun Du: https://github.com/jundu-dr-crater
.. _Dennise Valadez: https://github.com/dennvee
.. _Elizabeth Norman: https://github.com/enorman98
