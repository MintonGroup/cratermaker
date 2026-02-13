.. currentmodule:: cratermaker

What's New
==========

.. _whats-new.2026.2.0-alpha:

:release:`v2026.2.0-alpha`
---------------------------

- Subtantial API changes to the Simulation.run() and .populate(), as well as the Production.populate() methods. Time before present is now called "time" rather than "age." Instead of specifying "age" and "age_end" you now specify "time_start" and "time_end" in units of My before present. The "age" argument is kept for backwards compatability, though it now has a specific meaning of "time_start=age" and "time_end=0". This keeps the terminology consistent with the idea that time is counting down to the present rather than counting up from the past. `David Minton`_
- Added new export methods, including a Spatial Crater Count exporter that will export crater counts in `Craterstats <https://github.com/ggmichael/craterstats>`__ format. The default export method for Simulation is now "OpenCraterTool," which will produce an SCC output of crater counts along with a geotiff file that can be loaded up easily into QGIS with the `OpenCraterTool <https://github.com/thomasheyer/OpenCraterTool>`__ plugin. `David Minton`_
- Added a new "cleanup" function to clean out old simulation data. This is used in the examples to ensure that they can run without prompts about overwriting existing data. `David Minton`_
- Made a number of improvements to the Crater object so that craters are more efficiently processed for emplacement. `David Minton`_


.. _whats-new.2025.12.2-alpha:

:release:`v2025.12.2-alpha`
---------------------------

- Improved the camera focusing in the `show` method of Surface and HighResLocalSurface to better center on the specified location. `David Minton`_
- Fixed the formatting of the documentation so that the grid icons are clickable and also can resize properly on small screens. `David Minton`_

.. _whats-new.2025.12.1-alpha:

:release:`v2025.12.1-alpha`
---------------------------

- Added new DataSurface type that can initialize surface topography from DEM files. :pull:68 `David Minton`_
- Added new Counting class for crater counting and rim fitting. :pull:68 `David Minton`_


.. _whats-new.2025.11.3-alpha:

:release:`v2025.11.3-alpha`
---------------------------

- Ensure that the Simulation is saved before the components are exported.

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
