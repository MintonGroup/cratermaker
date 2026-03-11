.. currentmodule:: cratermaker

What's New
==========

.. _whats-new.2026.3.3-alpha:

:release:`v2026.3.3-alpha`
--------------------------

- Fixed bug that was causing cookie cutting to extend beyond the crater rim. :pull:`92` `David Minton`_

.. _whats-new.2026.3.2-alpha:

:release:`v2026.3.2-alpha`
--------------------------

- Major refactoring, fixes to plotting tools, and improved installation documentation. :pull:`91` `David Minton`_

   - Refactored the names of components to remove the word "simple" when not referring to simple crater morphology to avoid terminology confusion. This is a major API-breaking change, so be forewarned. `SimpleMoon` is now `BasicMoon`. `SimpleCount`is now `DepthCount`. `Du2025` is now `RealisticMoon`.
   - Global surface labels now have no line break by default and are placed above the main plot. 
   - Turned off scalebar by default in the default Simulation save_action.
   - Placed the plt.close() call after both show and save have been called so that the plot is active until the end.
   - Restructured installation guide page to put dependency install instructions at the beginning.
   

.. _whats-new.2026.3.1-alpha:

:release:`v2026.3.1-alpha`
--------------------------

- Improvements to tally calibration, core structure, and documentation. :pull:`89` `David Minton`_

   - Overhauled the crater tally system with a new set of function that correlate the measured crater depth-to-diameter with degradation state.
   - For long-running examples that use DataSurface, it will now just display a pre-rendered image instead of running the script so that the documentation can be built within the Readthdocs time limit.
   - Improved the way that components are imported in the top level of the project.
   - Improved the consistency of docstrings for class properties.
   - Added the missing Target API page.
   - Improved the organization of the API pages. 


.. _whats-new.2026.3.0-alpha:

:release:`v2026.3.0-alpha`
--------------------------

- Major improvements to documentation and recalibration of the degradation function. :pull:`87` `David Minton`_

   - Recalibrated the visibility and degradation functions based on how Cratermaker computes depth-to-diameter ratios, which is slightly different than how CTEM did it. 
   - Added a new example showing how to use the degradation state function to estimate the degradation state of a crater based on its depth-to-diameter ratio and how it compares with applied diffusion. 
   - Fixed issue causing cross references to not render as links in the documentation. 
   - Fixed bug in :py:func:`~cratermaker.utils.general_utils.format_large_units` that would cause it to give the wrong units for volume. Also added an area formatting option. 

.. _whats-new.2026.2.10-alpha:

:release:`v2026.2.10-alpha`
---------------------------

- Added new function `add_save_action()` that will append a single item to a component's `save_actions` list.
- Fixed a number of issues with saving and plotting HiResLocal surfaces.
- Added a new `relative_location` argument that can be used with the `MorphologyCrater` and its derivatives, and thus can be passed down via the emplace method. Updated the example that uses it to generate a spiral of craters.
- Now both the `Morphology` and `Counting` models have a Crater attribute that can be used to instantiate specialty craters by their chosen Morphology type.  When a Counting object is associated with a Morphology, its Crater attribute takes on the the class defined by its Morphology's Crater attribute.
- Added `ask_overwrite` flags to the export functions that can temporarily override the object's set value during an export operation.
- Changed `measure_degradation_state` to return the degradation state instead of a Crater. Dropped the number of faces down to for counting from 100 to 30 to capture 3 pix diameter craters. :pull:`86` `David Minton`_

.. _whats-new.2026.2.9-alpha:

:release:`v2026.2.9-alpha`
--------------------------

- Improved a number of issues with tallying and configuration file saving. :pull:`83` `David Minton`_

   - Improved counting stability by increasing the number of faces required for a crater to be countable to 100. 
   - the `.tally()` method in `Counting` now takes a `measure_rim` argument, which calls `fit_rim()`. Because of how expensive and unreliable rim fitting is currently, by default it is set to False.
   - Improved handling of config file saving and initializing simulations from config.
   - Fixed fixed problem that was causing neukum production functions to select the wrong version when the target body changed on a reload from config files.

.. _whats-new.2026.2.8-alpha:

:release:`v2026.2.8-alpha`
--------------------------

- Added a new "save_actions" parameter that can be set by any component class and is used to invoke specific method calls on the component when its .save() method is called. This is useful to trigger specific postprocesing actions, like plotting or exporting, at the end of each interval of a Simulation run. :pull:`82` `David Minton`_

.. _whats-new.2026.2.7-alpha:

:release:`v2026.2.7-alpha`
--------------------------

- Fixed bugs introduced due to hastily trying to get the last two releases out. `David Minton`_

.. _whats-new.2026.2.6-alpha:

:release:`v2026.2.6-alpha`
--------------------------

- Fixed minor issues with plotting and documentation. `David Minton`_

.. _whats-new.2026.2.5-alpha:

:release:`v2026.2.5-alpha`
--------------------------

- Fixed some issues with the tally system and added a new function for Surface that can compute the lon,lat location of a point given its distance and bearing from a known location, or the center of a local surface if called from one. This is useful to place a crater on a HiResLocalSurface at a specific location relative to the center without trying to figure out what lat,lon coordinates to use. `David Minton`_

.. _whats-new.2026.2.4-alpha:

:release:`v2026.2.4-alpha`
--------------------------

- Fixed a bug causing local grid files to be deleted on a Simulation reset. `David Minton`_

.. _whats-new.2026.2.3-alpha:

:release:`v2026.2.3-alpha`
--------------------------

- Increased performance by reducing the number of computations when forming craters and to hold off on sorting the tag layers on save and tally. This also improved memory usage by explicitly removing complex data from the crater lists after tallying is finished. :pull:`77` `David Minton`_ 

.. _whats-new.2026.2.2-alpha:

:release:`v2026.2.2-alpha`
--------------------------

- Minor changes to improve the reliability of the export tools when exporting multi-interval runs. `David Minton`_

.. _whats-new.2026.2.1-alpha:

:release:`v2026.2.1-alpha`
--------------------------

- Fixed a bug that was causing craters on slopes to be "stuck" in the tally. This involved a substantial overhaul of the way that craters are processed, including a complete redesign of the system for computing the mean plane that craters are emplaced on.  :pull:`73`: `David Minton`_
- Added new plotting tools for Counting and incorporated it into the Simulation with a "include_counting" option in .plot(). :pull:`74`: `David Minton`_

.. _whats-new.2026.2.0-alpha:

:release:`v2026.2.0-alpha`
--------------------------

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

- Added new DataSurface type that can initialize surface topography from DEM files. :pull:`68` `David Minton`_
- Added new Counting class for crater counting and rim fitting. :pull:`68` `David Minton`_


.. _whats-new.2025.11.3-alpha:

:release:`v2025.11.3-alpha`
---------------------------

- Ensure that the Simulation is saved before the components are exported.

.. _whats-new.2025.11.2-alpha:

:release:`v2025.11.2-alpha`
---------------------------

- Added export methods to Surface, LocalSurface, and Simulation. Surface data can now be exported to VTK, GeoPackage, ESRI Shapefile, and more. :pull:`67` `David Minton`_


.. _whats-new.2025.11.1-alpha:

:release:`v2025.11.1-alpha`
---------------------------

- Minor fix that allows example scripts to show interactive output in the Sphinx Gallery.



.. _whats-new.2025.11.0-alpha:

:release:`v2025.11.0-alpha`
---------------------------


- Major updates to documentation and GIS exporting tools. :pull:`66` `David Minton`_



.. _whats-new.2025.10.1-alpha:

:release:`v2025.10.1-alpha`
---------------------------


- Minor updates to exporting tools.


.. _whats-new.2025.10.0-alpha:

:release:`v2025.10.0-alpha`
---------------------------


- Added crater user guide by :pull:`55` `Dennise Valadez`_



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
   packages by `David Minton`_ :pull:`16`

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
