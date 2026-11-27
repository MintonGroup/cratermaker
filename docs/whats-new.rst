.. currentmodule:: cratermaker

What's New
==========

:release:`v2026.11.1-alpha`
---------------------------

.. _whats-new.2026.11.1-alpha:

- Replaced the experimental "realmoon" morphology model. A realmoon crater is now the "basicmoon" crater warped by an irregular rim outline, plus a random field of variation around it (terraces, slumps, wall and floor texture, and the texture of the ejecta), drawn from an empirical model of fresh lunar craters as a function of size and morphology type (the H-field, distributed as ``hfield_model.nc``). The earlier model based on the power spectra of profile parameters (``PSD1D``, ``psd1d_coeffs.nc`` and ``psd2d_coeffs.nc``) is removed. `David Minton`_

   - Transitional and complex realmoon craters get a central peak built on the surface: the faces where turbulence noise exceeds a threshold are raised to a Gaussian envelope (of standard deviation ``peak_width``, 0.08 crater radii; offset to ``peak_ring_radius`` for a peak ring), and the columns are relaxed by slope collapse at 34 degrees, with the amplitude iterated so the peak reaches ``peak_height``. Peak heights follow the peaks of fresh lunar craters, 0.26 (D / 60 km)^0.89 of the crater's relief with lognormal scatter. The peak's volume is taken from the crater floor, so crater formation still conserves volume.

- Updated the "basicmoon" crater morphology from the morphometry of fresh lunar craters. `David Minton`_

   - The size relations of the floor elevation, rim height, rim width, floor radius, wall curvature, and floor blend are new fits to fresh craters, with the scatter between craters modeled as normal or lognormal, whichever fits better. Outside the range of diameters of the data they are extrapolated by geometric similarity.
   - The structural uplift of the rim now decays exponentially with distance outside the rim (with a decay length of 0.25 crater radii).
   - Complex craters no longer have nested terrace rings. Multiring basins take the profile of Orientale, in proportion to their diameter: the main rim and two nested rings (the Outer and Inner Rook rings), with Orientale's fitted ejecta. Because a basin's final surface keeps little of its excavation, craters with nested rings are not held to volume conservation.
   - The ejecta thickness at the rim defaults to 0.3 of the rim height, and volume conservation now adjusts only the ejecta (``frac_ejrim`` and ``ejprofile``), never the rim height.

- Recalibrated the simple-to-complex transition and the morphology types from the curated fresh lunar craters. `David Minton`_

   - The transition diameter of silicate bodies (in the "montecarlo" and "ctem" scaling models) is now 14.4 km on the Moon, where fresh craters stop being mostly simple, consistent with the 15 km depth/diameter transition of Pike (1980) and Croft (1985). It was 12.2 km ("montecarlo") and 10 km ("ctem").
   - Morphology types now follow the staggered transitions of Pike (1980) and Schenk et al. (2021): transitional (immature complex) craters set in at the transition diameter and complex (central peak) craters at 1.8 times it, each over a logistic step in diameter. Previously, transitional craters could occur from half to twice the transition diameter, so many craters of 6-12 km on the Moon were transitional, where fresh craters are almost all simple.
   - In the "montecarlo" scaling model, the enlargement of complex craters relative to their equivalent simple crater is now 1.10 (D / D_transition)^0.19, measured from the cavity volumes of fresh lunar craters (within the range of Croft 1985), instead of 1.02 (D / D_transition)^0.079 (Holsapple 1993). Transient craters of large complex craters are smaller, by about 14% at 100 km on the Moon.

- Added a ``count_superdomain`` option to the counting models (``Counting.maker``, and so ``Simulation``) to count craters in the superdomain of a hireslocal surface as well as those in the local region. It can also be given to ``emplace`` for the craters of one call. Previously this was only possible by passing ``tag_superdomain`` to ``emplace``, which did not work because ``emplace`` did not pass its keyword arguments on to crater formation and counting (also fixed).
- Fixed the ``smallest_crater`` and ``largest_crater`` arguments of ``Simulation``, which were ignored when given to the constructor and lost when a simulation was resumed.
- ``show3d`` now draws the surface as a light grey, matte material by default, which shows the topography with more contrast against the dark background.
- Fixed the ``annotations`` and ``pickable`` arguments of ``show3d``, which were not passed on to PyVista.
- Fixed the face gradient on the mesh (used by ``compute_slope``, ``slope_collapse``, and the radial and azimuthal gradients), which returned half the true gradient. Slope collapse now limits slopes to the requested critical angle; previously a critical angle of 35 degrees allowed slopes of about 54 degrees. The critical angle of slope collapse after crater formation is now 40 degrees, from the measured interior slopes of fresh lunar craters, which sit at about 35 degrees on their walls but rarely exceed 40 degrees.
- Copying a "realmoon" crater now keeps the field it has already drawn, and no longer draws random numbers.


:release:`v2026.11.0-alpha`
---------------------------

.. _whats-new.2026.11.0-alpha:

- Improved performance by moving most of the crater formation into Rust. :issue:`150` :pull:`152` `David Minton`_

   - Crater formation (``form_crater``) is now a single compiled call per crater, using a Rust pipeline shared by the "basicmoon" and "realmoon" morphology models. Each batch of non-overlapping craters is formed in parallel.
   - Moved the reference surface fit, the degradation function, the search for the extent of a crater's region (``rmax``), and the projectile to transient crater scaling into Rust.
   - Faster region building for large craters, and the crater region is reused as the ejecta region when the two are the same.
   - Slope collapse now only revisits the faces near those that changed in the previous iteration, which greatly reduces its cost for large craters.
   - Faster crater creation, with fewer internal copies of each crater and faster sampling from truncated normal distributions.
   - The counting component no longer extracts a tagging region for each crater when ``save_observed`` is off.
   - Removed the Python implementations of ``crater_shape``, ``ejecta_shape``, ``compute_ejecta``, and ``ejecta_burial_degradation`` that are replaced by the Rust pipeline.
   - Results agree with the previous version to within rounding error, except where noted below.

- Each crater's projectile properties (velocity, impact angle, and direction) are now drawn once instead of twice. This changes the random numbers each crater receives, so a simulation run with the same ``rng_seed`` gives different, but statistically equivalent, results than earlier versions.
- ``Crater.maker(..., monte_carlo=False)`` now gives a vertical impact due north at the mean velocity, as documented, instead of reusing the angle and direction of an earlier crater.
- Fixed a bug where the "realmoon" ejecta blanket was always computed at a bearing of 0, so it did not vary with azimuth along with the crater rim. This also fixes a crash when using the "realmoon" model with ``dorays=True``.
- Fixed a bug where the ejecta fraction of the rim height (``frac_ejrim``) of the rings of complex and multiring craters was limited to be at least 1, rather than at most 1, after volume conservation.
- Fixed a bug where quasi-Monte Carlo craters read from a file were not reproducible with ``rng_seed``, because they were created without the simulation's random number generator. ``Crater.from_file`` and the related ``from_csv_file``, ``from_scc_file``, and ``from_xarray`` methods now pass their keyword arguments on to ``Crater.maker``.
- Copying a "realmoon" crater now keeps the PSD data it has already generated and its rings (as independent copies), and no longer draws random numbers.
- Added a ``morphology_type`` argument to ``Crater.maker`` to set the morphology type of a crater, and a new example showing the basic morphology types.
- Added the ``save_observed`` and ``save_emplaced`` arguments to separately control the saving of observed and emplaced craters.
- Emplaced craters are now saved even when counting is turned off. :issue:`142` :pull:`143`
- Fixed a grid aliasing problem and a resolution limitation of the HiResLocal surface. :issue:`144` :pull:`145`
- Increased the size of the crater id from uint32 to uint64 to prevent hash collisions, and made the crater id type consistent throughout. Crater ids are stored as int64 when exported to a GeoPackage, which has no unsigned 64-bit type. :issue:`146` :pull:`147`
- Fixed incorrect values returned by the Neukum production function. :issue:`148` :pull:`149`
- Fixed aliasing in ``Surface.to_raster``.
- Fixed a discontinuity in the crater profile.
- Fixed an artifact in "realmoon" profiles across the 0/360 degree bearing boundary, and improved the performance of the profile interpolation.
- Fixed ``LocalSurface.extract_subregion``, which was missing the buffer faces and nodes around the subregion.
- Fixed the value used when removing a crater tag from the surface.
- Fixed the velocity distribution in the random projectile velocity example.
- Other performance improvements, including building local surface arrays in Rust and caching the simulation directory.


:release:`v2026.10.3-alpha`
---------------------------

.. _whats-new.2026.10.3-alpha:

- Improved performance of experimental "realmoon" morphology model. :issue:`140` :pull:`141` `David Minton`_


:release:`v2026.10.2-alpha`
---------------------------

.. _whats-new.2026.10.2-alpha:

- Improved documentation and added a new save_action example. :issue:`138` :pull:`139` `David Minton`_
- Changed ambiguous arguments: `density` is now either `target_density` or `projectile_density`, likewise for `mass`, and other projectile parameters.

.. _whats-new.2026.10.1-alpha:

:release:`v2026.10.1-alpha`
---------------------------


- Added ability to pass callable functions as save_actions. :issue:`136` :pull:`137` `David Minton`_
- Updated ring generator for complex craters.


.. _whats-new.2026.10.0-alpha:

:release:`v2026.10.0-alpha`
---------------------------

- First draft of the RealisticMoonMorphology (accessible via morphology="realmoon" argument) is not available for preview and testing. :issue:`18` :pull:`134` `David Minton`_

   - Implemented all 1D PSD-based parameters: `radius`, `rim_height`, `rim_width`, and `floor_radius`. 
   - Made improvements to stability of crater profile function to avoid NaN values.
   - Improved aesthetics of the pyvista_plotter defaults, including adding a new PyVista theme called "cratermaker".
   - Increased the size of the region used to compute the reference surface for emplacing craters, but only use the portion inside the rim radius in order to prevent issues of small region producing anomolously large slopes.
   - Turn off volume conservation when reading in craters from file
   - Allow conserve_volume=False argument to be respected when a modified copy of a Crater is made
   - Switched to a linear rim uplift model with fewer free parameters. Fixed issue where rng was not being passed to the sampling functions when making craters.
   - Updated morphology of rims to use a gaussian rim uplift and now define ejecta thickness at rim as a fraction of the rim height. Switched to using relative elevation for rim as rim_height instead of absolute rim_elevation
   - Overhauled profile functions again and fixed issue that was causing old datasurface to not be recognized as the same when parameters were called the same. This was done by turning off Monte Carlo scaling when setting the superdomain values so that they are deterministicly computed.
   - Fixed problem that was causing reloads of datasurface to revert to the default resolution because it did not save the dem_file_list argument. The _user_defined argument list gets reset when calling super().__init__() functions inside of __init__, so any arguments set before calling the base class must be restored to this list.
   - Rearranged code for extracting pix values so that you don't have to fetch it from the source if you are reading in an old run where it's already defined.
   - Fixed inconsistency in the definition of pix when using DataComposer to retrieve cylindrical LOLA data. The returned value should be m/px not px/deg.
   - Fixed issue that was causing parameters to not be saved to file in DataSurface, so old surfaces were not being recognized even though they were the same.
   - Fixed issues fetching high resolution data and polar data above 87.5 degrees
   - Fixed issue with selecting PDS files for polar data when a feature straddles a +/-60 deg latitude line
   - Removed old macos-14 from build, test, and deploy actions
   - Removed ubuntu arm from testing due to dependency issue
   - Removed intel mac and arm64 ubuntu from action matrix
   - Improved the parallel efficiency of the apply_diffusion function
   - Removed the interpolated fill function as it causes issues near boundaries of NAC DTMs. Allow hillshade parameters azimuth, solar_angle, and blend_mode to be passied in as kwargs in the plot method.
   - Generalized the basic update_elevation and add_data methods to be able to use raster datasets, rather than needing a separate data_composer context manager. This allows any type of arbitrary dataset data to be added to any surface.
   - Fixed problem causing restarted runs to fail due to Counting not being able to read in all of the data because the Morphology component has to be initialized after Counting.
   - When applying a dataset to a surface, only update the surface with values that are not NaN
   - Fixed problem that occurs when passing overwrite=True to add_data when the dataset is new and there are NaN values in  the data. Dataset has to be created first before mask is made
   - Switched from positional arguments to keyword arguments for calls to add_data
   - Fixed problem where the input data source was altered when passed to the DataComposer.
   - Fixed issue that occurs when trying to use an SLDEM tile that goes past the SLDEM latitude boundaries of +/-60.
   - Allow the Counting plotter to accept custom surfaces
   - Fixed issue causing failure of datasurface for low resolution mid latitude locations because sldem is only available for 256 and 512 pix/deg
   - Fixed problem that was causing 429 errors when trying to download large datasets from the PDS line-by-line. Lines are now chunked into increments of 1000
   - Fixed problem fetching LOLA dems because unlike all of the other datasets, only sldem256 does not use zero-padded integers for the latitude range in the file name.
   - Restructured the code to move the elevation change step to the end of form_crater so that all of diffusion steps are applied to the pre-existing surface before drawing the main crater shape.
   - Added new azimuthal gradient calculator function
   - Simplified the hash variable for datasurface to avoid regridding when not necessary
   - Removed the fit_rim function and the dependence on OpenBLAS
   - Removed rim fitting example
   - Fixed issue where sometimes the ejecta_region was smaller than the crater_region. The ejecta_region radius is now adjusted to be no smaller than crater_region
   - Updated morphometric coefficients.
   - Fixed incompatabilities with newer versions of UxArray that don't allow calling .sel() using both n_face and n_node
   - Fixed issues getting correct indices for nodes due to changes in UxArray
   - Swapped out netcdf4 for h5netcdf engine for open_mfdataset
   - Added a reset method to the MorphologyCrater that is called after forming, which allows crater objects to be re-used in different contexts
   - Added a new load method that will load a particular interval of saved data into a Simulation
   - Reduced the time taken to get a volume estimation for small craters.
   - Fixed a bunch of issues related to getting an estimated initial volume
   - Made the monte_carlo flag more consistent throughout the components. Added ability to disable monte carlo sampling for a single Crater instance via an argument.
   - Lifted the restriction on Python 3.14.


.. _whats-new.2026.9.0-alpha:

:release:`v2026.9.0-alpha`
--------------------------

- Overhauled much of the morphology code in preparation for a set of much improved morphology models in the works. :issue:`130` :pull:`132` `David Minton`_
   - Built new profile function with controllable morphology inputs, and nesting for basins/complex craters/concentric craters.
   - Improved the models for crater morphometry initial conditions. 
   - Removed the `fit_rim` function. This is being re-written as part of the development of the improved morphology models. Removing this allows the code to be compiled without requiring an OpenBLAS library install.
   - Fixed a number of edge case bugs when fetching default datasets using the DataSurface.
   - Fixed problem that arose due to selecting both faces and nodes in a LocalSurface with newer versions of UxArray. Now only faces are selected.
- Added new functionality to the `Surface.add_data()` method that allows arbitrary raster data to be applied to any mesh as an associated UxArray DataArray. `David Minton`_ and `Leo Minton`_

.. _whats-new.2026.5.0-alpha:

:release:`v2026.5.0-alpha`
--------------------------

- Fixed :issue:`128` :pull:`129` by ensuring that the Scaling model gets assigned to the Morphology model, and that density can be passed as an argument to the `new_projectile` method. `David Minton`_


.. _whats-new.2026.4.4-alpha:

:release:`v2026.4.4-alpha`
--------------------------

- Fixed :issue:`124` :pull:`125` `David Minton`_
- Fixed a problem where the Simulation's rng was not being passed to the Crater wrapper, so that the simulation.Crater was not producing repeatable results. :pull:`125` `David Minton`_
- Removed "addopts = tests" option from the pytest.ini addopts so that individual tests can be run rather than always all of them. :pull:`125` `David Minton`_

.. _whats-new.2026.4.3-alpha:

:release:`v2026.4.3-alpha`
--------------------------

- Fixed :issue:`106` :pull:`123` `David Minton`_
- Fixed format_large_units to handle negative values.
- Added a colorbar and a minimum_plot_width option to the plot method of Surface. :pull:`123` `David Minton`_

.. _whats-new.2026.4.2-alpha:

:release:`v2026.4.2-alpha`
--------------------------

- Improved the documentation with more thorough discussion of how Quasi Monte Carlo works and how its behavior can be modified in the production user guide. :pull:`122` `David Minton`_
- Fixed a bug that was causing a failure when numpy>2.3. This was caused by sloppy conversions between 0-d Numpy arrays and floats. Conversions are now done explicitly for the affected root finding functions.
- Fixed :issue:`121` :pull:`122` `David Minton`_
- Fixed a bug causing the Neukum production functions to revert to the Moon if version was specified but not target, which caused the plots to show only the Moon versions instead of the Moon and Mars. `David Minton`_

.. _whats-new.2026.4.1-alpha:

:release:`v2026.4.1-alpha`
--------------------------

- Fixed minor issues with documentation to ensure all elements are included in the build. `David Minton`_
- Fixed problem where the Neukum production function was creating too many basins by adjusting the power law exponent of the upper extrapolation to be steeper by -5.0. :pull:`120` `David Minton`_

.. _whats-new.2026.4.0-alpha:

:release:`v2026.4.0-alpha`
--------------------------

- Fixed bug that was causing incorrect polar data files from being requested from the PDS for DataSurface :pull:`112` `David Minton`_ `Leo Minton`_
- Fixed bug that was causing `fill_value` and `dtype` arguments to be ignored in `add_data`. :pull:`114` `David Minton`_ `Leo Minton`_
- Created a new system to retrieve, merge, and apply data to surfaces :pull:`118` `Leo Minton`_
- Added Quasi-Monte Carlo :pull:`119` `David Minton`_ `Austin Blevins`_`

.. _whats-new.2026.3.7-alpha:

:release:`v2026.3.7-alpha`
--------------------------

- Include new argument pass-throughs for the show3d methods to better control PyVista output. :pull:`104` `David Minton`_
- Updated the README with better instructions. :pull:`105` `David Minton`_


.. _whats-new.2026.3.6-alpha:

:release:`v2026.3.6-alpha`
--------------------------

- More updates to the documentation, including the return of the long-running DataSurface examples thanks to an upgrade in our ReadtheDocs time limit, and a fancy new animated gif of the Moon being bombarded for the visualization user guide. :pull:`101` `David Minton`_

.. _whats-new.2026.3.5-alpha:

:release:`v2026.3.5-alpha`
--------------------------

- Gave the documentation pages a thorough makeover. :pull:`99` `David Minton`_

   - Updated the icons for the documentation pages to be more visually appealing and to better represent the content of each page. 
   - Improved the formatting of the documentation pages to make them more readable and easier to navigate. 
   - Added new sections to the user guide, including a new section on visualizing and exporting data. 

.. _whats-new.2026.3.4-alpha:

:release:`v2026.3.4-alpha`
--------------------------

- Fixed bug that was causing restarted runs to fail with the message "AttributeError: Crater has no attribute 'measured_depth_to_diameter'."  :pull:`96` `David Minton`_

.. _whats-new.2026.3.3-alpha:

:release:`v2026.3.3-alpha`
--------------------------

- Fixed bug that was causing cookie cutting to extend beyond the crater rim. :pull:`93` `David Minton`_

.. _whats-new.2026.3.2-alpha:

:release:`v2026.3.2-alpha`
--------------------------

- Major refactoring, fixes to plotting tools, and improved installation documentation. :pull:`91` `David Minton`_

   - Refactored the names of components to remove the word "simple" when not referring to simple crater morphology to avoid terminology confusion. This is a major API-breaking change, so be forewarned. `SimpleMoon` is now `BasicMoonMorphology`. `SimpleCount`is now `DepthCount`. `Du2025` is now `RealisticMoonMorphology`.
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
   - Added the missing earget API page.
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
- `Leo Minton`_

.. _David Minton: https://github.com/profminton
.. _Austin Blevins: https://github.com/austinblevins
.. _Jun Du: https://github.com/jundu-dr-crater
.. _Dennise Valadez: https://github.com/dennvee
.. _Elizabeth Norman: https://github.com/enorman98
.. _Leo Minton: https://github.com/leo-minton
