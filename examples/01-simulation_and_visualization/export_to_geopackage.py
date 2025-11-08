"""
Export Cratermaker data to GeoPackage and visualize with GeoPandas
==================================================================

.. rubric:: By David Minton

Cratermaker can export Simulation data to a variety of GIS vector formats, including GeoPackage, ESRI Shapefile, and more. We do this using the GeoPandas library, which provides a set of useful tools for working with geospatial data in Python. Using the "driver" argument of the `export` method, you can export to just about any format supported by GeoPandas (see `here <https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.to_file.html>`__). This is the same tool used to export data as a mesh in VTK format using the "VTK" driver, which does not use GeoPandas.

 One of the most versitile and efficient formats is GeoPackage (GPKG), which is a single file that can store multiple layers of vector data. This can be read by most GIS software, including QGIS and ArcGIS.  In this example, we will run a simple Cratermaker simulation, export the data to a GeoPackage file, and visualize it using GeoPandas. One of the disandvantages of the current implementation is that when exporting the global surface, faces that cross the antimeridian are not exported due an issue with the UxArray library used to handle the conversion. This will be fixed in a future release.

We will also show how to extract a local region of interest and plot it. The default CRS of local regions is a Azimuthal Equidistant projection centered on the local location, which is well-suited for small areas. Because this does not have the antimeridian issue, we can export and visualize the entire local region without any problems.


"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

import cratermaker

# Run a lunar simulation for 4 Gy and export the global surface to GeoPackage
sim = cratermaker.Simulation(gridlevel=6, ask_overwrite=False)
sim.run(age=4e3)
sim.export(driver="GPKG")

gdf = gpd.read_file(sim.surface.output_dir / "surface000001.gpkg")
gdf.plot("face_elevation", cmap="Greys_r")


# Extract a 1000km radius region at the equator/prime merdian and export to GeoPackage
local = sim.surface.extract_region(location=(0, 0), region_radius=1000e3)
# Pass the interval number explicitly so that the file number matches the last frame of the simulation
local.save(interval_number=sim.interval_number)
local.export(driver="GPKG")
gdf = gpd.read_file(sim.surface.output_dir / "local_surface000001.gpkg")
gdf.plot("face_elevation", cmap="Greys_r")
print("Done")
