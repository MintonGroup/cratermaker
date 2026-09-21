"""
Paint a surface mesh with data from the PDS
===========================================

.. rubric:: By David Minton

Cratermaker is not only able to initialize topography data using real-world data, it can also associate the faces of a mesh with other raster data. Cratemaker uses the powerful rasterio python package to load raster files with georeference metadata attached and associate them with a Surface objct. In this example, we will make a Surface centered on the fresh crater Bandfield. Then we will "paint" the local surface with bolometric temperature anomaly data from the LRO Diviner instrument, which will reveal that Bandfield is a spectacular "cold spot" crater.

"""

import numpy as np
import xarray as xr

from cratermaker import Simulation

simdir = "simdata-1_9"


location = (90.7652, -5.3949)
diameter = 1000.0

sim = Simulation(
    simdir=simdir, surface="datasurface", local_location=location, local_radius=20 * diameter, ask_overwrite=False, reset=True
)

# Add data only to the local surface to avoid downloading data to the low resolution superdomain.
sim.surface.local.add_data(
    data="https://pds-geosciences.wustl.edu/lro/urn-nasa-pds-lro_diviner_derived1/data_derived_ghrm/geotiff/dghrm_tbol_anom_70s70n_tif.xml",
    name="tbol_anom",
    long_name="Diviner Negative Bolometric Temperature Anomaly",
    units="K",
    resampling_order=3,  # By default, data will be applied without resampling (resampling_order=0, aka Nearest Neighbor). This will use cubic interpolation to smooth out the data
)

# If you alter data in-place, you do it on the global surface. Otherwise, you can just call add_data again
tbol_anom = sim.surface.uxds["tbol_anom"]

# Mask out the positive values. NaN values will show the shaded surface with no color
sim.surface.uxds["tbol_anom"] = xr.where(tbol_anom < -2.0, tbol_anom, np.nan)

sim.show3d(variable_name="tbol_anom", cmap="winter")
