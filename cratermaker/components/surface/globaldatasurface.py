from __future__ import annotations

from pathlib import Path
from typing import Any
from warnings import warn

import numpy as np
import xarray as xr
from numpy.typing import NDArray

from cratermaker.components.morphology import Morphology
from cratermaker.components.scaling import Scaling
from cratermaker.components.surface import Surface
from cratermaker.components.surface.datasurface import DataSurface
from cratermaker.components.surface.icosphere import IcosphereSurface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import (
    format_large_units,
    parameter,
)


@Surface.register("globaldatasurface")
class GlobalDataSurface(Surface):
    """
    A Surface subclass that generates an entire body using DEM data.

    Currently implements only lunar LOLA data from the PDS. Other data sources may be added in the future.

    Parameters
    ----------
    target : Target, optional
        The target body or name of a known target body for the impact simulation. If none provide, it will be either set to the default,
        or extracted from the scaling model if it is provied
    reset : bool, optional
        Flag to indicate whether to reset the surface. Default is True.
    regrid : bool, optional
        Flag to indicate whether to regrid the surface. Default is False.
    simdir : str | Path
        |simdir|
    pix : FloatLike | None, optional
        The approximate face size. This will be used to determine the target resolution of the DEM data to be used. The actual resolution may be different based on the available DEM data.
    dem_file :  str | Path, optional
        A global DEM files or urs link to a file to be used for the surface. If not provided, a DEM file will be selected based on pix.
    **kwargs : Any

    Returns
    -------
    GlobalDataSurface
        An instance of the GlobalDataSurface class initialized with faces and elevations drawn from DEM data

    Notes
    -----
    Currently this only uses LOLA DEM data for the Moon.

    """

    def __init__(
        self,
        target: Target | str | None = None,
        reset: bool = True,
        regrid: bool = False,
        simdir: str | Path | None = None,
        pix: FloatLike | None = None,
        dem_file: str | Path | None = None,
        **kwargs: Any,
    ):
        try:
            import rasterio  # noqa: F401
        except ImportError as e:
            # IMPORTANT: do not `return` from __init__.
            # Returning leaves a partially-initialized object that will later fail
            # with missing base-class attributes (e.g., _simdir, _output_dir_name).
            raise ImportError(
                "DataSurface requires the optional dependency 'rasterio'. "
                "Install it (and GDAL if needed) to use the 'datasurface' surface type."
            ) from e

        # Temporary storage for DEM data during initialization. This will be cleared after the elevation points are set.
        object.__setattr__(self, "_dem_data", None)

        # The location of the saved DEM data that can be retrieved later
        object.__setattr__(self, "_dem_output_file", None)
        object.__setattr__(self, "_dem_file", None)

        super().__init__(target=target, simdir=simdir, **kwargs)
        if dem_file is None and self.target.name != "Moon":
            raise ValueError("GlobalDataSurface currently only supports the Moon as a target if 'dem_file' is not provided.")
        self._dem_output_file = f"dem_data.{self._output_file_extension}"

        # Set the attributes directly to avoid triggering checks before the pix value is set
        self._pix = pix

        self.dem_file = dem_file

        super()._load_from_files(reset=reset, **kwargs)

        return

    def _get_lola_cylindrical_file(self, resolution) -> str:
        """
        Retrieve the appropriate cylindrically projected LOLA DEM file url for a given resolution from the PDS.

        Parameters
        ----------
        resolution : FloatLike
            Requested resolution in degrees per pixel. The closest available resolution will be used.
        """
        import rasterio
        from affine import Affine
        from rasterio.io import MemoryFile
        from rasterio.merge import merge
        from rasterio.vrt import WarpedVRT
        from rasterio.warp import Resampling
        from rasterio.windows import Window, from_bounds

        AVAILABLE_RESOLUTIONS = [4, 8, 16, 64]  #  pix / deg
        diffs = [abs(resolution - res) for res in AVAILABLE_RESOLUTIONS]
        pds_file_resolution = AVAILABLE_RESOLUTIONS[np.argmin(diffs)]

        return f"https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/float_img/ldem_{pds_file_resolution:d}_float.xml"

    def _get_dem_data(
        self,
    ):
        """
        Read and extract DEM data from one or more files for a specified region.

        Returns
        -------
        dem : dict or None
            Dictionary with x, y, z coordinates and elevation, or None if boundary case.

        """
        try:
            import rasterio
        except ImportError:
            warn("rasterio is not installed. Cannot use this feature.", stacklevel=2)
            return
        from affine import Affine
        from pyproj import Transformer
        from rasterio.io import MemoryFile
        from rasterio.merge import merge
        from rasterio.vrt import WarpedVRT
        from rasterio.warp import Resampling
        from rasterio.windows import Window, from_bounds

        if self.dem_file is None:
            raise ValueError("No DEM file provided to extract data from.")

        _EXPANSION_BUFFER = 1.05
        _NODATA = np.finfo(np.float32).min
        print("Reading DEM file:")
        print(f"  {self.dem_file}")
        with rasterio.open(self.dem_file) as dataset:
            self._pix = dataset.res[0]
            # Determine unit scaling
            units = getattr(dataset, "units", "") or ""
            if isinstance(units, (list, tuple)):
                units = units[0] if len(units) else ""
            scale_factor = 1000.0 if "KILOMETER" in str(units).upper() else 1.0

            # Read the data within the window
            elevation = dataset.read(1) * scale_factor

            # Generate row and column indices
            rows, cols = np.indices(elevation.shape)

            # Compute the x and y coordinates of each pixel in the CRS
            x_coords, y_coords = dataset.xy(rows, cols)
            x_coords = np.array(x_coords).flatten()
            y_coords = np.array(y_coords).flatten()
            elevation = elevation.flatten()

            # Handle nodata values
            if dataset.nodata is None or np.isnan(dataset.nodata):
                valid = np.isfinite(elevation)
            else:
                valid = np.isfinite(elevation) & (elevation != dataset.nodata)
            mean_elevation = np.mean(elevation[valid])
            elevation[~valid] = mean_elevation
            elevation = elevation.astype(np.float32)

            transformer_to_geodetic = Transformer.from_crs(dataset.crs, self.crs)

            # Transform to longitude and latitude
            longitudes, latitudes = transformer_to_geodetic.transform(x_coords, y_coords)
            dem_data = {
                "elevation": elevation,
                "latitudes": latitudes,
                "longitudes": longitudes,
            }

            self._dem_data = dem_data
        return

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.

        """
        print(f"Pixel size: {format_large_units(self.pix, quantity='length')}")

        self._get_dem_data()
        longitudes = np.radians(self._dem_data["longitudes"])
        latitudes = np.radians(self._dem_data["latitudes"])
        x = np.cos(latitudes) * np.cos(longitudes)
        y = np.cos(latitudes) * np.sin(longitudes)
        z = np.sin(latitudes)

        points = np.column_stack([x, y, z])
        print(f"Generated {len(points)} points.")

        decimals = max(6, -int(np.floor(np.log10(self.pix / self.radius))))

        points = np.array(points, dtype=np.float64)
        points = np.round(points, decimals=decimals)
        points = np.unique(points, axis=0).T

        return points

    def reset(self, **kwargs: Any) -> None:
        """
        Reset the surface to its initial state.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|

        """
        super().reset(**kwargs)
        if self._dem_data is not None:
            self._add_dem_elevation()
        elif (self.output_dir / self._dem_output_file).exists():
            with xr.open_dataset(self.output_dir / self._dem_output_file) as ds:
                elevation = np.concatenate([ds.isel(interval=0).face_elevation, ds.isel(interval=0).node_elevation])
                self.update_elevation(elevation)

        return

    def _regrid_if_needed(self, regrid: bool = False, **kwargs: Any) -> bool:
        """
        Check if the existing grid matches the desired parameters determine if regridding is necessary.

        This checks if the dem file exists and matches the current grid.

        Parameters
        ----------
        regrid: bool, optional
            Flag to force regridding even if the grid file exists. Default is False.

        Returns
        -------
        bool
            A boolean indicating whether the grid should be regenerated.
        """
        # DEM data can come from one of two sources. If we are in the middle of building the surface as part of _generate_face_distribution, then the DEM data should be stored in the _dem_data attribute.
        # If this is not set but we already have a grid generated and stored in the uxgrid attribute, then we need to check if a matching DEM file exists on disk. If neither of these sources are avaailable, then we need to regrid.
        if not regrid and self._dem_data is None:
            if not (self.output_dir / self._dem_output_file).exists():
                regrid = True
            elif self.uxgrid is not None:
                with xr.open_dataset(self.output_dir / self._dem_output_file) as ds:
                    if (
                        "face_elevation" not in ds
                        or "node_elevation" not in ds
                        or ds.face_elevation.shape[-1] != self.n_face
                        or ds.node_elevation.shape[-1] != self.n_node
                    ):
                        print(f"DEM datafile {self._dem_output_file} does not match current grid. Regridding needed.")
                        regrid = True

        return super()._regrid_if_needed(regrid=regrid, **kwargs)

    def _add_dem_elevation(self):
        """
        Sample the preserved DEM grid at the face and node centers using bilinear interpolation.

        """
        from scipy.interpolate import LinearNDInterpolator

        if self._dem_data is None:
            raise ValueError("DEM data is not available for interpolation.")
        lonlat = np.c_[self._dem_data["longitudes"], self._dem_data["latitudes"]]
        lut2 = LinearNDInterpolator(lonlat, self._dem_data["elevation"], fill_value=0.0)
        face_elevations = lut2(np.c_[self.face_lon, self.face_lat])
        node_elevations = lut2(np.c_[self.node_lon, self.node_lat])
        elevation = np.concatenate([face_elevations, node_elevations])
        self.update_elevation(elevation)
        self._dem_data = None  # Clear temporary DEM data storage
        return

    @parameter
    def dem_file(self) -> Path | str:
        """
        The file to use for the DEM data in the superdomain.

        """
        return self._dem_file

    @dem_file.setter
    def dem_file(self, value: int | None):
        if value is None:
            target_pds_resolution = np.pi / 180.0 * self.target.radius / self.pix
            self._dem_file = self._get_lola_cylindrical_file(target_pds_resolution)
            return
        if not isinstance(value, str | Path):
            raise TypeError("'dem_file' must be either a string, a Path object, or None.")

        self._dem_file = value

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.

        """
        return super()._hashvars + [self._dem_file]
