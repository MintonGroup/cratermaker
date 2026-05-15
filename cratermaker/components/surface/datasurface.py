from __future__ import annotations

from pathlib import Path
from typing import Any
from warnings import warn

import numpy as np
import rasterio
import xarray as xr
from numpy.typing import ArrayLike, NDArray

from cratermaker.components.surface import DataComposer, Surface
from cratermaker.components.surface.hireslocal import HiResLocalSurface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import format_large_units, parameter

_DEFAULT_N_FACES_LOCAL = 1e6


@Surface.register("datasurface")
class DataSurface(HiResLocalSurface):
    """
    A Surface subclass that generates a local region using DEM data.

    Currently implements only lunar LOLA data from the PDS. Other data sources may be added in the future.

    Parameters
    ----------
    local_radius : FloatLike
        The radius of the local region in meters.
    local_location : PairOfFloats
        The longitude and latitude of the location in degrees.
    superdomain_scale_factor : FloatLike, optional
        A factor defining relative size of the face at the antipode of the local region to the face size inside the local region.
        If not provided, construction of the surface will be deferred until the `set_superdomain` method is called.
        If a negative number is provided, it will be computed based on a provided (or default) scaling and morphology model, and will be set so that smallest craters that can be resolved on faces outside the local region could potentially deposit ejecta at the boundary of the local region.
        Default is -1 (which triggers automatic computation based on scaling and morphology models).
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
        The approximate face size inside the local region in meters. This will be used to determine the target resolution of the DEM data to be used. The actual resolution may be different based on the available DEM data. Note that if you provide a list of DEM files using the `dem_file_list` parameter, this value will be ignored. If None, is set, and no file(s) are provided, a default resolution that creates approximately 1e6 faces in the local region will be used.
    local_radius : FloatLike
    dem_file_list : list of str | Path, optional
        A list of one or more DEM files or urls links to files to be used for the local region. If not provided, a set of files will be determined based on the location and radius of the region. It will chose a set of files that fully covers the local region with approximately approximately 1e6 faces.
    superdomain_dem_file :  str | Path, optional
        A global DEM files or urs link to a file to be used for the superdomain region. If not provided, a default global DEM file will be used based on the target body.
    **kwargs : Any

    Returns
    -------
    DataSurface
        An instance of the DataSurface class initialized with faces and elevations drawn from DEM data

    Notes
    -----
    Currently this only uses LOLA DEM data for the Moon.

    """

    def __init__(
        self,
        local_radius: FloatLike,
        local_location: PairOfFloats,
        superdomain_scale_factor: FloatLike | None = -1,
        target: Target | str | None = None,
        reset: bool = True,
        regrid: bool = False,
        simdir: str | Path | None = None,
        pix: FloatLike | None = None,
        dem_file_list: list[str] | list[Path] | None = None,
        superdomain_dem_file: str | Path | None = None,
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
        object.__setattr__(self, "_local", None)

        # Temporary storage for DEM data during initialization. This will be cleared after the elevation points are set.
        object.__setattr__(self, "_local_dem_data", None)
        object.__setattr__(self, "_global_dem_data", None)

        # The location of the saved DEM data that can be retrieved later
        object.__setattr__(self, "_dem_output_file", None)
        object.__setattr__(self, "_dem_file_list", None)
        object.__setattr__(self, "_superdomain_dem_file", None)

        super(HiResLocalSurface, self).__init__(target=target, simdir=simdir, **kwargs)
        if local_radius > 1.49 * target.radius:
            raise ValueError(
                "The value of local_radius is too large. Consider using DataComposer on a global surface type, such as Icosphere."
            )
        if dem_file_list is None and self.target.name != "Moon":
            raise ValueError("DataSurface currently only supports the Moon as a target if 'dem_file_list' is not provided.")
        if superdomain_dem_file is None and self.target.name != "Moon":
            raise ValueError("DataSurface currently only supports the Moon as a target if 'superdomain_dem_file' is not provided.")
        if local_location is None:
            raise ValueError("local_location must be provided.")
        self._output_file_pattern += [f"local_{self._output_file_prefix}*.{self._output_file_extension}"]
        self._dem_output_file = f"dem_data.{self._output_file_extension}"

        # Set the attributes directly to avoid triggering checks before the pix value is set
        self._local_radius = local_radius
        self._local_location = local_location
        self._pix = pix
        self.dem_file_list = dem_file_list
        # Temporarily disable ask_overwrite to avoid prompts during initialization
        ask_overwrite = self.ask_overwrite
        self.ask_overwrite = False
        super_kwargs = {
            **kwargs,
            "pix": self.pix,
            "local_radius": local_radius,
            "local_location": local_location,
            "superdomain_scale_factor": superdomain_scale_factor,
            "target": self.target,
            "reset": reset,
            "regrid": regrid,
            "simdir": self.simdir,
            "ask_overwrite": self.ask_overwrite,
        }
        super().__init__(**super_kwargs)

        self._superdomain_dem_file = superdomain_dem_file
        self.ask_overwrite = ask_overwrite
        return

    def _get_location_extents(self):
        """
        Computes the longitude and latitude extents of the local region.

        This is used to determine which DEM files are needed to cover the region.

        Returns
        -------
        lon_min : float
            Minimum longitude in degrees.
        lon_max : float
            Maximum longitude in degrees.
        lat_min : float
            Minimum latitude in degrees.
        lat_max : float
            Maximum latitude in degrees.

        """
        R = self.target.radius
        lon0, lat0 = self.local_location
        alpha_deg = np.degrees(self.local_radius / R)

        lat_min = max(-90.0, lat0 - alpha_deg)
        lat_max = min(90.0, lat0 + alpha_deg)

        # If we hit a pole, longitudes span everything
        if abs(lat0) + alpha_deg >= 90.0:
            lon_min, lon_max = -180.0, 180.0
        else:
            half_lon_span = alpha_deg / np.cos(np.radians(lat0))
            half_lon_span = min(half_lon_span, 180.0)

            lon_min = lon0 - half_lon_span
            lon_max = lon0 + half_lon_span
        return float(lon_min), float(lon_max), float(lat_min), float(lat_max)

    def _get_local_dem_data(
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

        if self.dem_file_list is None:
            raise ValueError("No DEM files provided to extract data from.")

        _EXPANSION_BUFFER = 1.05
        _NODATA = np.finfo(np.float32).min
        region_radius = self.local_radius
        half_box_size = np.sqrt(2.0) * region_radius * _EXPANSION_BUFFER
        dst_crs = self.get_crs(
            radius=self.target.radius,
            name=self.target.name,
            location=self.local_location,
        )
        src_list = []
        print("Reading DEM files:")
        try:
            for f in self.dem_file_list:
                print(f"  {f}")
                src_list.append(rasterio.open(f))
        except Exception as e:
            raise RuntimeError(f"Error reading DEM file(s): {e}") from e
        nodata_val = src_list[0].nodata
        if nodata_val is None or np.isnan(nodata_val) or np.abs(nodata_val) < np.abs(_NODATA):
            nodata_val = _NODATA
        target_res = min(s.res[0] for s in src_list)
        dst_width = int(np.ceil(2 * half_box_size / target_res))
        dst_height = dst_width
        dst_transform = Affine(target_res, 0.0, -half_box_size, 0.0, -target_res, half_box_size)
        vrt_list = [
            WarpedVRT(
                src,
                crs=dst_crs,
                transform=dst_transform,
                width=dst_width,
                height=dst_height,
                resampling=Resampling.bilinear,
            )
            for src in src_list
        ]
        # Ensure nodata is respected; if missing, use NaN and float32
        mosaic, transform = merge(
            vrt_list,
            nodata=nodata_val,
            dtype="float32" if np.isnan(nodata_val) else None,
        )

        out_meta = src_list[0].meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": mosaic.shape[1],
                "width": mosaic.shape[2],
                "transform": transform,
                "crs": dst_crs,
                "nodata": nodata_val,
            }
        )

        with MemoryFile() as memfile, memfile.open(**out_meta) as src:
            src.units = src_list[0].units
            src.write(mosaic)

            window = Window(0, 0, src.width, src.height)
            if "KILOMETER" in src.units:
                scale_factor = 1000.0
            else:
                scale_factor = 1.0

            # Read the data within the window
            elevation = src.read(1, window=window) * scale_factor

            # Preserve the window grid and affine for later interpolation
            window_transform = src.window_transform(window)

            # Generate row and column indices
            rows, cols = np.indices(elevation.shape)

            # Compute the x and y coordinates of each pixel in the local CRS
            x_coords, y_coords = rasterio.transform.xy(window_transform, rows, cols, offset="center")
            x_coords = np.array(x_coords).flatten()
            y_coords = np.array(y_coords).flatten()
            elevation = elevation.flatten()

            # Handle nodata values
            mask_nodata = (elevation != src.nodata) & ~np.isnan(elevation)
            mean_elevation = np.mean(elevation[mask_nodata])
            elevation[~mask_nodata] = mean_elevation
            elevation = elevation.astype(np.float32)

            transformer_to_geodetic = Transformer.from_crs(dst_crs, self.crs, always_xy=True)

            # Transform to longitude and latitude
            longitudes, latitudes = transformer_to_geodetic.transform(x_coords, y_coords)
            r_vals = self.compute_distances(
                reference_location=self.local_location, locations=list(zip(longitudes, latitudes, strict=False))
            )
            local_dem_data = {
                "elevation": elevation,
                "mask": r_vals <= region_radius + self.pix / 4,
                "latitudes": latitudes,
                "longitudes": longitudes,
            }

        self._local_dem_data = local_dem_data
        return

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.

        """

        def _interior_distribution():
            self._get_local_dem_data()
            mask = self._local_dem_data["mask"]
            longitudes = np.radians(self._local_dem_data["longitudes"][mask])
            latitudes = np.radians(self._local_dem_data["latitudes"][mask])
            x = np.cos(latitudes) * np.cos(longitudes)
            y = np.cos(latitudes) * np.sin(longitudes)
            z = np.sin(latitudes)

            points = np.column_stack([x, y, z])

            return points

        def _exterior_distribution(theta):
            r = self.radius
            arc_distance = r * theta
            pix_local = self.superdomain_function(arc_distance)
            if (np.pi * r - arc_distance) < pix_local:
                return [], np.pi

            theta_next = theta + pix_local / r

            sin_theta = np.sin(theta)
            cos_theta = np.cos(theta)

            n_phi = max(1, int(round(2 * np.pi * r * sin_theta / pix_local)))
            phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)

            x = sin_theta * np.cos(phi)
            y = sin_theta * np.sin(phi)
            z = np.full_like(phi, cos_theta)

            points = np.column_stack((x, y, z))
            return points.tolist(), theta_next

        print(f"Center of local region: {self.local_location}")
        print(f"Radius of local region: {format_large_units(self.local_radius, quantity='length')}")
        print(f"Local region pixel size: {format_large_units(self.pix, quantity='length')}")

        interior_points = _interior_distribution()
        print(f"Generated {len(interior_points)} points in the local region.")

        theta = (self.local_radius + 0.5 * self.pix) / self.radius
        exterior_points = []
        while theta < np.pi:
            new_points, theta = _exterior_distribution(theta)
            exterior_points.extend(new_points)

        print(f"Generated {len(exterior_points)} points in the superdomain region.")
        exterior_points = self._rotate_point_cloud(exterior_points)

        points = np.vstack([interior_points, exterior_points])
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
        if self._local_dem_data is not None:
            self._add_local_dem_elevation()
            self._add_global_dem_elevation()
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
        # DEM data can come from one of two sources. If we are in the middle of building the surface as part of _generate_face_distribution, then the local DEM data should be stored in the _local_dem_data attribute.
        # If this is not set but we already have a grid generated and stored in the uxgrid attribute, then we need to check if a matching DEM file exists on disk. If neither of these sources are avaailable, then we need to regrid.
        if not regrid and self._local_dem_data is None:
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

    def _add_local_dem_elevation(self):
        """
        Sample the preserved DEM grid at the local face and node centers using bilinear interpolation.

        """
        from scipy.interpolate import LinearNDInterpolator

        if self._local_dem_data is None:
            raise ValueError("Local DEM data is not available for interpolation.")
        lonlat = np.c_[self._local_dem_data["longitudes"], self._local_dem_data["latitudes"]]
        lut2 = LinearNDInterpolator(lonlat, self._local_dem_data["elevation"], fill_value=0.0)
        face_elevations = lut2(np.c_[self.local.face_lon, self.local.face_lat])
        node_elevations = lut2(np.c_[self.local.node_lon, self.local.node_lat])
        elevation = np.concatenate([face_elevations, node_elevations])
        self.local.update_elevation(elevation)
        self._local_dem_data = None  # Clear temporary DEM data storage
        return

    def _add_global_dem_elevation(self):
        """
        Sample the preserved DEM grid at the global face and node centers using bilinear interpolation.

        """
        import rasterio
        from pyproj import Transformer
        from scipy.interpolate import LinearNDInterpolator

        ask_overwrite = self.ask_overwrite
        # If we haven't cached global DEM samples yet, build the exterior lon/lat lists
        if self._global_dem_data is None:
            # This will trigger setting the file automatically from the superdomain_scale_factor if not already set.
            self.superdomain_dem_file = self._superdomain_dem_file

        if self._global_dem_data is None:
            self.ask_overwrite = False
            self._global_dem_data = {}

            # Build index arrays for exterior faces/nodes (sorted for reproducibility)
            exterior_face_ind = np.array(sorted(set(range(self.n_face)) - set(self.local.face_indices)), dtype=int)
            exterior_node_ind = np.array(sorted(set(range(self.n_node)) - set(self.local.node_indices)), dtype=int)

            # Cache exterior point coordinates in geodetic lon/lat (degrees)
            self._global_dem_data["face_indices"] = exterior_face_ind
            self._global_dem_data["node_indices"] = exterior_node_ind
            self._global_dem_data["face_longitudes"] = np.asarray(self.face_lon)[exterior_face_ind]
            self._global_dem_data["face_latitudes"] = np.asarray(self.face_lat)[exterior_face_ind]
            self._global_dem_data["node_longitudes"] = np.asarray(self.node_lon)[exterior_node_ind]
            self._global_dem_data["node_latitudes"] = np.asarray(self.node_lat)[exterior_node_ind]

            # Read and sample the global DEM at those exterior points.
            print("Reading global DEM file:")
            print(f"  {self.superdomain_dem_file}")
            with rasterio.open(self.superdomain_dem_file) as dataset:
                if dataset.crs is None:
                    raise ValueError("Global DEM has no CRS; cannot sample elevations.")

                # Transform from geodetic lon/lat (self.crs) into the DEM's CRS
                to_dem = Transformer.from_crs(self.crs, dataset.crs, always_xy=True)

                # Determine unit scaling (match local DEM logic)
                units = getattr(dataset, "units", "") or ""
                if isinstance(units, (list, tuple, ArrayLike)):
                    units = units[0] if len(units) else ""
                scale_factor = 1000.0 if "KILOMETER" in str(units).upper() else 1.0

                nodata = dataset.nodata

                def _sample(lons_deg: NDArray, lats_deg: NDArray) -> NDArray:
                    xs, ys = to_dem.transform(lons_deg, lats_deg)
                    # rasterio.sample expects an iterable of (x, y)
                    samples = np.fromiter(
                        (v[0] for v in dataset.sample(zip(xs, ys, strict=False))),
                        dtype=np.float32,
                        count=len(lons_deg),
                    )
                    samples = samples.astype(np.float32) * scale_factor

                    # Replace nodata / NaN with mean of valid samples
                    if nodata is None or np.isnan(nodata):
                        valid = np.isfinite(samples)
                    else:
                        valid = np.isfinite(samples) & (samples != nodata)

                    if np.any(valid):
                        mean_val = float(np.mean(samples[valid]))
                    else:
                        mean_val = 0.0
                    samples[~valid] = mean_val
                    return samples

                face_samples = _sample(self._global_dem_data["face_longitudes"], self._global_dem_data["face_latitudes"])
                node_samples = _sample(self._global_dem_data["node_longitudes"], self._global_dem_data["node_latitudes"])

        # Create full-sized elevation arrays and insert the sampled values at the exterior indices
        face_elevations = np.zeros(self.n_face, dtype=np.float32)
        node_elevations = np.zeros(self.n_node, dtype=np.float32)
        face_elevations[self._global_dem_data["face_indices"]] = face_samples
        node_elevations[self._global_dem_data["node_indices"]] = node_samples

        # Apply elevations to the global (superdomain) surface
        elevation = np.concatenate([face_elevations, node_elevations])
        self.update_elevation(elevation)

        # Save so we can reload later without re-downloading / re-sampling the DEM
        self.save(filename=self._dem_output_file, skip_actions=True)
        self._global_dem_data = None  # Clear temporary DEM data storage
        self.ask_overwrite = ask_overwrite
        return

    @parameter
    def dem_file_list(self) -> list[str] | None:
        """
        The list of files to use for the DEM data in the high resolution local region.

        """
        return self._dem_file_list

    @dem_file_list.setter
    def dem_file_list(self, value: list[str] | None):
        if value is None:
            if self._pix is None:
                # Compute a reasonable default resolution that will contain approximately 1e6 faces based on the local radius
                self._pix = np.sqrt(np.pi * self.local_radius**2 / _DEFAULT_N_FACES_LOCAL)
            lon_min, lon_max, lat_min, lat_max = self.get_location_extents(self.local_location, self.local_radius)
            value, self._pix = DataComposer.get_lola_dem_file_list(
                pix=self._pix, lat_range=(lat_min, lat_max), lon_range=(lon_min, lon_max)
            )
        elif isinstance(value, list):
            if not all(isinstance(f, (str, Path)) for f in value):
                raise ValueError("All items in 'dem_file_list' must be strings or Path objects.")
        elif not isinstance(value, (str, Path)):
            raise TypeError("'dem_file_list' must be a list of strings or Path objects, or None.")

        self._dem_file_list = value

        # Set the pixel size based on the provided files. We take the highest resolution (smallest pixel size) among the files.

        pixvals = []
        for f in self._dem_file_list:
            with rasterio.open(f) as src:
                pixvals.append(src.res[0])
        self._pix = min(pixvals)

        return

    @parameter
    def superdomain_dem_file(self) -> str | Path | None:
        """
        The list of files to use for the DEM data in the superdomain.

        """
        return self._superdomain_dem_file

    @superdomain_dem_file.setter
    def superdomain_dem_file(self, value: str | Path | None):
        if value is None:
            # If the superdomain has not been set yet, we will defer setting this until later
            if self.superdomain_scale_factor is None:
                return
            min_global_pix = np.pi / 180.0 * self.target.radius / 128.0
            sdpix = self.superdomain_scale_factor * self.pix / 10.0
            if sdpix < min_global_pix:
                sdpix = min_global_pix
            self._superdomain_dem_file = DataComposer.get_lola_dem_file_list(pix=sdpix, lat_range=(-90, 90), lon_range=(-180, 180))[
                0
            ][0]
            return
        if not isinstance(value, (str, Path)):
            raise TypeError("'superdomain_dem_file' must be a strings or Path objects, or None.")

        self._superdomain_dem_file = value

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.

        """
        return super()._hashvars + [self._dem_file_list]
