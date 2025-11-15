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
from cratermaker.components.surface.hireslocal import HiResLocalSurface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import (
    format_large_units,
    parameter,
)


@Surface.register("datasurface")
class DataSurface(HiResLocalSurface):
    """
    A Surface subclass that generates a local region using DEM data.

    Currently implements only lunar LOLA data from the PDS. Other data sources may be added in the future.

    Parameters
    ----------
        The radius of the local region in meters.
    local_location : PairOfFloats
        The longitude and latitude of the location in degrees.
    superdomain_scale_factor : FloatLike, optional
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters
        that are modeled outside the local region are those whose ejecta could just reach the boundary. If not provide, then it will
        be computed based on a provided (or default) scaling and morphology model
    target : Target, optional
        The target body or name of a known target body for the impact simulation. If none provide, it will be either set to the default,
        or extracted from the scaling model if it is provied
    reset : bool, optional
        Flag to indicate whether to reset the surface. Default is True.
    regrid : bool, optional
        Flag to indicate whether to regrid the surface. Default is False.
    ask_overwrite : bool, optional
        If True, prompt the user for confirmation before deleting files. Default is False.
    simdir : str | Path
        The main project simulation directory. Default is the current working directory if None.
    pix : FloatLike | None, optional
        The approximate face size inside the local region in meters. This will be used to determine the target resolution of the DEM data to be used. The actual resolution may be different based on the available DEM data. Note that if you provide a list of DEM files using the `dem_file_list` parameter, this value will be ignored. If None, is set, and no file(s) are provided, a default resolution that creates approximately 1e6 faces in the local region will be used.
    local_radius : FloatLike
    dem_file_list : list of str | Path, optional
        A list of one or more DEM files or urls links to files to be used for the local region. If not provided, a set of files will be determined based on the location and radius of the region. It will chose a set of files that fully covers the local region with approximately approximately 1e6 faces.
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
        ask_overwrite: bool = False,
        simdir: str | Path | None = None,
        pix: FloatLike | None = None,
        dem_file_list: list[str] | list[Path] | None = None,
        **kwargs: Any,
    ):
        try:
            import rasterio
        except ImportError:
            warn(
                "rasterio is not installed. This is required for 'datasurface.' On some platforms, you may need to install GDAL first before installing rasterio.",
                stacklevel=2,
            )
            return
        object.__setattr__(self, "_local", None)

        # Temporary storage for DEM data during initialization. This will be cleared after the elevation points are set.
        object.__setattr__(self, "_dem_data", None)

        # The location of the saved DEM data that can be retrieved later
        object.__setattr__(self, "_dem_output_file", None)
        object.__setattr__(self, "_dem_file_list", None)

        super(HiResLocalSurface, self).__init__(target=target, simdir=simdir, **kwargs)
        if dem_file_list is None and self.target.name != "Moon":
            raise ValueError("DataSurface currently only supports the Moon as a target if 'dem_file_list' is not provided.")
        self._output_file_pattern += [f"local_{self._output_file_prefix}*.{self._output_file_extension}"]
        self._dem_output_file = f"dem_data.{self._output_file_extension}"

        # Set the attributes directly to avoid triggering checks before the pix value is set
        self._local_radius = local_radius
        self._local_location = local_location
        self._pix = pix
        self.dem_file_list = dem_file_list
        super().__init__(
            pix=self.pix,
            local_radius=local_radius,
            local_location=local_location,
            superdomain_scale_factor=superdomain_scale_factor,
            target=self.target,
            reset=reset,
            regrid=regrid,
            ask_overwrite=ask_overwrite,
            simdir=self.simdir,
            **kwargs,
        )

        return

    def _get_lola_dem_file_list(self):
        """
        Determine the set of LOLA DEM files needed to cover the local region if none are provided by the user.
        """
        lon_min, lon_max, lat_min, lat_max = self._get_location_extents()
        if lat_max < 60 and lat_min > -60:  # Use cylindrical files only for mid-latitudes
            target_pds_resolution = np.pi / 180.0 * self.target.radius / self.pix
            return self._get_lola_cylindrical_files_from_pds(resolution=target_pds_resolution)
        else:
            return self._get_lola_polar_files_from_pds(self.pix, latrange=(lat_min, lat_max))

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

    @staticmethod
    def _get_lola_cylindrical_url_from_pds(pds_file_resolution: int, location: PairOfFloats, boundary_offset=(0, 0)) -> Path:
        """
        Retrieve the appropriate LOLA DEM file url for a given location and resolution from the PDS.

        Parameters
        ----------
        pds_file_resolution : int
            DEM resolution in meters per pixel.
        location : PairOfFloats,
            The longitude and latitude of the location in degrees.
        boundary_offset : tuple of int, optional
            Offset to apply to tile index to access neighboring tiles. Default is (0, 0). This is used when the local region crosses one or more tile boundaries.

        Returns
        -------
        url: str
            Link to the DEM file

        Notes
        -----
        This is meant to be used for mid-latitude regions on the Moon between -60 and 60 degrees latitude.
        For resolutions of 128, 256, and 512 pix/deg, this will return the SLDEM2015 datasets.
        """
        lola_cylindrical_url = (
            "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/float_img/"
        )
        sldem_global_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/sldem2015/global/float_img/"
        sldem_tile_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/sldem2015/tiles/float_img/"

        def _get_pds_file_suffix(pds_file_resolution, location, boundary_offset):
            if location[0] < 0:
                location = (location[0] + 360, location[1])
            if pds_file_resolution == 512:
                dlon = 45
                dlat = 30
            elif pds_file_resolution == 256:
                dlon = 120
                dlat = 60

            latdir = "n" if location[1] + boundary_offset[1] * dlat > 0 else "s"
            latval = np.abs(location[1] / dlat) + boundary_offset[1]
            if latdir == "n":
                latlo = int(np.ceil(latval - 1) * dlat)
                lathi = int(latlo + dlat)
            else:
                lathi = int(np.floor(latval + 1) * dlat)
                latlo = int(lathi - dlat)

            lonval = np.abs(location[0] / dlon) + boundary_offset[0]

            # Pole crossing, so flip the longitudes by 180 degrees
            if lathi > 90 or latlo > 90:
                lathi = 90
                latlo = lathi - dlat
                lonval += 180 / dlon

            lonlo = int(np.floor(lonval) * dlon) % 360
            lonhi = int(lonlo + dlon)

            # Check to make sure the lo and hi values are in the right direction
            if (latdir.upper() == "N" and latlo > lathi) or (latdir.upper() == "S" and latlo < lathi):
                tmp = latlo
                latlo = lathi
                lathi = tmp

            if pds_file_resolution == 512:
                return f"{latlo:02d}{latdir.lower()}_{lathi:02d}{latdir.lower()}_{lonlo:03d}_{lonhi:03d}"
            elif pds_file_resolution == 256:
                return f"{latlo:d}{latdir.lower()}_{lathi:d}{latdir.lower()}_{lonlo:03d}_{lonhi:03d}"

        if pds_file_resolution >= 256:
            url = f"{sldem_tile_url}sldem2015_{pds_file_resolution:d}_{_get_pds_file_suffix(pds_file_resolution, location, boundary_offset)}_float.xml"
        elif pds_file_resolution == 128:
            url = f"{sldem_global_url}sldem2015_{pds_file_resolution:d}_60s_60n_000_360_float.xml"
        else:
            url = f"{lola_cylindrical_url}ldem_{pds_file_resolution:d}_float.xml"

        return url

    def _get_lola_cylindrical_files_from_pds(self, resolution) -> list[str] | str:
        """
        Retrieve the appropriate cylindrically projected LOLA DEM file url or list of urls for a given location and resolution from the PDS.

        This will determine which, if any, boundaries are crossed and will return a list of files to cover the entire local region.

        Parameters
        ----------
        resolution : FloatLike
            Requested resolution in degrees per pixel. The closest available resolution will be used.

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

        AVAILABLE_RESOLUTIONS = [4, 16, 64, 128, 256, 512]  #  pix / deg
        diffs = [abs(resolution - res) for res in AVAILABLE_RESOLUTIONS]
        pds_file_resolution = AVAILABLE_RESOLUTIONS[np.argmin(diffs)]

        # First, retrive the file for the centerpoint:
        filelist = [self._get_lola_cylindrical_url_from_pds(pds_file_resolution, self.local_location)]
        if pds_file_resolution < 256:
            return filelist  # These files cover the entire globe, no need to determine if boundaries are crossed

        lon_min, lon_max, lat_min, lat_max = self._get_location_extents()
        combo = [(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_min), (lon_max, lat_max)]

        for loc in combo:
            f = self._get_lola_cylindrical_url_from_pds(pds_file_resolution, loc)
            if f not in filelist:
                filelist.append(f)

        return filelist

    def _get_lola_polar_files_from_pds(self, resolution, latrange: tuple[float, float]) -> list[str] | str:
        """
        Retrieve the appropriate polar projected LOLA DEM file url or list of urls for a given location and resolution from the PDS.

        Parameters
        ----------
        resolution : FloatLike
            Requested resolution in meters per pixel. The closest available resolution will be used.
        latrange : tuple of float
            The (min_lat, max_lat) in degrees of the local region.
        """
        import rasterio

        src_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/polar/float_img/"

        AVAILABLE_RESOLUTIONS = [3000, 1000, 400, 240, 200, 120, 100, 80, 60, 40, 30, 20, 10, 5]
        MIN_LAT = [50, 50, 45, 60, 45, 60, 45, 80, 60, 80, 75, 80, 85, 87.5]

        if latrange[0] > 0:
            pole = "n"
            min_lat_index = 0  # in the northern hemisphere, the minimum latitude defines coverage
        else:
            pole = "s"
            min_lat_index = 1  # in the southern hemisphere, the maximum latitude defines coverage
        combo = [
            (res, minlat)
            for res, minlat in zip(AVAILABLE_RESOLUTIONS, MIN_LAT, strict=False)
            if abs(latrange[min_lat_index]) >= minlat
        ]
        if len(combo) == 0:
            raise ValueError("No available LOLA polar DEM files cover the requested latitude range.")

        valid_resolutions, valid_min_lat = zip(*combo, strict=False)

        diffs = [abs(resolution - res) for res in valid_resolutions]
        pds_file_resolution = valid_resolutions[np.argmin(diffs)]
        pds_lat_min = valid_min_lat[np.argmin(diffs)]
        filename = f"ldem_{pds_lat_min:d}{pole}_{pds_file_resolution:d}m"
        url = f"{src_url}{filename}_float.xml"

        self._pix = pds_file_resolution
        return [url]

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

        if self.dem_file_list is None:
            raise ValueError("No DEM files provided to extract data from.")

        _EXPANSION_BUFFER = 1.05
        _NODATA = -1.0e-9
        region_radius = self.local_radius
        box_size = 2 * np.sqrt(2.0) * region_radius * _EXPANSION_BUFFER
        half_box_size = box_size / 2
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
        if nodata_val is None or np.isnan(nodata_val):
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
                dst_nodata=nodata_val,
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
            r_vals = self.calculate_distance(self.local_location, list(zip(longitudes, latitudes, strict=False)))
            dem_data = {
                "elevation": elevation,
                "mask": r_vals <= region_radius + self.pix / 4,
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

        def _interior_distribution():
            self._get_dem_data()
            mask = self._dem_data["mask"]
            longitudes = np.radians(self._dem_data["longitudes"][mask])
            latitudes = np.radians(self._dem_data["latitudes"][mask])
            x = self.target.radius * np.cos(latitudes) * np.cos(longitudes)
            y = self.target.radius * np.cos(latitudes) * np.sin(longitudes)
            z = self.target.radius * np.sin(latitudes)

            points = np.column_stack([x, y, z]) / self.target.radius

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

    def reset(self, ask_overwrite: bool = False, **kwargs: Any) -> None:
        """
        Reset the surface to its initial state.

        Parameters
        ----------
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before deleting files. Default is False.
        **kwargs : Any
            Additional keyword arguments for subclasses.
        """
        super().reset(ask_overwrite=ask_overwrite, **kwargs)
        if self._dem_data is not None:
            self._add_dem_elevation()
        elif (self.output_dir / self._dem_output_file).exists():
            with xr.open_dataset(self.output_dir / self._dem_output_file) as ds:
                self.local.update_elevation(ds.isel(time=0).face_elevation)
                self.local.update_elevation(ds.isel(time=0).node_elevation)

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
                        or ds.face_elevation.shape[-1] != self.local.n_face
                        or ds.node_elevation.shape[-1] != self.local.n_node
                    ):
                        regrid = True

        return super()._regrid_if_needed(regrid=regrid, **kwargs)

    def _add_dem_elevation(self):
        """
        Sample the preserved DEM grid at the face and node centers using bilinear interpolation.

        """
        from scipy.interpolate import LinearNDInterpolator

        lonlat = np.c_[self._dem_data["longitudes"], self._dem_data["latitudes"]]
        lut2 = LinearNDInterpolator(lonlat, self._dem_data["elevation"], fill_value=0.0)
        face_elevations = lut2(np.c_[self.local.face_lon, self.local.face_lat])
        node_elevations = lut2(np.c_[self.local.node_lon, self.local.node_lat])

        self.local.update_elevation(face_elevations)
        self.local.update_elevation(node_elevations)

        # Now save the surface to a file that we can reload later if we want to avoid re-downloading the DEM data
        self.local.save(filename=self._dem_output_file)
        self._dem_data = None  # Clear temporary DEM data storage
        return

    @parameter
    def dem_file_list(self) -> int | None:
        """
        The list of files to use for the DEM data.
        """
        return self._dem_file_list

    @dem_file_list.setter
    def dem_file_list(self, value: int | None):
        if value is None:
            if self._pix is None:
                # Compute a reasonable default resolution based on the local radius
                self._pix = np.sqrt(np.pi * self.local_radius**2 / 1e6)
            self._dem_file_list = self._get_lola_dem_file_list()
            return
        if isinstance(value, list):
            if not all(isinstance(f, (str, Path)) for f in value):
                raise ValueError("All items in 'dem_file_list' must be strings or Path objects.")
        elif not isinstance(value, (str, Path)):
            raise TypeError("'dem_file_list' must be a list of strings or Path objects, or None.")

        self._dem_file_list = value

        return
