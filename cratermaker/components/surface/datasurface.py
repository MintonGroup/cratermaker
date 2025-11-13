from __future__ import annotations

from pathlib import Path
from typing import Any

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
    data_file : str | Path
        The path or url to the DEM data file to be used for the local region
    local_radius : FloatLike
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
    pds_file_resolution : int, optional
        The resolution of the DEM file to download from the PDS in meters per pixel. Valid values are [4,16,64,128,256,512]. If not provide, one will be chosen based on the local radius
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
        pds_file_resolution: int | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_local", None)
        object.__setattr__(self, "_demfile", None)
        object.__setattr__(
            self, "_dem_data", None
        )  # Temporary storage for DEM data during initialization. This will be cleared after the elevation points are set.
        object.__setattr__(self, "_valid_pds_file_resolutions", [4, 16, 64, 128, 256, 512])
        object.__setattr__(self, "_pds_file_resolution", None)

        super(HiResLocalSurface, self).__init__(target=target, simdir=simdir, **kwargs)
        if self.target.name != "Moon":
            raise ValueError("DataSurface currently only supports the Moon as a target.")
        self._output_file_pattern += [f"local_{self._output_file_prefix}*.{self._output_file_extension}"]
        self._demfile = f"dem_data.{self._output_file_extension}"

        # Set the attributes directly to avoid triggering checks before the pix value is set
        self._local_radius = local_radius
        self._local_location = local_location
        self.pds_file_resolution = pds_file_resolution
        pix = np.pi * self.target.radius / (180 * self.pds_file_resolution)
        super().__init__(
            pix=pix,
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

    def _get_lola_dem_file_from_pds(self, pds_file_resolution: int, location: PairOfFloats, boundary_offset=(0, 0)) -> Path:
        """
        Generate and download the appropriate LOLA DEM file for a given location and resolution.

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
        """
        if pds_file_resolution not in self._valid_pds_file_resolutions:
            raise ValueError("Invalid resolution.")

        src_url = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/float_img/"

        def _get_pds_file_suffix(pds_file_rsolution, location, boundary_offset):
            if pds_file_resolution == 1024:
                dlon = 30
                dlat = 15
            elif pds_file_resolution == 512:
                dlon = 90
                dlat = 45
            elif pds_file_resolution == 256:
                dlon = 180
                dlat = 90

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

            return f"{latlo:02d}{latdir.lower()}_{lathi:02d}{latdir.lower()}_{lonlo:03d}_{lonhi:03d}"

        if pds_file_resolution >= 256:
            filename = f"ldem_{pds_file_resolution:d}_{_get_pds_file_suffix(pds_file_resolution, location, boundary_offset)}"
        else:
            filename = f"ldem_{pds_file_resolution:d}"

        url = f"{src_url}{filename}_float.xml"
        return url

    def _get_dem_data(
        self,
        pds_file_resolution,
    ):
        """
        Retrieve DEM data for a specified crater location and extent.

        Parameters
        ----------
        pds_file_resolution : int
            Resolution of the NAC file in meters per pixel.

        Returns
        -------
        dem : dict
            Dictionary containing x, y, z coordinates and elevation.
        pix : float
            Average pixel size in meters.
        """
        isboundary = False

        # The high resolution DTM files need a slightly different location format
        location = self.local_location
        if pds_file_resolution < 256 and self.local_location[0] > 180:
            location = (self.local_location[0] - 360, self.local_location[1])
        elif pds_file_resolution >= 256 and self.local_location[0] < 0:
            location = (self.local_location[0] + 360, self.local_location[1])

        compute_boundary = pds_file_resolution >= 256
        filename = self._get_lola_dem_file_from_pds(pds_file_resolution, location)

        dem_data, isboundary, boundary_offsets = self._get_dem_from_file(
            location=location,
            filename=filename,
            compute_boundary=compute_boundary,
        )

        # Handle the case where the DEM is split over multiple files across a boundary
        if isboundary:
            filelist = [filename]
            for boundary_offset in boundary_offsets:
                filelist.append(self._get_lola_dem_file_from_pds(pds_file_resolution, location, boundary_offset))
            dem_data, _, _ = self._get_dem_from_file(
                location=location,
                filelist=filelist,
                compute_boundary=False,
            )
        self._dem_data = dem_data
        return

    def _get_dem_from_file(
        self,
        location,
        region_radius=None,
        filename=None,
        filelist=None,
        compute_boundary=True,
    ):
        """
        Read and extract DEM data from one or more files for a specified region.

        Parameters
        ----------
        location : tuple
            (longitude, latitude) in degrees.
        region_radius : float
            Radius of the region to extract in meters.
        filename : str, optional
            A single DEM file to read.
        filelist : list of str, optional
            List of DEM files to merge and read from.
        compute_boundary : bool, optional
            Whether to check for boundary overlap in the DEM. This is needed if the DEM straddles the boundary between two or more files

        Returns
        -------
        dem : dict or None
            Dictionary with x, y, z coordinates and elevation, or None if boundary case.
        isboundary : bool
            Whether the DEM window required boundary adjustment.
        boundary : list
            List indicating the directions of boundary overlap, if any.
        """
        import rasterio
        from affine import Affine
        from pyproj import Transformer
        from rasterio.io import MemoryFile
        from rasterio.merge import merge
        from rasterio.vrt import WarpedVRT
        from rasterio.warp import Resampling
        from rasterio.windows import from_bounds

        _EXPANSION_BUFFER = 1.05
        if region_radius is None:
            region_radius = self.local_radius
        box_size = 2 * np.sqrt(2.0) * region_radius * _EXPANSION_BUFFER
        isboundary = False
        boundary = [None, None]
        dst_crs = self.get_crs(
            radius=self.target.radius,
            name=self.target.name,
            location=location,
        )

        def _aeqd_window_bounds_in_src_crs(location, box_size_m, src_crs, dst_crs):
            """
            Compute (minx, miny, maxx, maxy) in the *source* CRS that tightly enclose the AEQD square window centered at `location` with side length `box_size_m`.

            Works for both geographic (lon/lat in degrees) and projected (linear units) source CRSs.
            """
            from pyproj import Transformer

            # Half the AEQD square (in meters, AEQD space)
            half = 0.5 * box_size_m
            left_m, right_m = -half, half
            bottom_m, top_m = -half, half

            # Sample the perimeter densely to capture extremal extents after inverse projection
            N = 64
            xs = np.concatenate(
                [
                    np.linspace(left_m, right_m, N),  # bottom edge x
                    np.full(N, right_m),  # right edge x
                    np.linspace(right_m, left_m, N),  # top edge x
                    np.full(N, left_m),  # left edge x
                ]
            )
            ys = np.concatenate(
                [
                    np.full(N, bottom_m),  # bottom edge y
                    np.linspace(bottom_m, top_m, N),  # right edge y
                    np.full(N, top_m),  # top edge y
                    np.linspace(top_m, bottom_m, N),  # left edge y
                ]
            )

            transformer = Transformer.from_crs(dst_crs, src_crs, always_xy=True)
            u, v = transformer.transform(xs, ys)

            # Determine whether the source CRS is geographic (degrees) or projected (linear)
            try:
                src_is_geographic = src_crs.is_geographic
            except Exception:
                # Fall back to assuming projected if CRS parsing fails
                src_is_geographic = False

            if src_is_geographic:
                # Normalize longitudes to [0, 360) and compute minimal covering arc on the circle
                lon_raw = np.asarray(u, dtype=float)
                lat = np.asarray(v, dtype=float)

                lon = np.mod(lon_raw, 360.0)
                lon = np.where(lon < 0.0, lon + 360.0, lon)

                # Sort and find the largest gap; the complement is the minimal arc covering all points
                lon_sorted = np.sort(lon)
                if lon_sorted.size == 0:
                    # Fallback: degenerate case
                    lon_min = lon_max = 0.0
                else:
                    diffs = np.diff(np.r_[lon_sorted, lon_sorted[0] + 360.0])
                    j = int(np.argmax(diffs))
                    # Minimal covering arc goes from next point after the largest gap to the point at the start of the gap
                    start = lon_sorted[(j + 1) % lon_sorted.size]
                    end = lon_sorted[j]
                    lon_min, lon_max = (
                        start,
                        end,
                    )  # Note: if start > end, this interval wraps 0°

                lat_min = float(np.min(lat))
                lat_max = float(np.max(lat))

                # Tiny pads for safety in angular units
                pad_lon = 1e-6 * max(1.0, (lon_max - lon_min) % 360.0)
                pad_lat = 1e-6 * max(1.0, abs(lat_max - lat_min))

                # Apply padding while keeping results in [0, 360)
                left = (lon_min - pad_lon) % 360.0
                right = (lon_max + pad_lon) % 360.0
                bottom = lat_min - pad_lat
                top = lat_max + pad_lat
            else:
                # Projected/linear units: use corner sampling and guard against seam/degenerate warps
                from pyproj import Transformer

                half = 0.5 * box_size_m
                corners_dst = np.array(
                    [
                        [-half, -half],
                        [half, -half],
                        [half, half],
                        [-half, half],
                    ],
                    dtype=float,
                )

                # Transform corners from AEQD (dst) -> source CRS
                to_src = Transformer.from_crs(dst_crs, src_crs, always_xy=True)
                x_c, y_c = to_src.transform(corners_dst[:, 0], corners_dst[:, 1])
                x_c = np.asarray(x_c, dtype=float)
                y_c = np.asarray(y_c, dtype=float)

                x_min = float(np.min(x_c))
                x_max = float(np.max(x_c))
                y_min = float(np.min(y_c))
                y_max = float(np.max(y_c))

                width = x_max - x_min
                height = y_max - y_min

                # Compute the expected side length in meters in source CRS near the center by
                # transforming two small AEQD steps and measuring the differential scaling.
                # Also compute the center in the source CRS for a safe fallback.
                to_src_small = to_src  # reuse
                # Small differential step (meters) to estimate local Jacobian scale
                d = min(1.0, half)
                x0_s, y0_s = to_src_small.transform(0.0, 0.0)
                x_dx, y_dx = to_src_small.transform(d, 0.0)
                x_dy, y_dy = to_src_small.transform(0.0, d)
                sx = max(1e-9, abs(x_dx - x0_s) / d)
                sy = max(1e-9, abs(y_dy - y0_s) / d)
                # Estimated expected box size mapped into source CRS
                exp_width = 2.0 * sx * half
                exp_height = 2.0 * sy * half

                # Heuristic: if the transformed box is unreasonably elongated or orders of magnitude
                # larger than the expected local mapping, fall back to a centered square around (0,0) in AEQD mapped to src.
                bad_aspect = (height > 0 and (width / height > 10.0)) or (width == 0 and height > 0)
                bad_scale = (width > 4.0 * max(exp_width, 1.0)) or (height > 4.0 * max(exp_height, 1.0))

                if bad_aspect or bad_scale:
                    # Fallback: construct a square in the source CRS centered at the mapped AEQD origin
                    # with side lengths matching the local linearized scale.
                    cx, cy = x0_s, y0_s
                    x_min = cx - exp_width * 0.5
                    x_max = cx + exp_width * 0.5
                    y_min = cy - exp_height * 0.5
                    y_max = cy + exp_height * 0.5

                pad_x = 1e-6 * max(1.0, abs(x_max - x_min))
                pad_y = 1e-6 * max(1.0, abs(y_max - y_min))
                left, bottom, right, top = (
                    x_min - pad_x,
                    y_min - pad_y,
                    x_max + pad_x,
                    y_max + pad_y,
                )

            return left, bottom, right, top

        def _merge_with_bounds_wrap(src_list, bounds, nodata_val=None, dtype=None):
            """Merge sources within `bounds`, handling 0/360 wrap only for geographic CRSs.

            If the source CRS is projected (linear units), perform a normal merge.
            """
            left, bottom, right, top = bounds

            # Detect if the source CRS is geographic (lon/lat)
            try:
                src_is_geographic = src_list[0].crs.is_geographic
            except Exception:
                src_is_geographic = False

            if not src_is_geographic:
                # No wrap logic in projected CRSs
                return merge(
                    src_list,
                    bounds=bounds,
                    nodata=nodata_val,
                    dtype=dtype,
                    res=src_list[0].res,
                )

            # Geographic: handle potential 0..360 crossing
            if left <= right:
                return merge(
                    src_list,
                    bounds=bounds,
                    nodata=nodata_val,
                    dtype=dtype,
                    res=src_list[0].res,
                )

            # Wrap case for 0..360 domain
            bounds1 = (left, bottom, 360.0, top)
            bounds2 = (0.0, bottom, right, top)

            mosa_parts = []
            trans_parts = []
            for b in (bounds1, bounds2):
                mosa, trans = merge(
                    src_list,
                    bounds=b,
                    nodata=nodata_val,
                    dtype=dtype,
                    res=src_list[0].res,
                )
                mosa_parts.append(mosa)
                trans_parts.append(trans)

            # Write partials to in-memory datasets and merge again to obtain a single mosaic
            mems = []
            dsets = []
            try:
                base = src_list[0]
                out_dtype = dtype or base.meta.get("dtype", mosa_parts[0].dtype.name)
                for idx, (mosa, trans) in enumerate(zip(mosa_parts, trans_parts, strict=False)):
                    # Shift the second chunk by +360° in x so both parts live in a continuous domain
                    trans_to_use = trans
                    if idx == 1:
                        trans_to_use = Affine(trans.a, trans.b, trans.c + 360.0, trans.d, trans.e, trans.f)

                    meta = base.meta.copy()
                    meta.update(
                        driver="GTiff",
                        count=mosa.shape[0],
                        height=mosa.shape[1],
                        width=mosa.shape[2],
                        transform=trans_to_use,
                        crs=base.crs,
                        dtype=out_dtype,
                        nodata=nodata_val,
                    )
                    mem = MemoryFile()
                    ds = mem.open(**meta)
                    ds.write(mosa.astype(out_dtype, copy=False))
                    mems.append(mem)
                    dsets.append(ds)

                mosa_full, trans_full = merge(dsets, nodata=nodata_val, dtype=out_dtype, res=base.res)
            finally:
                import contextlib

                for ds in dsets:
                    with contextlib.suppress(Exception):
                        ds.close()
                for mem in mems:
                    with contextlib.suppress(Exception):
                        mem.close()
            return mosa_full, trans_full

        with MemoryFile() as memfile:
            if filelist is not None:
                print(f"Merging files: {filelist}")
                src_list = []
                for filename in filelist:
                    src = rasterio.open(filename)
                    src_list.append(src)

                box_size = 2 * np.sqrt(2.0) * region_radius * _EXPANSION_BUFFER
                bounds = _aeqd_window_bounds_in_src_crs(location, box_size, src_list[0].crs, dst_crs)

                # Ensure nodata is respected; if missing, use NaN and float32
                nodata_val = src_list[0].nodata
                if nodata_val is None:
                    nodata_val = np.float32(np.nan)
                    mosaic, transform = _merge_with_bounds_wrap(src_list, bounds, nodata_val=nodata_val, dtype="float32")
                    mosaic = mosaic.astype("float32", copy=False)
                else:
                    mosaic, transform = _merge_with_bounds_wrap(src_list, bounds, nodata_val=nodata_val)

                out_meta = src.meta.copy()
                out_meta.update(
                    {
                        "driver": "GTiff",
                        "height": mosaic.shape[1],
                        "width": mosaic.shape[2],
                        "transform": transform,
                        "crs": src.crs,
                        "nodata": nodata_val,
                    }
                )
                with memfile.open(**out_meta) as dst:
                    dst.units = src_list[0].units
                    dst.write(mosaic)

                cm = memfile.open()
            else:
                print(f"Opening DEM file: {filename}")
                cm = rasterio.open(filename)
            with cm as src:
                # Desired box size in meters
                half_box_size = (
                    _EXPANSION_BUFFER * box_size / 2
                )  # Make the box slightly larger than the desired size, as it will be truncated into a square of the correct size later
                with WarpedVRT(src, crs=dst_crs, resampling=Resampling.cubic) as vrt:
                    # Compute the window in the destination CRS
                    x_min = -half_box_size
                    x_max = half_box_size
                    y_min = -half_box_size
                    y_max = half_box_size

                    # Compute the window in pixel coordinates
                    window_orig = from_bounds(x_min, y_min, x_max, y_max, vrt.transform)

                    if compute_boundary:
                        window = window_orig.intersection(rasterio.windows.Window(0, 0, vrt.width, vrt.height))
                        if window != window_orig:
                            isboundary = True
                            boundary_offsets = []
                            if window.col_off != window_orig.col_off:
                                boundary_offsets.append((-1, 0))
                            elif window.width != window_orig.width:
                                boundary_offsets.append((1, 0))
                            if window.row_off != window_orig.row_off:
                                boundary_offsets.append((0, 1))
                            elif window.height != window_orig.height:
                                boundary_offsets.append((0, -1))
                            if len(boundary_offsets) == 2:  # corner case
                                boundary_offsets.append(
                                    (
                                        boundary_offsets[0][0] + boundary_offsets[1][0],
                                        boundary_offsets[0][1] + boundary_offsets[1][1],
                                    )
                                )
                            return None, isboundary, boundary_offsets

                    else:
                        window = window_orig

                    if "KILOMETER" in src.units:
                        scale_factor = 1000.0
                    else:
                        scale_factor = 1.0
                    # Read the data within the window
                    elevation = vrt.read(1, window=window) * scale_factor

                    # Preserve the window grid and affine for later interpolation
                    window_transform = vrt.window_transform(window)

                    # Generate row and column indices
                    rows, cols = np.indices(elevation.shape)

                    # Compute the x and y coordinates of each pixel in the local CRS
                    x_coords, y_coords = rasterio.transform.xy(window_transform, rows, cols, offset="center")
                    x_coords = np.array(x_coords).flatten()
                    y_coords = np.array(y_coords).flatten()
                    r_vals = np.sqrt(x_coords**2 + y_coords**2)
                    elevation = elevation.flatten()

                    # Handle nodata values
                    mask_nodata = (elevation != vrt.nodata) & ~np.isnan(elevation)
                    mean_elevation = np.mean(elevation[mask_nodata])
                    elevation[~mask_nodata] = mean_elevation
                    elevation = elevation.astype(np.float32)

                    transformer_to_geodetic = Transformer.from_crs(dst_crs, self.crs, always_xy=True)

                    # Transform to longitude and latitude
                    longitudes, latitudes = transformer_to_geodetic.transform(x_coords, y_coords)
                    dem_data = {
                        "elevation": elevation,
                        "mask": r_vals <= region_radius + self.pix / 4,
                        "latitudes": latitudes,
                        "longitudes": longitudes,
                    }

        return dem_data, isboundary, boundary

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """

        def _interior_distribution():
            self._get_dem_data(self.pds_file_resolution)
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
        elif (self.output_dir / self._demfile).exists():
            with xr.open_dataset(self.output_dir / self._demfile) as ds:
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
            if not (self.output_dir / self._demfile).exists():
                regrid = True
            elif self.uxgrid is not None:
                with xr.open_dataset(self.output_dir / self._demfile) as ds:
                    if (
                        "face_elevation" not in ds
                        or "node_elevation" not in ds
                        or len(ds.face_elevation) != self.local.n_face
                        or len(ds.node_elevation) != self.local.n_node
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
        self.local.save(filename=self._demfile)
        self._dem_data = None  # Clear temporary DEM data storage
        return

    @parameter
    def pds_file_resolution(self) -> int | None:
        """The resolution of the DEM file to download from the PDS in meters per pixel."""
        return self._pds_file_resolution

    @pds_file_resolution.setter
    def pds_file_resolution(self, value: int | None):
        if value is None:
            # Valid resolution in pixels per degree in lat or lon
            # Choose one of the valid resolutions that gets you approximately 1000 pixels across the diameter of the local region
            local_size_deg = (2 * self.local_radius) / self.target.radius * (180.0 / np.pi)

            approx_pix_size = local_size_deg / 1000.0
            diffs = [abs(approx_pix_size - (360.0 / res)) for res in self._valid_pds_file_resolutions]
            value = self._valid_pds_file_resolutions[np.argmin(diffs)]
            print(f"Chosen PDS file resolution: {value} pix/deg")
        elif value not in self._valid_pds_file_resolutions:
            raise ValueError(f"Invalid pds_file_resolution {value}. Valid values are {self._valid_pds_file_resolutions}.")
        self._pds_file_resolution = value
        return
