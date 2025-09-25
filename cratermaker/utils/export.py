from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from numpy.typing import ArrayLike
from pyproj import CRS, Geod, Transformer
from tqdm import tqdm

from cratermaker.constants import (
    _CIRCLE_FILE_NAME,
    _COMBINED_DATA_FILE_NAME,
    _EXPORT_DIR,
    _GRID_FILE_NAME,
    _SURFACE_DIR,
    _VTK_FILE_EXTENSION,
    FloatLike,
)

if TYPE_CHECKING:
    from cratermaker.components.surface import Surface


def to_vtk(
    surface: Surface,
    interval_number: int = 0,
    time_variables: dict | None = None,
    save_geometry=True,
    **kwargs,
) -> None:
    """
    Export the surface mesh to a VTK file and stores it in the default export directory.
    """
    from vtk import (
        VTK_POLYGON,
        vtkPoints,
        vtkUnstructuredGrid,
        vtkWarpScalar,
        vtkXMLPolyDataWriter,
    )
    from vtkmodules.util.numpy_support import numpy_to_vtk
    from vtkmodules.vtkFiltersCore import vtkPolyDataNormals
    from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter

    from cratermaker import Surface

    if not isinstance(surface, Surface):
        raise TypeError("The surface argument must be an instance of the Surface class.")
    # Create the output directory if it doesn't exist
    out_dir = surface.simdir / _EXPORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    data_dir = surface.simdir / _SURFACE_DIR
    data_file_list = list(data_dir.glob("*.nc"))
    if surface.grid_file in data_file_list:
        data_file_list.remove(surface.grid_file)

    # Convert uxarray grid arrays to regular numpy arrays for vtk processing
    n_node = surface.n_node
    n_face = surface.n_face
    node_x = surface.node_x
    node_y = surface.node_y
    node_z = surface.node_z
    n_nodes_per_face = surface.n_nodes_per_face
    face_node_connectivity = surface.face_node_connectivity

    vtk_data = vtkUnstructuredGrid()
    nodes = vtkPoints()
    for i in range(n_node):
        nodes.InsertNextPoint(node_x[i], node_y[i], node_z[i])
    vtk_data.SetPoints(nodes)
    vtk_data.Allocate(n_face)
    for i, n in enumerate(n_nodes_per_face):
        point_ids = face_node_connectivity[i][0:n]
        vtk_data.InsertNextCell(VTK_POLYGON, n, point_ids)

    warp = vtkWarpScalar()
    warp.SetInputArrayToProcess(0, 0, 0, vtkUnstructuredGrid.FIELD_ASSOCIATION_POINTS, "node_elevation")

    writer = vtkXMLPolyDataWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()

    if save_geometry:
        # Saves the surface mesh and its geometry as a separate file
        geometry_variables = [
            "node_x",
            "node_y",
            "node_z",
            "node_lon",
            "node_lat",
            "face_x",
            "face_y",
            "face_z",
            "face_lon",
            "face_lat",
            "face_area",
            "face_size",
        ]
        current_grid = vtkUnstructuredGrid()
        current_grid.DeepCopy(vtk_data)

        for v in geometry_variables:
            # extract the attribute v from the surface object
            array = numpy_to_vtk(getattr(surface, v), deep=True)
            array.SetName(v)
            n = getattr(surface, v).size
            if n == surface.n_face:
                current_grid.GetCellData().AddArray(array)
            elif n == surface.n_node:
                current_grid.GetPointData().AddArray(array)

        geom_filter = vtkGeometryFilter()
        geom_filter.SetInputData(current_grid)
        geom_filter.Update()
        poly_data = geom_filter.GetOutput()

        normals_filter = vtkPolyDataNormals()
        normals_filter.SetInputData(poly_data)
        normals_filter.ComputeCellNormalsOn()
        normals_filter.ConsistencyOn()  # Tries to make normals consistent across shared edges
        normals_filter.AutoOrientNormalsOn()  # Attempt to orient normals consistently outward/inward
        normals_filter.SplittingOff()
        normals_filter.Update()
        poly_data_with_normals = normals_filter.GetOutput()

        output_filename = out_dir / _GRID_FILE_NAME.replace(".nc", f".{_VTK_FILE_EXTENSION}")
        writer.SetFileName(output_filename)
        writer.SetInputData(poly_data_with_normals)
        writer.Write()

    ds = surface.uxds.load()
    current_grid = vtkUnstructuredGrid()
    current_grid.DeepCopy(vtk_data)

    for v in ds.variables:
        array = numpy_to_vtk(ds[v].values, deep=True)
        array.SetName(v)
        n = ds[v].size
        if "n_face" in ds[v].dims:
            current_grid.GetCellData().AddArray(array)
        elif "n_node" in ds[v].dims:
            current_grid.GetPointData().AddArray(array)
            if v == "node_elevation":
                current_grid.GetPointData().SetActiveScalars(v)
        elif n == 1:
            current_grid.GetFieldData().AddArray(array)

    if time_variables is None:
        time_variables = {"elapsed_time": float(interval_number)}
    else:
        if not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")

    for k, v in time_variables.items():
        array = numpy_to_vtk(np.array([v]), deep=True)
        array.SetName(k)
        current_grid.GetFieldData().AddArray(array)

    geom_filter = vtkGeometryFilter()
    geom_filter.SetInputData(current_grid)
    geom_filter.Update()
    poly_data = geom_filter.GetOutput()

    normals_filter = vtkPolyDataNormals()
    normals_filter.SetInputData(poly_data)
    normals_filter.ComputeCellNormalsOn()
    normals_filter.ConsistencyOn()  # Tries to make normals consistent across shared edges
    normals_filter.AutoOrientNormalsOn()  # Attempt to orient normals consistently outward/inward
    normals_filter.SplittingOff()
    normals_filter.Update()
    poly_data_with_normals = normals_filter.GetOutput()

    warp.SetInputData(poly_data_with_normals)
    warp.Update()
    warped_output = warp.GetOutput()
    output_filename = out_dir / _COMBINED_DATA_FILE_NAME.replace(".nc", f"{interval_number:06d}.{_VTK_FILE_EXTENSION}")
    writer.SetFileName(output_filename)
    writer.SetInputData(warped_output)
    writer.Write()

    return


def to_gpkg(
    surface: Surface,
    interval_number: int = 0,
    **kwargs,
) -> None:
    """
    Export the surface data to a GeoPackage file and stores it in the default export directory.

    Parameters
    ----------
    surface : Surface
        The surface object containing the data to export.
    interval_number : int, optional
        The interval number to save, by default 0.
    **kwargs : Any
        Additional keyword arguments (not used).
    """
    out_dir = surface.simdir / _EXPORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    gpkg_path = out_dir / f"surface_{interval_number:06d}.gpkg"

    # load data and select the face-based variables
    ds = surface.uxds.load()
    variables = [v for v in ds.data_vars if any(dim == "n_face" for dim in ds[v].dims)]
    if not variables:
        raise ValueError("No face-based variables found to export to GeoPackage.")

    # Export each face-associated variable as its own layer in the GeoPackage file
    for var in variables:
        gdf = ds[var].to_geodataframe(engine="geopandas")

        if getattr(gdf, "crs", None) is None:
            gdf = gdf.set_crs(surface.crs)

        if var not in gdf.columns:
            value_col = var if var in gdf.columns else ("value" if "value" in gdf.columns else None)
            if value_col is not None and value_col != var:
                gdf[var] = gdf[value_col]

        gdf.to_file(gpkg_path, layer=var, driver="GPKG")

    return


def crater_layer(
    crater_data: dict,
    surface: Surface,
    interval_number: int = 0,
    layer_name: str = "craters",
    **kwargs,
) -> None:
    """
    Export the crater data to a GeoPackage file and stores it in the default export directory.
    """
    from shapely.geometry import Point, Polygon
    from shapely.ops import transform

    def _geodesic_ellipse_polygon(
        lon: float,
        lat: float,
        a: float,
        b: float,
        az_deg: float = 0.0,
        n: int = 150,
        R_pole: FloatLike = 1.0,
        R_equator: FloatLike = 1.0,
    ) -> Polygon:
        """
        Geodesic ellipse on a sphere: for each bearing theta from the center, we shoot a geodesic with distance r(theta) = (a*b)/sqrt((b*ct)^2 + (a*st)^2), then rotate all bearings by az_deg.

        Parameters
        ----------
        lon : float
            Longitude of the ellipse center in degrees.
        lat : float
            Latitude of the ellipse center in degrees.
        a : float
            Semi-major axis in meters.
        b : float
            Semi-minor axis in meters.
        az_deg : float, optional
            Azimuth rotation of the ellipse in degrees clockwise from north, by default 0.0.
        n : int, optional
            Number of points to use for the polygon, by default 150.
        R_pole : FloatLike, optional
            Planetary polar radius in units of meters, by default 1.0.
        R_equator : FloatLike, optional
            Planetary equatorial radius in units of meters, by default 1.0.

        Returns
        -------
        Returns a Shapely Polygon in lon/lat degrees.
        """
        geod = Geod(a=R_pole, b=R_equator)
        theta = np.linspace(0.0, 360.0, num=n, endpoint=False)

        # Polar radius of an axis-aligned ellipse in a Euclidean tangent plane
        ct = np.cos(np.deg2rad(theta))
        st = np.sin(np.deg2rad(theta))
        r = (a * b) / np.sqrt((b * ct) ** 2 + (a * st) ** 2)

        # Bearings (from east, CCW) rotated by azimuth
        bearings = (theta + az_deg) % 360.0

        # Forward geodesic for each bearing/distance
        lon, lat, _ = geod.fwd(lon * np.ones_like(bearings), lat * np.ones_like(bearings), bearings, r)

        # Normalize longitudes to [-180, 180] to play nice with GIS
        lon = ((np.asarray(lon) + 180.0) % 360.0) - 180.0

        return Polygon(zip(lon, lat, strict=False))

    out_dir = surface.simdir / _EXPORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    output_filename = out_dir / f"surface_{interval_number:06d}.gpkg"

    if "final_diameter" not in crater_data or "longitude" not in crater_data or "latitude" not in crater_data:
        raise ValueError("crater_data must contain 'final_diameter', 'longitude', and 'latitude' keys.")

    diameter = np.atleast_1d(crater_data["final_diameter"])
    longitude = np.atleast_1d(crater_data["longitude"])
    latitude = np.atleast_1d(crater_data["latitude"])

    attrs_df = pd.DataFrame({k: np.asarray(v) for k, v in crater_data.items()})

    # Check for length consistency
    if len(diameter) != len(longitude) or len(diameter) != len(latitude):
        raise ValueError("The diameter, latitude, and longitude, arguments must have the same length")

    # Validate non-negative values
    if np.any(diameter < 0):
        raise ValueError("All values in 'diameter' must be non-negative")

    geoms = []

    for lon, lat, diam in zip(longitude, latitude, diameter, strict=False):
        lon = float(lon)
        lat = float(lat)
        radius = float(diam) / 2.0
        poly = _geodesic_ellipse_polygon(lon, lat, a=radius, b=radius, R_pole=surface.radius, R_equator=surface.radius)
        geoms.append(poly)

    gdf = gpd.GeoDataFrame(attrs_df, geometry=geoms, crs=surface.crs)

    # Write to GeoPackage (layer name 'craters')
    gdf.to_file(output_filename, layer=layer_name, driver="GPKG")

    return


def to_geotiff(
    surface: Surface,
    interval_number: int = 0,
    bounds: tuple[float, float, float, float] | None = None,
    dtype: str = "float32",
    nodata: float | None = 0.0,
) -> None:
    """
    Rasterize a face-based elevation variable into a GeoTIFF using rasterio.

    Parameters
    ----------
    surface : Surface
        Source surface with an unstructured mesh and face-based data in UxArray.
    interval_number : int, optional
        Interval number to save, by default 0.
    bounds : tuple[float, float, float, float] | None, optional
        (minx, miny, maxx, maxy) bounds of the output raster in the surface CRS; if None, use the full extent of the data, by default None.
    dtype : str, optional
        Data type for the output raster, by default "float32".
    nodata : float | None, optional
        NoData value for the output raster; if None, no NoData value is set, by default np.nan.
    """
    import rasterio as rio
    from rasterio.features import rasterize

    out_dir = surface.simdir / _EXPORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = surface.uxds.load()
    variables = [v for v in ds.data_vars if any(dim == "n_face" for dim in ds[v].dims)]
    if not variables:
        raise ValueError("No face-based variables found to export to GeoTiff.")

    for var in variables:
        gdf = ds[var].to_geodataframe(engine="geopandas")

        # Ensure CRS is set
        if getattr(gdf, "crs", None) is None:
            gdf = gdf.set_crs(surface.crs)

        # Drop empty geometries
        gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notnull()].copy()

        # Choose bounds
        if bounds is None:
            minx, miny, maxx, maxy = gdf.total_bounds
        else:
            minx, miny, maxx, maxy = bounds

        # degrees per pixel (same for lat/lon if you want square pixels)
        deg_per_pix = 360.0 * surface.pix / (2 * np.pi * surface.radius)

        width = int(np.ceil((maxx - minx) / deg_per_pix))
        height = int(np.ceil((maxy - miny) / deg_per_pix))

        transform = rio.transform.from_origin(-180.0, 90.0, deg_per_pix, deg_per_pix)
        geo_axes = plt.axes(projection=cartopy.crs.PlateCarree())

        print(f"Rasterizing variable '{var}' to GeoTIFF...")
        raster = ds[var].to_raster(ax=geo_axes)
        profile = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": dtype,
            "crs": gdf.crs,
            "transform": transform,
        }
        if nodata is not None:
            profile["nodata"] = nodata
        output_file = out_dir / f"{var}_{interval_number:06d}.tiff"
        with rio.open(output_file, "w", **profile) as dst:
            dst.write(raster, 1)

    return
