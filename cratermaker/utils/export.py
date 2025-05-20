from pathlib import Path

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
from tqdm import tqdm

from cratermaker.components.surface import Surface
from cratermaker.constants import (
    _CIRCLE_FILE_NAME,
    _COMBINED_DATA_FILE_NAME,
    _EXPORT_DIR,
    _SURFACE_DIR,
    _VTK_FILE_EXTENSION,
    FloatLike,
)


def to_vtk(surface: Surface, *args, **kwargs) -> None:
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

    if not isinstance(surface, Surface):
        raise TypeError(
            "The surface argument must be an instance of the Surface class."
        )
    # Create the output directory if it doesn't exist
    out_dir = surface.simdir / _EXPORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    data_dir = surface.simdir / _SURFACE_DIR
    data_file_list = list(data_dir.glob("*.nc"))
    if surface.grid_file in data_file_list:
        data_file_list.remove(surface.grid_file)

    # Delete old export files if they exist
    old_export_files = list(out_dir.glob(f"*.{_VTK_FILE_EXTENSION}"))
    if len(old_export_files) > 0:
        for old_file in old_export_files:
            old_file.unlink()

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
    warp.SetInputArrayToProcess(
        0, 0, 0, vtkUnstructuredGrid.FIELD_ASSOCIATION_POINTS, "node_elevation"
    )

    writer = vtkXMLPolyDataWriter()
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()

    with xr.open_mfdataset(data_file_list) as ds:
        # Warp the surface based on node_elevation data so that the exported mesh is the true shape rather than a sphere
        for i in tqdm(
            range(len(ds.time)),
            desc="Exporting VTK files",
            unit="file",
            position=0,
            leave=True,
        ):
            ids = ds.isel(time=i).load()
            current_grid = vtkUnstructuredGrid()
            current_grid.DeepCopy(vtk_data)

            for v in ds.variables:
                array = numpy_to_vtk(ids[v].values, deep=True)
                array.SetName(v)
                n = ids[v].size
                if "n_face" in ids[v].dims:
                    current_grid.GetCellData().AddArray(array)
                elif "n_node" in ids[v].dims:
                    current_grid.GetPointData().AddArray(array)
                    if v == "node_elevation":
                        current_grid.GetPointData().SetActiveScalars(v)
                elif n == 1:
                    current_grid.GetFieldData().AddArray(array)

            geomFilter = vtkGeometryFilter()
            geomFilter.SetInputData(current_grid)
            geomFilter.Update()
            polyData = geomFilter.GetOutput()

            normalsFilter = vtkPolyDataNormals()
            normalsFilter.SetInputData(polyData)
            normalsFilter.ComputeCellNormalsOn()
            normalsFilter.ConsistencyOn()  # Tries to make normals consistent across shared edges
            normalsFilter.AutoOrientNormalsOn()  # Attempt to orient normals consistently outward/inward
            normalsFilter.SplittingOff()
            normalsFilter.Update()
            polyDataWithNormals = normalsFilter.GetOutput()

            warp.SetInputData(polyDataWithNormals)
            warp.Update()
            warped_output = warp.GetOutput()
            output_filename = out_dir / _COMBINED_DATA_FILE_NAME.replace(
                ".nc", f"{i:06d}.{_VTK_FILE_EXTENSION}"
            )
            writer.SetFileName(output_filename)
            writer.SetInputData(warped_output)
            writer.Write()

    return


def make_circle_file(
    surface,
    diameters: FloatLike | ArrayLike,
    longitudes: FloatLike | ArrayLike,
    latitudes: FloatLike | ArrayLike,
    output_filename: str | Path | None = None,
    *args,
    **kwargs,
) -> None:
    """
    Plot circles of diameter D centered at the given location.

    Parameters
    ----------
    diameters : FloatLike or ArrayLike
        Diameters of the circles in m.
    longitudes : FloatLike or ArrayLike of Floats
        Longitudes of the circle centers in degrees.
    latitudes : FloatLike or ArrayLike of Floats
        Latitudes of the circle centers in degrees.
    out_filename : str or Path, optional
        Name of the output file. If not provided, the default is "circle.vtp"
    """
    from vtk import (
        vtkCellArray,
        vtkPoints,
        vtkPolyData,
        vtkPolyLine,
        vtkXMLPolyDataWriter,
    )

    if output_filename is None:
        output_filename = surface.simdir / _EXPORT_DIR / _CIRCLE_FILE_NAME
    else:
        output_filename = surface.simdir / _EXPORT_DIR / output_filename

    diameters = np.atleast_1d(diameters)
    longitudes = np.atleast_1d(longitudes)
    latitudes = np.atleast_1d(latitudes)
    # Check for length consistency
    if len(diameters) != len(longitudes) or len(diameters) != len(latitudes):
        raise ValueError(
            "The diameters, latitudes, and longitudes, arguments must have the same length"
        )

    # Validate non-negative values
    if np.any(diameters < 0):
        raise ValueError("All values in 'diameters' must be non-negative")

    sphere_radius = surface.target.radius

    def create_circle(lon, lat, circle_radius, num_points=360):
        """
        Create a circle on the sphere's surface with a given radius and center.

        Parameters
        ----------
        lon : float
            Longitude of the circle's center in degrees.
        lat : float
            Latitude of the circle's center in degrees.
        circle_radius : float
            Radius of the circle in meters.
        num_points : int, optional
            Number of points to use to approximate the circle. The default is 360.
        """
        # Create an array of angle steps for the circle
        radians = np.linspace(0, 2 * np.pi, num_points)
        # Convert latitude and longitude to radians
        lat_rad = np.deg2rad(lat)
        lon_rad = np.deg2rad(lon)

        # Calculate the Cartesian coordinates for the circle's center
        center_x = sphere_radius * np.cos(lat_rad) * np.cos(lon_rad)
        center_y = sphere_radius * np.cos(lat_rad) * np.sin(lon_rad)
        center_z = sphere_radius * np.sin(lat_rad)

        # Calculate the vectors for the local east and north directions on the sphere's surface
        east = np.array([-np.sin(lon_rad), np.cos(lon_rad), 0])
        north = np.array(
            [
                -np.cos(lon_rad) * np.sin(lat_rad),
                -np.sin(lon_rad) * np.sin(lat_rad),
                np.cos(lat_rad),
            ]
        )

        # Initialize arrays to hold the circle points
        x = np.zeros_like(radians)
        y = np.zeros_like(radians)
        z = np.zeros_like(radians)

        # Calculate the points around the circle
        for i in range(num_points):
            x[i] = (
                center_x
                + circle_radius * np.cos(radians[i]) * east[0]
                + circle_radius * np.sin(radians[i]) * north[0]
            )
            y[i] = (
                center_y
                + circle_radius * np.cos(radians[i]) * east[1]
                + circle_radius * np.sin(radians[i]) * north[1]
            )
            z[i] = (
                center_z
                + circle_radius * np.cos(radians[i]) * east[2]
                + circle_radius * np.sin(radians[i]) * north[2]
            )

        return x, y, z

    points = vtkPoints()
    lines = vtkCellArray()
    point_id = 0  # Keep track of the point ID across all circles

    for lon, lat, diameter in zip(longitudes, latitudes, diameters):
        circle_radius = diameter / 2
        x, y, z = create_circle(lon, lat, circle_radius)

        for i in range(len(x)):
            points.InsertNextPoint(x[i], y[i], z[i])

        polyline = vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(x))
        for i in range(len(x)):
            polyline.GetPointIds().SetId(i, point_id)
            point_id += 1

        lines.InsertNextCell(polyline)

    # Create a polydata object and add points and lines to it
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)

    # Write the polydata to a VTK file
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(output_filename)
    writer.SetInputData(polydata)

    # Optional: set the data mode to binary to save disk space
    writer.SetDataModeToBinary()
    writer.Write()

    return
