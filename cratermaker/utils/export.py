from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
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
