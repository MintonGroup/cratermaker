import numpy as np
import os
import uxarray as uxr
import xarray as xr
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh

def make_target_mesh(mesh_file,data_file,target, pix):
    """
    Generate a tessellated mesh of a sphere using the jigsawpy library and convert it to GLB format.

    This function sets up Jigsaw mesh files and a mesh body, defines a basic sphere using an ellipsoid mesh model,
    sets mesh options, generates a tessellated mesh, and converts the mesh to a uxarray object. The generated mesh 
    is then saved in GLB format.

    Returns
    -------
    None
        The function does not return a value but updates the `self.mesh` attribute with the generated uxarray object.
    """
    def cellWidthVsLatLon(pix):
        """
        Create cell width array for this mesh on a regular latitude-longitude grid.
        Returns
        -------
        cellWidth : ndarray
            m x n array of cell width in km
        lon : ndarray
            longitude in degrees (length n and between -180 and 180)
        lat : ndarray
            longitude in degrees (length m and between -90 and 90)
        """
        dlat = 10
        dlon = 10
        constantCellWidth = np.sqrt(2.0) * pix * 1e-3

        nlat = int(180/dlat) + 1
        nlon = int(360/dlon) + 1

        lat = np.linspace(-90., 90., nlat)
        lon = np.linspace(-180., 180., nlon)

        cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
        return cellWidth, lon, lat
    
    
    cellWidth, lon, lat = cellWidthVsLatLon(pix)
    build_spherical_mesh(cellWidth, lon, lat, out_filename=mesh_file, earth_radius=target.radius,dir=".cache")
 
    primal_mesh = uxr.open_grid(mesh_file, use_dual=False) 
    elevation = uxr.UxDataArray(
        name="elevation",
        data=np.zeros(primal_mesh.n_node), 
        dims=['n_node'],
        uxgrid=primal_mesh
    )
    uxds = elevation.to_dataset()
    uxds.to_netcdf(data_file)
    
    return uxds

def load_target_mesh(mesh_file, data_file):
    """
    Load a target mesh from a file into the `self.mesh` attribute.

    This function uses the trimesh library to read a mesh file, which is expected to be in GLB format. The function
    handles the peculiarity of trimesh reading the mesh as a Scene instead of a Trimesh object and extracts the 
    actual mesh from the scene.

    Returns
    -------
    None
        The function does not return a value but updates the `self.mesh` attribute with the loaded trimesh object.
    """

    mesh = uxr.open_dataset(mesh_file, data_file)
    
    return mesh