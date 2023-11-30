import numpy as np
import xarray as xr
from mpas_tools.mesh.creation.build_mesh import build_spherical_mesh
import os
from pathlib import Path
from .general_utils import float_like
from typing import Tuple, List
from numpy.typing import NDArray

def generate(target_radius: float_like, 
             cell_size: float_like, 
             mesh_file: str | Path | os.PathLike = "surface_mesh.nc", 
             mesh_temp_dir:  str | Path | os.PathLike = ".mesh") -> xr.Dataset:
    """
    Generate a tessellated mesh of a sphere using the jigsaw-based mesh builder in MPAS-tools.

    This function generates temporary files in the `mesh_temp_dir` directory and saves the final mesh to `mesh_file`

    Returns
    -------
    xarray Dataset
       An MPAS-style dataset with the mesh data and variables 
    """
   
    cellWidth, lon, lat = make_uniform_cell_size(cell_size)
    orig_dir = os.getcwd()
    os.chdir(mesh_temp_dir)
    build_spherical_mesh(cellWidth, lon, lat, out_filename=str(mesh_file), earth_radius=target_radius, plot_cellWidth=False)
    os.chdir(orig_dir)

    mesh = xr.open_dataset(mesh_file)

    return mesh


def make_uniform_cell_size(cell_size: float_like) -> Tuple[NDArray,NDArray,NDArray]:
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
    constantCellWidth = cell_size * 1e-3 # build_spherical_mesh assumes units of km, so must be converted

    nlat = int(180/dlat) + 1
    nlon = int(360/dlon) + 1

    lat = np.linspace(-90., 90., nlat)
    lon = np.linspace(-180., 180., nlon)

    cellWidth = constantCellWidth * np.ones((lat.size, lon.size))
    return cellWidth, lon, lat


def set_elevation(nCells: int | None = None,
                  new_elev: NDArray[np.float64] | List[float_like] | None = None) -> xr.DataArray:
    """
    Set elevation data for the target's surface mesh.

    Parameters
    ----------
    nCells: int, optional
        The number of cells in the mesh. Must be passed if new_elev is not
    new_elev : array_like, optional
        New elevation data to be set. If None, the elevation is set to zero. Must be passed if nCells is not
    """        
    
    if new_elev is None:
        if nCells is None:
            raise ValueError("Must pass value for either nCells or new_elev.")
        new_elev = np.zeros(nCells,dtype=np.float64)
    else:
        if nCells is not None:
            raise ValueError("Cannot pass both nCells and new_elev")
            
    dem = xr.DataArray(
            data=new_elev,
            dims=["nCells"],
            attrs={"long_name":"elevation of cells"}
            )
    
    return dem


def get_cell_distance(mesh: xr.Dataset, 
                      location: Tuple[np.float64, np.float64],
                      target_radius: np.float64) -> xr.DataArray:
    """
    Computes the distances between cell centers and a given location.

    Parameters
    ----------
    location : Tuple[np.float64, np.float64]
        Tuple containing the longitude and latitude of the location in radians.

    Returns
    -------
    xarray.DataArray
        DataArray of distances for each cell in meters.
    """
    return calculate_haversine_distance(location[0],location[1],mesh.lonCell,mesh.latCell,target_radius)


def calculate_haversine_distance(lon1: float_like, 
                lat1: float_like, 
                lon2: float_like, 
                lat2: float_like,
                radius: float_like) -> np.float64:
    """
    Calculate the great circle distance between two points on a sphere.

    Parameters
    ----------
    lon1 : float_like
        Longitude of the first point in radians.
    lat1 : float_like
        Latitude of the first point in radians.
    lon2 : float_like
        Longitude of the second point in radians.
    lat2 : float_like
        Latitude of the second point in radians.
    radius : float_like
        Radius of the sphere in meters.

    Returns
    -------
    np.float64
        Great circle distance between the two points in meters.
    """
    # Calculate differences in coordinates
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    # Haversine formula
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return radius * c


def calculate_initial_bearing(lon1: float_like, 
                            lat1: float_like, 
                            lon2: float_like, 
                            lat2: float_like) -> np.float64:
    """
    Calculate the initial bearing from one point to another on the surface of a sphere.

    Parameters
    ----------
    lon1 : float_like
        Longitude of the first point in radians.
    lat1 : float_like
        Latitude of the first point in radians.
    lon2 : float_like
        Longitude of the second point in radians.
    lat2 : float_like
        Latitude of the second point in radians.

    Returns
    -------
    np.float64
        Initial bearing from the first point to the second point in radians.
    """
    # Calculate differences in coordinates
    dlon = np.mod(lon2 - lon1 + np.pi, 2 * np.pi) - np.pi

    # Haversine formula calculations
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    initial_bearing = np.arctan2(x, y)

    # Normalize bearing to 0 to 2*pi
    initial_bearing = (initial_bearing + 2 * np.pi) % (2 * np.pi)

    return initial_bearing


def get_cell_initial_bearing(mesh: xr.Dataset, 
                             location: Tuple[np.float64, np.float64]) -> xr.DataArray:
    """
    Computes the initial bearing between cell centers and a given location.

    Parameters
    ----------
    location : Tuple[np.float64, np.float64]
        Tuple containing the longitude and latitude of the location in radians.

    Returns
    -------
    xarray.DataArray
        DataArray of initial bearings for each cell in radians.
    """
    return calculate_initial_bearing(location[0], location[1], mesh.lonCell, mesh.latCell)


def get_average_surface(mesh: xr.Dataset, 
                        dem: xr.DataArray,
                        location: Tuple[np.float64, np.float64], 
                        radius: np.float64) -> Tuple[np.float64, np.float64]:
    """
    Calculate the orientation and radius of the hemispherical cap.

    Parameters
    ----------
    location : Tuple[float, float]
        Tuple containing the longitude and latitude of the reference location in radians.
    radius : float
        The reference radius of to compute the average over in meters.

    Returns
    -------
    cap_center_vector : ndarray
        The vector pointing to the center of the cap from the sphere's center.
    cap_radius : float
        The radius of the cap.
    """

    # Find cells within the crater radius
    cells_within_radius = mesh['crater_distance'] <= radius

    bearings = mesh['crater_bearing'].where(cells_within_radius, drop=True)
    distances = mesh['crater_distance'].where(cells_within_radius, drop=True)

    # Convert bearings to vector components
    # Bearing is angle from north, positive clockwise, but we need standard mathematical angle, positive counter-clockwise
    angles = np.deg2rad(90) - bearings  # Convert bearing to angle in radians
    x_components = np.cos(angles) * distances
    y_components = np.sin(angles) * distances

    # Calculate the weighted average vector components
    # Weight by the area of each cell to give more importance to larger cells
    cell_areas = mesh['areaCell'].where(cells_within_radius, drop=True)
    weighted_x = (x_components * cell_areas).sum() / cell_areas.sum()
    weighted_y = (y_components * cell_areas).sum() / cell_areas.sum()

    # Calculate the weighted mean elevation to get the z-component
    elevation_values = dem.where(cells_within_radius, drop=True)
    weighted_z = (elevation_values * cell_areas).sum() / cell_areas.sum()

    # Combine components to form the cap center vector
    cap_center_vector = np.array([weighted_x.item(), weighted_y.item(), weighted_z.item()])

    # The radius of the cap is the length of the cap center vector
    cap_radius = np.linalg.norm(cap_center_vector)

    return cap_center_vector, cap_radius

    
