use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray1,PyReadonlyArray2,PyArray1};

/// Computes the Haversine distance between a single point and an array of points on a sphere given their longitude and latitude in radians.
///
/// # Arguments
/// * `lon1`, `lat1` - Coordinates of the first point in radians.
/// * `lon2`, `lat2` - Array of coordinates of the second point in radians.
/// * `radius` - Radius of the sphere in meters.
///
/// # Returns
/// Distance in meters between the pairs of points along the surface of the sphere.
#[pyfunction]
pub fn compute_distances<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
    radius: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let lon2_v = lon2.as_array();
    let lat2_v = lat2.as_array();
    let result =  cratermaker_core::surface::compute_distances(
            lon1,
            lat1,
            lon2_v,
            lat2_v,
            radius
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes the initial bearing (forward azimuth) from a fixed point to each of a set of destination points.
///
/// The bearing is calculated on a spherical surface using great-circle paths and returned in radians,
/// normalized to the range [0, 2π).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `lon1` - Longitude of the reference point, in radians.
/// * `lat1` - Latitude of the reference point, in radians.
/// * `lon2` - Longitudes of destination points, in radians.
/// * `lat2` - Latitudes of destination points, in radians.
///
/// # Returns
///
/// * A NumPy array of initial bearing angles (radians), one for each (lon2, lat2) pair.
#[pyfunction]
pub fn compute_bearings<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let lon2_v = lon2.as_array();
    let lat2_v = lat2.as_array();
    let result =  cratermaker_core::surface::compute_bearings(
            lon1,
            lat1,
            lon2_v,
            lat2_v,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Applies one explicit diffusion update step over a surface mesh with variable diffusivity.
///
/// This function computes the change in elevation for each face on the mesh using a
/// face-centered finite-volume formulation of the operator:
///     ∂h/∂t = ∇ · (κ ∇h)
/// where κ varies per face. For each face, the flux with its neighbors is computed
/// using the expression (κ_f + κ_n)/2 * (h_n - h_f), summing over all neighbors n.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_kappa` - Topographic diffusivity (1D array of length n_face)
/// * `face_elevation` - Elevation value of the faces (1D array of length n_face)
/// * `face_area` - Area of the faces (1D array of length n_face)
/// * `edge_face_connectivity` - Indices of the faces saddle each edge. (2D array of shape n_edge x 2).
/// * `edge_face_distance` - Distances between the centers of the faces that saddle each edge in meters (1D array of length n_edge). 
///
/// # Returns
///
/// A NumPy array of shape (n_face,) giving the elevation change per face of the local mesh for this update step.
#[pyfunction]
pub fn apply_diffusion<'py>(
    py: Python<'py>,
    face_kappa: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_area: PyReadonlyArray1<'py, f64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_kappa_v = face_kappa.as_array();
    let face_elevation_v = face_elevation.as_array();
    let face_area_v = face_area.as_array();
    let edge_face_connectivity_v = edge_face_connectivity.as_array();
    let edge_face_distance_v = edge_face_distance.as_array();
    let edge_length_v = edge_length.as_array();
    let result =  cratermaker_core::surface::apply_diffusion(
            face_kappa_v,
            face_elevation_v,
            face_area_v,
            edge_face_connectivity_v,
            edge_face_distance_v,
            edge_length_v,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes the gradient vector in a radial direction defined by the face bearing at a face using the Green-Gauss method.
///
///
/// This function is designed to be parallel and returns a NumPy array of slopes
/// corresponding to the provided face indices.
///
/// # Arguments
/// * `py` - Python interpreter token.
/// * `variable` - The variable to compute the gradient for at each face (1D array).
/// * `face_lon` - Longitude at each face in radians (1D array).
/// * `face_lat` - Latitude at each face in radians (1D array).
/// * `face_bearing` - Bearing at each face in radians (1D array).
/// * `face_area` - Area of each face (1D array).
/// * `edge_face_connectivity` - Indices of the faces (global) that saddle each edge. (2D array of shape n_edge x 2).
/// * 'face_edge_connectivity` - Indices of the edges that surround each face (2D array of shape n_face x n_max_edges).
/// * `edge_face_distance` - Distances between the centers of the faces that saddle each edge in meters (1D array of length n_edge). 
///
/// # Returns
/// An arrays of radial gradient values same length as `face_indices`.
#[pyfunction]
pub fn compute_radial_gradient<'py>(
    py: Python<'py>,
    variable: PyReadonlyArray1<'py, f64>,
    face_lon: PyReadonlyArray1<'py, f64>,
    face_lat: PyReadonlyArray1<'py, f64>,
    face_bearing: PyReadonlyArray1<'py, f64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let variable_v = variable.as_array();
    let face_bearing_v = face_bearing.as_array();
    let face_lon_v = face_lon.as_array();
    let face_lat_v = face_lat.as_array();
    let edge_face_connectivity_v = edge_face_connectivity.as_array();
    let face_edge_connectivity_v = face_edge_connectivity.as_array();
    let edge_face_distance_v = edge_face_distance.as_array();
    let edge_length_v = edge_length.as_array();
    let result =  cratermaker_core::surface::compute_radial_gradient(
            variable_v,
            face_lon_v,
            face_lat_v,
            face_bearing_v,
            face_edge_connectivity_v,
            edge_face_connectivity_v,
            edge_face_distance_v,
            edge_length_v,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes the square root of the maximum squared slope at each face in a surface mesh.
///
/// For each face in the given `face_indices` subset, this function calculates the steepest
/// slope using pairs of neighboring faces, where slope is defined as the change in elevation
/// divided by the great-circle (haversine) distance. The result is the maximum root-sum-square
/// slope magnitude from adjacent neighbor pairs.
///
/// This function is designed to be parallel and returns a NumPy array of slopes
/// corresponding to the provided face indices.
///
/// # Arguments
/// * `py` - Python interpreter token.
/// * `face_elevation` - Elevation at each face (1D array).
/// * `face_area` - Area of each face (1D array).
/// * `edge_face_connectivity` - Indices of the faces (global) that saddle each edge. (2D array of shape n_edge x 2).
/// * 'face_edge_connectivity` - Indices of the edges that surround each face (2D array of shape n_face x n_max_edges).
/// * `edge_face_distance` - Distances between the centers of the faces that saddle each edge in meters (1D array of length n_edge). 
///
/// # Returns
/// A NumPy array of slope values (1D array), same length as `face_indices`.
#[pyfunction]
pub fn compute_slope<'py>(
    py: Python<'py>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_lon: PyReadonlyArray1<'py, f64>,
    face_lat: PyReadonlyArray1<'py, f64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_elevation = face_elevation.as_array();
    let face_lon = face_lon.as_array();
    let face_lat = face_lat.as_array();
    let edge_face_connectivity = edge_face_connectivity.as_array();
    let face_edge_connectivity = face_edge_connectivity.as_array();
    let edge_face_distance = edge_face_distance.as_array();
    let edge_length = edge_length.as_array();
    let result = cratermaker_core::surface::compute_slope(
            face_elevation,
            face_lon,
            face_lat,
            edge_face_connectivity,
            face_edge_connectivity,
            edge_face_distance,
            edge_length
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes the spatially varying diffusivity (`face_kappa`) for a slope collapse step.
///
/// Sets `kappa = diffmax` if any neighbor violates the critical slope threshold,
/// otherwise sets `kappa = 0.0`. `diffmax` is computed assuming a stable timestep of 1.0.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `critical_slope` - Maximum allowable slope (e.g., 0.7 for ~35 degrees).
/// * `face_elevation` - Elevation at each face (1D array).
/// * `face_area` - Area of each face (1D array).
/// * `edge_face_connectivity` - Indices of the faces that saddle each edge. (2D array of shape n_edge x 2).
/// * 'face_edge_connectivity` - Indices of the edges that surround each face (2D array of shape n_face x n_max_edges).
/// * `edge_face_distance` - Distances between the centers of the faces that saddle each edge in meters (1D array of length n_edge). 
///
/// # Returns
///
/// A NumPy array of `face_kappa` values.
#[pyfunction]
pub fn slope_collapse<'py>(
    py: Python<'py>,
    critical_slope: f64,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_lon: PyReadonlyArray1<'py, f64>,
    face_lat: PyReadonlyArray1<'py, f64>,
    face_area: PyReadonlyArray1<'py, f64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_elevation_v = face_elevation.as_array();
    let face_lon_v = face_lon.as_array();
    let face_lat_v = face_lat.as_array();
    let face_area_v = face_area.as_array();
    let edge_face_connectivity_v = edge_face_connectivity.as_array();
    let face_edge_connectivity_v = face_edge_connectivity.as_array();
    let edge_face_distance_v = edge_face_distance.as_array();
    let edge_length_v = edge_length.as_array();
    let result = cratermaker_core::surface::slope_collapse(
            critical_slope,
            face_elevation_v,
            face_lon_v,
            face_lat_v,
            face_area_v,
            edge_face_connectivity_v,
            face_edge_connectivity_v,
            edge_face_distance_v,
            edge_length_v
            )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes node elevations as area-weighted averages of adjacent face elevations.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_area` - Area of each face (1D array).
/// * `face_elevation` - Elevation at each face (1D array).
/// * `node_face_connectivity` - For each node, indices of connected faces (2D array).
///
/// # Returns
///
/// A NumPy array of node elevations (1D array).
#[pyfunction]
pub fn interpolate_node_elevation_from_faces<'py>(
    py: Python<'py>,
    face_area: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    node_face_connectivity: PyReadonlyArray2<'py, i64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_area_v = face_area.as_array();
    let face_elevation_v = face_elevation.as_array();
    let node_face_connectivity_v = node_face_connectivity.as_array();
    let result = cratermaker_core::surface::interpolate_node_elevation_from_faces(
            face_elevation_v, 
            face_area_v, 
            node_face_connectivity_v
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes 3D turbulence noise using multi-octave rotated simplex noise.
///
/// This function generates noise values over 3D positions defined by the arrays `x`, `y`, and `z`,
/// using a fractal sum of multiple noise octaves. Each octave is scaled by `freq^i` and `pers^i`,
/// and spatially rotated using the axis-angle values provided in `anchor`.
///
/// The result is normalized and scaled by `noise_height`.
///
/// # Arguments
/// * `py` - Python interpreter token.
/// * `x`, `y`, `z` - 1D arrays of the same length giving 3D positions.
/// * `noise_height` - Final amplitude scaling factor.
/// * `noise_width` - Base spatial scale (inverse frequency) of noise.
/// * `freq` - Frequency multiplier per octave (e.g., 2.0).
/// * `pers` - Amplitude multiplier per octave (e.g., 0.5).
/// * `anchor` - 2D array (num_octaves x 3) with rotation axis (x, y, z) per octave.
/// * `seed` - Seed for the noise generator.
///
/// # Returns
/// A 1D NumPy array of noise values, one per input coordinate triplet.
#[pyfunction]
pub fn turbulence_noise<'py>(
    py: Python<'py>,
    x: PyReadonlyArray1<'py, f64>,
    y: PyReadonlyArray1<'py, f64>,
    z: PyReadonlyArray1<'py, f64>,
    noise_height: f64,
    noise_width: f64,
    freq: f64,
    pers: f64,
    anchor: PyReadonlyArray2<'py, f64>,
    seed: u32,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let z_v = z.as_array();
    let anchor_v = anchor.as_array();
    let result = cratermaker_core::surface::turbulence_noise(
            x_v,
            y_v,
            z_v,
            noise_height,
            noise_width,
            freq,
            pers,
            anchor_v,
            seed
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Constructs the edge-other distances for a surface mesh
///
/// This will compute the haversine distance between the coordinates of the two components (faces, nodes, etc) that are associated with each edge.
///
/// # Arguments
/// * `py` - Python interpreter token.
/// * `edge_connectivity` - Indices of the elements associated each edge. (2D array of shape n_edge x 2).
/// * `lon1`, `lat1` - Coordinates of the reference point in radians.
///
/// # Returns
/// A 1D NumPy array of edge-other distances in meters, one for each edge.
#[pyfunction]
pub fn compute_edge_distances<'py>(
    py: Python<'py>,
    edge_connectivity: PyReadonlyArray2<'py, i64>,
    lon: PyReadonlyArray1<'py, f64>,
    lat: PyReadonlyArray1<'py, f64>,
    radius: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let edge_connectivity_v = edge_connectivity.as_array();
    let lon_v = lon.as_array();
    let lat_v = lat.as_array();
    let result = cratermaker_core::surface::compute_edge_distances(
            edge_connectivity_v,
            lon_v,
            lat_v,
            radius,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


