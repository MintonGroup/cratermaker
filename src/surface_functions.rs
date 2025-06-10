use std::f64::consts::PI;

use ndarray::Zip;
use noise::{NoiseFn, RotatePoint, ScalePoint, SuperSimplex};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Computes the positive modulus of `x` with respect to `m`.
/// Ensures the result is always in the range `[0, m)`.
#[inline]
fn positive_mod(x: f64, m: f64) -> f64 {
    ((x % m) + m) % m
}

/// Computes the Haversine distance between two points on a sphere given their longitude and latitude in radians.
///
/// # Arguments
/// * `lon1`, `lat1` - Coordinates of the first point in radians.
/// * `lon2`, `lat2` - Coordinates of the second point in radians.
/// * `radius` - Radius of the sphere in meters.
///
/// # Returns
/// Distance in meters between the two points along the surface of the sphere.
#[inline]
fn haversine_distance_scalar(lon1: f64, lat1: f64, lon2: f64, lat2: f64, radius: f64) -> f64 {
    let dlon = lon2 - lon1;
    let dlat = lat2 - lat1;
    let a = (dlat / 2.0).sin().powi(2)
        + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();
    radius * c
}


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
pub fn calculate_distance<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
    radius: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {

    let lon2 = lon2.as_array();
    let lat2 = lat2.as_array();
    let mut result = ndarray::Array1::<f64>::zeros(lon2.len());

    Zip::from(&mut result)
        .and(lon2)
        .and(lat2)
        .into_par_iter()
        .for_each(|(out, &lon2_i, &lat2_i)| {
            *out = haversine_distance_scalar(lon1, lat1, lon2_i, lat2_i, radius);
        });

    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes the initial bearing (forward azimuth) from point 1 to point 2 on a sphere.
#[inline]
fn compute_initial_bearing(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let dlon = positive_mod(lon2 - lon1 + PI, 2.0 * PI) - PI;
    let x = dlon.sin() * lat2.cos();
    let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
    x.atan2(y)
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
pub fn calculate_bearing<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let lon2 = lon2.as_array();
    let lat2 = lat2.as_array();
    let mut result = ndarray::Array1::<f64>::zeros(lon2.len());

    Zip::from(&mut result)
        .and(lon2)
        .and(lat2)
        .into_par_iter()
        .for_each(|(out, &lon2, &lat2)| {
            let initial_bearing = compute_initial_bearing(lon1, lat1, lon2, lat2);

            // Normalize bearing to 0 to 2*pi
            *out = (initial_bearing + 2.0 * PI) % (2.0 * PI);
        });

    Ok(PyArray1::from_owned_array(py, result))
}

/// Computes a conservative explicit timestep for stability using per-face aggregated flux contributions.
///
/// This method loops over all edges and accumulates inverse dt estimates for each face.
/// The final dt is the minimum over all faces of the reciprocal flux sum.
fn compute_dt_max(
    edge_face_distance: &ndarray::ArrayView1<'_, f64>,
    edge_face_connectivity: &ndarray::ArrayView2<'_, i64>,
    face_kappa: &ndarray::ArrayView1<'_, f64>,
    face_area: &ndarray::ArrayView1<'_, f64>,
    edge_length: &ndarray::ArrayView1<'_, f64>,
) -> f64 {
    let n_face = face_kappa.len();
    let mut inverse_dt_sum = ndarray::Array1::<f64>::zeros(n_face);

    for (e, faces) in edge_face_connectivity.outer_iter().enumerate() {
        let distance = edge_face_distance[e];
        if distance <= 0.0 {
            continue;
        }

        let [f1, f2] = match <[i64; 2]>::try_from(faces.as_slice().unwrap()) {
            Ok(pair) => pair,
            Err(_) => continue,
        };

        let f1 = f1 as usize;
        let f2 = f2 as usize;
        if f1 >= n_face || f2 >= n_face {
            continue;
        }

        let k1 = face_kappa[f1];
        let k2 = face_kappa[f2];
        let k_avg = 0.5 * (k1 + k2);
        let length = edge_length[e];

        if k_avg > 0.0 {
            let contrib_f1 = 2.0 * k_avg * length / (distance * face_area[f1]);
            let contrib_f2 = 2.0 * k_avg * length / (distance * face_area[f2]);
            inverse_dt_sum[f1] += contrib_f1;
            inverse_dt_sum[f2] += contrib_f2;
        }
    }

    inverse_dt_sum
        .iter()
        .filter(|&&x| x > 0.0)
        .map(|&x| 1.0 / x)
        .fold(f64::INFINITY, f64::min)
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
    let face_kappa = face_kappa.as_array();
    let face_elevation = face_elevation.as_array();
    let face_area = face_area.as_array();
    let edge_face_connectivity = edge_face_connectivity.as_array();
    let edge_face_distance = edge_face_distance.as_array();
    let edge_length = edge_length.as_array();
    let n_face = face_area.len();

    // Compute initial dt von neumann stability condition
    let dt_max = compute_dt_max(
        &edge_face_distance,
        &edge_face_connectivity,
        &face_kappa,
        &face_area,
        &edge_length,
    );

    if !dt_max.is_finite() || dt_max <= 0.0 {
        let result = ndarray::Array1::<f64>::zeros(n_face);
        return Ok(PyArray1::from_owned_array(py, result));
    }

    let nloops = (1.0 / dt_max).ceil() as usize;
    let dt = 1.0 / nloops as f64;
    let fac = dt / 2.0;

    let mut face_delta_elevation = ndarray::Array1::<f64>::zeros(n_face);

    for _ in 0..nloops {
        // Initialize dhdt as zeros with length n_face
        let mut dhdt = ndarray::Array1::<f64>::zeros(n_face);

        // Loop over edges, accumulate flux contributions to each face
        for (e, faces) in edge_face_connectivity.outer_iter().enumerate() {
            let [f1, f2] = <[i64; 2]>::try_from(faces.as_slice().unwrap()).unwrap();
            let f1 = f1 as usize;
            let f2 = f2 as usize;
            let distance = edge_face_distance[e];
            let length = edge_length[e];

            if distance <= 0.0 || f1 >= n_face || f2 >= n_face {
                continue;
            }
            if f1 == f2 {
                panic!("Edge {e} has identical face indices: {f1}");
            }
            if face_area[f1] == 0.0 || face_area[f2] == 0.0 {
                panic!("Zero area at f1={} or f2={} on edge {}", f1, f2, e);
            }

            let h1 = face_elevation[f1] + face_delta_elevation[f1];
            let h2 = face_elevation[f2] + face_delta_elevation[f2];
            let k1 = face_kappa[f1];
            let k2 = face_kappa[f2];

            let flux = (k1 + k2) * (h2 - h1) / distance * length;

            dhdt[f1] += flux / face_area[f1];
            dhdt[f2] -= flux / face_area[f2];
        }

        for f in 0..n_face {
            face_delta_elevation[f] += fac * dhdt[f];
        }
    }

    Ok(PyArray1::from_owned_array(py, face_delta_elevation))
}


/// Computes the slope squared at a face using the Green-Gauss method.
///
/// For each face, this function estimates the gradient vector (∂h/∂x, ∂h/∂y)
/// in the local tangent plane by summing contributions from all connected neighbors,
/// using the Green-Gauss theorem applied to the face's surrounding edges.
/// Returns the squared magnitude of this gradient as the slope squared.
///
/// For each neighbor, we take the vector (dx, dy) in the tangent plane and the elevation difference dh,
/// and estimate the gradient using the Green-Gauss formula.
/// Since we lack explicit (x, y) positions, we approximate dx, dy as unit vectors around the face.
/// This provides a reasonable local coordinate system for gradient estimation.
#[inline(always)]
fn compute_slope_squared(
    f: usize,
    face_elevation: &ndarray::Array1<f64>,
    connected_edges: &ndarray::ArrayView1<'_, i64>,
    edge_face_connectivity: &ndarray::ArrayView2<i64>,
    edge_face_distance: &ndarray::ArrayView1<f64>,
    edge_length: &ndarray::ArrayView1<f64>,
) -> f64 {
    let mut dh_dx = 0.0;
    let mut dh_dy = 0.0;
    let mut weight_sum = 0.0;

    let mut neighbors = Vec::new();
    for &edge_id in connected_edges.iter() {
        if edge_id < 0 {
            continue;
        }
        let e = edge_id as usize;
        let [f1, f2] = match <[i64; 2]>::try_from(edge_face_connectivity.row(e).as_slice().unwrap()) {
            Ok(pair) => pair,
            Err(_) => continue,
        };
        let f1 = f1 as usize;
        let f2 = f2 as usize;
        let other = if f == f1 { f2 } else if f == f2 { f1 } else { continue };
        let d = edge_face_distance[e];
        let l = edge_length[e];
        if d <= 0.0 || l <= 0.0 {
            continue;
        }
        neighbors.push((other, d, l));
    }

    let n = neighbors.len();
    if n < 2 {
        return 0.0;
    }

    for (i, &(other, distance, length)) in neighbors.iter().enumerate() {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        let nx = angle.cos();
        let ny = angle.sin();
        let dh = face_elevation[other] - face_elevation[f];
        let flux = dh * length / distance;

        dh_dx += flux * nx;
        dh_dy += flux * ny;
        weight_sum += length;
    }

    if weight_sum > 0.0 {
        dh_dx /= weight_sum;
        dh_dy /= weight_sum;
        return dh_dx * dh_dx + dh_dy * dh_dy;
    }

    0.0
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
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_elevation = face_elevation.as_array();
    let edge_face_connectivity = edge_face_connectivity.as_array();
    let face_edge_connectivity = face_edge_connectivity.as_array();
    let edge_face_distance = edge_face_distance.as_array();
    let edge_length = edge_length.as_array();
    let n_face = face_elevation.len();

    let face_elevation = face_elevation.to_owned();
    let slope_vec: Vec<_> = (0..n_face).into_par_iter()
        .map(|f| {
            let connected_edges = face_edge_connectivity.row(f);
            let slope_sq = compute_slope_squared(
                f,
                &face_elevation,
                &connected_edges,
                &edge_face_connectivity,
                &edge_face_distance,
                &edge_length,
            );
            slope_sq.sqrt()
        })
        .collect();
    let slope = ndarray::Array1::from_vec(slope_vec);

    Ok(PyArray1::from_owned_array(py, slope))
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
    face_area: PyReadonlyArray1<'py, f64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_area = face_area.as_array();
    let edge_face_connectivity = edge_face_connectivity.as_array();
    let face_edge_connectivity = face_edge_connectivity.as_array();
    let edge_face_distance = edge_face_distance.as_array();
    let edge_length = edge_length.as_array();
    let face_elevation_view = face_elevation.as_array();
    let n_face = face_edge_connectivity.nrows();
    let critical_slope_sq = critical_slope * critical_slope;

    // face_kappa_ones should be an array view
    let face_kappa = ndarray::Array1::<f64>::ones(n_face);
    let kappa_ref: &ndarray::ArrayBase<ndarray::ViewRepr<&f64>, ndarray::Dim<[usize; 1]>> = &face_kappa.view();

    let diffmax = compute_dt_max(
        &edge_face_distance,
        &edge_face_connectivity,
        &kappa_ref,
        &face_area,
        &edge_length,
    );
    let looplimit = 1000 as usize;

    let mut face_elevation = ndarray::Array1::<f64>::zeros(n_face);
    let mut face_delta_elevation = ndarray::Array1::<f64>::zeros(n_face);

    for _ in 0..looplimit {
        face_elevation.assign(&face_elevation_view);
        for f in 0..n_face {
            face_elevation[f] += face_delta_elevation[f];
        }

        // Compute (slope_sq, kappa) for each face in parallel
        let slope_kappa: Vec<_> = (0..n_face)
            .into_par_iter()
            .map(|f| {
                let connected_edges = face_edge_connectivity.row(f);
                let slope_sq = compute_slope_squared(
                    f,
                    &face_elevation,
                    &connected_edges,
                    &edge_face_connectivity,
                    &edge_face_distance,
                    &edge_length,
                );
                let kappa = if slope_sq > critical_slope_sq {
                    diffmax * (1.0 + slope_sq / critical_slope_sq) * 10.0
                } else {
                    0.0
                };
                kappa
            })
            .collect();

        let face_kappa: Vec<f64> = slope_kappa.into_iter().collect();

        let n_active = face_kappa.iter().filter(|&&k| k > 0.0).count();
        if n_active == 0 {
            break;
        }

        // Cast all of the arrays to the correct types for the Python bindings so that they can be passed to the apply_diffusion function
        let py_kappa = PyArray1::from_slice(py, &face_kappa).readonly();
        let py_face_elevation = PyArray1::from_array(py, &face_elevation).readonly();
        let py_face_area = PyArray1::from_array(py, &face_area).readonly();
        let py_edge_face_connectivity = PyArray2::from_array(py, &edge_face_connectivity).readonly();
        let py_edge_face_distance = PyArray1::from_array(py, &edge_face_distance).readonly();
        let py_edge_lengths = PyArray1::from_array(py, &edge_length).readonly();

        let delta = apply_diffusion(
            py,
            py_kappa,
            py_face_elevation,
            py_face_area,
            py_edge_face_connectivity,
            py_edge_face_distance,
            py_edge_lengths,
        )?;
        let delta_array = unsafe { delta.as_array() };
        face_delta_elevation += &delta_array;
    }
    Ok(PyArray1::from_owned_array(py, face_delta_elevation))
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
    let face_area = face_area.as_array();
    let face_elevation = face_elevation.as_array();
    let node_face_connectivity = node_face_connectivity.as_array();

    let n_nodes = node_face_connectivity.nrows();
    let mut result = ndarray::Array1::<f64>::zeros(n_nodes);

    for (node_id, row) in node_face_connectivity.outer_iter().enumerate() {
        let mut weighted_sum = 0.0;
        let mut area_sum = 0.0;

        for &face_id in row.iter() {
            if face_id < 0 {
                continue;
            }
            let f = face_id as usize;
            weighted_sum += face_elevation[f] * face_area[f];
            area_sum += face_area[f];
        }

        if area_sum > 0.0 {
            result[node_id] = weighted_sum / area_sum;
        }
    }

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
    let x = x.as_array();
    let y = y.as_array();
    let z = z.as_array();
    let anchor = anchor.as_array();

    // Get the maximum value in the x array as the scale
    let n_points = x.len();
    let mut result = ndarray::Array1::<f64>::zeros(n_points);
    let num_octaves = anchor.nrows();

    let mut norm = 0.5;
    for i in 0..num_octaves {
        let spatial_fac = freq.powi(i as i32) / noise_width;
        let noise_mag = pers.powi(i as i32);
        let rot_x = anchor[[i, 0]];
        let rot_y = anchor[[i, 1]];
        let rot_z = anchor[[i, 2]];
        norm += 0.5 * noise_mag;

        //let noise_source = Source::simplex(32345142).scale([spatial_fac, spatial_fac, spatial_fac]).rotate([rot_x,rot_y,rot_z]);
        let base = SuperSimplex::new(seed);
        let scaled = ScalePoint::new(base).set_scale(spatial_fac);
        let noise_source = RotatePoint::new(scaled).set_angles(rot_x, rot_y, rot_z, 0.0);

        Zip::from(&mut result)
            .and(&x)
            .and(&y)
            .and(&z)
            .into_par_iter()
            .for_each(|(r, &xv, &yv, &zv)| {
                *r += noise_source.get([xv, yv, zv]) * noise_mag;
            });
    }

    for val in result.iter_mut() {
        *val *= noise_height / norm;
    }

    Ok(PyArray1::from_owned_array(py, result))
}

/// Constructs the edge-other distances for a surface mesh
///
/// This will compute the haversine distance between the coordinates of the two components (faces, nodes, etc) that are associated with each edge.
///
/// # Arguments
/// * `py` - Python interpreter token.
/// * `edge_connectivity` - Indices of the elements associated each edge. (2D array of shape n_edge x 2).
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

    let edge_connectivity = edge_connectivity.as_array();
    let lon = lon.as_array();
    let lat = lat.as_array();
    let n_edge = edge_connectivity.nrows();

    let mut edge_distances = ndarray::Array1::<f64>::zeros(n_edge);

    Zip::from(&mut edge_distances)
        .and(edge_connectivity.outer_iter())
        .into_par_iter()
        .for_each(|(out, other)| {
            let o1 = other[0];
            let o2 = other[1];
            if o1 < 0 || o2 < 0 {
                *out = 0.0;
            } else {
                let o1 = o1 as usize;
                let o2 = o2 as usize;
                *out = haversine_distance_scalar(
                    lon[o1], lat[o1],
                    lon[o2], lat[o2],
                    radius,
                );
            }
        });

    Ok(PyArray1::from_owned_array(py, edge_distances))
}
