use std::f64::consts::PI;
use std::collections::HashMap;

use ndarray::Zip;
use noise::{NoiseFn, RotatePoint, ScalePoint, SuperSimplex};
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

#[inline]
fn positive_mod(x: f64, m: f64) -> f64 {
    ((x % m) + m) % m
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
pub fn calculate_initial_bearing<'py>(
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
            // Calculate differences in coordinates
            let dlon = positive_mod(lon2 - lon1 + PI, 2.0 * PI) - PI;

            // Haversine formula calculations
            let x = dlon.sin() * lat2.cos();
            let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
            let initial_bearing = f64::atan2(x, y);

            // Normalize bearing to 0 to 2*pi
            *out = (initial_bearing + 2.0 * PI) % (2.0 * PI);
        });

    Ok(PyArray1::from_owned_array(py, result))
}

/// View into a region of the surface mesh, consisting of node and face indices.
///
/// Used to localize crater effects to a subset of the full mesh.
#[derive(FromPyObject)]
pub struct SurfaceView<'py> {
    pub node_indices: PyReadonlyArray1<'py, i64>,
    pub face_indices: PyReadonlyArray1<'py, i64>,
}

/// Applies one explicit diffusion update step over a surface mesh with variable diffusivity.
///
/// This function computes the change in elevation for each face on the mesh using a
/// face-centered finite-volume formulation of the operator:
///     ∂h/∂t = ∇ · (κ ∇h)
/// where κ varies per face. For each face, the flux with its neighbors is computed
/// using the expression (κ_f + κ_n) * (h_n - h_f), summing over all neighbors n.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_areas` - Area of each face (1D array of length n_face).
/// * `face_kappa` - Topographic diffusivity at each face (1D array of length n_face).
/// * `face_elevation` - Elevation value at each face (1D array of length n_face).
/// * `face_face_connectivity` - For each face, the indices of its neighboring faces (2D array).
///
/// # Returns
///
/// A NumPy array of shape (n_face,) giving the elevation change per face for this update step.
#[pyfunction]
pub fn apply_diffusion<'py>(
    py: Python<'py>,
    face_areas: PyReadonlyArray1<'py, f64>,
    face_kappa: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_indices: PyReadonlyArray1<'py, i64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_areas = face_areas.as_array();
    let face_kappa = face_kappa.as_array();
    let face_elevation = face_elevation.as_array();
    let face_face_connectivity = face_face_connectivity.as_array();
    let n_max_face_faces = face_face_connectivity.ncols();
    let n_face = face_areas.len();
    let face_indices = face_indices.as_array();

    // This will map the face index to the local index in the result array
    let face_index_map: HashMap<usize, usize> = face_indices
        .iter()
        .enumerate()
        .map(|(i, &f)| (f as usize, i))
        .collect();

    // Compute max kappa for stability condition
    let max_kappa = face_kappa.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Compute initial dt using per-face stability condition
    let dt_initial = face_areas
        .iter()
        .map(|&a| a / (2.0 * max_kappa * n_max_face_faces as f64))
        .fold(f64::INFINITY, f64::min);

    if !dt_initial.is_finite() || dt_initial <= 0.0 {
        let result = ndarray::Array1::<f64>::zeros(face_indices.len());
        return Ok(PyArray1::from_owned_array(py, result));
    }

    let nloops = (1.0 / dt_initial).ceil() as usize;
    let dt = 1.0 / nloops as f64;

    let mut result = ndarray::Array1::<f64>::zeros(face_indices.len());

    for _ in 0..nloops {
        let dhdt: Vec<_> = Zip::from(&face_indices)
            .and(face_face_connectivity.outer_iter())
            .and(&result)
            .into_par_iter()
            .map(|(f, row, face_delta)| {
                let f = *f as usize;
                let h_f = face_elevation[f] + face_delta;
                let k_f = face_kappa[f];
                let mut dhdt_f = 0.0;
                for &f_n_raw in row {
                    if f_n_raw < 0 {
                        continue;
                    }
                    let f_n = f_n_raw as usize;
                    if f_n >= n_face {
                        continue;
                    }
                    let h_n = face_elevation[f_n]
                        + face_index_map
                            .get(&f_n)
                            .map_or(0.0, |&j| result[j]);
                    let k_n = face_kappa[f_n];
                    let flux = (k_f + k_n) * (h_n - h_f);

                    dhdt_f += flux;
                }
                dhdt_f
            })
            .collect();

        for (i, &f) in face_indices.iter().enumerate() {
            result[i] += dt * dhdt[i] / face_areas[f as usize];
        }
    }

    Ok(PyArray1::from_owned_array(py, result))
}

/// Computes node elevations as area-weighted averages of adjacent face elevations.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_areas` - Area of each face (1D array).
/// * `face_elevation` - Elevation at each face (1D array).
/// * `node_face_connectivity` - For each node, indices of connected faces (2D array).
///
/// # Returns
///
/// A NumPy array of node elevations (1D array).
#[pyfunction]
pub fn interpolate_node_elevation_from_faces<'py>(
    py: Python<'py>,
    face_areas: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    node_face_connectivity: PyReadonlyArray2<'py, i64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_areas = face_areas.as_array();
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
            weighted_sum += face_elevation[f] * face_areas[f];
            area_sum += face_areas[f];
        }

        if area_sum > 0.0 {
            result[node_id] = weighted_sum / area_sum;
        }
    }

    Ok(PyArray1::from_owned_array(py, result))
}

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
