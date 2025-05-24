use std::f64::consts::PI;
use std::collections::HashMap;

use ndarray::Zip;
use noise::{NoiseFn, RotatePoint, ScalePoint, SuperSimplex};
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyArrayMethods};
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
pub fn calculate_haversine_distance<'py>(
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


/// Builds a mapping from face indices to their corresponding positions in a result array.
/// Used to retrieve local indices efficiently from global face IDs.
fn compute_face_index_map(face_indices: &ndarray::ArrayView1<'_, i64>) -> HashMap<usize, usize> {
    face_indices
        .iter()
        .enumerate()
        .map(|(i, &f)| (f as usize, i))
        .collect()
}

/// Computes the stable time step for the diffusion equation based on face areas and maximum diffusivity.
///
/// # Arguments
/// * `face_areas` - Array of face areas.
/// * `max_kappa` - Maximum diffusivity value.
/// * `n_neighbors` - Number of neighboring faces per face.
///
/// # Returns
/// The minimum allowed timestep for stability.
fn compute_dt_initial(face_areas: &ndarray::ArrayView1<'_, f64>, max_kappa: f64, n_neighbors: usize) -> f64 {
    face_areas
        .iter()
        .map(|&a| a / (2.0 * max_kappa * n_neighbors as f64))
        .fold(f64::INFINITY, f64::min)
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
    let face_index_map = compute_face_index_map(&face_indices);

    // Compute max kappa for stability condition
    let max_kappa = face_kappa.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Compute initial dt using per-face stability condition
    let dt_initial = compute_dt_initial(&face_areas, max_kappa, n_max_face_faces);

    if !dt_initial.is_finite() || dt_initial <= 0.0 {
        let result = ndarray::Array1::<f64>::zeros(face_indices.len());
        return Ok(PyArray1::from_owned_array(py, result));
    }

    let nloops = (1.0 / dt_initial).ceil() as usize;
    let dt = 1.0 / nloops as f64;

    let mut face_delta_elevation = ndarray::Array1::<f64>::zeros(face_indices.len());

    for _ in 0..nloops {
        let dhdt: Vec<_> = Zip::from(&face_indices)
            .and(face_face_connectivity.outer_iter())
            .and(&face_delta_elevation)
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
                            .map_or(0.0, |&j| face_delta_elevation[j]);
                    let k_n = face_kappa[f_n];
                    let flux = (k_f + k_n) * (h_n - h_f);

                    dhdt_f += flux;
                }
                dhdt_f
            })
            .collect();

        for (i, &f) in face_indices.iter().enumerate() {
            face_delta_elevation[i] += dt * dhdt[i] / face_areas[f as usize];
        }
    }

    Ok(PyArray1::from_owned_array(py, face_delta_elevation))
}


/// Computes the maximum squared slope magnitude at a face using adjacent neighbor pairs.
///
/// This function assumes neighbors are ordered counterclockwise and uses Haversine distances
/// to determine gradient magnitudes.
///
/// # Arguments
/// * `f` - Index of the central face.
/// * `row` - Neighbor indices for the central face.
/// * `face_elevation` - Elevation array for all faces.
/// * `face_lon`, `face_lat` - Longitude and latitude arrays for all faces.
/// * `radius` - Radius of the sphere.
///
/// # Returns
/// Maximum squared slope from any adjacent neighbor pair around the face.
#[inline(always)]
fn compute_slope_squared(
    f: usize,
    row: &ndarray::ArrayView1<'_, i64>,
    face_elevation: &ndarray::Array1<f64>,
    face_lon: &ndarray::ArrayView1<'_, f64>,
    face_lat: &ndarray::ArrayView1<'_, f64>,
    radius: f64,
) -> f64 {
    let h_f = face_elevation[f];
    let mut max_slope_sq = 0.0;
    // Collect valid neighbor indices into a small array for pairing with wrapping.
    let mut valid_neighbors = [0usize; 24];
    let mut count = 0;
    for &fid in row.iter() {
        if fid >= 0 {
            let fid = fid as usize;
            if fid < face_elevation.len() {
                valid_neighbors[count] = fid;
                count += 1;
            } else {
                panic!("Neighbor index {} out of bounds (face_elevation.len() = {})", fid, face_elevation.len());
            }
        }
    }

    for i in 0..count {
        let j = (i + 1) % count;
        let f_j = valid_neighbors[i];
        let f_k = valid_neighbors[j];

        let h_j = face_elevation[f_j];
        let h_k = face_elevation[f_k];

        let slope_j = h_j - h_f;
        let slope_k = h_k - h_f;

        let d_j = haversine_distance_scalar(face_lon[f], face_lat[f], face_lon[f_j], face_lat[f_j], radius);
        let d_k = haversine_distance_scalar(face_lon[f], face_lat[f], face_lon[f_k], face_lat[f_k], radius);

        if d_j == 0.0 || d_k == 0.0 {
            continue;
        }

        let grad_j = slope_j / d_j;
        let grad_k = slope_k / d_k;
        let slope_sq = grad_j * grad_j + grad_k * grad_k;

        if slope_sq > max_slope_sq {
            max_slope_sq = slope_sq;
        }
    }

    max_slope_sq
}


/// Computes the spatially varying diffusivity (`face_kappa`) for a slope collapse step.
///
/// Sets `kappa = diffmax` if any neighbor violates the critical slope threshold,
/// otherwise sets `kappa = 0.0`. `diffmax` is computed assuming a stable timestep of 1.0.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_areas` - Area of each face (1D array).
/// * `face_elevation` - Elevation at each face (1D array).
/// * `face_face_connectivity` - Neighboring faces (2D array).
/// * `face_indices` - Subset of face indices (1D array).
/// * `critical_slope` - Maximum allowable slope (e.g., 0.7 for ~35 degrees).
/// * `face_lon` - Longitudes of faces (1D array).
/// * `face_lat` - Latitudes of faces (1D array).
/// * `radius` - Radius of the sphere.
///
/// # Returns
///
/// A NumPy array of `face_kappa` values.
#[pyfunction]
pub fn slope_collapse<'py>(
    py: Python<'py>,
    face_areas: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_face_connectivity: PyReadonlyArray2<'py, i64>,
    face_indices: PyReadonlyArray1<'py, i64>,
    face_lon: PyReadonlyArray1<'py, f64>,
    face_lat: PyReadonlyArray1<'py, f64>,
    radius: f64,
    critical_slope: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_areas_view = face_areas.as_array();
    let face_elevation_view = face_elevation.as_array();
    let face_face_connectivity_view = face_face_connectivity.as_array();
    let face_indices_view = face_indices.as_array();
    let face_lon_view = face_lon.as_array();
    let face_lat_view = face_lat.as_array();
    let n_max_face_faces = face_face_connectivity_view.ncols();
    let n_face = face_areas_view.len();


    let diffmax = compute_dt_initial(&face_areas_view, 1.0, n_max_face_faces);
    let looplimit = n_face as usize;

    let mut global_kappa = vec![0.0f64; n_face];
    let mut face_elevation = ndarray::Array1::<f64>::zeros(face_elevation_view.len());
    let mut face_delta_elevation = ndarray::Array1::<f64>::zeros(face_indices_view.len());

    for _ in (0..looplimit).rev() {
        face_elevation.assign(&face_elevation_view);
        for (i, &f) in face_indices_view.iter().enumerate() {
            face_elevation[f as usize] += face_delta_elevation[i];
        }
        let face_kappa: Vec<_> = Zip::from(&face_indices_view)
            .and(face_face_connectivity_view.outer_iter())
            .into_par_iter()
            .map(|(f, row)| {
                let f = *f as usize;
                let slope_sq = compute_slope_squared(
                    f,
                    &row,
                    &face_elevation,
                    &face_lon_view,
                    &face_lat_view,
                    radius,
                );
                if slope_sq > critical_slope {
                    diffmax
                } else {
                    0.0
                }
            })
            .collect();

        let n_active = face_kappa.iter().filter(|&&k| k > 0.0).count();

        if n_active == 0 {
            break;
        }
    
        for (i, &f) in face_indices_view.iter().enumerate() {
            let f = f as usize;
            assert!(f < global_kappa.len(), "f {} out of bounds for global_kappa (len = {})", f, global_kappa.len());
            global_kappa[f] = face_kappa[i];
        }
        assert_eq!(face_elevation.len(), n_face, "face_elevation length mismatch");
        let py_kappa_array = PyArray1::from_slice(py, &global_kappa);
        let py_kappa = py_kappa_array.readonly();
        let py_elev_array = PyArray1::from_array(py, &face_elevation); 
        let py_elev = py_elev_array.readonly();

        let delta = apply_diffusion(
            py,
            face_areas.clone(),
            py_kappa,
            py_elev,
            face_face_connectivity.clone(),
            face_indices.clone(),
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
