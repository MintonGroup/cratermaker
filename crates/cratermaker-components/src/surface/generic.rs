use std::f64::consts::PI;
use noise::{NoiseFn, RotatePoint, ScalePoint, SuperSimplex};
use numpy::ndarray::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use crate::ArrayResult;
use super::types::LocalSurfaceView;


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
/// * `face_kappa` - Topographic diffusivity (1D array of length n_face)
/// * `face_variable` - The variable to diffuse at each face (1D array of length n_face). Typically it will be region.face_elevation
/// * `region` - A LocalSurfaceView object representing the local mesh region.
///
/// # Returns
///
/// A NumPy array of shape (n_face,) giving the elevation change per face of the local mesh for this update step.
/// 
pub fn apply_diffusion(
    face_kappa: ArrayView1<'_, f64>,
    face_variable: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {
    // Compute initial dt von neumann stability condition
    let dt_max = compute_dt_max(
        face_kappa,
        region
    );

    let mut face_delta = Array1::<f64>::zeros(region.n_face);

    if !dt_max.is_finite() || dt_max <= 0.0 {
        return Ok(face_delta);
    }

    let nloops = (1.0 / dt_max).ceil() as usize;
    let dt = 1.0 / nloops as f64;
    let fac = dt / 2.0;


    for _ in 0..nloops {
        // Initialize dhdt as zeros with length n_face
        let mut dhdt = Array1::<f64>::zeros(region.n_face);

        // Loop over edges, accumulate flux contributions to each face
        for (e, faces) in region.edge_face_connectivity.outer_iter().enumerate() {

            let (f1, f2) = match extract_edge_faces(faces, e, region.n_face, "apply_diffusion") {
                Some(pair) => pair,
                None => continue,
            };

            let distance = region.edge_face_distance[e];
            let length = region.edge_length[e];

            debug_assert!(distance > 0.0, "non-positive distance at edge {}", e);


            let h1 = face_variable[f1] + face_delta[f1];
            let h2 = face_variable[f2] + face_delta[f2];
            let k1 = face_kappa[f1];
            let k2 = face_kappa[f2];

            let flux = (k1 + k2) * (h2 - h1) / distance * length;

            dhdt[f1] += flux / region.face_area[f1];
            dhdt[f2] -= flux / region.face_area[f2];
        }

        for f in 0..region.n_face {
            face_delta[f] += fac * dhdt[f];
        }
    }

    Ok(face_delta)
}

/// Computes the gradient vector in a radial direction defined by the face bearing at a face using the Green-Gauss method.
///
///
/// This function is designed to be parallel and returns a NumPy array of slopes
/// corresponding to the provided face indices.
///
/// # Arguments
/// * `variable` - The variable to compute the gradient for at each face (1D array).
/// * `region` - Reference to the local surface data structure containing mesh information. 
///
/// # Returns
/// An arrays of radial gradient values same length as `face_indices`.
/// 
pub fn compute_radial_gradient(
    variable: ArrayView1<'_,f64>,
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {
    let bearing = region.face_bearing.as_ref().ok_or("face_bearing required")?;
    let radgrad: Vec<f64> = (0..region.n_face).into_par_iter()
        .map(|f| {
            let (grad_zonal, grad_meridional) = compute_one_face_gradient(
                f,
                variable,
                region
            );
            grad_meridional * bearing[f].cos() + grad_zonal * bearing[f].sin()
        })
        .collect();
    Ok(Array1::from_vec(radgrad))
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
/// 
/// # Arguments
/// * `f` - Index of the face to compute the slope squared for.
/// * `face_elevation` - The variable to compute the gradient for at each face (1D array).
/// * `region` - Reference to the local surface data structure containing mesh information.
/// 
/// # Returns
/// 
/// The squared magnitude of the gradient vector at face `f`.
/// 
fn compute_face_slope_squared(
    f: usize,
    face_elevation: ArrayView1<'_,f64>,
    region: &LocalSurfaceView<'_>,
) -> f64 {

    let (dh_dx, dh_dy) = compute_one_face_gradient(
        f,
        face_elevation,
        region
    ); 
    return dh_dx * dh_dx + dh_dy * dh_dy;
    
}


/// Computes the gradient vector (∂h/∂x, ∂h/∂y) at a face using the Green-Gauss method.
///
/// For each face, this function estimates the gradient vector (∂h/∂x, ∂h/∂y)
/// in the local tangent plane by summing contributions from all connected neighbors,
/// using the Green-Gauss theorem applied to the face's surrounding edges.
/// Returns the squared magnitude of this gradient as the slope squared.
///
/// For each neighbor, we take the vector (dx, dy) in the tangent plane and the elevation difference dh,
/// and estimate the gradient using the Green-Gauss formula.
/// 
/// # Arguments
/// * `f` - Index of the face to compute the gradient for.
/// * `variable` - The variable to compute the gradient for at each face (1D array same length as number of faces).
/// * `region` - Reference to the local surface data structure containing mesh information.
/// 
/// # Returns
/// A tuple representing the zonal and meridional components of the gradient vector at face `f`.
/// 
#[inline(always)]
fn compute_one_face_gradient(
    f: usize,
    variable: ArrayView1<'_,f64>,
    region: &LocalSurfaceView<'_>,
) -> (f64, f64) {
    let connected_edges: ArrayView1<'_, i64> = region.face_edge_connectivity.row(f);
    let mut dh_zonal = 0.0;
    let mut dh_meridional = 0.0;
    let mut weight_sum = 0.0;
    let lon_f = region.face_lon[f];
    let lat_f = region.face_lat[f];

    // This will gather up the distance to each neighboring edge ands its length for the flux calculation
    let mut neighbors = Vec::new();
    for &edge_id in connected_edges.iter() {
        if edge_id < 0 {
            continue;
        }
        let e = edge_id as usize;

        let (f1, f2) = match extract_edge_faces(region.edge_face_connectivity.row(e), e, region.n_face, "compute_one_face_gradient") {
            Some(pair) => pair,
            None => continue,
        };

        let other = if f == f1 { f2 } else if f == f2 { f1 } else { continue };
        let d = region.edge_face_distance[e];
        let l = region.edge_length[e];
        if d <= 0.0 || l <= 0.0 {
            continue;
        }
        neighbors.push((other, d, l));
    }

    let n = neighbors.len();
    if n < 2 {
        return (0.0, 0.0);
    }

    // Loop over neighbors to compute gradient contributions
    for (_i, &(other, distance, length)) in neighbors.iter().enumerate() {
        let theta = compute_one_bearing(
            lon_f,
            lat_f,
            region.face_lon[other],
            region.face_lat[other],
        );  
        let dir_zonal = theta.sin();
        let dir_meridional = theta.cos();
        let dz = variable[other] - variable[f];
        let flux = dz * length / distance;

        dh_zonal += flux * dir_zonal;
        dh_meridional += flux * dir_meridional;
        weight_sum += length;
    }

    if weight_sum > 0.0 {
        dh_zonal /= weight_sum;
        dh_meridional /= weight_sum;
        return (dh_zonal, dh_meridional);
    }

    (0.0, 0.0)
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
/// * `region` - Reference to the local surface data structure containing mesh information.
///
/// # Returns
/// A NumPy array of slope values (1D array), same length as `face_indices`.
/// 
pub fn compute_slope(
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {

    let n_face = region.n_face;
    let slope_vec: Vec<_> = (0..n_face).into_par_iter()
        .map(|f| {
            let slope_sq = compute_face_slope_squared(
                f,
                region.face_elevation.view(),
                region
            );
            slope_sq.sqrt()
        })
        .collect();
    Ok(Array1::from_vec(slope_vec))
}


/// Computes the spatially varying diffusivity (`face_kappa`) for a slope collapse step.
///
/// Sets `kappa = diffmax` if any neighbor violates the critical slope threshold,
/// otherwise sets `kappa = 0.0`. `diffmax` is computed assuming a stable timestep of 1.0.
///
/// # Arguments
///
/// * `critical_slope` - Maximum allowable slope (e.g., 0.7 for ~35 degrees).
/// * `region` - Reference to the local surface data structure containing mesh information.
///
/// # Returns
///
/// A NumPy array of `face_kappa` values.
/// 
pub fn slope_collapse(
    critical_slope: f64,
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {
    let n_face = region.n_face;
    let critical_slope_sq = critical_slope * critical_slope;

    // face_kappa_ones should be an array view
    let face_kappa_ones = Array1::<f64>::ones(n_face);

    let diffmax = compute_dt_max(
        face_kappa_ones.view(),
        region,
    );
    let looplimit = 1000 as usize;

    let mut new_face_elevation = region.face_elevation.to_owned();
    let mut face_delta_elevation = Array1::<f64>::zeros(n_face);

    for _ in 0..looplimit {
        new_face_elevation.assign(&region.face_elevation);
        new_face_elevation += &face_delta_elevation;

        // Compute (slope_sq, kappa) for each face in parallel
        let slope_kappa: Vec<_> = (0..n_face)
            .into_par_iter()
            .map(|f| {
                let slope_sq = compute_face_slope_squared(
                    f,
                    new_face_elevation.view(),
                    region,
                );
                let kappa = if slope_sq > critical_slope_sq {
                    diffmax * (1.0 + slope_sq / critical_slope_sq) * 10.0
                } else {
                    0.0
                };
                kappa
            })
            .collect();

        let face_kappa=  Array1::from_vec(slope_kappa);

        let n_active = face_kappa.iter().filter(|&&k| k > 0.0).count();
        if n_active == 0 {
            break;
        }

        let delta = apply_diffusion(
            face_kappa.view(),
            new_face_elevation.view(),
            region,
        )?;
        face_delta_elevation += &delta;
    }
    Ok(face_delta_elevation)
}

/// Computes node elevations as area-weighted averages of adjacent face elevations.
///
/// # Arguments
///
/// * `region` - A LocalSurfaceView object representing the local mesh region.
///
/// # Returns
///
/// A NumPy array of node elevations (1D array).
/// 
pub fn interpolate_node_elevation_from_faces(
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {

    let mut result = Array1::<f64>::zeros(region.n_node);

    for (node_id, row) in region.node_face_connectivity.outer_iter().enumerate() {
        let mut weighted_sum = 0.0;
        let mut area_sum = 0.0;

        for &face_id in row.iter() {
            if face_id < 0 {
                continue;
            }
            let f = face_id as usize;
            weighted_sum += region.face_elevation[f] * region.face_area[f];
            area_sum += region.face_area[f];
        }

        if area_sum > 0.0 {
            result[node_id] = weighted_sum / area_sum;
        }
    }

    Ok(result)
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
pub fn turbulence_noise(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    z: ArrayView1<'_, f64>,
    noise_height: f64,
    noise_width: f64,
    freq: f64,
    pers: f64,
    anchor: ArrayView2<'_, f64>,
    seed: u32,
) -> ArrayResult {

    // Get the maximum value in the x array as the scale
    let n_points = x.len();
    let mut result = Array1::<f64>::zeros(n_points);
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


        let noise_values: Vec<f64> = (0..n_points)
            .into_par_iter()
            .map(|i| {
                let xv = x[i];
                let yv = y[i];
                let zv = z[i];
                noise_source.get([xv, yv, zv]) * noise_mag
            })
            .collect();
        for (r, &val) in result.iter_mut().zip(noise_values.iter()) {
            *r += val;
        }
    }

    for val in result.iter_mut() {
        *val *= noise_height / norm;
    }

    Ok(result)
}




/// Constructs the edge-other distances for a surface mesh
///
/// This will compute the haversine distance between the coordinates of the two components (faces, nodes, etc) that are associated with each edge.
///
/// # Arguments
/// * `edge_connectivity` - Indices of the elements associated each edge. (2D array of shape n_edge x 2).
/// * `lon1`, `lat1` - Coordinates of the reference point in radians.
///
/// # Returns
/// A 1D NumPy array of edge-other distances in meters, one for each edge.
pub fn compute_edge_distances(
    edge_connectivity: ArrayView2<'_, i64>,
    lon: ArrayView1<'_, f64>,
    lat: ArrayView1<'_, f64>,
    radius: f64,
) -> ArrayResult {

    let n_edge = edge_connectivity.nrows();
    let edge_distances_vec: Vec<f64> = (0..n_edge)
        .into_par_iter()
        .map(|e| {
            let other = edge_connectivity.row(e);
            let o1 = other[0];
            let o2 = other[1];
            if o1 < 0 || o2 < 0 {
                0.0
            } else {
                let o1u = o1 as usize;
                let o2u = o2 as usize;
                compute_one_distance(
                    lon[o1u], lat[o1u],
                    lon[o2u], lat[o2u],
                    radius,
                )
            }
        })
        .collect();
    Ok(Array1::from_vec(edge_distances_vec))
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
pub fn compute_distances(
    lon1: f64,
    lat1: f64,
    lon2: ArrayView1<'_, f64>,
    lat2: ArrayView1<'_, f64>,
    radius: f64,
) -> ArrayResult {

    let n = lon2.len();

    let result_vec: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let lon2_i = lon2[i];
            let lat2_i = lat2[i];
            compute_one_distance(lon1, lat1, lon2_i, lat2_i, radius)
        })
        .collect();

    Ok(Array1::from_vec(result_vec))
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
fn compute_one_distance(lon1: f64, lat1: f64, lon2: f64, lat2: f64, radius: f64) -> f64 {
    let dlon = lon2 - lon1;
    let dlat = lat2 - lat1;
    let a = (dlat / 2.0).sin().powi(2)
        + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();
    radius * c
}


/// Computes the initial bearing (forward azimuth) from a fixed point to each of a set of destination points.
///
/// The bearing is calculated on a spherical surface using great-circle paths and returned in radians,
/// normalized to the range [0, 2π).
///
/// # Arguments
///
/// * `lon1` - Longitude of the reference point, in radians.
/// * `lat1` - Latitude of the reference point, in radians.
/// * `lon2` - Longitudes of destination points, in radians.
/// * `lat2` - Latitudes of destination points, in radians.
///
/// # Returns
///
/// * A NumPy array of initial bearing angles (radians), one for each (lon2, lat2) pair.
pub fn compute_bearings(
    lon1: f64,
    lat1: f64,
    lon2: ArrayView1<'_, f64>,
    lat2: ArrayView1<'_, f64>,
) -> ArrayResult {
    let n = lon2.len();
    let result_vec: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let lon2_i = lon2[i];
            let lat2_i = lat2[i];
            let initial_bearing = compute_one_bearing(lon1, lat1, lon2_i, lat2_i);
            // Normalize bearing to 0 to 2*pi
            (initial_bearing + 2.0 * PI) % (2.0 * PI)
        })
        .collect();
    Ok(Array1::from_vec(result_vec))
}

/// Computes the initial bearing (forward azimuth) from point 1 to point 2 on a sphere.
#[inline]
fn compute_one_bearing(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let dlon = positive_mod(lon2 - lon1 + PI, 2.0 * PI) - PI;
    let x = dlon.sin() * lat2.cos();
    let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
    x.atan2(y)
}


/// Computes the positive modulus of `x` with respect to `m`.
/// Ensures the result is always in the range `[0, m)`.
#[inline]
fn positive_mod(x: f64, m: f64) -> f64 {
    ((x % m) + m) % m
}


/// Computes a conservative explicit timestep for stability using per-face aggregated flux contributions.
///
/// This method loops over all edges and accumulates inverse dt estimates for each face.
/// The final dt is the minimum over all faces of the reciprocal flux sum.
/// 
/// # Arguments
/// * `face_kappa` - Topographic diffusivity (1D array of length n_face)
/// * `region` - A LocalSurfaceView object representing the local mesh region.
/// 
/// # Returns
/// A single f64 value giving the maximum stable timestep for the diffusion operation.
/// 
fn compute_dt_max(
    face_kappa: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
) -> f64 {
    let mut inverse_dt_sum = Array1::<f64>::zeros(region.n_face);

    for (e, faces) in region.edge_face_connectivity.outer_iter().enumerate() {
        let distance = region.edge_face_distance[e];
        if distance <= 0.0 {
            continue;
        }

        let (f1, f2) = match extract_edge_faces(faces, e, region.n_face, "compute_dt_max") {
            Some(pair) => pair,
            None => continue,
        };

        let k1 = face_kappa[f1];
        let k2 = face_kappa[f2];
        let k_avg = 0.5 * (k1 + k2);
        let length = region.edge_length[e];
        if k_avg > 0.0 {
            let contrib_f1 = 2.0 * k_avg * length / (distance * region.face_area[f1]);
            let contrib_f2 = 2.0 * k_avg * length / (distance * region.face_area[f2]);
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


use numpy::ndarray::ArrayView1; // you already have prelude, so this is just for clarity

/// Safely extract a pair of face indices (f1, f2) from an edge-face row.
///
/// Returns `Some((f1, f2))` if:
/// - the row is contiguous,
/// - it has exactly two entries,
/// - both are within `[0, n_face)`.
/// Otherwise logs a message with `context` and returns `None`.
fn extract_edge_faces(
    row: ArrayView1<'_, i64>,
    edge_id: usize,
    n_face: usize,
    context: &str,
) -> Option<(usize, usize)> {
    let slice = match row.as_slice() {
        Some(s) => s,
        None => {
            eprintln!(
                "{}: non-contiguous edge_face row at edge {}, skipping",
                context, edge_id
            );
            return None;
        }
    };

    let [f1_raw, f2_raw] = match <[i64; 2]>::try_from(slice) {
        Ok(pair) => pair,
        Err(_) => {
            eprintln!(
                "{}: invalid face pair at edge {}, slice = {:?}",
                context, edge_id, slice
            );
            return None;
        }
    };

    let f1 = f1_raw as usize;
    let f2 = f2_raw as usize;
    if f1 >= n_face || f2 >= n_face {
        eprintln!(
            "{}: face index out of range at edge {}, f1={}, f2={}, n_face={}",
            context, edge_id, f1, f2, n_face
        );
        return None;
    }

    Some((f1, f2))
}