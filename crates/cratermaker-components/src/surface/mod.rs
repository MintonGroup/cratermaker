use crate::{ArrayResult, ArrayResult2D};
use itertools::izip;
use noise::{NoiseFn, RotatePoint, ScalePoint, SuperSimplex};
use numpy::ndarray::prelude::*;
use rayon::prelude::*;
use std::f64::consts::{PI, TAU};

/// Represents a local region of a surface mesh with various attributes accessible as array views.
pub struct LocalSurfaceView<'a> {
    pub n_face: usize,
    pub pix: f64,
    pub face_area: ArrayView1<'a, f64>,
    pub face_elevation: ArrayView1<'a, f64>,
    pub face_indices: ArrayView1<'a, i64>,
    pub face_lon: ArrayView1<'a, f64>,
    pub face_lat: ArrayView1<'a, f64>,
    pub face_x: ArrayView1<'a, f64>,
    pub face_y: ArrayView1<'a, f64>,
    pub face_z: ArrayView1<'a, f64>,

    pub n_node: usize,
    pub node_elevation: ArrayView1<'a, f64>,
    pub node_indices: ArrayView1<'a, i64>,
    pub node_lon: ArrayView1<'a, f64>,
    pub node_lat: ArrayView1<'a, f64>,
    pub node_x: ArrayView1<'a, f64>,
    pub node_y: ArrayView1<'a, f64>,
    pub node_z: ArrayView1<'a, f64>,

    pub n_edge: usize,
    pub edge_indices: ArrayView1<'a, i64>,
    pub edge_length: ArrayView1<'a, f64>,

    pub face_edge_connectivity: ArrayView2<'a, i64>,
    pub face_node_connectivity: ArrayView2<'a, i64>,
    pub face_face_connectivity: ArrayView2<'a, i64>,
    pub node_face_connectivity: ArrayView2<'a, i64>,
    pub edge_face_connectivity: ArrayView2<'a, i64>,
    pub edge_node_connectivity: ArrayView2<'a, i64>,
    pub edge_face_distance: ArrayView1<'a, f64>,

    // The following are optional, as they will not be present if the LocalSurfaceView represents a global Surface object
    pub face_proj_x: Option<ArrayView1<'a, f64>>,
    pub face_proj_y: Option<ArrayView1<'a, f64>>,
    pub node_proj_x: Option<ArrayView1<'a, f64>>,
    pub node_proj_y: Option<ArrayView1<'a, f64>>,
    pub face_distance: Option<ArrayView1<'a, f64>>,
    pub face_bearing: Option<ArrayView1<'a, f64>>,
    pub region_radius: Option<f64>,
    pub location: Option<(f64, f64)>,
    pub crater_id: Option<ArrayView2<'a, u64>>,
    pub desloped_face_elevation: Option<ArrayView1<'a, f64>>,
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
    apply_diffusion_with_topology(face_kappa, face_variable, region, &mut None)
}

/// Minimum number of items per rayon task in the diffusion loops. Small regions are processed on fewer threads, which avoids
/// fork/join overhead dominating the run time. This only affects scheduling, not results.
const DIFFUSION_MIN_PAR_LEN: usize = 2048;

/// Edge-to-face and face-to-edge incidence of a region, used by the diffusion update.
///
/// The face-to-edge lists are stored in CSR form, and each face's list is in increasing edge order (and, for the same edge,
/// the first face of the edge before the second). This is the same order in which the original edge-wise scatter added flux
/// contributions to each face, so gathering per face reproduces its floating-point summation order exactly.
pub struct DiffusionTopology {
    /// The two faces of each edge: [first, second]
    edge_faces: Vec<[usize; 2]>,
    /// incidence[offsets[f]..offsets[f + 1]] are the incidences of face f
    offsets: Vec<usize>,
    /// (edge index, true if the face is the second face of the edge)
    incidence: Vec<(usize, bool)>,
}

impl DiffusionTopology {
    /// The edges of face `f` as `(edge index, true if the face is the second face of the edge)`, in increasing edge order.
    fn incidence_of(&self, f: usize) -> &[(usize, bool)] {
        &self.incidence[self.offsets[f]..self.offsets[f + 1]]
    }

    pub fn new(region: &LocalSurfaceView<'_>) -> Self {
        let n_edges = region.edge_face_connectivity.nrows();
        let n_faces = region.n_face;
        let mut edge_faces = vec![[0usize; 2]; n_edges];
        let mut counts = vec![0usize; n_faces + 1];
        for e in 0..n_edges {
            let f1 = region.edge_face_connectivity[(e, 0)] as usize;
            let f2 = region.edge_face_connectivity[(e, 1)] as usize;
            assert!(f1 < n_faces && f2 < n_faces);
            edge_faces[e] = [f1, f2];
            counts[f1 + 1] += 1;
            counts[f2 + 1] += 1;
        }
        for f in 0..n_faces {
            counts[f + 1] += counts[f];
        }
        let offsets = counts;
        let mut cursor = offsets.clone();
        let mut incidence = vec![(0usize, false); 2 * n_edges];
        for (e, &[f1, f2]) in edge_faces.iter().enumerate() {
            incidence[cursor[f1]] = (e, false);
            cursor[f1] += 1;
            incidence[cursor[f2]] = (e, true);
            cursor[f2] += 1;
        }
        Self {
            edge_faces,
            offsets,
            incidence,
        }
    }
}

/// Same as `apply_diffusion`, but reuses (or lazily builds) the region's `DiffusionTopology`, so that repeated calls on the same
/// region (e.g. from `slope_collapse`) do not rebuild it.
fn apply_diffusion_with_topology(
    face_kappa: ArrayView1<'_, f64>,
    face_variable: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
    topology: &mut Option<DiffusionTopology>,
) -> ArrayResult {
    // Compute initial dt von neumann stability condition
    let dt_max = compute_dt_max(face_kappa, region);

    let mut face_delta = Array1::<f64>::zeros(region.n_face);

    if !dt_max.is_finite() || dt_max <= 0.0 {
        return Ok(face_delta);
    }

    let nloops = (1.0 / dt_max).ceil() as usize;
    let dt = 1.0 / nloops as f64;
    let fac = dt / 2.0;

    let topology = topology.get_or_insert_with(|| DiffusionTopology::new(region));
    let edge_faces = &topology.edge_faces;
    let offsets = &topology.offsets;
    let incidence = &topology.incidence;
    let mut flux = vec![0.0f64; edge_faces.len()];

    for _ in 0..nloops {
        // Flux across each edge, using the face_delta from the previous step
        flux.par_iter_mut()
            .with_min_len(DIFFUSION_MIN_PAR_LEN)
            .enumerate()
            .for_each(|(e, flux_e)| {
                let [f1, f2] = edge_faces[e];
                let distance = region.edge_face_distance[e];
                let length = region.edge_length[e];

                debug_assert!(distance > 0.0, "non-positive distance at edge {}", e);

                let h1 = face_variable[f1] + face_delta[f1];
                let h2 = face_variable[f2] + face_delta[f2];
                let k1 = face_kappa[f1];
                let k2 = face_kappa[f2];

                *flux_e = (k1 + k2) * (h2 - h1) / distance * length;
            });

        // Gather each face's flux contributions in edge order and update it
        face_delta
            .as_slice_mut()
            .expect("face_delta is contiguous")
            .par_iter_mut()
            .with_min_len(DIFFUSION_MIN_PAR_LEN)
            .enumerate()
            .for_each(|(f, delta_f)| {
                let area = region.face_area[f];
                let mut dhdt = 0.0f64;
                for &(e, is_second) in &incidence[offsets[f]..offsets[f + 1]] {
                    if is_second {
                        dhdt += -flux[e] / area;
                    } else {
                        dhdt += flux[e] / area;
                    }
                }
                *delta_f += fac * dhdt;
            });
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
    variable: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {
    let bearings = region
        .face_bearing
        .as_ref()
        .ok_or("face_bearing required")?;
    let bearings_rad = bearings.mapv(|b| b.to_radians());
    let radgrad: Vec<f64> = (0..region.n_face)
        .into_par_iter()
        .map(|f| {
            let (grad_zonal, grad_meridional) = compute_one_face_gradient(f, variable, region);
            let bearing = bearings_rad[f];
            grad_meridional * bearing.cos() + grad_zonal * bearing.sin()
        })
        .collect();
    Ok(Array1::from_vec(radgrad))
}

/// Computes the gradient vector in a azimuthal direction defined by the face bearing at a face using the Green-Gauss method.
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
pub fn compute_azimuthal_gradient(
    variable: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
) -> ArrayResult {
    let bearings = region
        .face_bearing
        .as_ref()
        .ok_or("face_bearing required")?;
    let bearings_rad = bearings.mapv(|b| b.to_radians());
    let azgrad: Vec<f64> = (0..region.n_face)
        .into_par_iter()
        .map(|f| {
            let (grad_zonal, grad_meridional) = compute_one_face_gradient(f, variable, region);
            let bearing = bearings_rad[f];
            grad_meridional * bearing.sin() - grad_zonal * bearing.cos()
        })
        .collect();
    Ok(Array1::from_vec(azgrad))
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
    face_elevation: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
) -> f64 {
    let (dh_dx, dh_dy) = compute_one_face_gradient(f, face_elevation, region);
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
///
/// A tuple representing the zonal and meridional components of the gradient vector at face `f`.
///
#[inline(always)]
fn compute_one_face_gradient(
    f: usize,
    variable: ArrayView1<'_, f64>,
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

        let (f1, f2) = match extract_edge_faces(
            region.edge_face_connectivity.row(e),
            e,
            region.n_face,
            "compute_one_face_gradient",
        ) {
            Some(pair) => pair,
            None => continue,
        };

        let other = if f == f1 {
            f2
        } else if f == f2 {
            f1
        } else {
            continue;
        };
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
            lon_f.to_radians(),
            lat_f.to_radians(),
            region.face_lon[other].to_radians(),
            region.face_lat[other].to_radians(),
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
pub fn compute_slope(region: &LocalSurfaceView<'_>) -> ArrayResult {
    let n_face = region.n_face;
    let slope_vec: Vec<_> = (0..n_face)
        .into_par_iter()
        .map(|f| {
            let slope_sq = compute_face_slope_squared(f, region.face_elevation.view(), region);
            slope_sq.sqrt()
        })
        .collect();
    Ok(Array1::from_vec(slope_vec))
}

/// Collapses slopes steeper than the critical slope by diffusion, and returns the resulting face elevation change.
///
/// Each iteration gives every face that is steeper than the critical slope a diffusivity that grows with its steepness, and
/// applies one diffusion step, until no face is too steep (or for at most 1000 iterations).
///
/// A diffusion step only changes the faces on the edges of steep faces, and a face's slope depends only on its own elevation and
/// those of its edge neighbors. So after the first iteration, only the faces that changed and their neighbors have their slopes
/// re-evaluated, and each diffusion step is applied only across the edges of steep faces. This gives the same result as
/// re-evaluating and diffusing the whole region every iteration, at a cost proportional to the size of the collapsing area.
///
/// # Arguments
///
/// * `critical_slope` - Maximum allowable slope (rise over run, e.g. 0.7 for about 35 degrees).
/// * `region` - Reference to the local surface data structure containing mesh information.
///
/// # Returns
///
/// The change in elevation of each face.
pub fn slope_collapse(critical_slope: f64, region: &LocalSurfaceView<'_>) -> ArrayResult {
    const LOOP_LIMIT: usize = 1000;
    let n_face = region.n_face;
    let critical_slope_sq = critical_slope * critical_slope;
    // The diffusivity that is stable with a time step of 1
    let diffmax = compute_dt_max(Array1::<f64>::ones(n_face).view(), region);

    let base = region.face_elevation;
    let mut total_delta = Array1::<f64>::zeros(n_face);
    let mut elevation = base.to_owned();
    let mut slope_sq: Vec<f64> = (0..n_face)
        .into_par_iter()
        .with_min_len(DIFFUSION_MIN_PAR_LEN)
        .map(|f| compute_face_slope_squared(f, elevation.view(), region))
        .collect();
    let mut steep: Vec<usize> = (0..n_face)
        .filter(|&f| slope_sq[f] > critical_slope_sq)
        .collect();
    if steep.is_empty() {
        return Ok(total_delta);
    }

    let topology = DiffusionTopology::new(region);
    let mut scratch = ActiveDiffusionScratch::new(n_face, topology.edge_faces.len());
    let mut kappa = vec![0.0f64; n_face];
    let mut face_mark = vec![false; n_face];

    for _ in 0..LOOP_LIMIT {
        if steep.is_empty() {
            break;
        }
        for &f in &steep {
            kappa[f] = diffmax * (1.0 + slope_sq[f] / critical_slope_sq) * 10.0;
        }
        let touched = diffuse_active_edges(&steep, &kappa, elevation.view(), region, &topology, &mut scratch);
        for &(f, step) in &touched {
            total_delta[f] += step;
            elevation[f] = base[f] + total_delta[f];
        }
        for &f in &steep {
            kappa[f] = 0.0;
        }

        // The faces whose slope may have changed: those that changed, and their neighbors
        let mut candidates: Vec<usize> = Vec::with_capacity(touched.len() * 2);
        for &(f, _) in &touched {
            if !face_mark[f] {
                face_mark[f] = true;
                candidates.push(f);
            }
        }
        for &(f, _) in &touched {
            for &(e, is_second) in topology.incidence_of(f) {
                let other = topology.edge_faces[e][usize::from(!is_second)];
                if !face_mark[other] {
                    face_mark[other] = true;
                    candidates.push(other);
                }
            }
        }
        for &f in &candidates {
            face_mark[f] = false;
        }
        let new_slopes: Vec<f64> = candidates
            .par_iter()
            .with_min_len(DIFFUSION_MIN_PAR_LEN)
            .map(|&f| compute_face_slope_squared(f, elevation.view(), region))
            .collect();
        steep.clear();
        for (&f, &s) in candidates.iter().zip(new_slopes.iter()) {
            slope_sq[f] = s;
            if s > critical_slope_sq {
                steep.push(f);
            }
        }
    }
    Ok(total_delta)
}

/// Reusable work arrays for [`diffuse_active_edges`], kept all-zero (or unset) between calls.
struct ActiveDiffusionScratch {
    /// Per face: the change in elevation during the current step
    step: Vec<f64>,
    /// Per face: the sum of the inverse stable time steps of its active edges
    inverse_dt_sum: Vec<f64>,
    /// Per face: whether it is in the list of touched faces
    touched: Vec<bool>,
    /// Per edge: its position in the list of active edges, or `usize::MAX` if it is not active
    edge_position: Vec<usize>,
}

impl ActiveDiffusionScratch {
    fn new(n_face: usize, n_edge: usize) -> Self {
        Self {
            step: vec![0.0; n_face],
            inverse_dt_sum: vec![0.0; n_face],
            touched: vec![false; n_face],
            edge_position: vec![usize::MAX; n_edge],
        }
    }
}

/// Applies one diffusion step (as `apply_diffusion` does) in which only the faces in `active` have a nonzero diffusivity.
///
/// Only the edges of active faces carry any flux, so the step is computed on those edges and the faces they touch, summing each
/// face's fluxes in the same order as the full diffusion step.
///
/// # Returns
///
/// The faces the step changed, with their change in elevation.
fn diffuse_active_edges(
    active: &[usize],
    kappa: &[f64],
    elevation: ArrayView1<'_, f64>,
    region: &LocalSurfaceView<'_>,
    topology: &DiffusionTopology,
    scratch: &mut ActiveDiffusionScratch,
) -> Vec<(usize, f64)> {
    // The edges of the active faces, in increasing order, and the faces on them
    let mut edges: Vec<usize> = Vec::new();
    for &f in active {
        for &(e, _) in topology.incidence_of(f) {
            if scratch.edge_position[e] == usize::MAX {
                scratch.edge_position[e] = 0;
                edges.push(e);
            }
        }
    }
    edges.sort_unstable();
    let mut faces: Vec<usize> = Vec::new();
    for (i, &e) in edges.iter().enumerate() {
        scratch.edge_position[e] = i;
        for &f in &topology.edge_faces[e] {
            if !scratch.touched[f] {
                scratch.touched[f] = true;
                faces.push(f);
            }
        }
    }

    // The largest stable time step, as `compute_dt_max` computes it (only edges with a positive diffusivity contribute)
    for &e in &edges {
        let distance = region.edge_face_distance[e];
        if distance <= 0.0 {
            continue;
        }
        let [f1, f2] = topology.edge_faces[e];
        let k_avg = 0.5 * (kappa[f1] + kappa[f2]);
        let length = region.edge_length[e];
        if k_avg > 0.0 {
            scratch.inverse_dt_sum[f1] += 2.0 * k_avg * length / (distance * region.face_area[f1]);
            scratch.inverse_dt_sum[f2] += 2.0 * k_avg * length / (distance * region.face_area[f2]);
        }
    }
    let dt_max = faces
        .iter()
        .map(|&f| scratch.inverse_dt_sum[f])
        .filter(|&x| x > 0.0)
        .map(|x| 1.0 / x)
        .fold(f64::INFINITY, f64::min);
    for &f in &faces {
        scratch.inverse_dt_sum[f] = 0.0;
    }

    if dt_max.is_finite() && dt_max > 0.0 {
        let nloops = (1.0 / dt_max).ceil() as usize;
        let fac = (1.0 / nloops as f64) / 2.0;
        for _ in 0..nloops {
            let step = &scratch.step;
            let flux: Vec<f64> = edges
                .par_iter()
                .with_min_len(DIFFUSION_MIN_PAR_LEN)
                .map(|&e| {
                    let [f1, f2] = topology.edge_faces[e];
                    let h1 = elevation[f1] + step[f1];
                    let h2 = elevation[f2] + step[f2];
                    (kappa[f1] + kappa[f2]) * (h2 - h1) / region.edge_face_distance[e] * region.edge_length[e]
                })
                .collect();
            let edge_position = &scratch.edge_position;
            let dhdt: Vec<f64> = faces
                .par_iter()
                .with_min_len(DIFFUSION_MIN_PAR_LEN)
                .map(|&f| {
                    let area = region.face_area[f];
                    let mut dhdt = 0.0f64;
                    for &(e, is_second) in topology.incidence_of(f) {
                        let i = edge_position[e];
                        if i == usize::MAX {
                            continue;
                        }
                        if is_second {
                            dhdt += -flux[i] / area;
                        } else {
                            dhdt += flux[i] / area;
                        }
                    }
                    dhdt
                })
                .collect();
            for (&f, &d) in faces.iter().zip(dhdt.iter()) {
                scratch.step[f] += fac * d;
            }
        }
    }

    let touched: Vec<(usize, f64)> = faces.iter().map(|&f| (f, scratch.step[f])).collect();
    for &f in &faces {
        scratch.step[f] = 0.0;
        scratch.touched[f] = false;
    }
    for &e in &edges {
        scratch.edge_position[e] = usize::MAX;
    }
    touched
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
pub fn interpolate_node_elevation_from_faces(region: &LocalSurfaceView<'_>) -> ArrayResult {
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
                compute_one_distance(lon[o1u], lat[o1u], lon[o2u], lat[o2u], radius)
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

/// Computes the (lon, lat) coordinates of a point on a sphere given a distance and bearing.
///
/// # Arguments
/// * `lon1`, `lat` - Coordinates of the origin point in radians
/// * `distance` - Distance from origin in meters.
/// * `bearing` - Initial bearing from the origin in radians.
/// * `radius` - Radius of the sphere in meters
///
/// # Returns
///
/// The (lon, lat) coordinates in radians.
#[inline]
fn compute_one_location_from_distance_bearing(
    lon1: f64,
    lat1: f64,
    distance: f64,
    bearing: f64,
    radius: f64,
) -> (f64, f64) {
    let lat2 = (lat1.sin() * (distance / radius).cos()
        + lat1.cos() * (distance / radius).sin() * bearing.cos())
    .asin();
    let lon2 = lon1
        + (bearing.sin() * (distance / radius).sin() * lat1.cos())
            .atan2((distance / radius).cos() - lat1.sin() * lat2.sin());
    (lon2, lat2)
}

pub fn compute_location_from_distance_bearing(
    lon1: f64,
    lat1: f64,
    distances: ArrayView1<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    radius: f64,
) -> ArrayResult2D {
    let n = distances.len();
    let mut result = Array2::<f64>::zeros((n, 2));
    for i in 0..n {
        let (lon2, lat2) = compute_one_location_from_distance_bearing(
            lon1,
            lat1,
            distances[i],
            bearings[i],
            radius,
        );
        result[[i, 0]] = positive_mod(lon2 + PI, TAU) - PI;
        result[[i, 1]] = positive_mod(lat2 + PI, TAU) - PI;
    }
    Ok(result)
}

/// Computes the Haversine distance between two points on a sphere given their longitude and latitude in radians.
///
/// # Arguments
/// * `lon1`, `lat1` - Coordinates of the first point in radians.
/// * `lon2`, `lat2` - Coordinates of the second point in radians.
/// * `radius` - Radius of the sphere in meters.
///
/// # Returns
///
/// Distance in meters between the two points along the surface of the sphere.
#[inline]
fn compute_one_distance(lon1: f64, lat1: f64, lon2: f64, lat2: f64, radius: f64) -> f64 {
    // Make sure to deal with edge cases near the dateline
    let dlon = positive_mod(lon2 - lon1 + PI, TAU) - PI;
    let dlat = lat2 - lat1;
    let a = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();
    radius * c
}

/// Degrees to radians, computed exactly as NumPy's `np.radians` does (`x * (π / 180)`), so results are bit-identical.
#[inline(always)]
fn deg_to_rad(x: f64) -> f64 {
    x * (PI / 180.0)
}

/// Radians to degrees, computed exactly as NumPy's `np.degrees` does (`x * (180 / π)`), so results are bit-identical.
#[inline(always)]
fn rad_to_deg(x: f64) -> f64 {
    x * (180.0 / PI)
}

/// Minimum number of points per rayon task for the distance/bearing calculations, so small regions are not dominated by
/// fork/join overhead. This only affects scheduling, not results.
const DISTANCE_MIN_PAR_LEN: usize = 4096;

/// Computes great-circle distances and/or initial bearings from origin points to a subset of points given in degrees.
///
/// The points are `(lon[i], lat[i])` for each `i` in `indices` (or every point if `indices` is `None`). The distance to each
/// point is computed from `distance_origin`, and the bearing from `bearing_origin`; either output is skipped if its origin is
/// `None`. Each value is identical to what `compute_distances` and `compute_bearings` return for the same points after
/// converting their coordinates to radians with NumPy, and the bearings are converted back to degrees in the same way.
///
/// # Arguments
///
/// * `lon`, `lat` - Longitude and latitude of all points, in degrees.
/// * `indices` - Indices of the points to use, or `None` for all of them.
/// * `distance_origin` - `(lon, lat)` in degrees of the point to measure distances from, or `None` to skip distances.
/// * `bearing_origin` - `(lon, lat)` in degrees of the point to measure bearings from, or `None` to skip bearings.
/// * `radius` - Radius of the sphere in meters.
/// * `with_projection` - If true, also return the projected coordinates `x = distance * sin(bearing)` and
///   `y = distance * cos(bearing)` (an azimuthal equidistant projection about the origin). Requires both origins.
///
/// # Returns
///
/// * `(distances in meters, bearings in degrees in [0, 360), projected x in meters, projected y in meters)`, each `None` if it
///   was not requested.
#[allow(clippy::type_complexity)]
pub fn compute_distances_and_bearings(
    lon: ArrayView1<'_, f64>,
    lat: ArrayView1<'_, f64>,
    indices: Option<ArrayView1<'_, i64>>,
    distance_origin: Option<(f64, f64)>,
    bearing_origin: Option<(f64, f64)>,
    radius: f64,
    with_projection: bool,
) -> Result<(Option<Array1<f64>>, Option<Array1<f64>>, Option<Array1<f64>>, Option<Array1<f64>>), String> {
    if with_projection && (distance_origin.is_none() || bearing_origin.is_none()) {
        return Err("with_projection requires both a distance origin and a bearing origin".to_string());
    }
    if lon.len() != lat.len() {
        return Err(format!("lon and lat have different lengths ({} and {})", lon.len(), lat.len()));
    }
    let n_total = lon.len();
    let idx: Option<Vec<usize>> = match indices {
        Some(ind) => Some(
            ind.iter()
                .map(|&i| {
                    if i >= 0 && (i as usize) < n_total {
                        Ok(i as usize)
                    } else {
                        Err(format!("index {i} is out of range [0, {n_total})"))
                    }
                })
                .collect::<Result<_, _>>()?,
        ),
        None => None,
    };
    let n = idx.as_ref().map_or(n_total, |v| v.len());
    let point = |k: usize| -> (f64, f64) {
        let i = idx.as_ref().map_or(k, |v| v[k]);
        (deg_to_rad(lon[i]), deg_to_rad(lat[i]))
    };
    // Per-origin constants (cos and sin of the origin latitude) are computed once, and per-point terms (cos and sin of the
    // point latitude, and the longitude difference when both outputs share an origin longitude) are shared between the distance
    // and bearing calculations. The formulas are the same as compute_one_distance and compute_one_bearing, but the results can
    // differ from those functions in the last bits because the compiler is free to evaluate the shared sin/cos differently.
    let origin = |(lo, la): (f64, f64)| {
        let (lon1, lat1) = (deg_to_rad(lo), deg_to_rad(la));
        let (sin_lat1, cos_lat1) = lat1.sin_cos();
        (lon1, lat1, cos_lat1, sin_lat1)
    };
    let d_origin = distance_origin.map(origin);
    let b_origin = bearing_origin.map(origin);
    let shared_dlon = matches!((d_origin, b_origin), (Some(d), Some(b)) if d.0.to_bits() == b.0.to_bits());

    let values: Vec<[f64; 4]> = (0..n)
        .into_par_iter()
        .with_min_len(DISTANCE_MIN_PAR_LEN)
        .map(|k| {
            let (lon2, lat2) = point(k);
            let (sin_lat2, cos_lat2) = lat2.sin_cos();
            let mut dlon_d = f64::NAN;
            let distance = match d_origin {
                Some((lon1, lat1, cos_lat1, _)) => {
                    // Haversine distance, as in compute_one_distance
                    dlon_d = positive_mod(lon2 - lon1 + PI, TAU) - PI;
                    let dlat = lat2 - lat1;
                    let a = (dlat / 2.0).sin().powi(2) + cos_lat1 * cos_lat2 * (dlon_d / 2.0).sin().powi(2);
                    radius * (2.0 * a.sqrt().asin())
                }
                None => f64::NAN,
            };
            let bearing_rad = match b_origin {
                Some((lon1, _, cos_lat1, sin_lat1)) => {
                    // Initial bearing, as in compute_one_bearing, normalized to [0, 2π)
                    let dlon = if shared_dlon { dlon_d } else { positive_mod(lon2 - lon1 + PI, TAU) - PI };
                    let (sin_dlon, cos_dlon) = dlon.sin_cos();
                    let x = sin_dlon * cos_lat2;
                    let y = cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon;
                    positive_mod(x.atan2(y), TAU)
                }
                None => f64::NAN,
            };
            let (proj_x, proj_y) = if with_projection {
                let (sin_b, cos_b) = bearing_rad.sin_cos();
                (distance * sin_b, distance * cos_b)
            } else {
                (f64::NAN, f64::NAN)
            };
            [distance, rad_to_deg(bearing_rad), proj_x, proj_y]
        })
        .collect();

    let column = |j: usize| Array1::from_iter(values.iter().map(|v| v[j]));
    Ok((
        distance_origin.map(|_| column(0)),
        bearing_origin.map(|_| column(1)),
        with_projection.then(|| column(2)),
        with_projection.then(|| column(3)),
    ))
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
/// * A NumPy array of initial bearing angles (degrees), one for each (lon2, lat2) pair.
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
            // Normalize bearing to 0 to 2*PI radians
            positive_mod(initial_bearing, TAU)
        })
        .collect();
    Ok(Array1::from_vec(result_vec))
}

/// Computes the initial bearing (forward azimuth) from point 1 to point 2 on a sphere in radians.
#[inline]
fn compute_one_bearing(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let dlon = positive_mod(lon2 - lon1 + PI, TAU) - PI;
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
fn compute_dt_max(face_kappa: ArrayView1<'_, f64>, region: &LocalSurfaceView<'_>) -> f64 {
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

// This function takes arrays of x, y, z coordinates and a target radius, and returns new arrays where each point is projected radially to the specified radius from the origin.
//
// # Arguments
// * `x`, `y`, `z` - 1D arrays of the same length giving 3D positions.
// * `radius` - The target radius to project points onto.
//
// # Returns
// A tuple of three 1D NumPy arrays (x_out, y_out, z_out) where each point has been scaled to lie on the sphere of the given radius.
pub fn reset_radial_distances(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    z: ArrayView1<'_, f64>,
    r: ArrayView1<'_, f64>,
) -> Result<
    (
        numpy::ndarray::Array1<f64>,
        numpy::ndarray::Array1<f64>,
        numpy::ndarray::Array1<f64>,
    ),
    String,
> {
    let n = x.len();
    if y.len() != n || z.len() != n || r.len() != n {
        return Err("Input arrays must have the same length".to_string());
    }

    let mut x_out = Array1::<f64>::zeros(n);
    let mut y_out = Array1::<f64>::zeros(n);
    let mut z_out = Array1::<f64>::zeros(n);
    for (i, (xh, yh, zh, rh)) in izip!(&x, &y, &z, &r).enumerate() {
        let mag = (xh * xh + yh * yh + zh * zh).sqrt();
        if mag > 0.0 {
            x_out[i] = rh * xh / mag;
            y_out[i] = rh * yh / mag;
            z_out[i] = rh * zh / mag;
        } else {
            x_out[i] = 0.0;
            y_out[i] = 0.0;
            z_out[i] = 0.0;
        }
    }

    Ok((x_out, y_out, z_out))
}

thread_local! {
    // Reusable global-to-local index map for `remap_connectivity_to_local`. It is sized to the largest global index domain seen so
    // far and every entry is kept at -1 between calls, so each call costs O(number of entries touched) rather than O(mesh size).
    static REMAP_SCRATCH: std::cell::RefCell<Vec<i64>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Extracts rows of a global connectivity array and remaps their values from global to local indices.
///
/// Row `i` of the output is row `row_indices[i]` of `connectivity`, with each value `v` replaced by the position of `v` in
/// `value_indices` (if `v` appears more than once, its last position is used). Values equal to `fill_value`, values outside the
/// range `[0, total_value_size)`, and values not present in `value_indices` are set to `fill_value`.
///
/// # Arguments
///
/// * `connectivity` - The global 2D connectivity array of shape (N, K).
/// * `row_indices` - Global indices of the rows to retain.
/// * `value_indices` - Global indices of the value domain, whose positions become the local indices (e.g., the face indices).
/// * `total_value_size` - Total number of items in the global value domain (e.g., the total number of faces).
/// * `fill_value` - Fill value marking invalid entries.
///
/// # Returns
///
/// * A 2D array of shape (len(row_indices), K) of local indices, or an error if a row index is out of range.
pub fn remap_connectivity_to_local(
    connectivity: ArrayView2<'_, i64>,
    row_indices: ArrayView1<'_, i64>,
    value_indices: ArrayView1<'_, i64>,
    total_value_size: usize,
    fill_value: i64,
) -> Result<Array2<i64>, String> {
    let rows = to_row_indices(row_indices, connectivity.nrows())?;
    REMAP_SCRATCH.with(|scratch| {
        let mut map = scratch.borrow_mut();
        Ok(with_value_map(&mut map, value_indices, total_value_size, |map| {
            remap_rows(connectivity, &rows, map, total_value_size, fill_value)
        }))
    })
}

/// Converts row indices to `usize`, accepting negative indices that count from the end as in NumPy.
fn to_row_indices(row_indices: ArrayView1<'_, i64>, n_rows: usize) -> Result<Vec<usize>, String> {
    row_indices
        .iter()
        .map(|&r| {
            if r >= 0 && (r as usize) < n_rows {
                Ok(r as usize)
            } else if r < 0 && ((-r) as usize) <= n_rows {
                Ok((n_rows as i64 + r) as usize) // Negative indices count from the end, as in NumPy
            } else {
                Err(format!("row index {r} is out of bounds for connectivity with {n_rows} rows"))
            }
        })
        .collect()
}

/// Fills the global-to-local map for `value_indices`, runs `f` with it, and restores the map to all -1 afterwards.
///
/// The map must be all -1 on entry. Each value `v` in `[0, total_value_size)` maps to its position in `value_indices` (the last
/// position, if it appears more than once).
fn with_value_map<T>(
    map: &mut Vec<i64>,
    value_indices: ArrayView1<'_, i64>,
    total_value_size: usize,
    f: impl FnOnce(&[i64]) -> T,
) -> T {
    if map.len() < total_value_size {
        map.resize(total_value_size, -1);
    }
    let in_domain = |v: i64| v >= 0 && (v as usize) < total_value_size;
    for (local, &v) in value_indices.iter().enumerate() {
        if in_domain(v) {
            map[v as usize] = local as i64;
        }
    }
    let result = f(map);
    for &v in value_indices.iter() {
        if in_domain(v) {
            map[v as usize] = -1;
        }
    }
    result
}

/// Extracts `rows` of `connectivity`, mapping each value through `map` (global to local index, -1 if absent). Values equal to
/// `fill_value`, outside `[0, total_value_size)`, or absent from the map become `fill_value`.
fn remap_rows(
    connectivity: ArrayView2<'_, i64>,
    rows: &[usize],
    map: &[i64],
    total_value_size: usize,
    fill_value: i64,
) -> Array2<i64> {
    let in_domain = |v: i64| v >= 0 && (v as usize) < total_value_size;
    let ncols = connectivity.ncols();
    let mut data = vec![fill_value; rows.len() * ncols];
    if ncols > 0 {
        // Each output row depends only on its own input row, so rows are filled in parallel. The minimum chunk length keeps
        // small regions on a single thread.
        data.par_chunks_mut(ncols)
            .zip(rows.par_iter())
            .with_min_len(REMAP_MIN_PAR_ROWS)
            .for_each(|(out_row, &r)| {
                for (o, &v) in out_row.iter_mut().zip(connectivity.row(r).iter()) {
                    if v != fill_value && in_domain(v) {
                        let local = map[v as usize];
                        if local >= 0 {
                            *o = local;
                        }
                    }
                }
            });
    }
    Array2::from_shape_vec((rows.len(), ncols), data).expect("the shape matches the data length")
}

/// Minimum number of rows per rayon task in `remap_rows`.
const REMAP_MIN_PAR_ROWS: usize = 4096;

/// The global connectivity arrays of a surface mesh.
pub struct GlobalConnectivity<'a> {
    pub face_edge: ArrayView2<'a, i64>,
    pub face_node: ArrayView2<'a, i64>,
    pub face_face: ArrayView2<'a, i64>,
    pub node_face: ArrayView2<'a, i64>,
    pub edge_face: ArrayView2<'a, i64>,
    pub edge_node: ArrayView2<'a, i64>,
    pub n_face: usize,
    pub n_node: usize,
    pub n_edge: usize,
    pub fill_value: i64,
}

/// The six connectivity arrays of a local region, in local indices.
pub struct LocalConnectivity {
    pub face_edge: Array2<i64>,
    pub face_node: Array2<i64>,
    pub face_face: Array2<i64>,
    pub node_face: Array2<i64>,
    pub edge_face: Array2<i64>,
    pub edge_node: Array2<i64>,
}

/// Builds all six local connectivity arrays of a region in one pass.
///
/// Each output is identical to calling `remap_connectivity_to_local` with the corresponding rows and value domain (for example,
/// `face_edge` uses the face indices as rows and the edge indices as values). Doing them together means the global-to-local map
/// of each value domain (faces, nodes, edges) is built once rather than once per array.
pub fn build_local_connectivity(
    face_indices: ArrayView1<'_, i64>,
    node_indices: ArrayView1<'_, i64>,
    edge_indices: ArrayView1<'_, i64>,
    global: &GlobalConnectivity<'_>,
) -> Result<LocalConnectivity, String> {
    let fill = global.fill_value;
    let face_rows = to_row_indices(face_indices, global.n_face)?;
    let node_rows = to_row_indices(node_indices, global.n_node)?;
    let edge_rows = to_row_indices(edge_indices, global.n_edge)?;
    for (name, conn, n) in [
        ("face_edge", &global.face_edge, global.n_face),
        ("face_node", &global.face_node, global.n_face),
        ("face_face", &global.face_face, global.n_face),
        ("node_face", &global.node_face, global.n_node),
        ("edge_face", &global.edge_face, global.n_edge),
        ("edge_node", &global.edge_node, global.n_edge),
    ] {
        if conn.nrows() != n {
            return Err(format!("{name} connectivity has {} rows, expected {n}", conn.nrows()));
        }
    }

    REMAP_SCRATCH.with(|scratch| {
        let mut map = scratch.borrow_mut();
        // Values that are faces
        let (face_face, node_face, edge_face) = with_value_map(&mut map, face_indices, global.n_face, |map| {
            (
                remap_rows(global.face_face, &face_rows, map, global.n_face, fill),
                remap_rows(global.node_face, &node_rows, map, global.n_face, fill),
                remap_rows(global.edge_face, &edge_rows, map, global.n_face, fill),
            )
        });
        // Values that are nodes
        let (face_node, edge_node) = with_value_map(&mut map, node_indices, global.n_node, |map| {
            (
                remap_rows(global.face_node, &face_rows, map, global.n_node, fill),
                remap_rows(global.edge_node, &edge_rows, map, global.n_node, fill),
            )
        });
        // Values that are edges
        let face_edge = with_value_map(&mut map, edge_indices, global.n_edge, |map| {
            remap_rows(global.face_edge, &face_rows, map, global.n_edge, fill)
        });
        Ok(LocalConnectivity {
            face_edge,
            face_node,
            face_face,
            node_face,
            edge_face,
            edge_node,
        })
    })
}

/// A [`UniqueCollector`] that has collected at least `1 / DENSE_COLLECTION_DIVISOR` of its index domain produces its sorted
/// output by scanning the marks rather than sorting. Around this point the two cost about the same.
const DENSE_COLLECTION_DIVISOR: usize = 32;

thread_local! {
    // Reusable membership marks for `build_region_indices`. Sized to the largest index domain seen so far and kept all-false
    // between uses, so each use costs O(number of entries touched) rather than O(mesh size).
    static UNIQUE_MARKS: std::cell::RefCell<Vec<bool>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Accumulates the sorted, unique values of a set of global indices, equivalent to `np.unique` with the fill value removed.
///
/// Membership is tracked in a reusable mark array (all-false on entry, restored to all-false by `finish`), so the cost is
/// proportional to the number of values added, not to the size of the mesh.
struct UniqueCollector<'m> {
    marks: &'m mut Vec<bool>,
    out: Vec<i64>,
    n_total: usize,
    fill_value: i64,
    error: Option<String>,
    what: &'static str,
}

impl<'m> UniqueCollector<'m> {
    fn new(marks: &'m mut Vec<bool>, n_total: usize, fill_value: i64, what: &'static str) -> Self {
        if marks.len() < n_total {
            marks.resize(n_total, false);
        }
        Self {
            marks,
            out: Vec::new(),
            n_total,
            fill_value,
            error: None,
            what,
        }
    }

    #[inline]
    fn add(&mut self, v: i64) {
        if v == self.fill_value {
            return;
        }
        if v < 0 || v as usize >= self.n_total {
            if self.error.is_none() {
                self.error = Some(format!("{} index {v} is out of range [0, {})", self.what, self.n_total));
            }
            return;
        }
        let m = &mut self.marks[v as usize];
        if !*m {
            *m = true;
            self.out.push(v);
        }
    }

    /// Adds every value in the given rows of a row-major connectivity array with `ncols` columns.
    fn add_rows(&mut self, conn: &[i64], ncols: usize, rows: &[i64]) {
        for &r in rows {
            let start = r as usize * ncols;
            for &v in &conn[start..start + ncols] {
                self.add(v);
            }
        }
    }

    /// Returns the unique values in the order they were first added, and clears the marks.
    fn finish_unsorted(mut self) -> Result<Vec<i64>, String> {
        let out = std::mem::take(&mut self.out);
        for &v in &out {
            self.marks[v as usize] = false;
        }
        match self.error.take() {
            Some(e) => Err(e),
            None => Ok(out),
        }
    }

    fn finish(mut self) -> Result<Vec<i64>, String> {
        let mut out = std::mem::take(&mut self.out);
        if out.len() >= self.n_total / DENSE_COLLECTION_DIVISOR {
            // A large fraction of the domain was collected (as for regions that cover much of the globe), so reading the marks
            // in index order is cheaper than sorting. It gives the same sorted, unique list, and clears the marks as it goes.
            out.clear();
            for (i, m) in self.marks[..self.n_total].iter_mut().enumerate() {
                if *m {
                    *m = false;
                    out.push(i as i64);
                }
            }
        } else {
            for &v in &out {
                self.marks[v as usize] = false;
            }
            out.sort_unstable();
        }
        if let Some(e) = self.error.take() {
            return Err(e);
        }
        Ok(out)
    }
}

/// Determines the global face, node, and edge indices of a local region from its seed faces.
///
/// This is the region construction used by `Surface.extract_region` and `LocalSurface.extract_subregion`:
/// 1. `edges` are the edges of the seed faces.
/// 2. `faces` are the seed faces plus a buffer of every face that shares a node with them. (Faces that share an edge with a
///    seed face also share that edge's nodes, so this includes every face that shares a face, node, or edge with them.)
/// 3. `nodes` are the nodes of all of `faces`, including the buffer.
///
/// Each list is sorted and unique, with fill values removed, exactly as `np.unique` would produce.
///
/// # Returns
///
/// * `(faces, nodes, edges)`, or an error if any index is out of range.
pub fn build_region_indices(
    seed_faces: ArrayView1<'_, i64>,
    global: &GlobalConnectivity<'_>,
) -> Result<(Vec<i64>, Vec<i64>, Vec<i64>), String> {
    let fill = global.fill_value;
    let seeds: Vec<i64> = seed_faces.iter().copied().collect();
    if let Some(&bad) = seeds.iter().find(|&&f| f < 0 || f as usize >= global.n_face) {
        return Err(format!("seed face index {bad} is out of range [0, {})", global.n_face));
    }
    let check_rows = |name: &str, conn: &ArrayView2<'_, i64>, n: usize| {
        if conn.nrows() == n {
            Ok(())
        } else {
            Err(format!("{name} connectivity has {} rows, expected {n}", conn.nrows()))
        }
    };
    check_rows("face_edge", &global.face_edge, global.n_face)?;
    check_rows("face_node", &global.face_node, global.n_face)?;
    check_rows("node_face", &global.node_face, global.n_node)?;

    // Row-major views (no copy when the arrays are already C-contiguous, as they are coming from uxarray)
    let face_edge = global.face_edge.as_standard_layout();
    let face_node = global.face_node.as_standard_layout();
    let node_face = global.node_face.as_standard_layout();
    let fe = face_edge.as_slice().expect("standard layout");
    let fnn = face_node.as_slice().expect("standard layout");
    let nf = node_face.as_slice().expect("standard layout");

    UNIQUE_MARKS.with(|scratch| {
        let mut marks = scratch.borrow_mut();

        let mut c = UniqueCollector::new(&mut marks, global.n_edge, fill, "edge");
        c.add_rows(fe, face_edge.ncols(), &seeds);
        let edges = c.finish()?;

        // The seed nodes are only used to find the faces around them, so they do not need to be sorted
        let mut c = UniqueCollector::new(&mut marks, global.n_node, fill, "node");
        c.add_rows(fnn, face_node.ncols(), &seeds);
        let seed_nodes = c.finish_unsorted()?;

        let mut c = UniqueCollector::new(&mut marks, global.n_face, fill, "face");
        for &f in &seeds {
            c.add(f);
        }
        c.add_rows(nf, node_face.ncols(), &seed_nodes);
        let faces = c.finish()?;

        let mut c = UniqueCollector::new(&mut marks, global.n_node, fill, "node");
        c.add_rows(fnn, face_node.ncols(), &faces);
        let nodes = c.finish()?;
        Ok((faces, nodes, edges))
    })
}

/// Computes the elevation of a mean reference plane fitted to the points within a reference radius of a region's center.
///
/// The plane `z = a*x + b*y + c` is fitted by least squares to the points whose radial distance is at most `reference_radius`,
/// using coordinates projected onto the local tangent plane. Points inside the radius are replaced by the plane's elevation;
/// points outside it keep their original elevation. If fewer than 5 points lie within the radius, the original elevations are
/// returned unchanged.
///
/// The fit is done on coordinates centered on the mean of the fitted points, which keeps the 2x2 normal equations well
/// conditioned even though the projected coordinates can be many orders of magnitude larger than the elevations.
///
/// # Arguments
///
/// * `x` - 1D array of projected x coordinates (meters).
/// * `y` - 1D array of projected y coordinates (meters).
/// * `elevation` - 1D array of elevations (meters).
/// * `distance` - 1D array of radial distances from the region center (meters).
/// * `reference_radius` - Radius of the region used to fit the reference plane (meters).
///
/// # Returns
///
/// A 1D array of reference elevations, one per input point.
///
/// # Errors
///
/// Returns an error if the input arrays have mismatched lengths.
pub fn reference_surface(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    elevation: ArrayView1<'_, f64>,
    distance: ArrayView1<'_, f64>,
    reference_radius: f64,
) -> ArrayResult {
    let n = elevation.len();
    if x.len() != n || y.len() != n || distance.len() != n {
        return Err("x, y, elevation, and distance must all have the same length".to_string());
    }

    // Means of the fitted points
    let (mut count, mut xm, mut ym, mut zm) = (0usize, 0.0, 0.0, 0.0);
    for i in 0..n {
        if distance[i] <= reference_radius {
            count += 1;
            xm += x[i];
            ym += y[i];
            zm += elevation[i];
        }
    }
    if count < 5 {
        return Ok(elevation.to_owned());
    }
    let inv_count = 1.0 / count as f64;
    xm *= inv_count;
    ym *= inv_count;
    zm *= inv_count;

    // Centered second moments
    let (mut sxx, mut sxy, mut syy, mut sxz, mut syz) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        if distance[i] <= reference_radius {
            let dx = x[i] - xm;
            let dy = y[i] - ym;
            let dz = elevation[i] - zm;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
            sxz += dx * dz;
            syz += dy * dz;
        }
    }

    let (a, b) = solve_plane_slopes(sxx, sxy, syy, sxz, syz);
    let c = zm - a * xm - b * ym;

    let result = Array1::from_shape_fn(n, |i| {
        if distance[i] <= reference_radius {
            a * x[i] + b * y[i] + c
        } else {
            elevation[i]
        }
    });
    Ok(result)
}

/// Solves the centered 2x2 normal equations `[[sxx, sxy], [sxy, syy]] [a, b] = [sxz, syz]` for the plane slopes.
///
/// When the system is singular (all fitted points collinear or coincident), returns the minimum-norm solution, matching the
/// behavior of a least-squares solver.
fn solve_plane_slopes(sxx: f64, sxy: f64, syy: f64, sxz: f64, syz: f64) -> (f64, f64) {
    let trace = sxx + syy;
    if trace <= 0.0 {
        return (0.0, 0.0);
    }
    let det = sxx * syy - sxy * sxy;
    if det > 1e-12 * trace * trace {
        return ((sxz * syy - syz * sxy) / det, (syz * sxx - sxz * sxy) / det);
    }
    // Rank-1: the matrix is trace * v v^T for the unit vector v along the points' line, and the minimum-norm solution is
    // (v . rhs) / trace * v.
    let (vx, vy) = if sxx >= syy { (sxx, sxy) } else { (sxy, syy) };
    let norm = (vx * vx + vy * vy).sqrt();
    let (vx, vy) = (vx / norm, vy / norm);
    let s = (vx * sxz + vy * syz) / trace;
    (s * vx, s * vy)
}

#[cfg(test)]
mod reference_surface_tests {
    use super::*;

    #[test]
    fn recovers_exact_plane_and_leaves_outside_points() {
        let n = 50;
        let x = Array1::from_shape_fn(n, |i| 1.0e5 * ((i as f64) * 0.37).cos());
        let y = Array1::from_shape_fn(n, |i| 1.0e5 * ((i as f64) * 0.61).sin());
        let d = Array1::from_shape_fn(n, |i| (x[i] * x[i] + y[i] * y[i]).sqrt());
        let z = Array1::from_shape_fn(n, |i| 1.0e-3 * x[i] - 2.0e-3 * y[i] + 1234.5 + if i % 2 == 0 { 0.0 } else { 7.0 });
        let r = 8.0e4;
        let out = reference_surface(x.view(), y.view(), z.view(), d.view(), r).unwrap();
        for i in 0..n {
            if d[i] > r {
                assert_eq!(out[i], z[i]);
            }
        }
        // Fit should pass through the mean offset of the alternating +0/+7 noise, with the exact slopes
        let inside: Vec<usize> = (0..n).filter(|&i| d[i] <= r).collect();
        let resid_mean: f64 = inside.iter().map(|&i| z[i] - out[i]).sum::<f64>() / inside.len() as f64;
        assert!(resid_mean.abs() < 1e-9);
    }

    #[test]
    fn too_few_points_returns_input() {
        let x = Array1::from(vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        let d = x.clone();
        let z = Array1::from(vec![1.0, 2.0, 0.0, 5.0, 3.0, 9.0]);
        let out = reference_surface(x.view(), x.view(), z.view(), d.view(), 3.0).unwrap();
        assert_eq!(out, z);
    }

    #[test]
    fn collinear_points_give_min_norm_slopes() {
        // Points on the line y = x: only the slope along the line is determined
        let x = Array1::from(vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        let z = Array1::from_shape_fn(6, |i| 2.0 * x[i]);
        let d = Array1::zeros(6);
        let out = reference_surface(x.view(), x.view(), z.view(), d.view(), 1.0).unwrap();
        for i in 0..6 {
            assert!((out[i] - z[i]).abs() < 1e-12);
        }
    }
}
