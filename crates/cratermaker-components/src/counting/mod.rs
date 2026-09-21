use crate::ArrayResult;
use crate::crater::Crater;
use crate::surface::LocalSurfaceView;
use numpy::ndarray::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::f64::consts::{FRAC_PI_2, TAU};
const _RIM_RADIUS_MAX: f64 = 1.20;
const _RIM_RADIUS_MIN: f64 = 0.95;
const _FLOOR_RADIUS: f64 = 0.2;

pub fn mask_crater_faces(
    region: &LocalSurfaceView<'_>,
    crater: &Crater,
) -> Result<Array1<bool>, &'static str> {
    let ids = region.crater_id.as_ref().ok_or("crater_id required")?;

    // `ids` is shaped (n_face, n_layer). A face is part of the crater if ANY layer id matches.
    let mask_vec: Vec<bool> = ids
        .axis_iter(Axis(0))
        .map(|row| row.iter().any(|&id| id == crater.id))
        .collect();

    Ok(Array1::from(mask_vec))
}

fn filter_crater_faces(
    arr: &Array1<f64>,
    mask: &Array1<bool>,
) -> Result<Array1<f64>, &'static str> {
    if arr.len() != mask.len() {
        return Err("filter_crater_faces: arr and mask must have the same length");
    }

    let vals: Vec<f64> = arr
        .iter()
        .zip(mask.iter())
        .filter_map(|(&e, &m)| if m && !e.is_nan() { Some(e) } else { None })
        .collect();

    Ok(Array1::from(vals))
}

pub fn measure_floor_elevation(
    region: &LocalSurfaceView<'_>,
    crater: &Crater,
) -> Result<f64, &'static str> {
    let mut elev_sum: f64 = 0.0;
    let mut n_floor: usize = 0;
    let mask_crater = mask_crater_faces(region, crater)?;
    let elevation = &region
        .desloped_face_elevation
        .unwrap_or(region.face_elevation)
        .to_owned();
    let elevation = filter_crater_faces(&elevation, &mask_crater)?;
    let n_face = elevation.len();
    let floor_distance = _FLOOR_RADIUS * crater.measured_radius;
    let face_distance = region.face_distance.ok_or("face_distance required")?;
    let face_distance = filter_crater_faces(&face_distance.to_owned(), &mask_crater)?;

    let mask_floor: Vec<bool> = (0..n_face)
        .into_par_iter()
        .map(|f| face_distance[f] <= floor_distance)
        .collect();

    for (&elev, &mask) in elevation.iter().zip(mask_floor.iter()) {
        if mask && !elev.is_nan() {
            elev_sum += elev;
            n_floor += 1;
        }
    }
    if n_floor == 0 {
        return Err("no valid faces inside floor boundary");
    }
    Ok(elev_sum / (n_floor as f64))
}

pub fn measure_rim_height(
    region: &LocalSurfaceView<'_>,
    crater: &Crater,
) -> Result<f64, &'static str> {
    let mut elev_sum: f64 = 0.0;
    let mut n_rim: usize = 0;
    let mask_crater = mask_crater_faces(region, crater)?;
    let elevation = &region
        .desloped_face_elevation
        .unwrap_or(region.face_elevation)
        .to_owned();
    let elevation = filter_crater_faces(&elevation, &mask_crater)?;
    let n_face = elevation.len();
    let rim_distance_min = _RIM_RADIUS_MIN * crater.measured_radius;
    let rim_distance_max = _RIM_RADIUS_MAX * crater.measured_radius;
    let face_distance = region.face_distance.ok_or("face_distance required")?;
    let face_distance = filter_crater_faces(&face_distance.to_owned(), &mask_crater)?;

    let mask_rim: Vec<bool> = (0..n_face)
        .into_par_iter()
        .map(|f| face_distance[f] >= rim_distance_min && face_distance[f] <= rim_distance_max)
        .collect();

    for (&elev, &mask) in elevation.iter().zip(mask_rim.iter()) {
        if mask && !elev.is_nan() {
            elev_sum += elev;
            n_rim += 1;
        }
    }
    if n_rim == 0 {
        return Err("no valid faces inside rim boundary");
    }
    Ok(elev_sum / (n_rim as f64))
}


/// Computes the signed radial distance from points to an ellipse.
///
/// For each input point `(x[i], y[i])`, this function:
///
/// 1. Computes the polar coordinates `(r, θ)` of the point relative to the
///    ellipse center `(x0, y0)`.
/// 2. Computes the radius of the ellipse along the same direction, using the
///    semi-major axis `a`, semi-minor axis `b`, and orientation.
/// 3. Returns the signed difference `r - r_ellipse(θ)`.
///
/// By this convention:
///
/// * Values **> 0** indicate points outside the ellipse.
/// * Values **< 0** indicate points inside the ellipse.
/// * Values **≈ 0** lie close to the ellipse boundary.
///
/// The orientation is given in radians.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `a` - Semi-major axis of the ellipse.
/// * `b` - Semi-minor axis of the ellipse.
/// * `orientation` - Orientation of the ellipse in radians (see above).
/// * `x0` - x-coordinate of the ellipse center.
/// * `y0` - y-coordinate of the ellipse center.
///
/// # Returns
///
/// On success, returns an array of signed radial distances, one value
/// per input point.
///
/// # Errors
///
/// Currently this function always returns `Ok(...)`. The `ArrayResult`
/// return type allows future extensions to return descriptive errors
/// (for example, if the input arrays have inconsistent lengths).
pub fn radial_distance_to_crater_rim(
    x: &ArrayView1<'_, f64>,
    y: &ArrayView1<'_, f64>,
    crater: &Crater,
    x0: f64,
    y0: f64,
) -> ArrayResult {
    let phi = TAU - FRAC_PI_2 - crater.measured_orientation.to_radians();
    let a = crater.measured_semimajor_axis;
    let b = crater.measured_semiminor_axis;

    let result_vec: Vec<f64> = (0..x.len())
        .into_par_iter()
        .map(|i| {
            let xi = x[i];
            let yi = y[i];
            let dx = xi - x0;
            let dy = yi - y0;
            let r = (dx * dx + dy * dy).sqrt();
            let theta = dy.atan2(dx);
            let alpha = theta - phi;
            let ca = alpha.cos();
            let sa = alpha.sin();
            r - (a * b) / ((b * ca).powi(2) + (a * sa).powi(2)).sqrt()
        })
        .collect();

    Ok(Array1::from(result_vec))
}

