use std::f64::{
    self,
    consts::{PI, SQRT_2},
};

use itertools::Itertools;
use ndarray::ArrayView1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};
use rand::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{VSMALL, EJPROFILE};

const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.5;


/// Computes the ejecta profile scaling at a given radial distance.
///
/// This function returns a value based on a power-law decay that modifies the
/// ejecta distribution intensity outside the crater rim. Values inside the crater return zero.
///
/// # Arguments
///
/// * `r_actual` - Radial distance from the crater center (in meters).
/// * `crater_radius` - Radius of the crater (in meters).
/// * `ejrim` - Rim elevation parameter used to scale the profile.
///
/// # Returns
///
/// * Scaled profile value representing the ejecta contribution at distance `r_actual`.
#[inline]
pub fn profile_function(r_actual: f64, crater_radius: f64, ejrim: f64) -> f64 {
    if r_actual >= crater_radius {
        let r = r_actual / crater_radius;
        ejrim * r.powf(-EJPROFILE)
    } else {
        0.0
    }
}


/// Computes only the radial ejecta profile without ray modulation.
///
/// This is a simple power-law decay of ejecta intensity with radial distance.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distance` - 1D array of radial distances from crater center.
/// * `crater_diameter` - Diameter of the crater (meters).
/// * `ejrim` - Profile scaling factor.
///
/// # Returns
///
/// * A NumPy array of ejecta profile values.
#[pyfunction]
pub fn profile<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    Ok(PyArray1::from_vec(
        py,
        radial_distance
            .as_array()
            .iter()
            .map(|&r| profile_function(r, crater_diameter / 2.0, ejrim))
            .collect(),
    ))
}

/// Computes the ray intensity contribution for a single ejecta point.
///
/// # Arguments
///
/// * `r` - Radial distance from the crater center (in meters).
/// * `theta0` - Initial bearing angle (radians).
/// * `crater_radius` - Radius of the crater (in meters).
/// * `rmin` - Minimum normalized radial distance.
/// * `rmax` - Maximum normalized radial distance.
/// * `thetari` - Precomputed ray azimuth angles.
/// * `random_numbers` - Set of random numbers used for pattern variation.
///
/// # Returns
///
/// * Computed ray intensity value for this point.
fn ray_intensity_point(
    r: f64,
    theta0: f64,
    crater_radius: f64,
    rmin: f64,
    rmax: f64,
    thetari: &[f64],
    random_numbers: &[f64],
) -> f64 {
    if r < crater_radius {
        return 0.0;
    }
    (0..NPATT as usize)
        .into_par_iter()
        .map(|j| {
            let rn = random_numbers[j];
            let theta = (theta0 + rn * 2.0 * PI) % (2.0 * PI);
            let r_pattern = r / crater_radius - rn;
            FRAYREDUCTION.powi(j as i32 - 1)
                * ray_intensity_func(r_pattern, theta, rmin, rmax, thetari, 2.348 * crater_radius.powf(0.006))
        })
        .sum()
}


/// Computes the ray-modulated ejecta intensity values for a set of input points.
///
/// Each intensity value is the sum of contributions from multiple rays, modulated by random
/// angular perturbations and decaying with distance.
///
/// # Arguments
///
/// * `radial_distance` - 1D array of radial distances (meters).
/// * `initial_bearing` - 1D array of initial bearing angles (radians).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A vector of normalized ray-modulated intensity values.
pub fn ray_intensity_internal<'py>(
    radial_distance: ArrayView1<'py, f64>,
    initial_bearing: ArrayView1<'py, f64>,
    crater_diameter: f64,
) -> PyResult<Vec<f64>> {
    if radial_distance.len() != initial_bearing.len() {
        return Err(PyValueError::new_err(
            "initial_bearing and radial_distance arrays must be the same size.",
        ));
    }
    let crater_radius = crater_diameter / 2.0;
    let rmax = 100.0;
    let rmin = 1.0;

    let mut rng = rand::rng();

    // Distribute ray patterns evenly around the crater
    let mut thetari = (0..NRAYMAX)
        .map(|i| PI * 2.0 * (i + 1) as f64 / NRAYMAX as f64)
        .collect_vec();
    thetari.shuffle(&mut rng);

    let random_numbers: Vec<f64> = rng.random_iter().take(NPATT as usize).collect_vec();

    let mut intensity: Vec<f64> = (0..radial_distance.len())
        .into_par_iter()
        .map(|i| {
            ray_intensity_point(
                radial_distance[i],
                initial_bearing[i],
                crater_radius,
                rmin,
                rmax,
                &thetari,
                &random_numbers,
            )
        })
        .collect();
    let max_val = intensity
        .iter()
        .zip(radial_distance.iter())
        .filter_map(|(&val, &r)| if r >= crater_radius { Some(val) } else { None })
        .fold(f64::MIN, |a, b| a.max(b));
    intensity = intensity
        .into_iter()
        .zip(radial_distance)
        .map(|(intensity, &radial_distance)| {
            if radial_distance >= crater_radius {
                intensity / max_val
            } else {
                intensity
            }
        })
        .collect_vec();
    Ok(intensity)
}

/// Python wrapper for computing ray-modulated ejecta intensity field.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distance` - 1D array of radial distances from crater center.
/// * `initial_bearing` - 1D array of bearing angles (radians).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A NumPy array of ray-modulated intensity values.
#[pyfunction]
pub fn ray_intensity<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let intensity = ray_intensity_internal(
        radial_distance.as_array(),
        initial_bearing.as_array(),
        crater_diameter,
    )?;
    Ok(intensity.into_pyarray(py))
}

/// Computes the intensity contribution of all rays for a single radial/angular location.
///
/// Applies a distance-dependent decay of ray count, and for each ray evaluates its contribution
/// based on length and angular width using a modified Gaussian model.
///
/// # Arguments
///
/// * `r` - Normalized radial distance (unitless).
/// * `theta` - Angle from crater center (radians).
/// * `rmin` - Minimum normalized radius.
/// * `rmax` - Maximum normalized radius (truncation).
/// * `thetari` - Ray azimuth directions.
/// * `minray` - Minimum ejecta blanket radius (scaling for ray length).
///
/// # Returns
///
/// * Total ray intensity contribution at the point.
fn ray_intensity_func(
    r: f64,
    theta: f64,
    rmin: f64,
    rmax: f64,
    thetari: &[f64],
    minray: f64,
) -> f64 {
    const RAYP: f64 = 4.0;
    if !r.is_finite() || r <= 0.0 || r > rmax {
        return 0.0;
    } else if r < 1.0 {
        1.0
    } else {
        let tmp = (NRAYMAX as f64).powf(RAYP)
            - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln();
        let n = if tmp < 0.0 {
            NRAYMAX // "Nrays" in Minton et al. (2019)
        } else {
            // Exponential decay of ray number with distance
            ((((NRAYMAX as f64).powf(RAYP)
                - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln())
            .powf(1.0 / RAYP))
            .floor() as i32)
                .clamp(1, NRAYMAX)
        };
        (0..NRAYMAX)
            .into_par_iter()
            .map(|i| {
                let length = minray
                    * ((rmax / minray).ln()
                        * ((((NRAYMAX - i + 2) as f64).powf(RAYP) - 1.0)
                            / ((NRAYMAX as f64).powf(RAYP) - 1.0)))
                        .exp();
                if r > length {
                    0.0 // Skip this ray contribution
                } else {
                    let w = (rmax / length).powf(1.0);
                    let rw = PI / (w * NRAYMAX as f64)
                        * (rmin / r)
                        * (1.0 - (1.0 - w / rmin) * (1.0 - (r / rmin).powi(2)).exp());
                    ejecta_ray_func(theta, thetari[i as usize], r, n, rw)
                }
            })
            .sum()
    }
}

/// Computes the contribution of a single ray to the ejecta field using a Gaussian profile.
///
/// The ray is centered at `thetar` with angular width `w`, evaluated at angle `theta`.
/// The function ensures normalization across rays and returns zero for values below machine epsilon.
///
/// # Arguments
///
/// * `theta` - Sample angle (radians).
/// * `thetar` - Ray center angle (radians).
/// * `r` - Radial distance (meters).
/// * `n` - Number of rays contributing at this distance.
/// * `w` - Angular width of the ray.
///
/// # Returns
///
/// * Ray intensity value at the specified angle and distance.
fn ejecta_ray_func(theta: f64, thetar: f64, r: f64, n: i32, w: f64) -> f64 {
    let c = w / r;
    let b = thetar;
    let dtheta = f64::min(2.0 * PI - (theta - b).abs(), (theta - b).abs());
    let logval = -dtheta.powi(2) / (2.0 * c.powi(2));
    if logval < VSMALL.ln() {
        0.0
    } else {
        let a = (2.0 * PI).sqrt() / (n as f64 * c * (PI / (2.0 * SQRT_2 * c)).erf());
        a * logval.exp()
    }
}



/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject)]
pub struct Crater {
    pub final_diameter: f64,
}

/// Morphological parameters for generating and modifying lunar surface craters.
///
/// Includes floor geometry, rim height, and whether to apply ray modulation.
#[derive(FromPyObject)]
pub struct SimpleMoonMorphology {
    pub floor_depth: f64,
    pub floor_diameter: f64,
    pub rim_height: f64,
    pub ejrim: f64,
    pub crater: Crater,
}

/// View into a region of the surface mesh, consisting of node and face indices.
///
/// Used to localize crater effects to a subset of the full mesh.
#[derive(FromPyObject)]
pub struct SurfaceView<'py> {
    pub node_indices: PyReadonlyArray1<'py, i64>,
    pub face_indices: PyReadonlyArray1<'py, i64>,
}

