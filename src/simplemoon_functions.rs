use itertools::Itertools;
use ndarray::ArrayView1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};
use rand::prelude::*;
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use std::f64::{
    self,
    consts::{PI, SQRT_2},
};

use crate::{EJPROFILE, RIMDROP, VSMALL};

const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.90;

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

/// Calculates the elevation of a crater as a function of distance from the center.
///
/// This function applies a polynomial profile for the crater interior (r < 1.0) and a rim dropoff
/// function for the exterior (r ≥ 1.0). It is based on the crater profile model described in
/// Fassett and Thomson (2014).
///
/// Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and
/// the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271.
/// https://doi.org/10.1002/2014JE004698
///
/// This function is split off from the `profile` function for clarity, and is not intended to be
/// called directly.
///
/// # Arguments
///
/// * `r` - Normalized radial distance (unitless, where 1.0 corresponds to the crater rim).
/// * `elevation` - Baseline elevation before crater modification.
/// * `c0`, `c1`, `c2`, `c3` - Polynomial coefficients for the crater profile interior.
/// * `rim_height` - Height of the crater rim.
/// * `ejrim` - Rim dropoff parameter.
///
/// # Returns
///
/// * Adjusted elevation according to crater profile at distance `r`.
#[inline]
fn crater_profile_function(
    r: f64,
    elevation: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    rim_height: f64,
    ejrim: f64,
) -> f64 {
    if r >= 1.0 {
        elevation + (rim_height - ejrim) * r.powf(-RIMDROP)
    } else {
        elevation + c0 + c1 * r + c2 * r.powi(2) + c3 * r.powi(3)
    }
}

/// Computes a crater profile elevation array from input radial distances and reference elevations.
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `diameter` - Total diameter of the crater (in meters).
/// * `floor_depth` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_height` - Height of the crater rim above mean surface level (in meters).
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
#[pyfunction]
pub fn crater_profile<'py>(
    py: Python<'py>,
    r_array: PyReadonlyArray1<'py, f64>,
    reference_elevation_array: PyReadonlyArray1<'py, f64>,
    diameter: f64,
    floor_depth: f64,
    floor_diameter: f64,
    rim_height: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances = r_array.as_array();
    let reference_elevations = reference_elevation_array.as_array();
    if radial_distances.len() != reference_elevations.len() {
        return Err(PyValueError::new_err(
            "input arrays must have the same length",
        ));
    }
    const A: f64 = 4.0 / 11.0;
    const B: f64 = -32.0 / 187.0;
    
    // Calculate the floor radius relative to the final crater radius
    let flrad = floor_diameter / diameter;
    let radius = diameter / 2.0;

    // Use polynomial crater profile similar to that of Fassett and Thomson (2014), but the parameters are set by the crater dimensions
    let c1 = (-floor_depth - rim_height)
        / (flrad - 1.0 + A * (flrad.powi(2) - 1.0) + B * (flrad.powi(3) - 1.0));
    let c0 = rim_height - c1 * (1.0 + A + B);
    let c2 = A * c1;
    let c3 = B * c1;

    let ninc = radial_distances.iter().filter(|&&x| x <= radius).count();
    let meanref = if ninc == 0 {
        *radial_distances
            .iter()
            .zip(reference_elevations)
            .min_by(|(&radius_a, _), (&radius_b, _)| radius_a.partial_cmp(&radius_b).unwrap())
            .unwrap()
            .1
    } else {
        radial_distances
            .iter()
            .zip(reference_elevations)
            .filter(|(&r, _)| r <= radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let min_elevation = meanref - floor_depth;

    Ok(PyArray1::from_iter(
        py,
        reference_elevations
            .iter()
            .zip(radial_distances)
            .map(|(&elevation, &radial_distance)| {
                let r = radial_distance / radius;
                (
                    crater_profile_function(r, elevation, c0, c1, c2, c3, rim_height, ejrim),
                    radial_distance,
                )
            })
            .map(|(elevation, radial_distance)| {
                if radial_distance <= radius {
                    elevation.max(min_elevation)
                } else {
                    elevation
                }
            }),
    ))
}

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
pub fn ejecta_profile_function(r_actual: f64, crater_radius: f64, ejrim: f64) -> f64 {
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
pub fn ejecta_profile<'py>(
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
            .map(|&r| ejecta_profile_function(r, crater_diameter / 2.0, ejrim))
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
    let mut sum = 0.0;
    for j in 0..NPATT as usize {
        let rn = random_numbers[j];
        let theta = (theta0 + rn * 2.0 * PI) % (2.0 * PI);
        let r_pattern = r / crater_radius - rn;
        sum += FRAYREDUCTION.powi(j as i32 - 1)
            * ray_intensity_func(
                r_pattern,
                theta,
                rmin,
                rmax,
                thetari,
                2.348 * crater_radius.powf(0.006),
            );
    }
    sum
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
    seed: u64,
) -> PyResult<Vec<f64>> {
    if radial_distance.len() != initial_bearing.len() {
        return Err(PyValueError::new_err(
            "initial_bearing and radial_distance arrays must be the same size.",
        ));
    }
    let crater_radius = crater_diameter / 2.0;
    let rmax = 100.0;
    let rmin = 1.0;

    let mut rng = ChaCha12Rng::seed_from_u64(seed);

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
    seed: u64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let intensity = ray_intensity_internal(
        radial_distance.as_array(),
        initial_bearing.as_array(),
        crater_diameter,
        seed,
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
    const RAYP: f64 = 4.0; // Power law exponent for ray number decay
    const RAY_WIDTH_EXPONENT: f64 = 1.25; // Exponent for angular width decay
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
                        * (1.0 - (1.0 - w / rmin) * (1.0 - (r / rmin).powf(RAY_WIDTH_EXPONENT)).exp());
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
