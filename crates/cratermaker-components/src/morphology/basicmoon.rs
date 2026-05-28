use crate::ArrayResult;
use itertools::Itertools;
use libm::erf;
use numpy::ndarray::prelude::*;
use pyo3::FromPyObject;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::f64::{
    self,
    consts::{PI, SQRT_2, TAU},
};

const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.90;

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject, Clone, Debug)]
pub struct BasicMoonCrater {
    pub id: u32,
    pub diameter: f64,
    pub radius: f64,
    pub semimajor_axis: f64,
    pub semiminor_axis: f64,
    pub orientation: f64,
    pub transient_diameter: f64,
    pub projectile_diameter: f64,
    pub projectile_velocity: f64,
    pub projectile_angle: f64,
    pub projectile_density: f64,
    pub location: (f64, f64),
    pub morphology_type: String,
    pub measured_semimajor_axis: f64,
    pub measured_semiminor_axis: f64,
    pub measured_orientation: f64,
    pub measured_diameter: f64,
    pub measured_radius: f64,
    pub measured_location: (f64, f64),
    pub time: Option<f64>,
    pub floor_elevation: f64,
    pub floor_radius: f64,
    pub wall_curvature: f64,
    pub rim_width: f64,
    pub rim_elevation: f64,
    pub rimdrop: f64,
    pub rim_flank_radius: f64,
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_offset: f64,
}

// Computes a either a crater and/or ejecta 1D profile array from input radial distances and reference elevations using the model of
// Minton et al. (2026)
//
//
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_radius` - Radius of the crater's rim (in meters).
/// * `floor_elevation` - Elevation of the crater floor relative to the reference surface (in meters and should be negative).
/// * `floor_radius` - Radius of the flat portion of the crater floor (in meters).
/// * `wall_curvature` - Parameter controlling the curvature of the crater wall (>1 for more curvature)
/// * `rim_width` - Width of the crater rim (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters).
/// * `rimdrop` - Exponent for the rim dropoff function (typically -4 to -6)
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
/// * `eprofile` - Exponent for the ejecta dropoff function (typically -3)
/// * `peak_height` - Height of the central peak above the crater floor (in meters).
/// * `peak_width` - Width of the central peak (in meters).
/// * `peak_offset` - Radial offset of the central peak from the crater center (in meters). "concentric").
/// * `include_crater` - Whether to include the crater profile in the output (true/false).
/// * `include_ejecta` - Whether to include the ejecta profile in the output (true/false).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
pub fn basicmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &BasicMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> ArrayResult {
    assert_eq!(radial_distances.len(), reference_elevations.len());

    // Compute the weighted elevation profile relative to the reference plane
    let ninc = radial_distances
        .iter()
        .filter(|&&x| x <= crater.radius)
        .count();
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
            .filter(|(&r, _)| r <= crater.radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let min_elevation = meanref + crater.floor_elevation;

    Ok(Array1::from_iter(
        reference_elevations
            .iter()
            .zip(radial_distances.iter().copied())
            .map(|(href, r)| {
                let mut hcrat = crater_profile_function(r, crater.radius, crater.floor_elevation, crater.floor_radius, crater.wall_curvature, crater.rim_width, crater.rim_elevation, crater.rimdrop, crater.peak_height, crater.peak_width, crater.peak_offset);
                let mut hej = ejecta_profile_function(r, crater.radius, crater.ejrim, crater.ejprofile);
                if r < crater.radius && r > crater.floor_radius {
                    hej += (hcrat - (crater.rim_elevation - crater.ejrim)).max(0.0);
                }

                if include_crater {
                    if r > crater.radius || hcrat > 0.0 {
                        hcrat = (hcrat - hej).max(0.0);
                    } 
                } else {
                    hcrat = 0.0;
                }
                if !include_ejecta {
                    hej = 0.0;
                }

                let h = href + hcrat + hej;

                if r <= crater.radius {
                    h.max(min_elevation)
                } else {
                    h
                }
            }),
    ))
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
/// * `r` - Radial distance from crater center
/// * `radius` - Radius of the crater rim (in meters).
/// * `hf` - Elevation of the crater floor relative to the reference plane.
/// * `rf` - Radius of the crater floor.
/// * `beta` - Wall curvature parameter (1.0 means straight, > 1.0 means curved)
/// * `rw` - Width of the crater rim
/// * `hr` - Height of the crater rim above the reference plane.
/// * `prd` - Exponent for the rim dropoff function.
/// * `hc` - Height of the central peak above the floor.
/// * `rc` - Radius of the central peak.
/// * `ro` - Radial offset of the central peak from the crater center .
///
/// # Returns
///
/// * Elevation values at distance `r`.
#[inline]
pub fn crater_profile_function(r: f64, radius: f64, hf: f64, rf: f64, beta: f64, rw: f64, hr: f64, prd: f64, hc: f64, rc: f64, ro: f64) -> f64 {
    let rw_half = rw / 2.0;
    let fc = hc * (-((r - ro) / rc).powi(2)).exp(); // Central peak contribution. Compute this separately to avoid sharp discontinuities
    if r <= rf {
        fc + hf
    } else {
        let fe = hr * (r / radius).powf(prd);
        if r >= radius + rw_half {
            fe
        } else {
            let r0 = (r - rf) / (radius - rf);
            let c = (hr - hf) * ((-beta / 2.0).exp() + 1.0) / (beta.exp() - 1.0);
            let t = (r - (radius - rw_half)) / rw;
            let phi = 6.0 * t.powi(5) - 15.0 * t.powi(4) + 10.0 * t.powi(3);
            let fw = c * ((beta * r0).exp() - beta.exp()) / (1.0 + (beta * (r0 - 0.5)).exp()) + hr;
            if r <= radius - rw_half {
                fw + fc
            } else {
                (1.0 - phi) * fw + phi * fe + fc
            }
        }
    }
}

/// Computes the ejecta profile scaling at a given radial distance.
///
/// This function returns a value based on a power-law decay that modifies the
/// ejecta distribution intensity outside the crater rim. Values inside the crater return zero.
///
/// # Arguments
///
/// * `r` - Radial distance from the crater center (in meters).
/// * `radius` - Radius of the crater (in meters).
/// * `ejrim` - Rim elevation parameter used to scale the profile.
/// *  `ejprofile` - Exponent for the power-law decay of the ejecta profile.
///
/// # Returns
///
/// * Scaled profile value representing the ejecta contribution at distance `r_actual`.
#[inline]
pub fn ejecta_profile_function(r: f64, radius: f64, ejrim: f64, ejprofile: f64) -> f64 {
    if r >= radius {
        let rej = r / radius;
        ejrim * rej.powf(ejprofile)
    } else {
        0.0
    }
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
/// * `thetari` - Precomputed ray azimuth angles (in radians).
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
        let theta = (theta0 + rn * TAU) % TAU;
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
/// * `radial_distances` - 1D array of radial distances (meters).
/// * `initial_bearing` - 1D array of initial bearing angles (degrees).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A vector of normalized ray-modulated intensity values.
pub fn ray_intensity(
    radial_distances: ArrayView1<'_, f64>,
    initial_bearings: ArrayView1<'_, f64>,
    crater_diameter: f64,
    seed: u64,
) -> ArrayResult {
    if radial_distances.len() != initial_bearings.len() {
        return Err("radial_distances and reference_elevations must have same length".into());
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

    let mut intensity: Vec<f64> = (0..radial_distances.len())
        .into_par_iter()
        .map(|i| {
            ray_intensity_point(
                radial_distances[i],
                initial_bearings[i],
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
        .zip(radial_distances.iter())
        .filter_map(|(&val, &r)| if r >= crater_radius { Some(val) } else { None })
        .fold(f64::MIN, |a, b| a.max(b));
    intensity = intensity
        .into_iter()
        .zip(radial_distances.iter())
        .map(|(intensity, &r)| {
            if r >= crater_radius {
                intensity / max_val
            } else {
                intensity
            }
        })
        .collect_vec();
    Ok(Array1::from_vec(intensity))
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
                        * (1.0
                            - (1.0 - w / rmin) * (1.0 - (r / rmin).powf(RAY_WIDTH_EXPONENT)).exp());
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
    if logval < crate::VSMALL.ln() {
        0.0
    } else {
        let denom = erf(PI / (2.0 * SQRT_2 * c));
        if denom.abs() <= crate::VSMALL {
            return 0.0;
        }
        let a = (2.0 * PI).sqrt() / (n as f64 * c * denom);
        a * logval.exp()
    }
}
