use crate::ArrayResult;
use itertools::Itertools;
use libm::{erf, exp};
use numpy::ndarray::prelude::*;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::f64::{
    self,
    consts::{PI, SQRT_2, TAU},
};

const RIMDROP: f64 = -6.00; // The exponent for the uplifted rim dropoff.
const EJPROFILE: f64 = -3.0; // The exponent for the ejecta profile
const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.90;

pub fn crater_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater_diameter: f64,
    floor_elevation: f64,
    floor_diameter: f64,
    rim_elevation: f64,
    ejrim: f64,
    fassett_yang_fraction: f64,
    morphology_subtype: &str,
) -> ArrayResult {
    let p1 = fassett2020_profile(
        radial_distances,
        reference_elevations,
        crater_diameter,
        floor_elevation,
        floor_diameter,
        rim_elevation,
        ejrim,
    )?;
    if fassett_yang_fraction < 1.0 {
        let p2 = yang2021_profile(
            radial_distances,
            reference_elevations,
            crater_diameter,
            floor_elevation,
            floor_diameter,
            rim_elevation,
            ejrim,
            morphology_subtype,
        )?;
        Ok(p1 * fassett_yang_fraction + p2 * (1.0 - fassett_yang_fraction))
    } else {
        Ok(p1)
    }
}

/// Computes a crater profile elevation array from input radial distances and reference elevations using modification of the model
/// given in Fassett and Thomson (2014) and Fassett et al. (2020).
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
///
/// # Arguments
///
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_diameter` - Total diameter of the crater (in meters).
/// * `floor_elevation` - Depth of the crater floor below mean surface level (in meters, relative to datum).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters, relative to datum).
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
///
/// # References
///
/// Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and the rate of erosion on the
///     Moon. J. Geophys. Res. 119, 2014JE004698-2271. https://doi.org/10.1002/2014JE004698
///
/// Fassett, C.I., Beyer, R.A., Deutsch, A.N., Hirabayashi, M., Leight, C., Mahanti, P., Nypaver, C.A., Thomson, B.J.,
///     Minton, D.A., 2022. Topographic Diffusion Revisited: Small Crater Lifetime on the Moon and Implications for Volatile
///     Exploration. Journal of Geophysical Research: Planets 127, e2022JE007510. https://doi.org/10.1029/2022JE007510

pub fn fassett2020_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater_diameter: f64,
    floor_elevation: f64,
    floor_diameter: f64,
    rim_elevation: f64,
    ejrim: f64,
) -> ArrayResult {
    assert_eq!(radial_distances.len(), reference_elevations.len());
    const A: f64 = 4.0 / 11.0;
    const B: f64 = -32.0 / 187.0;

    // Calculate the floor radius relative to the final crater radius
    let flrad = floor_diameter / crater_diameter;
    let crater_radius = crater_diameter / 2.0;

    // Use polynomial crater profile similar to that of Fassett and Thomson (2014), but the parameters are set by the crater dimensions
    let c1 = (floor_elevation - rim_elevation)
        / (flrad - 1.0 + A * (flrad.powi(2) - 1.0) + B * (flrad.powi(3) - 1.0));
    let c0 = rim_elevation - c1 * (1.0 + A + B);
    let c2 = A * c1;
    let c3 = B * c1;

    let ninc = radial_distances
        .iter()
        .filter(|&&x| x <= crater_radius)
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
            .filter(|(&r, _)| r <= crater_radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let min_elevation = meanref + floor_elevation;

    Ok(Array1::from_iter(
        reference_elevations
            .iter()
            .zip(radial_distances)
            .map(|(&elevation, &radial_distances)| {
                let r = radial_distances / crater_radius;
                (
                    fassett2020_profile_function(
                        r,
                        elevation,
                        c0,
                        c1,
                        c2,
                        c3,
                        rim_elevation,
                        ejrim,
                    ),
                    radial_distances,
                )
            })
            .map(|(elevation, radial_distances)| {
                if radial_distances <= crater_radius {
                    elevation.max(min_elevation)
                } else {
                    elevation
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
/// * `r` - Normalized radial distance (unitless, where 1.0 corresponds to the crater rim).
/// * `elevation` - Baseline elevation before crater modification.
/// * `c0`, `c1`, `c2`, `c3` - Polynomial coefficients for the crater profile interior.
/// * `rim_elevation` - Height of the crater rim.
/// * `ejrim` - Rim dropoff parameter.
///
/// # Returns
///
/// * Adjusted elevation according to crater profile at distance `r`.
#[inline]
fn fassett2020_profile_function(
    r: f64,
    elevation: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    rim_elevation: f64,
    ejrim: f64,
) -> f64 {
    if r >= 1.0 {
        elevation + (rim_elevation - ejrim) * r.powf(RIMDROP)
    } else {
        elevation + c0 + c1 * r + c2 * r.powi(2) + c3 * r.powi(3)
    }
}

/// Computes a crater profile elevation array from input radial distances and reference elevations using modification of the model
/// given in Yang et al. (2021) for small craters.
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
///
/// # Arguments
///
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_diameter` - Total diameter of the crater (in meters).
/// * `floor_elevation` - Depth of the crater floor below mean surface level (in meters, relative to datum).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters, relative to datum).
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
///
/// # References
///
/// Yang, X., Fa, W., Du, J., Xie, M., Liu, T., 2021. Effect of Topographic Degradation on Small Lunar Craters: Implications for
///     Regolith Thickness Estimation. Geophysical Research Letters 48, e2021GL095537. https://doi.org/10.1029/2021GL095537
pub fn yang2021_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater_diameter: f64,
    floor_elevation: f64,
    floor_diameter: f64,
    rim_elevation: f64,
    ejrim: f64,
    morphology_subtype: &str,
) -> ArrayResult {
    assert_eq!(radial_distances.len(), reference_elevations.len());
    let hr = rim_elevation / crater_diameter;
    let d0 = -floor_elevation / crater_diameter + hr;
    let alpha = -3.1906;
    let crater_radius = crater_diameter * 0.5;
    let rb = floor_diameter / crater_diameter;
    let he = ejrim / crater_diameter;

    Ok(Array1::from_iter(
        reference_elevations
            .iter()
            .zip(radial_distances)
            .map(|(&elevation, &radial_distances)| {
                let r = radial_distances / crater_radius;
                let h = elevation / crater_diameter;
                if morphology_subtype == "normal" {
                    (yang2021_normal_profile(r, alpha, d0, hr, he) + hr) * crater_diameter
                } else if morphology_subtype == "central mound" {
                    let rm = 0.293 * crater_diameter.powf(-0.086);
                    let hm = 0.23e-3 * crater_diameter.powf(0.64);
                    (yang2021_centralmound_profile(r, alpha, d0, hr, he, rb, rm, hm) + hr + h)
                        * crater_diameter
                } else if morphology_subtype == "flat-bottomed" {
                    (yang2021_flatbottom_profile(r, alpha, d0, hr, he, rb) + hr + h)
                        * crater_diameter
                } else if morphology_subtype == "concentric" {
                    let c3 = 0.0155 * crater_diameter.powf(0.343);
                    let ri = 0.383 * crater_diameter.powf(0.053);
                    let ro = 0.421 * crater_diameter.powf(0.102);
                    (yang2021_concentric_profile(r, alpha, d0, hr, he, ri, ro, c3) + hr + h)
                        * crater_diameter
                } else {
                    panic!("Unknown morphology_subtype: {}", morphology_subtype);
                }
            }),
    ))
}

fn yang2021_normal_profile(r: f64, alpha: f64, d0: f64, hr: f64, he: f64) -> f64 {
    if r >= 1.0 {
        hr * (r.powf(alpha) - 1.0) - ejecta_profile_function(r, 1.0, he)
    } else {
        let a = -2.8567;
        let b = 5.8270;
        let c = d0 * (exp(a) + 1.0) / (exp(b) - 1.0);

        c * (exp(b * r) - exp(b)) / (1.0 + exp(a + b * r))
    }
}

fn yang2021_centralmound_profile(
    r: f64,
    alpha: f64,
    d0: f64,
    hr: f64,
    he: f64,
    rb: f64,
    rm: f64,
    hm: f64,
) -> f64 {
    if r <= rm {
        (1.0 - r / rm) * hm - d0
    } else if r <= rb {
        -d0
    } else if r <= 1.0 {
        let a = -2.6921;
        let b = 6.1678;
        let c = d0 * (exp(a) + 1.0) / (exp(b) - 1.0);
        let r0 = (r - rb) / (1.0 - rb);

        c * (exp(b * r0) - exp(b)) / (1.0 + exp(a + b * r0))
    } else {
        hr * (r.powf(alpha) - 1.0) - ejecta_profile_function(r, 1.0, he)
    }
}

fn yang2021_flatbottom_profile(r: f64, alpha: f64, d0: f64, hr: f64, he: f64, rb: f64) -> f64 {
    if r <= rb {
        -d0
    } else if r <= 1.0 {
        let a = -2.6003;
        let b = 5.8783;
        let c = d0 * (exp(a) + 1.0) / (exp(b) - 1.0);
        let r0 = (r - rb) / (1.0 - rb);

        c * (exp(b * r0) - exp(b)) / (1.0 + exp(a + b * r0))
    } else {
        hr * (r.powf(alpha) - 1.0) - ejecta_profile_function(r, 1.0, he)
    }
}

fn yang2021_concentric_profile(
    r: f64,
    alpha: f64,
    d0: f64,
    hr: f64,
    he: f64,
    ri: f64,
    ro: f64,
    c3: f64,
) -> f64 {
    let c1 = 0.192;
    let c2 = 0.01;
    let r0 = (r - ro) / (1.0 - ro);
    let a = -1.6536;
    let b = 4.7626;
    let f0 = c1 * r.powi(2) + c2 * r - d0;
    let h1 = c1 * ri.powi(2) + c2 * ri - d0;
    let f1 = c3 * (r - ri) + h1;
    let h2 = c3 * (ro - ri) + h1;
    let c = -h2 * (exp(a) + 1.0) / (exp(b) - 1.0);
    let f2 = c * (exp(b * r0) - exp(b)) / (1.0 + exp(a + b * r0));

    if r <= ri {
        f0
    } else if r <= ro {
        f1
    } else if r <= 1.0 {
        f2
    } else {
        hr * (r.powf(alpha) - 1.0) - ejecta_profile_function(r, 1.0, he)
    }
}

/// Computes only the radial ejecta profile without ray modulation.
///
/// This is a simple power-law decay of ejecta intensity with radial distance.
///
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from crater center.
/// * `crater_diameter` - Diameter of the crater (meters).
/// * `ejrim` - Profile scaling factor.
///
/// # Returns
///
/// * A NumPy array of ejecta profile values.
pub fn ejecta_profile(
    radial_distances: ArrayView1<'_, f64>,
    crater_diameter: f64,
    ejrim: f64,
) -> ArrayResult {
    Ok(Array1::from_vec(
        radial_distances
            .iter()
            .map(|&r| ejecta_profile_function(r, crater_diameter / 2.0, ejrim))
            .collect(),
    ))
}

/// Computes the ejecta profile scaling at a given radial distance.
///
/// This function returns a value based on a power-law decay that modifies the
/// ejecta distribution intensity outside the crater rim. Values inside the crater return zero.
///
/// # Arguments
///
/// * `r` - Radial distance from the crater center (in meters).
/// * `crater_radius` - Radius of the crater (in meters).
/// * `ejrim` - Rim elevation parameter used to scale the profile.
///
/// # Returns
///
/// * Scaled profile value representing the ejecta contribution at distance `r_actual`.
#[inline]
pub fn ejecta_profile_function(r: f64, crater_radius: f64, ejrim: f64) -> f64 {
    if r >= crater_radius {
        let rej = r / crater_radius;
        ejrim * rej.powf(EJPROFILE)
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
