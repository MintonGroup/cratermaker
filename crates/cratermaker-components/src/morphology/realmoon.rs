use crate::{ArrayResult,ArrayResult2D};
use crate::morphology::basicmoon::{crater_profile_function, ejecta_profile_function};
use std::collections::HashMap;
use std::f64::consts::TAU;
use rand::prelude::*;
use rand::SeedableRng;
use rand_distr::{Normal,Uniform};
use rand_chacha::ChaCha12Rng;
use numpy::ndarray::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
pub struct RealMoonCrater<'a> {
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
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_ring_radius: f64,
    pub peak_center_distance: f64,
    pub peak_center_bearing: f64,
    pub elevation_offset: f64,
    pub rim_radius_rng_seed: u64,
    pub rim_elevation_rng_seed: u64,
    pub floor_radius_rng_seed: u64,
    pub wall_texture_rng_seed: u64,
    pub ejecta_texture_rng_seed: u64,
    pub floor_texture_rng_seed: u64,
    pub rim_radius_psd: ArrayView2<'a, f64>,
    pub floor_radius_psd: ArrayView2<'a, f64>,
}

// Creates a profile of the crater
//
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `bearings` - 1D array of bearings corresponding to each radius.
/// * `reference_elevations` - 1D array of reference elevations at each radius (in meters).
/// * `crater` - A `RealMoonCrater` struct containing the crater parameters and PSDs.
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
pub fn realmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &RealMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> ArrayResult {
    assert_eq!(radial_distances.len(), reference_elevations.len(), "Input arrays must have the same length");
    assert_eq!(radial_distances.len(), bearings.len(), "Input arrays must have the same length");
    let n_points = radial_distances.len();
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
    let rim_elevation = crater.rim_elevation - crater.elevation_offset;
    let floor_elevation = crater.floor_elevation - crater.elevation_offset;
    let min_elevation = meanref + crater.floor_elevation;

    let (rim_radius_profile, floor_radius_profile) =
        crossbeam::thread::scope(|s| {
            let h1 = s.spawn(|_| {
                profile_from_psd(
                    crater.radius,
                    crater.radius,
                    crater.rim_radius_psd,
                    bearings,
                    None,
                    crater.rim_radius_rng_seed,
                )
            });

            let h2 = s.spawn(|_| {
                if include_crater {
                    profile_from_psd(
                        crater.radius,
                        crater.floor_radius,
                        crater.floor_radius_psd,
                        bearings,
                        None,
                        crater.floor_radius_rng_seed,
                    )
                } else {
                    Ok(Array1::<f64>::from_elem(bearings.len(), crater.floor_radius))
                }
            });

            (h1.join().unwrap(), h2.join().unwrap())
        })
        .map_err(|_| "crossbeam scope panicked".to_string())?;

    let rim_radius_profile = rim_radius_profile?;
    let floor_radius_profile = floor_radius_profile?;

    let out: Vec<f64> = (0..n_points)
        .into_par_iter() 
        .map(|i| {
            let r = radial_distances[i];
            let href = reference_elevations[i];
            let rim_r = rim_radius_profile[i];
            let floor_r = floor_radius_profile[i];
            let rim_elev = rim_elevation * ((rim_r - floor_r) / (crater.radius - crater.floor_radius));
            let mut hcrat = crater_profile_function(
                r,
                rim_r,                 // per-angle rim radius
                floor_elevation,
                floor_r,               // per-angle floor radius
                crater.wall_curvature,
                crater.rim_width,
                rim_elev,              // per-angle rim elevation
                crater.rimdrop,
                crater.peak_height,
                crater.peak_width,
                crater.peak_ring_radius,
            );

            let mut hej = ejecta_profile_function(r, rim_r, crater.ejrim, crater.ejprofile);

            if r < rim_r && r > floor_r {
                hej += (hcrat - (rim_elev - crater.ejrim)).max(0.0);
            }

            if include_crater {
                if r > rim_r || hcrat > 0.0 {
                    hcrat = (hcrat - hej).max(0.0);
                }
                hcrat += crater.elevation_offset;
            } else {
                hcrat = 0.0;
            }

            if !include_ejecta {
                hej = 0.0;
            }

            let h = href + hcrat + hej; 
            if r <= rim_r { h.max(min_elevation) } else { h }
        }).collect();

    Ok(Array1::from_vec(out))
}
///
///
/// Computes a target 1D power spectral density distribution based on control points defining a piecewise linear function in log-log space, with optional Gaussian noise added to the log power values.
///
/// # Arguments
/// * `control_points` - A dictionary containing the control points for the piecewise linear function. The expected keys are:
/// - "s12": The slope of the first segment (y2-y1)/(x2-x1) 
/// - "y1": The y-coordinate of the first breakpoint in ln(power).
/// - "x2": The x-coordinate of the second breakpoint in ln(wavelength). 
/// - "y2": The y-coordinate of the second breakpoint in ln(power).
/// - "y3": The y-coordinate at the third breakpoint (2nd highest wavelength) in log(power).
/// - "y4": The y-coordinate of the highest wavelength in ln(power).
///  * `npoints` - The number of points in the output PSD, which determines the wavelength resolution and the Nyquist frequency.
///  * `add_noise` - Whether to add Gaussian noise to the ln(power) values to simulate natural variability in the PSD.
///  * `seed` - The random seed for reproducibility of the noise if `add_noise` is true.
///
pub fn get_1d_psd_from_control_points(
    control_points: &HashMap<String, f64>,
    npoints: usize,
    add_noise: bool,
    rng_seed: u64,
) -> ArrayResult2D {
    let s12 = control_points["s12"];
    let x2 = control_points["x2"];
    let y2 = control_points["y2"];
    let y3 = control_points["y3"];
    let y4 = control_points["y4"];

    let x4 = TAU.ln();
    let x3 = (TAU * 0.5).ln();

    // Same spacing logic as Python: interval = exp(bp4_)x / npoints
    let interval = (x4).exp() / npoints as f64;

    // Equivalent to rfft sizing in Python:
    // dfft.size = npoints/2 + 1, iend = dfft.size - 1
    let iend = npoints / 2;
    let nrows = iend.saturating_sub(1); // wavelength from bins 1..iend-1

    let mut psd = Array2::<f64>::zeros((nrows, 2));

    // wavelength[k-1] = 1 / freq[k], freq[k] = k / (npoints * interval)
    // => wavelength = (npoints * interval) / k
    let base = npoints as f64 * interval;
    for k in 1..iend {
        let row = k - 1;
        psd[[row, 0]] = base / k as f64;
    }

    // Find x2_index (matching Python loop behavior)
    let threshold = (x2).exp();
    let mut x2_index = nrows.saturating_sub(1);
    for i in 0..nrows {
        if psd[[i, 0]] < threshold {
            x2_index = i.saturating_sub(1);
            break;
        }
    }

    // Piecewise lines in log-log space
    let s23 = (y3 - y2) / (x3 - x2);
    let b23 = y3 - s23 * x3;
    let b12 = y2 - s12 * x2;

    if nrows > 0 {
        psd[[0, 1]] = y4.exp(); 
    }
    if nrows > 1 {
        psd[[1, 1]] = y3.exp(); 
    }

    // psd[2 : x2_index + 1, 1]
    if x2_index >= 2 {
        for i in 2..=x2_index {
            let log_w = psd[[i, 0]].ln();
            psd[[i, 1]] = (s23 * log_w + b23).exp();
        }
    }

    // psd[x2_index + 1 :, 1]
    for i in (x2_index + 1)..nrows {
        let log_w = psd[[i, 0]].ln();
        psd[[i, 1]] = (s12 * log_w + b12).exp();
    }

    // flipud
    let mut flipped = Array2::<f64>::zeros((nrows, 2));
    for i in 0..nrows {
        flipped.row_mut(i).assign(&psd.row(nrows - 1 - i));
    }

    // Optional Gaussian noise in log power
    if add_noise {
        let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
        let normal = Normal::new(0.0, 0.55).expect("valid normal distribution");
        for i in 0..nrows {
            let log_power = flipped[[i, 1]].ln();
            let noisy_log_power = log_power + normal.sample(&mut rng);
            flipped[[i, 1]] = (noisy_log_power).exp();
        }
    }

    Ok(flipped)
}

///
/// Generates a surface profile based on a 1D power spectral density (PSD) and optional phase information, simulating a crater surface with specified roughness characteristics.
///
/// # Arguments
/// * `crater_radius` - The radius of the crater (in meters), which scales the amplitude
/// * `ymean` - The mean elevation of the surface (in meters), which serves as a baseline for the profile.
/// * `psd` - A 2D array where the first column contains wavelengths and the second column contains power values, defining the roughness characteristics of the surface.
/// * `theta` - A 1D array of angular positions (in radians) at which to compute the profile, typically ranging from 0 to 2π.
///    - A 1D array of polar angle (in radians) at which to compute the profile, typically ranging from 0 to 2π.
/// * `phases` - An optional 1D array of phase values (in radians) corresponding to each frequency in the PSD. If not provided, random phases will be generated.
/// * `rng_seed` - The random seed for reproducibility when generating random phases if `phases` is not provided.
/// # Returns
/// * A 1D array of values corresponding to the input angles, representing the linear profile generated from the PSD and phase information.
///
pub fn profile_from_psd(
    crater_radius: f64,
    ymean: f64,
    psd: ArrayView2<'_, f64>,
    theta: ArrayView1<'_, f64>,
    phases: Option<ArrayView1<'_, f64>>,
    rng_seed: u64,
) -> ArrayResult {
    let nfreq = psd.nrows();
    let ntheta = theta.len();

    let phase_values: Array1<f64> = if let Some(p) = phases {
        p.to_owned()
    } else {
        let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
        let uniform = Uniform::new(0.0, TAU).expect("valid uniform distribution");
        Array1::from_iter((0..nfreq).map(|_| uniform.sample(&mut rng)))
    };

    let amplitude: Array1<f64> = psd.column(1).mapv(|p| (p * TAU).sqrt());

    let out: Vec<f64> = (0..ntheta)
        .into_par_iter()
        .map(|j| {
            let t = theta[j];
            let mut dy = 0.0;
            for i in 0..nfreq {
                let freq = 1.0 / psd[[i, 0]];
                dy += amplitude[i] * (TAU * freq * t + phase_values[i]).cos();
            }
            dy * crater_radius + ymean 
        })
        .collect();

    Ok(Array1::from_vec(out))
}