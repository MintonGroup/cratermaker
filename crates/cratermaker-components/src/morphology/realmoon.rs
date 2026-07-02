use crate::morphology::basicmoon::{crater_profile_function, ejecta_profile_function};
use crate::{ArrayResult, ArrayResult2D};
use interp::{InterpMode, interp};
use numpy::ndarray::prelude::*;
use rand::SeedableRng;
use rand::prelude::*;
use rand_chacha::ChaCha12Rng;
use rand_distr::{Normal, Uniform};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rustfft::{FftPlanner, num_complex::Complex};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
pub struct RealMoonCrater<'a> {
    pub diameter: f64,
    pub radius: f64,
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
    assert_eq!(
        radial_distances.len(),
        reference_elevations.len(),
        "Input arrays must have the same length"
    );
    assert_eq!(
        radial_distances.len(),
        bearings.len(),
        "Input arrays must have the same length"
    );
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
            .min_by(|&(&radius_a, _), &(&radius_b, _)| radius_a.partial_cmp(&radius_b).unwrap())
            .unwrap()
            .1
    } else {
        radial_distances
            .iter()
            .zip(reference_elevations)
            .filter(|&(&r, _)| r <= crater.radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let rim_elevation = crater.rim_elevation - crater.elevation_offset;
    let floor_elevation = crater.floor_elevation - crater.elevation_offset;
    let min_elevation = meanref + crater.floor_elevation;

    // Create profile functions that will be interpolated later
    let rim_radius_profile =
        compute_profile_from_psd(crater.radius, crater.radius, crater.rim_radius_psd);
    let floor_radius_profile =
        compute_profile_from_psd(crater.radius, crater.floor_radius, crater.floor_radius_psd);

    let mut rimtheta: Vec<f64> = Vec::new();
    let mut floortheta: Vec<f64> = Vec::new();
    let nrim = rim_radius_profile.len();
    let nfloor = floor_radius_profile.len();
    for i in 0..nrim {
        rimtheta.push(TAU * ((i - 1) as f64 / (nrim - 2) as f64));
    }
    for i in 0..nfloor {
        floortheta.push(TAU * ((i - 1) as f64 / (nfloor - 2) as f64));
    }

    let out: Vec<f64> = (0..n_points)
        .into_par_iter()
        .map(|i| {
            let r = radial_distances[i];
            let href = reference_elevations[i];
            let theta = bearings[i];
            let rim_r = interp(
                &rimtheta,
                &rim_radius_profile,
                theta,
                &InterpMode::default(),
            );
            let floor_r = interp(
                &floortheta,
                &floor_radius_profile,
                theta,
                &InterpMode::default(),
            );
            let rim_elev =
                rim_elevation * (rim_r / crater.radius) * (crater.floor_radius / floor_r);
            let mut hcrat = crater_profile_function(
                r,
                rim_r, // per-angle rim radius
                floor_elevation,
                floor_r, // per-angle floor radius
                crater.wall_curvature,
                crater.rim_width,
                rim_elev, // per-angle rim elevation
                crater.rimdrop,
                crater.peak_height,
                crater.peak_width,
                crater.peak_ring_radius,
            );

            let mut hej = ejecta_profile_function(r, rim_r, crater.ejrim, crater.ejprofile);

            if r < rim_r && r > floor_r {
                hej += hcrat - rim_elev + crater.ejrim;
                hej = hej.clamp(0.0, crater.ejrim);
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
        })
        .collect();

    Ok(Array1::from_vec(out))
}
///
///
/// Computes a target 1D power spectral density distribution based on control points defining a piecewise linear function in log-log space, with optional Gaussian noise added to the log power values.
///
/// # Arguments
/// * `control_points` - A dictionary containing the control points for the piecewise linear function. The expected keys are:
/// - "sn": The slope of the first segment (y2-y1)/(x2-x1)
/// - "yn": The y-coordinate of the first breakpoint in ln(power).
/// - "y1": The x-coordinate of the second breakpoint in ln(wavelength).
/// - "y2": The y-coordinate of the second breakpoint in ln(power).
/// - "y3": The y-coordinate at the third breakpoint (2nd highest wavelength) in log(power).
/// - "y4": The y-coordinate of the highest wavelength in ln(power).
/// - "y5": The y-coordinate of the highest wavelength in ln(power).
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
    let sn = control_points["sn"];
    let yn = control_points["yn"];
    let y1 = control_points["y1"];
    let y2 = control_points["y2"];
    let y3 = control_points["y3"];
    let y4 = control_points["y4"];
    let y5 = control_points["y5"];

    let x1 = TAU.ln();

    let interval = x1.exp() / npoints as f64;

    // Equivalent to rfft sizing in Python:
    let iend = npoints / 2;
    let nrows = iend.saturating_sub(1); // wavelength from bins 1..iend-1
    let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
    let mut psd = Array2::<f64>::zeros((nrows, 3));
    let base = npoints as f64 * interval;
    let uniform = Uniform::new(0.0, TAU).expect("valid uniform distribution");
    for k in 1..iend {
        let row = k - 1;
        psd[[row, 0]] = base / k as f64; // Wavelength
        psd[[row, 2]] = uniform.sample(&mut rng); // Phase (randomized)
    }

    psd[[0, 1]] = y1.exp();
    psd[[1, 1]] = y2.exp();
    psd[[2, 1]] = y3.exp();
    psd[[3, 1]] = y4.exp();
    psd[[4, 1]] = y5.exp();

    let xn = psd[[5, 0]].ln();
    for i in 5..nrows {
        let log_x = psd[[i, 0]].ln();
        psd[[i, 1]] = (yn + sn * (log_x - xn)).exp();
    }

    let mut flipped = Array2::<f64>::zeros((nrows, 3));
    for i in 0..nrows {
        flipped.row_mut(i).assign(&psd.row(nrows - 1 - i));
    }

    // Optional Gaussian noise in log power
    if add_noise {
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
/// * `psd` - A (nfreq,3) array where the first column contains wavelengths, the second column contains power values, and the third column contains the phases, which define the variability of the profile.
/// # Returns
/// * A 1D array of values corresponding to the input signal psd, representing the linear profile generated from the PSD and phase information.
///
pub fn compute_profile_from_psd(
    crater_radius: f64,
    ymean: f64,
    psd: ArrayView2<'_, f64>,
) -> Vec<f64> {
    let nfreq = psd.nrows();
    let num_points = 2 * nfreq;
    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(num_points);
    let mut buffer: Vec<Complex<f64>> = vec![Complex { re: 0.0, im: 0.0 }; num_points];

    for i in 0..nfreq {
        let freq_f64 = TAU / psd[[i, 0]];
        let k = freq_f64.round() as usize;

        // Safety check to prevent out-of-bounds if a wavenumber exceeds the Nyquist limit
        if k >= num_points / 2 {
            continue;
        }

        let amplitude = (psd[[i, 1]] * TAU).sqrt();
        let phase = psd[[i, 2]];

        if k == 0 {
            // DC Component (wavenumber 0)
            buffer[0] = Complex::from_polar(amplitude, phase);
        } else {
            // Positive frequency (index k)
            buffer[k] = Complex::from_polar(amplitude / 2.0, phase);
            // Negative frequency (index N - k) uses the negative phase (conjugate)
            buffer[num_points - k] = buffer[k].conj();
        }
    }

    ifft.process(&mut buffer);

    let mut out: Vec<f64> = Vec::with_capacity(num_points);
    for val in buffer.iter() {
        let z = val.re;
        out.push(crater_radius * z + ymean);
    }

    let first_val = out[0];
    let last_val = out[out.len() - 1];
    out.insert(0, last_val);
    out.push(first_val);

    out
}
