use crate::morphology::basicmoon::{BasicMoonCrater, basicmoon_profile_one};
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
    pub floor_blend: f64,
    pub rim_width: f64,
    pub rim_height: f64,
    pub frac_ejrim: f64,
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
    rings: &Option<Vec<RealMoonCrater>>,
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
    let min_elevation = meanref + crater.floor_elevation;

    // Create profile functions that will be interpolated later
    let (rimtheta, rim_profile) =
        compute_profile_from_psd(crater.radius, crater.radius, crater.rim_radius_psd);
    let (floortheta, floor_profile) =
        compute_profile_from_psd(crater.radius, crater.floor_radius, crater.floor_radius_psd);
    let rings_rim_profile: Option<Vec<(Vec<f64>, Vec<f64>)>> = rings.as_ref().map(|rings_vec| {
        rings_vec
            .iter()
            .map(|ring| compute_profile_from_psd(ring.radius, ring.radius, ring.rim_radius_psd))
            .collect()
    });
    let rings_floor_profile: Option<Vec<(Vec<f64>, Vec<f64>)>> = rings.as_ref().map(|rings_vec| {
        rings_vec
            .iter()
            .map(|ring| {
                compute_profile_from_psd(ring.radius, ring.floor_radius, ring.floor_radius_psd)
            })
            .collect()
    });

    let out: Vec<f64> = (0..n_points)
        .into_par_iter()
        .map(|i| {
            let r = radial_distances[i];
            let href = reference_elevations[i];
            let theta = bearings[i];
            let rim_r = interp(&rimtheta, &rim_profile, theta, &InterpMode::default());
            let floor_r = interp(&floortheta, &floor_profile, theta, &InterpMode::default());
            let rim_elev =
                crater.rim_height * (rim_r / crater.radius) * (crater.floor_radius / floor_r);
            let icrater = realtobasic(crater, rim_r, floor_r, rim_elev);
            let irings: Option<Vec<BasicMoonCrater>> =
                if let (Some(rings_vec), Some(rim_profiles), Some(floor_profiles)) =
                    (&rings, &rings_rim_profile, &rings_floor_profile)
                {
                    let new_rings = rings_vec
                        .iter()
                        .zip(rim_profiles.iter())
                        .zip(floor_profiles.iter())
                        .map(
                            |((ring, (irimtheta, irimprofile)), (ifloortheta, ifloorprofile))| {
                                let iring_rim_r =
                                    interp(irimtheta, irimprofile, theta, &InterpMode::default());
                                let iring_floor_r = interp(
                                    ifloortheta,
                                    ifloorprofile,
                                    theta,
                                    &InterpMode::default(),
                                );

                                let irim_elev = ring.rim_height
                                    * (iring_rim_r / ring.radius)
                                    * (ring.floor_radius / iring_floor_r);

                                realtobasic(&ring, iring_rim_r, iring_floor_r, irim_elev)
                            },
                        )
                        .collect();
                    Some(new_rings)
                } else {
                    None
                };

            let h =
                basicmoon_profile_one(r, href, &icrater, &irings, include_crater, include_ejecta);
            if r <= rim_r { h.max(min_elevation) } else { h }
        })
        .collect();

    Ok(Array1::from_vec(out))
}

fn realtobasic(
    crater: &RealMoonCrater,
    radius: f64,
    floor_radius: f64,
    rim_height: f64,
) -> BasicMoonCrater {
    BasicMoonCrater {
        diameter: 2.0 * radius,
        radius: radius,
        floor_radius: floor_radius,
        rim_height: rim_height,
        floor_elevation: crater.floor_elevation,
        wall_curvature: crater.wall_curvature,
        floor_blend: crater.floor_blend,
        rim_width: crater.rim_width,
        frac_ejrim: crater.frac_ejrim,
        ejprofile: crater.ejprofile,
        peak_height: crater.peak_height,
        peak_width: crater.peak_width,
        peak_ring_radius: crater.peak_ring_radius,
        peak_center_distance: crater.peak_center_distance,
        peak_center_bearing: crater.peak_center_bearing,
        elevation_offset: crater.elevation_offset,
    }
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
) -> (Vec<f64>, Vec<f64>) {
    let nfreq = psd.nrows();
    let mut num_points = 2 * nfreq;
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
    num_points = out.len();

    let mut theta: Vec<f64> = Vec::with_capacity(num_points);
    for i in 0..num_points {
        theta.push(TAU * (i as f64 - 1.0) / (num_points as f64 - 2.0));
    }

    (theta, out)
}
