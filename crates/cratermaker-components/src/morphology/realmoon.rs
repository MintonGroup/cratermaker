use crate::ArrayResult;
use crate::morphology::basicmoon::{BasicMoonCrater, basicmoon_profile_one};
use crate::morphology::blend_band;
use numpy::ndarray::prelude::*;
use paste::paste;
use rand::SeedableRng;
use rand::prelude::*;
use rand_chacha::ChaCha12Rng;
use rand_distr::{Normal, Uniform};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rustfft::{FftPlanner, num_complex::Complex};
use std::collections::HashMap;
use std::f64::consts::TAU;

const BLEND_WIDTH_FRACTION: f64 = 0.05; // Relative width of the blending region between PSD-derived parameter values and their means

/// Represents a local region of a surface mesh with various attributes accessible as array views.
pub struct PSD1DView<'a> {
    pub nprofile: usize,
    pub nfreq: usize,
    pub normalization_length: f64,
    pub mean: f64,
    pub pix: f64,
    pub wavelength: ArrayView1<'a, f64>,
    pub amplitude: ArrayView1<'a, f64>,
    pub phase: ArrayView1<'a, f64>,
}

// This will build a boilerplate struct that contains the PSD objects for all of the crater object properties with the pattern [name]_psd.
macro_rules! define_crater_structs {
    (
        // List the base float fields
        floats: [ $($float_field:ident),* ],
        // List the PSD features (e.g., radius, floor_radius, wall_texture)
        psds: [ $($psd_field:ident),* ]
    ) => {
        // 1. Generate PyReadonlyRealMoonCrater
        paste! {
            pub struct RealMoonCrater<'a> {
                $( pub $float_field: f64, )*
                $( pub [<$psd_field _psd>]: PSD1DView<'a>, )*
            }
        }
    };
}

define_crater_structs! {
    floats: [
        diameter, radius, floor_elevation, floor_radius, wall_curvature,
        floor_blend, rim_width, rim_height, frac_ejrim, ejprofile,
        peak_height, peak_width, peak_ring_radius, peak_center_distance,
        peak_center_bearing
    ],
    psds: [ radius, floor_radius, rim_width, rim_height]
}

/// Evaluates a periodic bearing profile generated from a PSD by piecewise-linear interpolation.
///
/// The PSD profile is sampled at `nprofile` uniformly spaced bearings over [0, 2π), so the last sample sits one interval short of
/// 2π. The profile is periodic, so the grid is padded with one wrapped sample on each side (the last sample shifted by -2π in
/// front and the first sample shifted by +2π at the end). Together with wrapping each query bearing into [0, 2π), this makes the
/// interpolation continuous across the 0/2π seam instead of extrapolating past the last sample.
///
/// The segment slopes are computed once at construction, so each evaluation is a binary search plus one multiply-add.
pub struct ProfileEvaluator {
    /// Padded bearing grid: [x[n-1] - 2π, x[0], ..., x[n-1], x[0] + 2π]
    bearing_grid: Vec<f64>,
    /// Profile values matching `bearing_grid`: [y[n-1], y[0], ..., y[n-1], y[0]]
    profile_values: Vec<f64>,
    /// slopes[i] is the slope of the segment between bearing_grid[i] and bearing_grid[i + 1]
    slopes: Vec<f64>,
}

impl ProfileEvaluator {
    pub fn new(psd: &PSD1DView) -> Self {
        let (bearing, profile) = profile_from_psd(psd);
        Self::from_samples(&bearing, &profile)
    }

    /// Builds an evaluator from one period of uniformly ordered samples over [0, 2π).
    pub fn from_samples(bearing: &[f64], profile: &[f64]) -> Self {
        let n = bearing.len().min(profile.len());
        if n == 0 {
            return Self {
                bearing_grid: Vec::new(),
                profile_values: Vec::new(),
                slopes: Vec::new(),
            };
        }
        let mut bearing_grid = Vec::with_capacity(n + 2);
        let mut profile_values = Vec::with_capacity(n + 2);
        bearing_grid.push(bearing[n - 1] - TAU);
        profile_values.push(profile[n - 1]);
        bearing_grid.extend_from_slice(&bearing[..n]);
        profile_values.extend_from_slice(&profile[..n]);
        bearing_grid.push(bearing[0] + TAU);
        profile_values.push(profile[0]);

        let slopes = bearing_grid
            .windows(2)
            .zip(profile_values.windows(2))
            .map(|(x, y)| (y[1] - y[0]) / (x[1] - x[0]))
            .collect();
        Self {
            bearing_grid,
            profile_values,
            slopes,
        }
    }

    #[inline(always)]
    pub fn eval(&self, bearing: f64) -> f64 {
        if self.slopes.is_empty() {
            return self.profile_values.first().copied().unwrap_or(0.0);
        }
        // Wrap into [0, 2π]. rem_euclid can round up to exactly 2π for tiny negative inputs, which the padded grid also covers.
        let b = bearing.rem_euclid(TAU);
        // Index of the segment containing b, clamped to the valid range of segments
        let i = self
            .bearing_grid
            .partition_point(|&probe| probe <= b)
            .saturating_sub(1)
            .min(self.slopes.len() - 1);
        self.profile_values[i] + (b - self.bearing_grid[i]) * self.slopes[i]
    }
}

macro_rules! define_evaluator_bundle {
    (
        psds: [ $($psd_field:ident),* ]
    ) => {
        paste! {
            /// Bundles the interpolated spatial profiles for all dynamic PSD fields
            pub struct CraterProfileEvaluators {
                $( pub [<$psd_field _eval>]: ProfileEvaluator, )*
            }

            impl CraterProfileEvaluators {
                pub fn new(crater: &RealMoonCrater) -> Self {
                    Self {
                        $( [<$psd_field _eval>]: ProfileEvaluator::new(&crater.[<$psd_field _psd>]), )*
                    }
                }

                /// Evaluates all dynamic properties at a given bearing, returning
                /// a customized BasicMoonCrater mapping to this slice of the profile.
                #[inline]
                pub fn eval_to_basic(&self, base_crater: &RealMoonCrater, r: f64, bearing: f64) -> BasicMoonCrater {
                    // Evaluate all dynamic PSD profiles
                    $(
                        let [<$psd_field _raw>] = self.[<$psd_field _eval>].eval(bearing);
                    )*

                    // Blend PSD-derived parameters with their means when the radius value is far from its boundary
                    // This smooths out the sharp radial features in the full crater profile
                    let rtrans_rim = radius_raw;
                    let rwidth_rim = BLEND_WIDTH_FRACTION * rtrans_rim;
                    let rtrans_floor = floor_radius_raw;
                    let rwidth_floor = 2.0 * BLEND_WIDTH_FRACTION * rtrans_floor;

                    let radius_val = blend_band(r, radius_raw, base_crater.radius, rtrans_rim, rwidth_rim);
                    let rim_width_val = blend_band(r, rim_width_raw, base_crater.rim_width, rtrans_rim, rwidth_rim);
                    let rim_height_val = blend_band(r, rim_height_raw, base_crater.rim_height, rtrans_rim, rwidth_rim);
                    let floor_radius_val = blend_band(r, floor_radius_raw, base_crater.floor_radius, rtrans_floor, rwidth_floor);

                    let diameter_val = 2.0 * radius_val;
                    BasicMoonCrater {
                        diameter: diameter_val,
                        radius: radius_val,
                        floor_radius: floor_radius_val,
                        rim_height: rim_height_val,
                        rim_width: rim_width_val,
                        // Static non-PSD fields
                        frac_ejrim: base_crater.frac_ejrim,
                        ejprofile: base_crater.ejprofile,
                        floor_elevation: base_crater.floor_elevation,
                        peak_height: base_crater.peak_height,
                        peak_width: base_crater.peak_width,
                        peak_ring_radius: base_crater.peak_ring_radius,
                        peak_center_distance: base_crater.peak_center_distance,
                        peak_center_bearing: base_crater.peak_center_bearing,
                        wall_curvature: base_crater.wall_curvature,
                        floor_blend: base_crater.floor_blend,
                    }
                }
            }
        }
    };
}

// Generate the dynamic bundle automatically for all psd-controlled parameters
define_evaluator_bundle! {
    psds: [ radius, floor_radius, rim_width, rim_height]
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
        "Input arrays must have same length"
    );
    assert_eq!(
        radial_distances.len(),
        bearings.len(),
        "Input arrays must have same length"
    );

    let n_points = radial_distances.len();

    // 0. Compute boundary thresholds
    let ninc = radial_distances
        .iter()
        .filter(|&&x| x <= crater.radius)
        .count();
    let meanref = if ninc == 0 {
        *radial_distances
            .iter()
            .zip(reference_elevations)
            .min_by(|&(&ra, _), &(&rb, _)| ra.partial_cmp(&rb).unwrap())
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

    let main_evaluators = CraterProfileEvaluators::new(crater);

    let rings_evaluators: Option<Vec<(CraterProfileEvaluators, &RealMoonCrater)>> =
        rings.as_ref().map(|rings_vec| {
            rings_vec
                .iter()
                .map(|ring| (CraterProfileEvaluators::new(ring), ring))
                .collect()
        });

    let out: Vec<f64> = (0..n_points)
        .into_par_iter()
        .map(|i| {
            let r = radial_distances[i];
            let href = reference_elevations[i];
            let bearing = bearings[i];

            // Generate the Basic representation dynamically at this local bearing
            let icrater = main_evaluators.eval_to_basic(crater, r, bearing);

            // Generate Basic representations for all sub-rings at this bearing
            let irings: Option<Vec<BasicMoonCrater>> = rings_evaluators.as_ref().map(|evals_vec| {
                evals_vec
                    .iter()
                    .map(|(evaluator, ring)| evaluator.eval_to_basic(ring, r, bearing))
                    .collect()
            });

            // Compute the profile using the basic profile function using the values at this bearing value
            let h =
                basicmoon_profile_one(r, href, &icrater, &irings, include_crater, include_ejecta);

            // Enforce minimum crater floor constraint
            if r <= icrater.radius {
                h.max(min_elevation)
            } else {
                h
            }
        })
        .collect();

    Ok(Array1::from_vec(out))
}

///
///
/// Computes a target 1D power spectral density distribution based on control points defining a piecewise linear function in log-log space, with optional Gaussian noise added to the log amplitude values.
///
/// # Arguments
/// * `control_points` - A dictionary containing the control points for the piecewise linear function. The expected keys are:
/// - "sn": The slope of the first segment (y2-y1)/(x2-x1)
/// - "yn": The y-coordinate of the first breakpoint in ln(amplitude).
/// - "y1": The x-coordinate of the second breakpoint in ln(wavelength).
/// - "y2": The y-coordinate of the second breakpoint in ln(amplitude).
/// - "y3": The y-coordinate at the third breakpoint (2nd highest wavelength) in log(amplitude).
/// - "y4": The y-coordinate of the highest wavelength in ln(amplitude).
/// - "y5": The y-coordinate of the highest wavelength in ln(amplitude).
///  * `nprofile` - The number of points in the output PSD, which determines the wavelength resolution and the Nyquist frequency.
///  * `add_noise` - Whether to add Gaussian noise to the ln(amplitude) values to simulate natural variability in the PSD.
///  * `seed` - The random seed for reproducibility of the noise if `add_noise` is true.
///
///
/// # Returns
///
/// * A tuple of Arrays containing the wavelength, amplitude, and phase angles of the AC components of the signal. The DC component is left unfilled.
///
pub fn get_1d_psd_from_control_points(
    control_points: &HashMap<String, f64>,
    mean: f64,
    nprofile: usize,
    add_noise: bool,
    rng_seed: u64,
) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>), String> {
    let sn = control_points["sn"];
    let yn = control_points["yn"];
    let y1 = control_points["y1"];
    let y2 = control_points["y2"];
    let y3 = control_points["y3"];
    let y4 = control_points["y4"];
    let y5 = control_points["y5"];

    let interval = TAU / nprofile as f64;

    // Equivalent to rfft sizing in Python:
    let nfreq = nprofile / 2 + 1;

    let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
    let mut wavelength = Array1::<f64>::zeros(nfreq);
    let mut amplitude = Array1::<f64>::zeros(nfreq);
    let mut phase = Array1::<f64>::zeros(nfreq);
    let base = nprofile as f64 * interval;
    let uniform = Uniform::new(0.0, TAU).expect("valid uniform distribution");
    wavelength[0] = f64::NAN;
    phase[0] = 0.0;
    for i in 1..nfreq {
        wavelength[i] = base / i as f64;
        phase[i] = uniform.sample(&mut rng); // randomized phases
    }
    amplitude[0] = mean;
    amplitude[1] = y1.exp();
    amplitude[2] = y2.exp();
    amplitude[3] = y3.exp();
    amplitude[4] = y4.exp();
    amplitude[5] = y5.exp();

    let xn = wavelength[6].ln();
    for i in 6..nfreq {
        let log_x = wavelength[i].ln();
        amplitude[i] = (yn + sn * (log_x - xn)).exp();
    }

    // Optional Gaussian noise in log amplitude
    if add_noise {
        let normal = Normal::new(0.0, 0.55).expect("valid normal distribution");
        for i in 1..nfreq {
            let log_amplitude = amplitude[i].ln();
            let noisy_log_amplitude = log_amplitude + normal.sample(&mut rng);
            amplitude[i] = (noisy_log_amplitude).exp();
        }
    }
    Ok((wavelength, amplitude, phase))
}

///
/// Generates a surface profile based on a 1D power spectral density (PSD) and optional phase information, simulating a crater surface with specified roughness characteristics.
///
/// The AC components (index 1 and up) describe the spectrum of the natural log of the profile, and the DC component
/// (`amplitude[0]`) is the arithmetic mean of the profile in units of `normalization_length`. The log profile is synthesized
/// with zero mean, exponentiated, and rescaled so that its arithmetic mean equals `amplitude[0]`. This keeps the profile
/// positive for any spectrum, and gives it the right-skewed distribution seen in measured crater profiles.
///
/// # Arguments
/// * `psd` - The struct containing wavelength, amplitude, and phase, which define the variability of the profile.
/// # Returns
/// * A 1D array of values corresponding to the input signal psd, representing the linear profile generated from the PSD and phase information.
///
pub fn profile_from_psd(psd: &PSD1DView) -> (Vec<f64>, Vec<f64>) {
    let wavelength = psd.wavelength;
    let amplitude = psd.amplitude;
    let phase = psd.phase;

    // nfreq is the number of AC components. Total components = nfreq + 1 (for DC).
    let nfreq = wavelength.len();
    let nprofile = 2 * (nfreq - 1);
    let nyquist_limit = nprofile / 2;

    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(nprofile);
    // The log profile has zero mean. The mean of the profile itself is applied after exponentiating.
    let mut buffer: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nprofile];

    // The AC components
    for k in 1..nyquist_limit {
        let amp = amplitude[k] * nprofile as f64 / 2.0;
        let complex_val = Complex::from_polar(amp, phase[k]);
        buffer[k] = complex_val;
        buffer[nprofile - k] = complex_val.conj();
    }

    // The Nyquist component is at index nfreq - 1
    if nprofile % 2 == 0 {
        buffer[nyquist_limit] = Complex::from_polar(
            amplitude[nyquist_limit] * nprofile as f64,
            phase[nyquist_limit],
        );
    }
    ifft.process(&mut buffer);

    // Exponentiate the log profile, then rescale so the arithmetic mean is amplitude[0] in units of normalization_length
    let mut out: Vec<f64> = buffer
        .iter()
        .map(|val| (val.re / nprofile as f64).exp())
        .collect();
    let exp_mean = out.iter().sum::<f64>() / nprofile as f64;
    let scale_factor = amplitude[0] * psd.normalization_length / exp_mean;
    out.iter_mut().for_each(|v| *v *= scale_factor);

    // Generate the corresponding bearing values for this profile
    let mut bearing: Vec<f64> = Vec::with_capacity(nprofile);
    for i in 0..nprofile {
        bearing.push(TAU * (i as f64 / nprofile as f64));
    }
    (bearing, out)
}

pub fn psd_from_profile(profile: &ArrayView1<'_, f64>) -> (Vec<f64>, Vec<f64>) {
    let nprofile = profile.len();
    let nyquist_limit = nprofile / 2;

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(nprofile);
    let mut buffer: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); nprofile];

    // Reconstruct the full spectrum. The loop now handles the DC component (k=0) naturally.
    for k in 0..nprofile {
        buffer[k] = profile[k].into();
    }

    fft.process(&mut buffer);

    // Slice only up to the Nyquist limit (N/2 + 1) to get the standard one-sided PSD representation
    let one_sided_len = (nprofile / 2) + 1;
    let mut amplitude = Vec::with_capacity(one_sided_len);
    let mut phase = Vec::with_capacity(one_sided_len);

    for k in 0..one_sided_len {
        let (amp, ph) = buffer[k].to_polar();
        let mut scaled_amp = amp / nprofile as f64;
        if k > 0 && k < nyquist_limit {
            scaled_amp *= 2.0;
        }
        amplitude.push(scaled_amp);
        phase.push(ph);
    }

    (amplitude, phase)
}

#[cfg(test)]
mod tests {
    use super::*;
    use interp::{InterpMode, interp};

    /// One period of a random PSD profile, sampled on the uniform [0, 2π) grid used by profile_from_psd.
    fn random_profile(nfreq: usize, rng: &mut ChaCha12Rng) -> (Vec<f64>, Vec<f64>) {
        let normal = Normal::<f64>::new(0.0, 1.0).unwrap();
        let nprofile = 2 * (nfreq - 1);
        let wavelength = Array1::from_iter((0..nfreq).map(|i| TAU / i.max(1) as f64));
        let amplitude = Array1::from_iter((0..nfreq).map(|_| normal.sample(rng).abs()));
        let phase = Array1::from_iter((0..nfreq).map(|_| normal.sample(rng)));
        let psd = PSD1DView {
            nprofile,
            nfreq,
            normalization_length: 1234.5,
            mean: amplitude[0],
            pix: 1.0,
            wavelength: wavelength.view(),
            amplitude: amplitude.view(),
            phase: phase.view(),
        };
        profile_from_psd(&psd)
    }

    #[test]
    fn profile_evaluator_is_exact_at_samples_and_periodic() {
        let mut rng = ChaCha12Rng::seed_from_u64(42);
        for &nfreq in &[2usize, 3, 9, 33, 257] {
            let (x, y) = random_profile(nfreq, &mut rng);
            let n = x.len();
            let evaluator = ProfileEvaluator::from_samples(&x, &y);

            // Exact at every sample, including one period either side
            for k in 0..n {
                assert_eq!(evaluator.eval(x[k]), y[k], "nfreq={nfreq}, sample {k}");
            }

            // Continuous across the seam: just below 2π approaches y[0], and the last interval interpolates toward y[0]
            let dx = TAU / n as f64;
            let tol = 1e-9 * y.iter().fold(1.0f64, |a, &v| a.max(v.abs()));
            assert!((evaluator.eval(TAU.next_down()) - y[0]).abs() < tol, "nfreq={nfreq}: no continuity at 2π");
            assert!((evaluator.eval(TAU) - y[0]).abs() < tol);
            let mid_last = x[n - 1] + 0.5 * dx;
            assert!((evaluator.eval(mid_last) - 0.5 * (y[n - 1] + y[0])).abs() < tol, "nfreq={nfreq}: seam interval");

            // Matches a periodically padded linear-interpolation reference, and is periodic outside [0, 2π)
            let mut xp = vec![x[n - 1] - TAU];
            xp.extend_from_slice(&x);
            xp.push(x[0] + TAU);
            let mut yp = vec![y[n - 1]];
            yp.extend_from_slice(&y);
            yp.push(y[0]);
            for _ in 0..5000 {
                let b: f64 = rng.random_range(0.0..TAU);
                let got = evaluator.eval(b);
                let expected = interp(&xp, &yp, b, &InterpMode::Extrapolate);
                assert!((got - expected).abs() < tol, "nfreq={nfreq}, bearing={b}: {got} != {expected}");
                for shift in [-2.0 * TAU, -TAU, TAU, 3.0 * TAU] {
                    assert!((evaluator.eval(b + shift) - got).abs() < tol, "nfreq={nfreq}: not periodic at {b}+{shift}");
                }
            }
            // Tiny negative bearings wrap to just below 2π, which must be close to y[0]
            assert!((evaluator.eval(-1e-300) - y[0]).abs() < tol);
            assert!((evaluator.eval(-0.0) - y[0]).abs() < tol);
        }
    }

    #[test]
    fn profile_from_psd_is_positive_and_keeps_mean() {
        let mut rng = ChaCha12Rng::seed_from_u64(7);
        for &nfreq in &[3usize, 9, 33, 257, 2049] {
            let mean = 0.3;
            let normalization_length = 1234.5;
            // Large log amplitudes that would drive a linear synthesis well below zero
            let wavelength = Array1::from_iter((0..nfreq).map(|i| TAU / i.max(1) as f64));
            let mut amplitude = Array1::from_iter((0..nfreq).map(|_| 0.5 + rng.random_range(0.0..0.5)));
            amplitude[0] = mean;
            let phase = Array1::from_iter((0..nfreq).map(|_| rng.random_range(0.0..TAU)));
            let psd = PSD1DView {
                nprofile: 2 * (nfreq - 1),
                nfreq,
                normalization_length,
                mean,
                pix: 1.0,
                wavelength: wavelength.view(),
                amplitude: amplitude.view(),
                phase: phase.view(),
            };
            let (_, y) = profile_from_psd(&psd);
            assert!(y.iter().all(|&v| v > 0.0), "nfreq={nfreq}: profile has non-positive values");
            let ymean = y.iter().sum::<f64>() / y.len() as f64;
            assert!((ymean / (mean * normalization_length) - 1.0).abs() < 1e-12, "nfreq={nfreq}: mean not preserved");
        }
    }
}
