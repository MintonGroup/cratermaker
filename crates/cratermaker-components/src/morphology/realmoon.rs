use crate::ArrayResult;
use crate::morphology::basicmoon::{BasicMoonCrater, basicmoon_profile_one_local, crater_elevation};
use crate::morphology::smootherstep;
use numpy::ndarray::prelude::*;
use paste::paste;
use rand::SeedableRng;
use rand::prelude::*;
use rand_chacha::ChaCha12Rng;
use rand_distr::Uniform;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rustfft::{FftPlanner, num_complex::Complex};
use std::collections::HashMap;
use std::f64::consts::TAU;


/// Represents a local region of a surface mesh with various attributes accessible as array views.
pub struct PSD1DView<'a> {
    pub nprofile: usize,
    pub nfreq: usize,
    pub normalization_length: f64,
    pub mean: f64,
    pub pix: f64,
    /// Target skewness of the log profile. Zero leaves the phases as they are.
    pub skewness: f64,
    /// Factor applied to the log profile (to all of its non-mean Fourier amplitudes)
    pub log_scale: f64,
    /// Bounds on ln(profile / mean profile); values beyond them are clamped
    pub log_ratio_min: f64,
    pub log_ratio_max: f64,
    /// Largest value of the profile as a fraction of the crater radius at the same bearing
    pub max_radius_fraction: f64,
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

                /// The parameters of the crater (or ring) along one bearing, from its PSD profiles.
                #[inline]
                pub fn at_bearing(&self, base_crater: &RealMoonCrater, bearing: f64) -> BearingFeature {
                    $(
                        let [<$psd_field _raw>] = self.[<$psd_field _eval>].eval(bearing);
                    )*
                    // Geometric guards: no value a larger fraction of the radius than in any measured profile
                    let floor_radius = floor_radius_raw.min(base_crater.floor_radius_psd.max_radius_fraction * radius_raw);
                    let rim_width = rim_width_raw.min(base_crater.rim_width_psd.max_radius_fraction * radius_raw);
                    let rim_height = rim_height_raw.min(base_crater.rim_height_psd.max_radius_fraction * radius_raw);
                    BearingFeature::new(base_crater, radius_raw, floor_radius, rim_width, rim_height)
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

            // The parameters of the crater and its rings along this bearing
            let main = main_evaluators.at_bearing(crater, bearing);
            let ring_features: Vec<BearingFeature> = rings_evaluators
                .as_ref()
                .map(|evals_vec| evals_vec.iter().map(|(evaluator, ring)| evaluator.at_bearing(ring, bearing)).collect())
                .unwrap_or_default();

            // Each feature is its mean profile evaluated at its own warped radial coordinate
            let (icrater, r_main) = main.at(crater, r);
            let (irings, ring_r): (Vec<BasicMoonCrater>, Vec<f64>) = rings
                .as_ref()
                .map(|rings_vec| ring_features.iter().zip(rings_vec).map(|(f, ring)| f.at(ring, r)).unzip())
                .unwrap_or_default();

            // Each ring sits at the elevation of the feature around it (the crater for the outermost ring, otherwise the next
            // ring out) where the ring's rim is along this bearing
            let offsets: Vec<f64> = match rings {
                Some(rings_vec) => (0..rings_vec.len())
                    .map(|j| {
                        let rj = ring_features[j].radius;
                        let (c, rc) = if j == 0 { main.at(crater, rj) } else { ring_features[j - 1].at(&rings_vec[j - 1], rj) };
                        crater_elevation(&c, rc)
                    })
                    .collect(),
                None => Vec::new(),
            };

            let irings = if irings.is_empty() { None } else { Some(irings) };
            let h = basicmoon_profile_one_local(r_main, &ring_r, href, &icrater, &irings, &offsets, include_crater, include_ejecta);

            // Enforce minimum crater floor constraint
            if r_main <= icrater.radius {
                h.max(min_elevation)
            } else {
                h
            }
        })
        .collect();

    Ok(Array1::from_vec(out))
}

/// Distance, in mean rim widths of a crater (or ring), over which the effect of its parameters along a bearing fades out:
/// beyond its rim zone the rim height and rim width return to their means, and away from its rim and floor edge the radial
/// warp that places them returns to the identity.
const LOCALITY_RIM_WIDTHS: f64 = 2.0;

/// The parameters of a crater (or ring) along one bearing, and how they enter its profile.
///
/// The profile along a bearing is the mean profile of the crater evaluated at a warped radial coordinate. The warp moves the
/// mean floor edge and rim to where they are along this bearing, and is the identity away from them, so the position of a
/// feature affects the profile only near that feature, and each feature appears exactly once. The rim height and rim width
/// take their values along this bearing across the rim zone and return to their means away from it, so they too affect the
/// profile only near the rim.
#[derive(Clone, Copy, Debug)]
pub struct BearingFeature {
    pub radius: f64,
    pub floor_radius: f64,
    pub rim_width: f64,
    pub rim_height: f64,
    warp: RadialWarp,
}

impl BearingFeature {
    pub fn new(base: &RealMoonCrater, radius: f64, floor_radius: f64, rim_width: f64, rim_height: f64) -> Self {
        let scale = (LOCALITY_RIM_WIDTHS * base.rim_width).max(1e-6 * base.radius);
        Self {
            radius,
            floor_radius,
            rim_width,
            rim_height,
            warp: RadialWarp::new(floor_radius, base.floor_radius, radius, base.radius, scale),
        }
    }

    /// The crater to evaluate at radial distance `r` along this bearing, and the radial coordinate to evaluate it at.
    #[inline]
    pub fn at(&self, base: &RealMoonCrater, r: f64) -> (BasicMoonCrater, f64) {
        let rw = self.warp.eval(r);
        let scale = (LOCALITY_RIM_WIDTHS * base.rim_width).max(1e-6 * base.radius);
        let beyond = ((rw - base.radius).abs() - base.rim_width).max(0.0);
        let phi = smootherstep(beyond / scale);
        let crater = BasicMoonCrater {
            diameter: 2.0 * base.radius,
            radius: base.radius,
            floor_radius: base.floor_radius,
            rim_height: (1.0 - phi) * self.rim_height + phi * base.rim_height,
            rim_width: (1.0 - phi) * self.rim_width + phi * base.rim_width,
            frac_ejrim: base.frac_ejrim,
            ejprofile: base.ejprofile,
            floor_elevation: base.floor_elevation,
            peak_height: base.peak_height,
            peak_width: base.peak_width,
            peak_ring_radius: base.peak_ring_radius,
            peak_center_distance: base.peak_center_distance,
            peak_center_bearing: base.peak_center_bearing,
            wall_curvature: base.wall_curvature,
            floor_blend: base.floor_blend,
        };
        (crater, rw)
    }
}

/// A smooth, monotone map of the radial distance along a bearing to the radial coordinate of the mean crater profile.
///
/// It maps the floor radius along the bearing to the mean floor radius and the rim radius along the bearing to the mean rim
/// radius, and is the identity inside `floor_radius - l_floor` and beyond `radius + l_rim`, where each length is `scale`, or
/// twice the shift of its feature if that is larger (which keeps the map monotone). Between the knots it is a monotone cubic
/// Hermite interpolant (Fritsch-Butland slopes, with slope 1 where it joins the identity).
#[derive(Clone, Copy, Debug)]
struct RadialWarp {
    x: [f64; 5],
    y: [f64; 5],
    d: [f64; 5],
    n: usize,
}

impl RadialWarp {
    fn new(floor_radius: f64, mean_floor_radius: f64, radius: f64, mean_radius: f64, scale: f64) -> Self {
        let l_floor = scale.max(2.0 * (floor_radius - mean_floor_radius).abs());
        let l_rim = scale.max(2.0 * (radius - mean_radius).abs());
        let mut x = [0.0; 5];
        let mut y = [0.0; 5];
        let mut fixed = [false; 5];
        let mut n = 0;
        let inner = floor_radius - l_floor;
        if inner > 0.0 {
            x[n] = inner;
            y[n] = inner;
            fixed[n] = true;
        } else {
            x[n] = 0.0;
            y[n] = 0.0;
        }
        n += 1;
        for (xi, yi) in [(floor_radius, mean_floor_radius), (radius, mean_radius)] {
            // Keep the knots strictly increasing (the floor edge can only coincide with the rim in degenerate profiles)
            if xi > x[n - 1] && yi > y[n - 1] {
                x[n] = xi;
                y[n] = yi;
                n += 1;
            }
        }
        x[n] = radius + l_rim;
        y[n] = radius + l_rim;
        fixed[n] = true;
        n += 1;

        let mut sec = [0.0; 4];
        for k in 0..n - 1 {
            sec[k] = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
        }
        let mut d = [0.0; 5];
        for k in 0..n {
            d[k] = if fixed[k] {
                1.0
            } else if k == 0 {
                sec[0]
            } else if k == n - 1 {
                sec[n - 2]
            } else {
                // Weighted harmonic mean of the neighboring secants (Fritsch-Butland), which keeps the interpolant monotone
                let (h0, h1) = (x[k] - x[k - 1], x[k + 1] - x[k]);
                let (w0, w1) = (2.0 * h1 + h0, h1 + 2.0 * h0);
                (w0 + w1) / (w0 / sec[k - 1] + w1 / sec[k])
            };
        }
        Self { x, y, d, n }
    }

    #[inline]
    fn eval(&self, r: f64) -> f64 {
        if r <= self.x[0] {
            // Inside the first knot the map is the identity, or a scaling if the first knot is at the center
            return if self.x[0] > 0.0 { r } else { r * self.d[0] };
        }
        if r >= self.x[self.n - 1] {
            return r;
        }
        let mut k = 0;
        while r > self.x[k + 1] {
            k += 1;
        }
        let h = self.x[k + 1] - self.x[k];
        let t = (r - self.x[k]) / h;
        let (t2, t3) = (t * t, t * t * t);
        (2.0 * t3 - 3.0 * t2 + 1.0) * self.y[k]
            + (t3 - 2.0 * t2 + t) * h * self.d[k]
            + (-2.0 * t3 + 3.0 * t2) * self.y[k + 1]
            + (t3 - t2) * h * self.d[k + 1]
    }
}

/// Highest wavenumber constrained by the measured profiles, which are sampled at 360 bearings (Nyquist at k = 180)
const PSD_DATA_KMAX: usize = 180;

/// Minimum log-log slope of ln(amplitude) vs ln(wavelength) used to extend the spectrum beyond PSD_DATA_KMAX.
///
/// The variance contributed by wavenumbers above the data limit scales as Σ k^(-2s), which only converges for s > 0.5. Many
/// measured slopes are shallower than that, so extending them unchanged would make the profile variance grow without bound
/// as the surface resolution increases. Steepening the tail to at least this slope keeps the structure finer than the data
/// resolution while making the profile statistics independent of resolution once it is fine enough.
const PSD_MIN_EXTRAPOLATION_SLOPE: f64 = 2.0;

/// Euler-Mascheroni constant
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Mean of ln(R) for a Rayleigh variate R with unit scale: E[ln R] = ln(2)/2 - γ/2
const LN_RAYLEIGH_MEAN: f64 = 0.5 * std::f64::consts::LN_2 - 0.5 * EULER_GAMMA;

/// Converts a uniform variate u in [0, 1) into a Rayleigh variate scaled so that the mean of its log is zero.
///
/// The log of the result has standard deviation π/√24 ≈ 0.641 and skewness ≈ -1.14.
#[inline]
fn rayleigh_log_centered(u: f64) -> f64 {
    let r = (-2.0 * (1.0 - u).ln()).sqrt();
    r * (-LN_RAYLEIGH_MEAN).exp()
}

///
///
/// Computes a target 1D power spectral density distribution based on control points defining a piecewise linear function in log-log space, with optional Rayleigh noise applied to the amplitudes of the power-law part of the spectrum.
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
///  * `add_noise` - Whether to multiply the amplitudes from index 6 up by Rayleigh variates (scaled to zero mean log) to simulate the natural variability of the PSD.
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

    // Power law with slope sn up to the data limit, continued beyond it with a slope of at least PSD_MIN_EXTRAPOLATION_SLOPE
    let xn = wavelength[6].ln();
    let x_data = (base / PSD_DATA_KMAX as f64).ln();
    let y_data = yn + sn * (x_data - xn);
    let s_ext = sn.max(PSD_MIN_EXTRAPOLATION_SLOPE);
    for i in 6..nfreq {
        let log_x = wavelength[i].ln();
        let log_amplitude = if i <= PSD_DATA_KMAX {
            yn + sn * (log_x - xn)
        } else {
            y_data + s_ext * (log_x - x_data)
        };
        amplitude[i] = log_amplitude.exp();
    }

    // Optional noise in the amplitudes. The Fourier amplitudes of a random profile are Rayleigh distributed, and the
    // residuals of the measured log amplitudes about each crater's power law match this (sd ≈ 0.61, skew ≈ -1.1, with no
    // dependence on feature, crater size, or wavenumber). Each amplitude is multiplied by a Rayleigh variate scaled so that
    // the mean of its log is zero, so the power law remains the mean of ln(amplitude). The first five components are
    // control points themselves, whose noise already carries their measured scatter, so the noise is only added to the
    // power-law part of the spectrum.
    if add_noise {
        let unit = Uniform::new(0.0, 1.0).expect("valid uniform distribution");
        for i in 6..nfreq {
            amplitude[i] *= rayleigh_log_centered(unit.sample(&mut rng));
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
        let amp = psd.log_scale * amplitude[k] * nprofile as f64 / 2.0;
        let complex_val = Complex::from_polar(amp, phase[k]);
        buffer[k] = complex_val;
        buffer[nprofile - k] = complex_val.conj();
    }

    // The Nyquist component is at index nfreq - 1
    if nprofile % 2 == 0 {
        buffer[nyquist_limit] = Complex::from_polar(
            psd.log_scale * amplitude[nyquist_limit] * nprofile as f64,
            phase[nyquist_limit],
        );
    }
    // The target spectrum, kept for imposing skewness
    let target_magnitude: Vec<f64> = buffer.iter().map(|val| val.norm()).collect();
    ifft.process(&mut buffer);
    let mut log_profile: Vec<f64> = buffer.iter().map(|val| val.re / nprofile as f64).collect();
    if psd.skewness != 0.0 && nprofile >= 4 {
        log_profile = impose_skewness(&log_profile, &target_magnitude, psd.skewness, &mut planner);
    }

    // Exponentiate the log profile, then rescale so the arithmetic mean is amplitude[0] in units of normalization_length
    let mut out: Vec<f64> = log_profile.iter().map(|val| val.exp()).collect();
    let exp_mean = out.iter().sum::<f64>() / nprofile as f64;
    let target_mean = amplitude[0] * psd.normalization_length;
    let scale_factor = target_mean / exp_mean;
    // Guard against extremes: no value further from the mean than the bounds allow
    let (vmin, vmax) = (target_mean * psd.log_ratio_min.exp(), target_mean * psd.log_ratio_max.exp());
    // max/min rather than clamp, which would panic on a NaN bound (e.g. 0 * inf for a zero mean)
    out.iter_mut().for_each(|v| *v = (*v * scale_factor).max(vmin).min(vmax));

    // Generate the corresponding bearing values for this profile
    let mut bearing: Vec<f64> = Vec::with_capacity(nprofile);
    for i in 0..nprofile {
        bearing.push(TAU * (i as f64 / nprofile as f64));
    }
    (bearing, out)
}

/// Number of iterations used to impose a target skewness on a log profile while keeping its spectrum
const SKEWNESS_ITERATIONS: usize = 20;

/// Sample skewness (biased, m3 / m2^1.5) of a set of values
fn sample_skewness(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let (m2, m3) = values.iter().fold((0.0, 0.0), |(m2, m3), &v| {
        let d = v - mean;
        (m2 + d * d, m3 + d * d * d)
    });
    let (m2, m3) = (m2 / n, m3 / n);
    if m2 > 0.0 { m3 / m2.powf(1.5) } else { 0.0 }
}

/// Curvature of the skewing transform on the compressed side of the median, as a fraction of its curvature on the stretched
/// side. With equal curvatures (1) the transform is (e^(λz) - 1)/λ, which compresses one side so strongly that profiles
/// dominated by a few long-wavelength components develop flat plateaus, which measured profiles do not have. With 0 the
/// compressed side is left unchanged, which gives profiles that are too peaked. 0.25 matches the measured profiles in both
/// respects (tested on the 852 measured profiles with 360 bearings, each with its own spectrum and skewness).
const SKEW_COMPRESSED_CURVATURE: f64 = 0.25;

/// Largest exponent used on the stretched side of the skewing transform, which keeps it from overflowing. Realistic
/// skewness targets never come close to it.
const SKEW_MAX_EXPONENT: f64 = 50.0;

/// Values of the skewing transform of z (standardized about its median).
///
/// On the side of the median being stretched (below it for λ < 0, which gives deep dips; above it for λ > 0, which gives
/// tall peaks) the transform is (e^(λz) - 1)/λ. On the other side it is the same function with curvature
/// SKEW_COMPRESSED_CURVATURE * λ, which compresses it more gently. The transform is increasing and continuous with a
/// continuous slope (1 at the median) for either sign of λ.
fn skew_transform(z: &[f64], lambda: f64) -> Vec<f64> {
    if lambda.abs() < 1e-12 {
        return z.to_vec();
    }
    let lambda_compressed = SKEW_COMPRESSED_CURVATURE * lambda;
    z.iter()
        .map(|&v| {
            let stretched = (lambda < 0.0 && v < 0.0) || (lambda > 0.0 && v > 0.0);
            if stretched {
                (lambda * v).min(SKEW_MAX_EXPONENT).exp_m1() / lambda
            } else if lambda_compressed.abs() < 1e-12 {
                v
            } else {
                (lambda_compressed * v).exp_m1() / lambda_compressed
            }
        })
        .collect()
}

/// Finds λ such that the skewness of the transformed values is the target, by bisection. The skewness increases
/// monotonically with λ. If the target is beyond what the values can reach, the closest reachable skewness is used.
fn solve_skew_lambda(z: &[f64], target: f64) -> f64 {
    const LAMBDA_MAX: f64 = 64.0;
    let skew_at = |lambda: f64| sample_skewness(&skew_transform(z, lambda));
    let (mut lo, mut hi) = if target > sample_skewness(z) { (0.0, 1.0) } else { (-1.0, 0.0) };
    while skew_at(hi) < target && hi < LAMBDA_MAX {
        lo = hi;
        hi *= 2.0;
    }
    while skew_at(lo) > target && lo > -LAMBDA_MAX {
        hi = lo;
        lo *= 2.0;
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if skew_at(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Imposes a target skewness on a zero-mean log profile while keeping its spectrum.
///
/// Profiles synthesized with random phases have symmetric distributions of values, but measured crater profiles are skewed
/// (rim heights have dips but not comparable peaks, crater radii have outward bulges). The values are passed through the
/// skewing transform (see `skew_transform`), with λ chosen to give the target skewness, which defines the target distribution.
/// The profile is then iterated between the target spectrum and the target distribution (the IAAFT method): impose the
/// target Fourier magnitudes while keeping the phases, then replace the values by the target values of the same rank. The
/// result has exactly the target distribution, and so the target skewness, and a spectrum within a few percent of the
/// target. The phase coupling this produces gives the measured relation between skewness and kurtosis.
///
/// # Arguments
/// * `log_profile` - The zero-mean log profile synthesized from the PSD
/// * `target_magnitude` - The magnitudes of the full (two-sided) unnormalized DFT of the log profile
/// * `skewness` - The target skewness
/// * `planner` - FFT planner
fn impose_skewness(
    log_profile: &[f64],
    target_magnitude: &[f64],
    skewness: f64,
    planner: &mut FftPlanner<f64>,
) -> Vec<f64> {
    let n = log_profile.len();
    let mean = log_profile.iter().sum::<f64>() / n as f64;
    let sd = (log_profile.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
    if sd == 0.0 {
        return log_profile.to_vec();
    }
    // Standardize about the median, which divides the profile into the stretched and compressed sides of the transform
    let mut sorted = log_profile.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if n % 2 == 0 { 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]) } else { sorted[n / 2] };
    let z: Vec<f64> = log_profile.iter().map(|v| (v - median) / sd).collect();

    // Target values: the transformed values, standardized and rescaled to the spread of the profile
    let lambda = solve_skew_lambda(&z, skewness);
    let mut target = skew_transform(&z, lambda);
    let tmean = target.iter().sum::<f64>() / n as f64;
    let tsd = (target.iter().map(|v| (v - tmean).powi(2)).sum::<f64>() / n as f64).sqrt();
    target.iter_mut().for_each(|v| *v = (*v - tmean) / tsd * sd);
    target.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);
    let mut order: Vec<usize> = (0..n).collect();
    let rank_remap = |values: &mut Vec<f64>, order: &mut Vec<usize>| {
        order.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap());
        for (rank, &i) in order.iter().enumerate() {
            values[i] = target[rank];
        }
    };

    let mut values = log_profile.to_vec();
    rank_remap(&mut values, &mut order);
    let mut buffer: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); n];
    for _ in 0..SKEWNESS_ITERATIONS {
        for (b, &v) in buffer.iter_mut().zip(values.iter()) {
            *b = Complex::new(v, 0.0);
        }
        fft.process(&mut buffer);
        for (b, &m) in buffer.iter_mut().zip(target_magnitude.iter()) {
            let norm = b.norm();
            *b = if norm > 0.0 { *b * (m / norm) } else { Complex::new(m, 0.0) };
        }
        ifft.process(&mut buffer);
        for (v, b) in values.iter_mut().zip(buffer.iter()) {
            *v = b.re / n as f64;
        }
        rank_remap(&mut values, &mut order);
    }
    values
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
    use rand_distr::Normal;

    /// The radial warp maps the floor edge and rim of a bearing onto those of the mean profile, is the identity away from them,
    /// and is strictly increasing, for small and large shifts in either direction.
    #[test]
    fn radial_warp_is_local_and_monotone() {
        let (mean_rf, mean_r, scale) = (300.0, 1000.0, 200.0);
        for (rf, r) in [(300.0, 1000.0), (250.0, 1150.0), (380.0, 820.0), (100.0, 1400.0), (600.0, 700.0)] {
            let w = RadialWarp::new(rf, mean_rf, r, mean_r, scale);
            assert!((w.eval(rf) - mean_rf).abs() < 1e-9, "floor edge {rf} maps to {}", w.eval(rf));
            assert!((w.eval(r) - mean_r).abs() < 1e-9, "rim {r} maps to {}", w.eval(r));
            let far = r + scale.max(2.0 * (r - mean_r).abs()) + 1.0;
            assert!((w.eval(far) - far).abs() < 1e-12);
            let mut prev = w.eval(0.0);
            for i in 1..=4000 {
                let x = 2.0 * mean_r * i as f64 / 4000.0;
                let y = w.eval(x);
                assert!(y > prev, "not increasing at {x} for rf {rf}, r {r}");
                prev = y;
            }
        }
    }

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
            skewness: 0.0,
            log_scale: 1.0,
            log_ratio_min: f64::NEG_INFINITY,
            log_ratio_max: f64::INFINITY,
            max_radius_fraction: f64::INFINITY,
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
    fn extrapolated_variance_converges_with_resolution() {
        // A shallow slope whose variance would diverge if the power law were extended unchanged
        let control_points: HashMap<String, f64> = [("sn", 0.2), ("yn", -3.0), ("y1", -3.0), ("y2", -3.0), ("y3", -3.0), ("y4", -3.0), ("y5", -3.0)]
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        let variance = |nprofile: usize| {
            let (_, amplitude, _) = get_1d_psd_from_control_points(&control_points, 1.0, nprofile, false, 0).unwrap();
            amplitude.iter().skip(1).map(|a| 0.5 * a * a).sum::<f64>()
        };
        let (v_data, v_fine, v_finer) = (variance(2 * PSD_DATA_KMAX), variance(20_000), variance(200_000));
        assert!(v_fine > v_data, "the tail beyond the data limit should add variance");
        assert!((v_finer / v_fine - 1.0).abs() < 0.01, "variance not converged: {v_fine} vs {v_finer}");
        // Unchanged up to the data limit
        let (_, a_data, _) = get_1d_psd_from_control_points(&control_points, 1.0, 2 * PSD_DATA_KMAX, false, 0).unwrap();
        let (_, a_fine, _) = get_1d_psd_from_control_points(&control_points, 1.0, 20_000, false, 0).unwrap();
        for k in 1..=PSD_DATA_KMAX {
            assert!((a_data[k] / a_fine[k] - 1.0).abs() < 1e-12, "k={k}");
        }
    }

    #[test]
    fn imposed_skewness_is_exact_and_keeps_spectrum() {
        let mut rng = ChaCha12Rng::seed_from_u64(11);
        let nfreq = 181;
        let wavelength = Array1::from_iter((0..nfreq).map(|i| TAU / i.max(1) as f64));
        let mut amplitude = Array1::from_iter((0..nfreq).map(|k| 0.1 * (k.max(1) as f64).powf(-0.8)));
        amplitude[0] = 0.2;
        let phase = Array1::from_iter((0..nfreq).map(|_| rng.random_range(0.0..TAU)));
        let nprofile = 2 * (nfreq - 1);
        for &target in &[-2.0, -0.5, 0.7] {
            let psd = PSD1DView {
                nprofile,
                nfreq,
                normalization_length: 1000.0,
                mean: 0.2,
                pix: 1.0,
                skewness: target,
                log_scale: 1.0,
                log_ratio_min: f64::NEG_INFINITY,
                log_ratio_max: f64::INFINITY,
                max_radius_fraction: f64::INFINITY,
                wavelength: wavelength.view(),
                amplitude: amplitude.view(),
                phase: phase.view(),
            };
            let (_, y) = profile_from_psd(&psd);
            let ymean = y.iter().sum::<f64>() / y.len() as f64;
            assert!((ymean / 200.0 - 1.0).abs() < 1e-12, "mean not preserved");
            let s: Vec<f64> = y.iter().map(|v| v.ln()).collect();
            assert!((sample_skewness(&s) - target).abs() < 1e-6, "skewness {} != {target}", sample_skewness(&s));
            // The spectrum of the log profile stays close to the target
            let mut buffer: Vec<Complex<f64>> = s.iter().map(|&v| Complex::new(v, 0.0)).collect();
            FftPlanner::new().plan_fft_forward(nprofile).process(&mut buffer);
            let (mut err, mut norm) = (0.0, 0.0);
            for k in 1..nfreq - 1 {
                let a = 2.0 * buffer[k].norm() / nprofile as f64;
                err += (a - amplitude[k]).powi(2);
                norm += amplitude[k].powi(2);
            }
            let rel = (err / norm).sqrt();
            assert!(rel < 0.05, "spectrum error {rel} for skewness {target}");
        }
    }

    #[test]
    fn profile_respects_ratio_bounds_and_log_scale() {
        let mut rng = ChaCha12Rng::seed_from_u64(5);
        let nfreq = 181;
        let wavelength = Array1::from_iter((0..nfreq).map(|i| TAU / i.max(1) as f64));
        let mut amplitude = Array1::from_iter((0..nfreq).map(|k| 0.5 * (k.max(1) as f64).powf(-0.3)));
        amplitude[0] = 0.1;
        let phase = Array1::from_iter((0..nfreq).map(|_| rng.random_range(0.0..TAU)));
        let make = |log_scale: f64, lo: f64, hi: f64| {
            let psd = PSD1DView {
                nprofile: 2 * (nfreq - 1),
                nfreq,
                normalization_length: 1000.0,
                mean: 0.1,
                pix: 1.0,
                skewness: 0.0,
                log_scale,
                log_ratio_min: lo,
                log_ratio_max: hi,
                max_radius_fraction: f64::INFINITY,
                wavelength: wavelength.view(),
                amplitude: amplitude.view(),
                phase: phase.view(),
            };
            profile_from_psd(&psd).1
        };
        let sd = |y: &[f64]| {
            let s: Vec<f64> = y.iter().map(|v| v.ln()).collect();
            let m = s.iter().sum::<f64>() / s.len() as f64;
            (s.iter().map(|v| (v - m).powi(2)).sum::<f64>() / s.len() as f64).sqrt()
        };
        // Scaling the log profile scales its spread
        let (y1, y2) = (make(1.0, f64::NEG_INFINITY, f64::INFINITY), make(0.5, f64::NEG_INFINITY, f64::INFINITY));
        assert!((sd(&y2) / sd(&y1) - 0.5).abs() < 1e-9);
        // Bounds hold, relative to the target mean of 100
        let y = make(1.0, -0.5, 0.4);
        assert!(y.iter().all(|&v| v >= 100.0 * (-0.5f64).exp() * (1.0 - 1e-12) && v <= 100.0 * 0.4f64.exp() * (1.0 + 1e-12)));
        assert!(y1.iter().any(|&v| v < 100.0 * (-0.5f64).exp()), "test should exercise the lower bound");
    }

    #[test]
    fn rayleigh_noise_has_zero_mean_log() {
        let mut rng = ChaCha12Rng::seed_from_u64(3);
        let unit = Uniform::new(0.0, 1.0).unwrap();
        let n = 1_000_000;
        let logs: Vec<f64> = (0..n).map(|_| rayleigh_log_centered(unit.sample(&mut rng)).ln()).collect();
        let mean = logs.iter().sum::<f64>() / n as f64;
        let sd = (logs.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
        assert!(mean.abs() < 3e-3, "mean of log = {mean}");
        assert!((sd - std::f64::consts::PI / 24f64.sqrt()).abs() < 3e-3, "sd of log = {sd}");
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
            skewness: 0.0,
            log_scale: 1.0,
            log_ratio_min: f64::NEG_INFINITY,
            log_ratio_max: f64::INFINITY,
            max_radius_fraction: f64::INFINITY,
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
