use crate::{ArrayResult,ArrayResult2D};
use std::collections::HashMap;
use rand::prelude::*;
use rand::SeedableRng;
use rand_distr::{Normal,Uniform};
use rand_chacha::ChaCha12Rng;
use numpy::ndarray::prelude::*;
use pyo3::FromPyObject;
use crate::morphology::basicmoon::{crater_profile_function, ejecta_profile_function};

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject, Clone, Debug)]
pub struct RealMoonCrater {
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
    pub rim_flank_radius: f64,
    pub rimdrop: f64,
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_offset: f64,
    pub rim_radius_control: HashMap<String, f64>,
    pub rim_flank_radius_control: HashMap<String, f64>,
    pub rim_elevation_control: HashMap<String, f64>,
    pub floor_elevation_control: HashMap<String, f64>,
    pub wall_texture_control: HashMap<String, f64>,
    pub ejecta_texture_control: HashMap<String, f64>,
    pub floor_texture_control: HashMap<String, f64>,
}

// Creates a profile of the crater
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
pub fn realmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &RealMoonCrater,
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
///
///
/// Computes a target 1D power spectral density distribution based on control points defining a piecewise linear function in log-log space, with optional Gaussian noise added to the log power values.
///
/// # Arguments
/// * `control_points` - A dictionary containing the control points for the piecewise linear function. The expected keys are:
/// - "Slope_12": The slope of the first segment (largest wavelengths).
/// - "Breakpoint_2_x": The x-coordinate of the breakpoint between the first and second segments (in log10(wavelength)).
/// - "Breakpoint_2_y": The y-coordinate of the breakpoint between the first and second
///   segments (in log10(power)).
/// - "Breakpoint_3_y": The y-coordinate of the breakpoint between the second and third
///   segments (in log10(power)).
/// - "Breakpoint_4_y": The y-coordinate of the breakpoint for the smallest wavelengths (in log10(power)).
///  * `npoints` - The number of points in the output PSD, which determines the wavelength resolution and the Nyquist frequency.
///  * `add_noise` - Whether to add Gaussian noise to the log10(power) values to simulate natural variability in the PSD.
///  * `seed` - The random seed for reproducibility of the noise if `add_noise` is true.
///
pub fn get_1d_psd_from_control_points(
    control_points: &HashMap<String, f64>,
    npoints: usize,
    add_noise: bool,
    rng_seed: u64,
) -> ArrayResult2D {
    let slope_12 = control_points["Slope_12"];
    let bp2_x = control_points["Breakpoint_2_x"];
    let bp2_y = control_points["Breakpoint_2_y"];
    let bp3_y = control_points["Breakpoint_3_y"];
    let bp4_y = control_points["Breakpoint_4_y"];

    let bp4_x = (2.0 * std::f64::consts::PI).log10();
    let bp3_x = (10_f64.powf(bp4_x) / 2.0).log10();

    // Same spacing logic as Python: interval = 10**bp4_x / npoints
    let interval = 10_f64.powf(bp4_x) / npoints as f64;

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

    // Find bp2_x_index (matching Python loop behavior)
    let threshold = 10_f64.powf(bp2_x);
    let mut bp2_x_index = nrows.saturating_sub(1);
    for i in 0..nrows {
        if psd[[i, 0]] < threshold {
            bp2_x_index = i.saturating_sub(1);
            break;
        }
    }

    // Piecewise lines in log10-log10 space
    let k_23 = (bp3_y - bp2_y) / (bp3_x - bp2_x);
    let b_23 = bp3_y - k_23 * bp3_x;
    let k_12 = slope_12;
    let b_12 = bp2_y - k_12 * bp2_x;

    if nrows > 0 {
        psd[[0, 1]] = 10_f64.powf(bp4_y);
    }
    if nrows > 1 {
        psd[[1, 1]] = 10_f64.powf(bp3_y);
    }

    // psd[2 : bp2_x_index + 1, 1]
    if bp2_x_index >= 2 {
        for i in 2..=bp2_x_index {
            let log_w = psd[[i, 0]].log10();
            psd[[i, 1]] = 10_f64.powf(k_23 * log_w + b_23);
        }
    }

    // psd[bp2_x_index + 1 :, 1]
    for i in (bp2_x_index + 1)..nrows {
        let log_w = psd[[i, 0]].log10();
        psd[[i, 1]] = 10_f64.powf(k_12 * log_w + b_12);
    }

    // flipud
    let mut flipped = Array2::<f64>::zeros((nrows, 2));
    for i in 0..nrows {
        flipped.row_mut(i).assign(&psd.row(nrows - 1 - i));
    }

    // Optional Gaussian noise in log10 power
    if add_noise {
        let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
        let normal = Normal::new(0.0, 0.55).expect("valid normal distribution");
        for i in 0..nrows {
            let log_power = flipped[[i, 1]].log10();
            let noisy_log_power = log_power + normal.sample(&mut rng);
            flipped[[i, 1]] = 10_f64.powf(noisy_log_power);
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
/// * `bearings
///    - A 1D array of angular bearings (in radians) at which to compute the profile, typically ranging from 0 tto 2π.
/// * `phases` - An optional 1D array of phase values (in radians) corresponding to each frequency in the PSD. If not provided, random phases will be generated.
/// * `rng_seed` - The random seed for reproducibility when generating random phases if `phases` is not provided.
/// # Returns
/// * A 1D array of values corresponding to the input bearings, representing the linear profile generated from the PSD and phase information.
///
pub fn profile_from_psd(
    crater_radius: f64,
    ymean: f64,
    psd: ArrayView2<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    phases: Option<ArrayView1<'_, f64>>,
    rng_seed: u64,
) -> ArrayResult {
    let period_total: f64 = psd[[psd.nrows() - 1, 0]];
    let nfreq: usize = psd.nrows();
    let nbearings: usize = bearings.len();

    // Generate or use provided phases
    let phase_values: Array1<f64> = if let Some(p) = phases {
        p.to_owned()
    } else {
        let mut rng = ChaCha12Rng::seed_from_u64(rng_seed);
        let uniform = Uniform::new(1.0, period_total).expect("valid uniform distribution");
        Array1::from_iter((0..nfreq).map(|_| uniform.sample(&mut rng)))
    };

    // Compute amplitudes: sqrt(psd[:, 1] * period_total / (nfreq^2))
    let amplitude: Array1<f64> = psd
        .column(1)
        .mapv(|p| (p * period_total / (nbearings as f64 * nbearings as f64)).sqrt());

    // Compute y_ind: amplitude[i] * sin(2π * (1/psd[i,0]) * (theta + phase[i]))
    let mut delta_y = Array1::<f64>::zeros(nbearings);

    for i in 0..nfreq {
        let wavelength = psd[[i, 0]];
        let freq = 1.0 / wavelength;
        let phase = phase_values[i];
        let amp = amplitude[i];

        for j in 0..nbearings {
            let theta = bearings[j];
            let y = amp * (2.0 * std::f64::consts::PI * freq * (theta + phase)).sin();
            delta_y[j] += y;
        }
    }

    // Scale by crater_radius and add mean
    Ok(Array1::from_iter(delta_y.iter().map(|dy| dy * crater_radius + ymean)))
}