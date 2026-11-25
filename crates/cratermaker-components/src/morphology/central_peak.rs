//! Central peaks and peak rings: a patchy population of tall columns relaxed by slope collapse.
//!
//! Turbulence noise evaluated at each face of the crater (in the plane about the crater center) is used as an indicator: where
//! it is more than [`PEAK_THRESHOLD`] standard deviations above its mean, the face is raised to the height of a Gaussian
//! envelope (centered, or offset as a ring), and elsewhere it is not raised. Slope collapse at [`PEAK_CRITICAL_ANGLE`] then
//! relaxes the columns into patchy mounds. The amplitude of the envelope is iterated so that the highest point of the relaxed
//! peak has the peak's height. The relaxation is done on a flat copy of the region, so that it changes only the peak and not
//! the crater around it, and the volume of the peak is removed from the crater floor as a uniform layer, so that crater
//! formation still conserves volume.
//!
//! The constants are those that reproduce the typical central peak of the fresh lunar craters of 35-140 km measured at the
//! native SLDEM2015 spacing (59 m): an envelope of standard deviation 0.08 R, columns where the noise is 0.5 standard
//! deviations above its mean, a noise wavelength of 1.5 envelope widths, and a critical angle of 34 degrees.

use crate::ArrayResult;
use crate::morphology::form_crater::CraterRegionView;
use crate::surface::{slope_collapse, turbulence_noise};
use numpy::ndarray::prelude::*;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

/// Noise level, in standard deviations above its mean, above which a face is a column.
pub const PEAK_THRESHOLD: f64 = 0.5;
/// Wavelength of the largest noise octave, as a multiple of the envelope width.
pub const PEAK_NOISE_WAVELENGTH: f64 = 1.5;
/// Number of noise octaves.
pub const PEAK_NOISE_OCTAVES: usize = 3;
/// Amplitude ratio between successive noise octaves.
pub const PEAK_NOISE_PERSISTENCE: f64 = 0.5;
/// Critical slope angle of the collapse of the columns (degrees).
pub const PEAK_CRITICAL_ANGLE: f64 = 34.0;
/// Extent of the columns from the envelope's center or ring, in envelope widths.
const PEAK_SUPPORT: f64 = 2.5;
/// Relative tolerance on the height of the relaxed peak, and the most collapses used to reach it.
const HEIGHT_TOLERANCE: f64 = 0.01;
const MAX_ITERATIONS: usize = 10;

/// The parameters of a central peak or peak ring.
#[derive(Clone, Copy, Debug)]
pub struct CentralPeak {
    /// Height of the highest point of the relaxed peak above the floor (m).
    pub height: f64,
    /// Standard deviation of the Gaussian envelope (m).
    pub width: f64,
    /// Radius of the envelope's ring (m); 0 for a central peak.
    pub ring_radius: f64,
    /// Radius of the crater floor (m), from which the peak's volume is removed.
    pub floor_radius: f64,
}

impl CentralPeak {
    /// Whether the peak has any relief.
    pub fn is_present(&self) -> bool {
        self.height > 0.0 && self.width > 0.0
    }
}

/// The face elevation change of a central peak or peak ring on a crater region.
///
/// # Arguments
///
/// * `region` - The geometry of the crater's region; distances and planar coordinates are about the crater center.
/// * `peak` - The peak's parameters.
/// * `seed` - Seed of the noise.
///
/// # Returns
///
/// The face elevation change: the relaxed peak, less a uniform layer over the floor of the same volume.
pub fn central_peak(region: &CraterRegionView<'_>, peak: &CentralPeak, seed: u64) -> ArrayResult {
    let nf = region.n_face();
    let mut change = Array1::<f64>::zeros(nf);
    if !peak.is_present() {
        return Ok(change);
    }

    // The faces the columns can occupy
    let support: Vec<usize> = (0..nf)
        .filter(|&f| (region.face_distance[f] - peak.ring_radius).abs() <= PEAK_SUPPORT * peak.width)
        .collect();
    if support.is_empty() {
        return Ok(change);
    }

    // Noise in the plane about the crater center, normalized over the support
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let anchor = Array2::from_shape_fn((PEAK_NOISE_OCTAVES, 3), |_| {
        rng.random_range(0.0..std::f64::consts::TAU)
    });
    let x = Array1::from_iter(support.iter().map(|&f| region.face_proj_x[f]));
    let y = Array1::from_iter(support.iter().map(|&f| region.face_proj_y[f]));
    let z = Array1::<f64>::zeros(support.len());
    let noise = turbulence_noise(
        x.view(),
        y.view(),
        z.view(),
        1.0,
        PEAK_NOISE_WAVELENGTH * peak.width,
        2.0,
        PEAK_NOISE_PERSISTENCE,
        anchor.view(),
        rng.random::<u32>(),
    )?;
    let n = noise.len() as f64;
    let mean = noise.sum() / n;
    let sd = (noise.mapv(|v| (v - mean).powi(2)).sum() / n).sqrt();
    if !(sd > 0.0) {
        return Ok(change);
    }

    // The columns, in units of the envelope's amplitude
    let mut shape = Array1::<f64>::zeros(nf);
    for (i, &f) in support.iter().enumerate() {
        if (noise[i] - mean) / sd > PEAK_THRESHOLD {
            let u = (region.face_distance[f] - peak.ring_radius) / peak.width;
            shape[f] = (-0.5 * u * u).exp();
        }
    }
    if !shape.iter().any(|&s| s > 0.0) {
        return Ok(change);
    }

    // Collapse the columns on a flat region, iterating the amplitude (secant steps in log amplitude and log height) until the
    // relaxed peak has the required height
    let critical_slope = PEAK_CRITICAL_ANGLE.to_radians().tan();
    let zero_nodes = Array1::<f64>::zeros(region.n_node());
    let relax = |amplitude: f64| -> ArrayResult {
        let columns = &shape * amplitude;
        let delta = {
            let view = region.kernel_view(columns.view(), zero_nodes.view());
            slope_collapse(critical_slope, &view)?
        };
        Ok(columns + delta)
    };
    let mut amplitude = 2.0 * peak.height;
    let mut relaxed = relax(amplitude)?;
    let mut history: Vec<(f64, f64)> = Vec::with_capacity(MAX_ITERATIONS);
    for _ in 0..MAX_ITERATIONS {
        let height = relaxed.fold(f64::NEG_INFINITY, |m, &h| m.max(h));
        if !(height > 0.0) {
            return Ok(change);
        }
        history.push((amplitude.ln(), height.ln()));
        if (height / peak.height - 1.0).abs() < HEIGHT_TOLERANCE {
            break;
        }
        let gain = match history.as_slice() {
            [.., (a0, h0), (a1, h1)] if (a1 - a0).abs() > 1e-12 => ((h1 - h0) / (a1 - a0)).max(0.1),
            _ => 1.0,
        };
        let (la, lh) = *history.last().expect("history is not empty");
        amplitude = (la + (peak.height.ln() - lh) / gain).exp();
        relaxed = relax(amplitude)?;
    }

    // Remove the peak's volume from the floor as a uniform layer
    let volume: f64 = relaxed.iter().zip(region.face_area.iter()).map(|(&h, &a)| h * a).sum();
    let floor: Vec<usize> = (0..nf)
        .filter(|&f| region.face_distance[f] <= peak.floor_radius)
        .collect();
    let floor_area: f64 = floor.iter().map(|&f| region.face_area[f]).sum();
    change.assign(&relaxed);
    if floor_area > 0.0 {
        let layer = volume / floor_area;
        for &f in &floor {
            change[f] -= layer;
        }
    }
    Ok(change)
}
