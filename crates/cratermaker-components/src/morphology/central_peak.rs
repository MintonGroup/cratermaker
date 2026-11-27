//! Central peaks and peak rings: a patchy population of tall columns relaxed by slope collapse.
//!
//! Turbulence noise evaluated at each face of the crater (in the plane about the crater center) is used as an indicator: where
//! it is more than [`PEAK_THRESHOLD`] standard deviations above its mean, the face is raised to the height of a Gaussian
//! envelope (centered, or offset as a ring), scaled up by [`PEAK_TOP_VARIATION`] times the excess of the noise over the
//! threshold so that the columns have uneven tops, and elsewhere it is not raised. Slope collapse then relaxes the columns
//! into patchy mounds. The critical angle of the collapse varies from face to face around [`PEAK_CRITICAL_ANGLE`] with a
//! second, independent noise field ([`PEAK_ANGLE_SPREAD`], [`PEAK_ANGLE_NOISE_WAVELENGTH`]), so that the flanks are not
//! planar. The amplitude of the envelope is iterated so that the highest point of the relaxed
//! peak has the peak's height. The relaxation is done on a flat copy of the region, so that it changes only the peak and not
//! the crater around it, and the volume of the peak is removed from the crater floor as a layer that thins to zero at the edge
//! of the floor (in proportion to 1 - (r / floor radius)^2, so that it leaves no step there), so that crater formation still
//! conserves volume.
//!
//! The constants are those that reproduce the typical central peak of the fresh lunar craters of 35-140 km measured at the
//! native SLDEM2015 spacing (59 m): an envelope of standard deviation 0.08 R, columns where the noise is 0.5 standard
//! deviations above its mean, a noise wavelength of 1.5 envelope widths, and a critical angle of 34 degrees on average.

use crate::ArrayResult;
use crate::morphology::form_crater::CraterRegionView;
use crate::surface::{slope_collapse_field, turbulence_noise};
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
/// Variability of the column heights: a column is raised to the envelope times 1 + PEAK_TOP_VARIATION x (noise - threshold),
/// with the noise in standard deviations.
pub const PEAK_TOP_VARIATION: f64 = 1.0;
/// Mean critical slope angle of the collapse of the columns (degrees).
pub const PEAK_CRITICAL_ANGLE: f64 = 34.0;
/// Standard deviation of the critical slope angle (degrees), which varies with an independent noise field (clipped at
/// [`PEAK_ANGLE_CLIP`] standard deviations).
pub const PEAK_ANGLE_SPREAD: f64 = 6.0;
const PEAK_ANGLE_CLIP: f64 = 2.5;
/// Wavelength of the largest octave of the critical angle noise, as a multiple of the envelope width.
pub const PEAK_ANGLE_NOISE_WAVELENGTH: f64 = 0.5;
/// Extent of the columns from the envelope's center or ring, in envelope widths.
const PEAK_SUPPORT: f64 = 2.5;
/// Relative tolerance on the height of the relaxed peak, and the most collapses used to reach it.
const HEIGHT_TOLERANCE: f64 = 0.01;
const MAX_ITERATIONS: usize = 10;
/// The largest factor by which one step of the iteration changes the amplitude.
const MAX_AMPLITUDE_STEP: f64 = 4.0;
/// The largest amplitude of the columns, as a multiple of the peak's height. A peak that needs taller columns is narrower than
/// the mesh can resolve: its columns stand on a few faces, and slope collapse spreads taller columns further, so that the
/// relaxed peak does not grow with them. (Without a bound, the iteration then drives the amplitude up without limit, and the
/// slope collapse of the columns, whose number of diffusion steps grows with the square of their slopes, does not finish.)
const MAX_AMPLITUDE: f64 = 10.0;

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
/// The face elevation change: the relaxed peak, less a layer of the same volume over the floor that thins to zero at its edge.
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

    // Two independent noise fields in the plane about the crater center, normalized over the support: one for the columns, and
    // one for the critical angle of their collapse
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let Some(noise) = normalized_noise(region, &support, PEAK_NOISE_WAVELENGTH * peak.width, &mut rng)? else {
        return Ok(change);
    };
    let Some(angle_noise) = normalized_noise(region, &support, PEAK_ANGLE_NOISE_WAVELENGTH * peak.width, &mut rng)? else {
        return Ok(change);
    };

    // The columns, in units of the envelope's amplitude
    let mut shape = Array1::<f64>::zeros(nf);
    for (i, &f) in support.iter().enumerate() {
        if noise[i] > PEAK_THRESHOLD {
            let u = (region.face_distance[f] - peak.ring_radius) / peak.width;
            shape[f] = (-0.5 * u * u).exp() * (1.0 + PEAK_TOP_VARIATION * (noise[i] - PEAK_THRESHOLD));
        }
    }
    if !shape.iter().any(|&s| s > 0.0) {
        return Ok(change);
    }

    // Collapse the columns on a flat region, iterating the amplitude (secant steps in log amplitude and log height) until the
    // relaxed peak has the required height
    // The critical slope of each face: the mean angle outside the support, varying with the angle noise inside it
    let mut critical_slope = Array1::<f64>::from_elem(nf, PEAK_CRITICAL_ANGLE.to_radians().tan());
    for (i, &f) in support.iter().enumerate() {
        let z = angle_noise[i].clamp(-PEAK_ANGLE_CLIP, PEAK_ANGLE_CLIP);
        critical_slope[f] = (PEAK_CRITICAL_ANGLE + PEAK_ANGLE_SPREAD * z).to_radians().tan();
    }
    let zero_nodes = Array1::<f64>::zeros(region.n_node());
    let relax = |amplitude: f64| -> ArrayResult {
        let columns = &shape * amplitude;
        let delta = {
            let view = region.kernel_view(columns.view(), zero_nodes.view());
            slope_collapse_field(critical_slope.view(), &view)?
        };
        Ok(columns + delta)
    };
    // The steps are limited to a factor MAX_AMPLITUDE_STEP and the amplitude to MAX_AMPLITUDE times the height. If the height
    // stops growing with the amplitude, the mesh cannot hold a taller peak, and the iteration stops. The relaxed peak closest
    // to the required height is kept.
    let ln_target = peak.height.ln();
    let ln_max_amplitude = (MAX_AMPLITUDE * peak.height).ln();
    let ln_max_step = MAX_AMPLITUDE_STEP.ln();
    let mut amplitude = 2.0 * peak.height;
    let mut relaxed = relax(amplitude)?;
    let mut history: Vec<(f64, f64)> = Vec::with_capacity(MAX_ITERATIONS);
    let mut best: Option<(f64, Array1<f64>)> = None;
    for _ in 0..MAX_ITERATIONS {
        let height = relaxed.fold(f64::NEG_INFINITY, |m, &h| m.max(h));
        if !(height > 0.0) {
            break;
        }
        let (la, lh) = (amplitude.ln(), height.ln());
        let miss = (lh - ln_target).abs();
        if best.as_ref().is_none_or(|(m, _)| miss < *m) {
            best = Some((miss, relaxed.clone()));
        }
        if (height / peak.height - 1.0).abs() < HEIGHT_TOLERANCE {
            break;
        }
        let gain = match history.last() {
            Some(&(a0, h0)) if (la - a0).abs() > 1e-12 => (lh - h0) / (la - a0),
            _ => 1.0,
        };
        history.push((la, lh));
        if !(gain > 0.0) {
            break; // the height no longer grows with the amplitude
        }
        let step = ((ln_target - lh) / gain).clamp(-ln_max_step, ln_max_step);
        let next = (la + step).min(ln_max_amplitude);
        if (next - la).abs() < 1e-12 {
            break; // at the largest amplitude
        }
        amplitude = next.exp();
        relaxed = relax(amplitude)?;
    }
    let Some((_, relaxed)) = best else {
        return Ok(change);
    };

    // Remove the peak's volume from the floor as a layer that thins to zero at its edge
    let volume: f64 = relaxed.iter().zip(region.face_area.iter()).map(|(&h, &a)| h * a).sum();
    let taper: Vec<(usize, f64)> = (0..nf)
        .filter_map(|f| {
            let u = region.face_distance[f] / peak.floor_radius;
            (u < 1.0).then(|| (f, 1.0 - u * u))
        })
        .collect();
    let taper_volume: f64 = taper.iter().map(|&(f, w)| w * region.face_area[f]).sum();
    change.assign(&relaxed);
    if taper_volume > 0.0 {
        let layer = volume / taper_volume;
        for &(f, w) in &taper {
            change[f] -= layer * w;
        }
    }
    Ok(change)
}

/// Turbulence noise at the given faces (in the plane about the crater center), normalized to zero mean and unit standard
/// deviation over them, or `None` if it does not vary.
fn normalized_noise(
    region: &CraterRegionView<'_>,
    faces: &[usize],
    wavelength: f64,
    rng: &mut ChaCha8Rng,
) -> Result<Option<Array1<f64>>, String> {
    let anchor = Array2::from_shape_fn((PEAK_NOISE_OCTAVES, 3), |_| rng.random_range(0.0..std::f64::consts::TAU));
    let x = Array1::from_iter(faces.iter().map(|&f| region.face_proj_x[f]));
    let y = Array1::from_iter(faces.iter().map(|&f| region.face_proj_y[f]));
    let z = Array1::<f64>::zeros(faces.len());
    let noise = turbulence_noise(
        x.view(),
        y.view(),
        z.view(),
        1.0,
        wavelength,
        2.0,
        PEAK_NOISE_PERSISTENCE,
        anchor.view(),
        rng.random::<u32>(),
    )?;
    let n = noise.len() as f64;
    let mean = noise.sum() / n;
    let sd = (noise.mapv(|v| (v - mean).powi(2)).sum() / n).sqrt();
    if !(sd > 0.0) {
        return Ok(None);
    }
    Ok(Some(noise.mapv(|v| (v - mean) / sd)))
}
