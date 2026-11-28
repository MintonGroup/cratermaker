//! Central peaks and peak rings: a patchy population of columns with flanks at the critical slope.
//!
//! Turbulence noise evaluated at each face of the crater (in the plane about the crater center) is used as an indicator: where
//! it is more than `threshold` standard deviations above its mean, the face is a column, raised to the height of a Gaussian
//! envelope (centered, or offset as a ring) times (noise - threshold)^`ramp_exponent`, and elsewhere it is not raised. The
//! columns so rise from zero at the threshold, and the massifs grow out of the floor rather than standing on it as mesas with
//! scarps as tall as the envelope. Each column is then given flanks at the critical slope:
//! the surface is raised to the upper envelope of cones falling from the columns at the critical slope of each face it falls
//! through, so that every massif, large or small, keeps a sharp summit and flanks at the angle of repose. The critical angle
//! varies from face to face around `critical_angle` with a second, independent noise field (`angle_spread`,
//! `angle_noise_wavelength`), so that the flanks are not planar. The highest point of the peak is the highest column, so the
//! amplitude of the envelope follows directly from the peak's height. The peak is built on a flat copy of the region, so that it
//! changes only the peak and not the crater around it, and its volume is removed from the crater floor as a layer that thins to
//! zero at the edge of the floor (in proportion to (1 - (r / floor radius)^2)^2, so that it leaves neither a step nor a break in
//! slope there), so that crater
//! formation still conserves volume.
//!
//! The peak is offset from the crater center (by the crater's peak center distance and bearing, or, if none is given, by an
//! offset drawn from its seed) and elongated along a random direction drawn from its seed: the envelope is
//! evaluated in coordinates centered on the offset point and stretched along the elongation direction (by a factor whose
//! logarithm is `elongation_sigma` times a half-normal deviate, preserving area). The noise of the columns and of their
//! critical angle is not stretched, so that the peak is elongated as a whole while its texture stays rugged in every
//! direction. Peak rings are not elongated (the elongation fades out as the peak opens into a ring), but their center is
//! offset. The offset is a two-dimensional Gaussian of standard deviation `offset_sigma` envelope widths in each
//! direction, unrelated to the impact direction: the central peaks of oblique impact craters on Venus are offset randomly, by
//! 0.031 crater diameters on average, indistinguishably from those of near-vertical impacts (Ekholm and Melosh 2001), as the
//! peaks of fresh lunar complex craters are (0.033 D on average; Lunar_Morphometry_Revisited, peak_measure.py).
//!
//! The parameters of the model ([`PeakModel`]) are calibrated against the central peaks of fresh lunar craters measured at the
//! native SLDEM2015 spacing, and the peak rings of Schrodinger and Orientale (Lunar_Morphometry_Revisited,
//! peak_model_calibrate.py).

use crate::ArrayResult;
use crate::morphology::form_crater::CraterRegionView;
use crate::morphology::smootherstep;
use crate::surface::{slope_limited_envelope, turbulence_noise};
use numpy::ndarray::prelude::*;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

/// The clipping of the critical angle noise, in standard deviations.
const PEAK_ANGLE_CLIP: f64 = 2.5;

/// The parameters of the central peak model (see the module documentation).
#[derive(Clone, Copy, Debug)]
pub struct PeakModel {
    /// Noise level, in standard deviations above its mean, above which a face is a column.
    pub threshold: f64,
    /// Wavelength of the largest octave of the column noise, as a multiple of the envelope width.
    pub noise_wavelength: f64,
    /// Number of noise octaves.
    pub octaves: usize,
    /// Amplitude ratio between successive noise octaves.
    pub persistence: f64,
    /// Exponent of the rise of the column heights with the excess of the noise over the threshold: a column is raised to the
    /// envelope times (noise - threshold)^ramp_exponent, with the noise in standard deviations.
    pub ramp_exponent: f64,
    /// Mean critical slope angle of the flanks (degrees).
    pub critical_angle: f64,
    /// Standard deviation of the critical slope angle (degrees), which varies with an independent noise field.
    pub angle_spread: f64,
    /// Wavelength of the largest octave of the critical angle noise, as a multiple of the envelope width.
    pub angle_noise_wavelength: f64,
    /// Extent of the columns from the envelope's center or ring, in envelope widths.
    pub support: f64,
    /// Standard deviation of each component of the offset of the envelope's center from the crater center, in envelope widths.
    pub offset_sigma: f64,
    /// Standard deviation of the half-normal distribution of the logarithm of the elongation (the ratio of the long to the
    /// short axis) of the envelope (the noise is not stretched); 0 for round peaks.
    pub elongation_sigma: f64,
}

impl Default for PeakModel {
    /// The parameters calibrated against the central peaks of the 14 fresh complex craters larger than 39 km measured at the
    /// native SLDEM2015 spacing (their volume, footprint, slopes, flank angle, massifs, elongation, solidity, offset and texture)
    /// and the peak rings of Schrodinger and Orientale (Lunar_Morphometry_Revisited, peak_model_calibrate.py). The columns extend
    /// to where the envelope falls below 1e-4 of its amplitude, so that the peak tails off into its surroundings.
    fn default() -> Self {
        Self {
            threshold: -0.68,
            noise_wavelength: 1.39,
            octaves: 3,
            persistence: 0.5,
            ramp_exponent: 0.70,
            critical_angle: 29.2,
            angle_spread: 9.4,
            angle_noise_wavelength: 0.5,
            support: 4.3,
            offset_sigma: 0.65,
            elongation_sigma: 0.65,
        }
    }
}

/// The parameters of a central peak or peak ring.
#[derive(Clone, Copy, Debug)]
pub struct CentralPeak {
    /// Height of the highest point of the peak above the floor (m).
    pub height: f64,
    /// Standard deviation of the Gaussian envelope (m).
    pub width: f64,
    /// Radius of the envelope's ring (m); 0 for a central peak.
    pub ring_radius: f64,
    /// Radius of the crater floor (m), from which the peak's volume is removed.
    pub floor_radius: f64,
    /// The offset of the envelope's center from the crater center: its distance (m) and bearing (degrees, clockwise from
    /// north). None to draw it from the seed (a two-dimensional Gaussian of standard deviation `PeakModel::offset_sigma`
    /// envelope widths in each direction).
    pub center: Option<(f64, f64)>,
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
/// * `model` - The parameters of the peak model.
/// * `seed` - Seed of the noise.
///
/// # Returns
///
/// The face elevation change: the peak, less a layer of the same volume over the floor that thins to zero at its edge.
pub fn central_peak(region: &CraterRegionView<'_>, peak: &CentralPeak, model: &PeakModel, seed: u64) -> ArrayResult {
    let nf = region.n_face();
    let mut change = Array1::<f64>::zeros(nf);
    if !peak.is_present() {
        return Ok(change);
    }
    let lowest_angle = model.critical_angle - PEAK_ANGLE_CLIP * model.angle_spread;
    if !(lowest_angle > 0.0 && model.critical_angle < 90.0) {
        return Err(format!(
            "the critical angle must stay between 0 and 90 degrees: {} +/- {} x {} degrees",
            model.critical_angle, PEAK_ANGLE_CLIP, model.angle_spread
        ));
    }

    // The offset of the peak and its elongation (see the module documentation), and the coordinates of each face in the peak's
    // frame: centered on the offset point, rotated to the elongation direction, and stretched along it
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let (dx, dy) = match peak.center {
        Some((distance, bearing)) => {
            let (sb, cb) = bearing.to_radians().sin_cos();
            (distance * sb, distance * cb)
        }
        None => (model.offset_sigma * peak.width * normal(&mut rng), model.offset_sigma * peak.width * normal(&mut rng)),
    };
    // Peak rings stay circular: the elongation fades out as the peak opens into a ring (its offset moves the ring's center)
    let round = smootherstep(peak.ring_radius / (2.0 * peak.width));
    let stretch = ((1.0 - round) * model.elongation_sigma * normal(&mut rng).abs()).exp().sqrt();
    let angle = rng.random_range(0.0..std::f64::consts::PI);
    let (sin, cos) = angle.sin_cos();
    let frame: Vec<(f64, f64)> = (0..nf)
        .map(|f| {
            let (px, py) = (region.face_proj_x[f] - dx, region.face_proj_y[f] - dy);
            ((cos * px + sin * py) / stretch, (-sin * px + cos * py) * stretch)
        })
        .collect();
    let distance: Vec<f64> = frame.iter().map(|&(u, v)| u.hypot(v)).collect();
    // The unstretched coordinates of each face, for the noise
    let plain: Vec<(f64, f64)> = (0..nf).map(|f| (region.face_proj_x[f], region.face_proj_y[f])).collect();

    // The faces the columns can occupy
    let support: Vec<usize> = (0..nf)
        .filter(|&f| (distance[f] - peak.ring_radius).abs() <= model.support * peak.width)
        .collect();
    if support.is_empty() {
        return Ok(change);
    }

    // Two independent noise fields in the plane about the crater center, normalized over the support: one for the columns, and
    // one for the critical angle of their flanks
    let Some(noise) = normalized_noise(&plain, &support, model.noise_wavelength * peak.width, model, &mut rng)? else {
        return Ok(change);
    };
    let Some(angle_noise) = normalized_noise(&plain, &support, model.angle_noise_wavelength * peak.width, model, &mut rng)? else {
        return Ok(change);
    };

    // The columns, in units of the envelope's amplitude
    let mut shape = Array1::<f64>::zeros(nf);
    for (i, &f) in support.iter().enumerate() {
        if noise[i] > model.threshold {
            let u = (distance[f] - peak.ring_radius) / peak.width;
            shape[f] = (-0.5 * u * u).exp() * (noise[i] - model.threshold).powf(model.ramp_exponent);
        }
    }
    let top = shape.fold(0.0f64, |m, &v| m.max(v));
    if !(top > 0.0) {
        return Ok(change);
    }

    // The critical slope of each face: the mean angle outside the support, varying with the angle noise inside it
    let mut critical_slope = Array1::<f64>::from_elem(nf, model.critical_angle.to_radians().tan());
    for (i, &f) in support.iter().enumerate() {
        let z = angle_noise[i].clamp(-PEAK_ANGLE_CLIP, PEAK_ANGLE_CLIP);
        critical_slope[f] = (model.critical_angle + model.angle_spread * z).to_radians().tan();
    }

    // The columns, scaled so that the highest is the peak's height, with flanks at the critical slope
    let columns = shape * (peak.height / top);
    let zero_nodes = Array1::<f64>::zeros(region.n_node());
    let peak_relief = {
        let view = region.kernel_view(columns.view(), zero_nodes.view());
        slope_limited_envelope(columns.view(), critical_slope.view(), &view)?
    };

    // Remove the peak's volume from the floor as a layer that thins to zero at its edge
    let volume: f64 = peak_relief.iter().zip(region.face_area.iter()).map(|(&h, &a)| h * a).sum();
    let taper: Vec<(usize, f64)> = (0..nf)
        .filter_map(|f| {
            let u = region.face_distance[f] / peak.floor_radius;
            (u < 1.0).then(|| (f, (1.0 - u * u).powi(2)))
        })
        .collect();
    let taper_volume: f64 = taper.iter().map(|&(f, w)| w * region.face_area[f]).sum();
    change.assign(&peak_relief);
    if taper_volume > 0.0 {
        let layer = volume / taper_volume;
        for &(f, w) in &taper {
            change[f] -= layer * w;
        }
    }
    Ok(change)
}

/// A standard normal deviate (Box-Muller).
fn normal(rng: &mut ChaCha8Rng) -> f64 {
    let u1: f64 = rng.random_range(f64::MIN_POSITIVE..1.0);
    let u2: f64 = rng.random_range(0.0..1.0);
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

/// Turbulence noise at the given faces (at the given planar coordinates), normalized to zero mean and unit standard
/// deviation over them, or `None` if it does not vary.
fn normalized_noise(
    frame: &[(f64, f64)],
    faces: &[usize],
    wavelength: f64,
    model: &PeakModel,
    rng: &mut ChaCha8Rng,
) -> Result<Option<Array1<f64>>, String> {
    let anchor = Array2::from_shape_fn((model.octaves, 3), |_| rng.random_range(0.0..std::f64::consts::TAU));
    let x = Array1::from_iter(faces.iter().map(|&f| frame[f].0));
    let y = Array1::from_iter(faces.iter().map(|&f| frame[f].1));
    let z = Array1::<f64>::zeros(faces.len());
    let noise = turbulence_noise(
        x.view(),
        y.view(),
        z.view(),
        1.0,
        wavelength,
        2.0,
        model.persistence,
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
