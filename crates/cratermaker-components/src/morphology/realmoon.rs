//! The realistic lunar crater model: the basic lunar crater profile, warped by an irregular rim outline, plus a random field of
//! variation around the crater (terraces, slumps, wall and floor texture, and the texture of the ejecta), the H-field.
//!
//! A crater is described on a polar grid of bearing theta and warped radius rho = r / R_rim(theta), where R_rim(theta) is the
//! radius of its rim crest at each bearing. The elevation is
//!
//! ```text
//! h(r, theta) = basic profile(rho R) + D(theta, rho)
//! ```
//!
//! where the basic profile is that of `basicmoon_profile` evaluated at the warped radius (so its rim crest follows the outline),
//! and D is the field, given on a regular grid of `n_theta` bearings (theta_i = (i + 1/2) 360 / n_theta degrees) and `n_rho`
//! radii from 0 to `rho_max` (inclusive), in meters. The field is drawn in Python from the H-field model of fresh lunar
//! craters (Lunar_Morphometry_Revisited). Beyond `rho_max` it is zero.
//!
//! The field is the variation of the whole topography about the mean profile, so where the basic profile stands above the
//! reference surface (the rim and the ejecta blanket) it is shared between the crater (structural uplift) and the ejecta in
//! proportion to their parts of the basic profile's height there: the ejecta thickness e becomes e + s D, where s = e / h is
//! the ejecta's share of the height h of the basic profile (at most 1: on the upper inner wall the basic ejecta can be thicker
//! than the profile is high, and an unbounded share would multiply the field there), or zero where that is negative, and the
//! crater takes the rest of the field. The ejecta is so never negative, and crater formation never clips it.
//! (Given to the ejecta alone, the variation of the rim, which is mostly uplift, would make the ejecta negative just outside the
//! rim, and clipped to zero there it would leave flat patches with a sharp edge at the rim.) Below the reference surface (the
//! crater's bowl) the field belongs to the crater.
//!
//! The walls of complex craters and basins are stacked levels (see `WallLevels`): the crest and the terrace scarps below it, each a
//! continuous, nearly round curve with lobes of its own, whose scarps merge where they collide. For these craters the field holds
//! only the fine texture, and the rim outline is the outermost level.

use crate::ArrayResult;
use crate::morphology::basicmoon::{BasicMoonCrater, basicmoon_profile, scarp_function};
use crate::morphology::smootherstep;
use noise::{NoiseFn, SuperSimplex};
use numpy::ndarray::prelude::*;
use rayon::prelude::*;
use std::sync::OnceLock;

/// The terrace between two levels of the wall (see `WallLevels`): going inward from a level, a scarp drops over the outer
/// `TERRACE_SCARP_FRACTION` of the annulus, then a bench slopes gently to the next level, keeping `TERRACE_BENCH_SHARE` of the
/// annulus's drop.
pub const TERRACE_SCARP_FRACTION: f64 = 0.35;
pub const TERRACE_BENCH_SHARE: f64 = 0.3;
/// The elevation offset of the crest fades out over this distance (in rim radii) outside it.
pub const WALL_FADE: f64 = 0.2;
/// The wall turbulence is tapered over this distance (in rim radii) at the crest and at the lowest level.
pub const WALL_TURBULENCE_TAPER: f64 = 0.05;
/// Octaves of the wall turbulence, and the amplitude ratio between successive octaves.
pub const WALL_TURBULENCE_OCTAVES: u32 = 6;
pub const WALL_TURBULENCE_PERSISTENCE: f64 = 0.5;

/// The stacked levels of the wall of a complex crater: the rim crest and the terrace scarps below it, each a continuous,
/// nearly round closed curve with lobes of its own, with a scarp whose strength and an elevation offset vary along it. The
/// levels wobble into each other: at each bearing they are taken in order of radius, and where two collide their scarps merge.
/// Between consecutive levels the crater part of the profile is a terrace (see `TERRACE_SCARP_FRACTION`) whose heights are those
/// of the basic profile at the levels' nominal radii plus their offsets, blended (by the scarp strength) with the basic profile
/// between them; inside the lowest level the basic profile follows that level's outline. Over the wall (from the lowest level to
/// the crest) a turbulence field is added: octaves of the absolute value of simplex noise, normalized to zero mean and unit
/// standard deviation, times `turbulence`.
pub struct WallLevels<'a> {
    /// The nominal radii of the levels (rim radii), from the crest (1) inward.
    pub nominal: ArrayView1<'a, f64>,
    /// For each level (n_levels x n_bearings, at the bearings of the field): its radius (rim radii), the strength of its scarp
    /// (0 to 1) and its elevation offset (m), stacked along the first axis (3 x n_levels x n_bearings).
    pub levels: ArrayView3<'a, f64>,
    /// The standard deviation of the wall turbulence (m).
    pub turbulence: f64,
    /// The wavelength of the largest octave of the wall turbulence (rim radii).
    pub wavelength: f64,
    /// The seed of the wall turbulence.
    pub seed: u32,
}

/// The mean and standard deviation of the sum of octaves of |simplex noise| (see `WallLevels`), sampled once.
fn turbulence_moments() -> (f64, f64) {
    static MOMENTS: OnceLock<(f64, f64)> = OnceLock::new();
    *MOMENTS.get_or_init(|| {
        let gens: Vec<SuperSimplex> = (0..WALL_TURBULENCE_OCTAVES).map(SuperSimplex::new).collect();
        let n = 256;
        let values: Vec<f64> = (0..n * n).map(|i| raw_turbulence(&gens, (i % n) as f64 * 0.137, (i / n) as f64 * 0.137)).collect();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let sd = (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64).sqrt();
        (mean, sd)
    })
}

#[inline]
fn raw_turbulence(gens: &[SuperSimplex], x: f64, y: f64) -> f64 {
    let mut value = 0.0;
    let mut amplitude = 1.0;
    let mut scale = 1.0;
    for g in gens {
        value += amplitude * g.get([x * scale, y * scale]).abs();
        amplitude *= WALL_TURBULENCE_PERSISTENCE;
        scale *= 2.0;
    }
    value
}

#[inline]
fn terrace(t: f64) -> f64 {
    let edge = 1.0 - TERRACE_SCARP_FRACTION;
    if t < edge {
        TERRACE_BENCH_SHARE * t / edge
    } else {
        let s = ((t - edge) / TERRACE_SCARP_FRACTION).clamp(0.0, 1.0);
        TERRACE_BENCH_SHARE + (1.0 - TERRACE_BENCH_SHARE) * s * s * (3.0 - 2.0 * s)
    }
}

/// The rim outline and field of variation of a crater.
pub struct RealMoonField<'a> {
    /// R_rim(theta) / R at the `n_theta` bearings of the field.
    pub rim: ArrayView1<'a, f64>,
    /// The field D (m), `n_theta` x `n_rho`.
    pub field: ArrayView2<'a, f64>,
    /// The largest warped radius of the field.
    pub rho_max: f64,
    /// The radius and height of the ring scarp of a multiring basin around the basin, as factors of `scarp_radius` and
    /// `scarp_height`, at bearings (i + 1/2) 360 / n degrees; empty for a crater without a scarp (or a uniform scarp).
    pub scarp_radius: ArrayView1<'a, f64>,
    pub scarp_height: ArrayView1<'a, f64>,
    /// The stacked levels of the wall of a complex crater (see `WallLevels`), or None.
    pub wall: Option<WallLevels<'a>>,
}

/// The periodic Catmull-Rom stencil at a bearing (radians) for a function given at `n` bearings (i + 1/2) 360 / n degrees: the
/// indices of the four bearings around it and their weights.
#[inline]
fn bearing_stencil(n: usize, bearing: f64) -> ([usize; 4], [f64; 4]) {
    let n = n as i64;
    let x = (bearing.to_degrees().rem_euclid(360.0)) * n as f64 / 360.0 - 0.5;
    let i = x.floor();
    let w = cubic_weights(x - i);
    let i = i as i64;
    ([-1, 0, 1, 2].map(|d| (i + d).rem_euclid(n) as usize), w)
}

/// A periodic function given at its bearings (see `bearing_stencil`), interpolated with a stencil.
#[inline]
fn apply_stencil(values: ArrayView1<'_, f64>, stencil: &([usize; 4], [f64; 4])) -> f64 {
    (0..4).map(|a| stencil.1[a] * values[stencil.0[a]]).sum()
}

/// A periodic function given at `n` bearings (i + 1/2) 360 / n degrees, interpolated at a bearing (radians) with a periodic
/// Catmull-Rom cubic.
#[inline]
fn periodic_at(values: ArrayView1<'_, f64>, bearing: f64) -> f64 {
    apply_stencil(values, &bearing_stencil(values.len(), bearing))
}

/// The warped radius inside the ring scarp of a multiring basin: the peak-ring basin inside the scarp (its floor, inner wall and
/// the scarp itself) follows the outline of the scarp, so inside the scarp's crest the warped radius is divided by the scarp's
/// radius factor too, a factor that fades out (a smootherstep in the warped radius) between the crest and the main rim.
#[inline]
fn scarp_warp(rho: f64, bearing: f64, scarp_rho: f64, scarp_radius: ArrayView1<'_, f64>) -> f64 {
    let weight = 1.0 - smootherstep((rho - scarp_rho) / (1.0 - scarp_rho));
    rho / periodic_at(scarp_radius, bearing).powf(weight)
}

impl RealMoonField<'_> {
    fn n_theta(&self) -> usize {
        self.field.nrows()
    }

    /// The rim outline R_rim / R at a bearing (radians), interpolated with a periodic cubic.
    #[inline]
    pub fn rim_at(&self, bearing: f64) -> f64 {
        periodic_at(self.rim, bearing)
    }

    /// The field (m) at a bearing (radians) and warped radius, interpolated bicubically (Catmull-Rom; periodic in bearing,
    /// with the edge values held in radius), so that it is smooth between the points of the field; zero beyond `rho_max`.
    #[inline]
    pub fn field_at(&self, bearing: f64, rho: f64) -> f64 {
        if !(rho >= 0.0) || rho > self.rho_max {
            return 0.0;
        }
        let n_rho = self.field.ncols() as i64;
        let y = rho / self.rho_max * (n_rho - 1) as f64;
        let j = (y.floor() as i64).min(n_rho - 2);
        let v = y - j as f64;
        let jdx = [-1, 0, 1, 2].map(|d| (j + d).clamp(0, n_rho - 1) as usize);
        let (idx, wt) = bearing_stencil(self.n_theta(), bearing);
        let wr = cubic_weights(v);
        let f = &self.field;
        let mut value = 0.0;
        for a in 0..4 {
            let row: f64 = (0..4).map(|b| wr[b] * f[[idx[a], jdx[b]]]).sum();
            value += wt[a] * row;
        }
        value
    }
}

/// The weights of the four points of a Catmull-Rom cubic for a position `t` (0 to 1) between the middle two. They sum to one,
/// and give the middle points exactly at `t` = 0 and 1.
#[inline]
fn cubic_weights(t: f64) -> [f64; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    [
        0.5 * (-t3 + 2.0 * t2 - t),
        0.5 * (3.0 * t3 - 5.0 * t2 + 2.0),
        0.5 * (-3.0 * t3 + 4.0 * t2 + t),
        0.5 * (t3 - t2),
    ]
}

/// Elevations of a realistic lunar crater relative to `reference_elevations` (see the module documentation).
///
/// # Arguments
///
/// * `radial_distances` - Distances from the crater center (m).
/// * `bearings` - Bearings from the crater center (radians, clockwise from north).
/// * `reference_elevations` - The surface the profile is drawn relative to (m).
/// * `crater` - The basic crater.
/// * `data` - The rim outline and field.
/// * `include_crater`, `include_ejecta` - Which parts of the profile to include. The field is shared between them in proportion
///   to their parts of the height of the basic profile (see the module documentation).
pub fn realmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &BasicMoonCrater,
    data: &RealMoonField,
    include_crater: bool,
    include_ejecta: bool,
) -> ArrayResult {
    if bearings.len() != radial_distances.len()
        || reference_elevations.len() != radial_distances.len()
    {
        return Err(
            "radial_distances, bearings and reference_elevations must have the same length".into(),
        );
    }
    if data.rim.len() != data.n_theta() || data.field.ncols() < 2 || data.n_theta() < 2 {
        return Err("the rim outline must have one value per bearing of the field, and the field at least 2 x 2 values".into());
    }
    let radius = crater.radius;
    // The warped radius of each point, and the distance at which the basic profile is evaluated (for a multiring basin, with
    // the warp of its ring scarp, see `scarp_warp`)
    let varying_scarp = crater.has_scarp() && !data.scarp_radius.is_empty() && data.scarp_radius.len() == data.scarp_height.len();
    let scarp_rho = crater.scarp_radius / radius;
    let rho: Vec<f64> = radial_distances
        .iter()
        .zip(bearings.iter())
        .map(|(&r, &b)| {
            let rho0 = r / (radius * data.rim_at(b));
            if varying_scarp { scarp_warp(rho0, b, scarp_rho, data.scarp_radius) } else { rho0 }
        })
        .collect();
    let warped = Array1::from_iter(rho.iter().map(|&x| x * radius));
    // The central peak is built on the surface by the central peak model (see `central_peak`), not by the basic profile's
    // Gaussian peak. A ring scarp that varies around the basin is added below, so it is left out of the basic profile.
    let mut basic = *crater;
    basic.peak_height = 0.0;
    if varying_scarp {
        basic.scarp_height = 0.0;
    }
    let mut h = basicmoon_profile(
        warped.view(),
        reference_elevations,
        &basic,
        include_crater,
        include_ejecta,
    )?;
    // The basic profile's height and its ejecta part, which set the ejecta's share of the field
    let (total, ejecta) = if include_crater && include_ejecta {
        (Array1::<f64>::zeros(0), Array1::<f64>::zeros(0)) // the whole field is included, so its split does not matter
    } else {
        let zeros = Array1::<f64>::zeros(warped.len());
        (
            basicmoon_profile(warped.view(), zeros.view(), &basic, true, true)?,
            basicmoon_profile(warped.view(), zeros.view(), &basic, false, true)?,
        )
    };
    let add: Vec<f64> = (0..rho.len())
        .into_par_iter()
        .with_min_len(4096)
        .map(|i| {
            if !include_crater && !include_ejecta {
                return 0.0;
            }
            let field = data.field_at(bearings[i], rho[i]);
            if include_crater && include_ejecta {
                return field;
            }
            // The change of the ejecta: its thickness becomes e (1 + D / h), or zero where that is negative
            let (t, e) = (total[i], ejecta[i]);
            let ejecta_change = if t > 0.0 && e > 0.0 {
                let share = (e / t).min(1.0);
                (e + share * field).max(0.0) - e
            } else {
                0.0
            };
            if include_ejecta { ejecta_change } else { field - ejecta_change }
        })
        .collect();
    h.iter_mut().zip(add.iter()).for_each(|(v, a)| *v += a);
    // The ring scarp, as the rim of the peak-ring basin inside the multiring basin: its height varies with bearing, and its
    // radius through the warped radius (above). It is part of the crater.
    if varying_scarp && include_crater {
        h.iter_mut().enumerate().for_each(|(i, v)| {
            let mut local = *crater;
            local.scarp_height = crater.scarp_height * periodic_at(data.scarp_height, bearings[i]).max(0.0);
            *v += scarp_function(warped[i], &local);
        });
    }
    if include_crater && let Some(wall) = &data.wall {
        let scarp = if varying_scarp { Some(data.scarp_radius) } else { None };
        let delta = wall_delta(radial_distances, bearings, &basic, wall, &warped, scarp, data.n_theta())?;
        h.iter_mut().zip(delta.iter()).for_each(|(v, d)| *v += d);
    }
    Ok(h)
}

/// The change of the crater part of the profile made by the wall levels (see `WallLevels`): at each point, the terraced wall
/// (or, inside the lowest level, the basic profile following that level's outline, and outside the crest, the crest's fading
/// elevation offset) less the basic profile at the warped radius, plus the wall turbulence.
fn wall_delta(
    radial_distances: ArrayView1<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    basic: &BasicMoonCrater,
    wall: &WallLevels<'_>,
    warped: &Array1<f64>,
    scarp_radius: Option<ArrayView1<'_, f64>>,
    n_theta: usize,
) -> ArrayResult {
    let n_levels = wall.nominal.len();
    let radius = basic.radius;
    if n_levels < 2 || wall.levels.shape() != [3, n_levels, n_theta] {
        return Err("the wall levels must be 3 x n_levels x n_theta, with at least two levels".into());
    }
    let n = radial_distances.len();
    let zeros = Array1::<f64>::zeros(n);
    // The basic crater part at the warped radius, and at the nominal radii of the levels
    let base = basicmoon_profile(warped.view(), zeros.view(), basic, true, false)?;
    let knots_r = Array1::from_iter(wall.nominal.iter().map(|&x| x * radius));
    let knots = basicmoon_profile(knots_r.view(), Array1::<f64>::zeros(n_levels).view(), basic, true, false)?;
    // At each point: the zone, the terrace parameters and the radius at which the basic profile is needed
    let scarp_rho = basic.scarp_radius / radius;
    let gens: Vec<SuperSimplex> = (0..WALL_TURBULENCE_OCTAVES).map(|o| SuperSimplex::new(wall.seed.wrapping_add(o))).collect();
    let (t_mean, t_sd) = turbulence_moments();
    struct Point {
        query: f64,
        value: f64,
        kind: u8, // 0: no change, 1: terrace (lo, hi, ...; query: the smooth profile), 2: inside (query), 3: outside (value: offset)
        hi: f64,
        lo: f64,
        k_hi: f64,
        k_lo: f64,
        f_terrace: f64,
        strength: f64,
        turb: f64,
    }
    let points: Vec<Point> = (0..n)
        .into_par_iter()
        .with_min_len(2048)
        .map(|i| {
            let b = bearings[i];
            let x = radial_distances[i] / radius;
            let stencil = bearing_stencil(n_theta, b);
            let mut lv: Vec<(f64, f64, f64)> = (0..n_levels)
                .map(|j| {
                    (
                        apply_stencil(wall.levels.slice(s![0, j, ..]), &stencil),
                        apply_stencil(wall.levels.slice(s![1, j, ..]), &stencil).clamp(0.0, 1.0),
                        apply_stencil(wall.levels.slice(s![2, j, ..]), &stencil),
                    )
                })
                .collect();
            lv.sort_by(|a, c| c.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            let crest = lv[0].0;
            let lowest = lv[n_levels - 1].0;
            let mut p = Point { query: 0.0, value: 0.0, kind: 0, hi: 0.0, lo: 0.0, k_hi: 0.0, k_lo: 0.0, f_terrace: 0.0, strength: 0.0, turb: 0.0 };
            if x >= crest {
                if x < crest + WALL_FADE {
                    let u = 1.0 - (x - crest) / WALL_FADE;
                    p.kind = 3;
                    p.value = lv[0].2 * u * u * (3.0 - 2.0 * u);
                }
                return p;
            }
            if wall.turbulence > 0.0 {
                let w = ((crest - x).min(x - lowest) / WALL_TURBULENCE_TAPER + 0.5).clamp(0.0, 1.0);
                if w > 0.0 {
                    let (sb, cb) = b.sin_cos();
                    let raw = raw_turbulence(&gens, x * sb / wall.wavelength, x * cb / wall.wavelength);
                    p.turb = wall.turbulence * w * (raw - t_mean) / t_sd;
                }
            }
            if x < lowest {
                // The basic profile follows the lowest level's outline (and, in a multiring basin, the scarp's)
                let mut rho_in = x * wall.nominal[n_levels - 1] / lowest;
                if let Some(scarp) = scarp_radius {
                    rho_in = scarp_warp(rho_in, b, scarp_rho, scarp);
                }
                p.kind = 2;
                p.query = rho_in * radius;
                return p;
            }
            let j = (0..n_levels - 1).find(|&j| x < lv[j].0 && x >= lv[j + 1].0).unwrap_or(n_levels - 2);
            let t = (x - lv[j + 1].0) / (lv[j].0 - lv[j + 1].0).max(1e-12);
            p.kind = 1;
            p.k_hi = knots[j];
            p.k_lo = knots[j + 1];
            p.hi = knots[j] + lv[j].2;
            p.lo = knots[j + 1] + if j + 1 == n_levels - 1 { 0.0 } else { lv[j + 1].2 };
            p.f_terrace = terrace(t);
            p.strength = lv[j].1;
            p.query = (wall.nominal[j + 1] + t * (wall.nominal[j] - wall.nominal[j + 1])) * radius;
            p
        })
        .collect();
    let queries = Array1::from_iter(points.iter().map(|p| p.query));
    let at_query = basicmoon_profile(queries.view(), zeros.view(), basic, true, false)?;
    Ok(Array1::from_iter(points.iter().enumerate().map(|(i, p)| match p.kind {
        1 => {
            let span = p.k_hi - p.k_lo;
            let smooth = if span.abs() > 0.0 { (at_query[i] - p.k_lo) / span } else { p.f_terrace };
            let f = p.strength * p.f_terrace + (1.0 - p.strength) * smooth;
            p.lo + (p.hi - p.lo) * f - base[i] + p.turb
        }
        2 => at_query[i] - base[i] + p.turb,
        3 => p.value,
        _ => 0.0,
    })))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A crater without peak or scarp, with the given radius, floor elevation, floor radius, rim width and rim height.
    fn test_crater(radius: f64, floor_elevation: f64, floor_radius: f64, rim_width: f64, rim_height: f64) -> BasicMoonCrater {
        BasicMoonCrater {
            diameter: 2.0 * radius,
            radius,
            floor_elevation,
            floor_radius,
            wall_curvature: 0.2,
            floor_blend: 0.5,
            rim_width,
            rim_height,
            frac_ejrim: 0.3,
            ejprofile: -3.0,
            peak_height: 0.0,
            peak_width: 0.0,
            peak_ring_radius: 0.0,
            peak_center_distance: 0.0,
            peak_center_bearing: 0.0,
            peak_elongation: 1.0,
            peak_elongation_bearing: 0.0,
            peak_seed: 0,
            scarp_radius: 0.0,
            scarp_height: 0.0,
            scarp_width: 0.0,
        }
    }

    fn data(n_theta: usize, n_rho: usize, value: f64) -> (Array1<f64>, Array2<f64>) {
        (
            Array1::from_elem(n_theta, 1.0),
            Array2::from_elem((n_theta, n_rho), value),
        )
    }

    #[test]
    fn round_levels_without_scarps_or_offsets_give_the_basic_profile() {
        // With round levels at their nominal radii, no scarps, no offsets and no turbulence, the terraced wall is the basic profile
        let crater = test_crater(50e3, -3000.0, 30e3, 8e3, 1000.0);
        let n_theta = 16;
        let (rim, field) = data(n_theta, 11, 0.0);
        let nominal = Array1::from(vec![1.0, 0.9, 0.78, 0.65]);
        let mut levels = Array3::<f64>::zeros((3, 4, n_theta));
        for j in 0..4 {
            levels.slice_mut(s![0, j, ..]).fill(nominal[j]);
        }
        let empty: [f64; 0] = [];
        let with_wall = RealMoonField {
            rim: rim.view(),
            field: field.view(),
            rho_max: 2.0,
            scarp_radius: ArrayView1::from(&empty[..]),
            scarp_height: ArrayView1::from(&empty[..]),
            wall: Some(WallLevels { nominal: nominal.view(), levels: levels.view(), turbulence: 0.0, wavelength: 0.25, seed: 1 }),
        };
        let r = Array1::linspace(0.0, 80e3, 161);
        let b = Array1::linspace(0.0, 6.0, 161);
        let z = Array1::zeros(161);
        let a = realmoon_profile(r.view(), b.view(), z.view(), &crater, &with_wall, true, true).unwrap();
        let e = basicmoon_profile(r.view(), z.view(), &crater, true, true).unwrap();
        for (x, y) in a.iter().zip(e.iter()) {
            assert!((x - y).abs() < 1e-6, "{x} != {y}");
        }
    }

    #[test]
    fn field_interpolation_is_periodic_and_zero_outside() {
        let (rim, mut field) = data(8, 5, 0.0);
        for i in 0..8 {
            field.row_mut(i).fill(i as f64);
        }
        let d = RealMoonField {
            rim: rim.view(),
            field: field.view(),
            rho_max: 2.0,
            scarp_radius: ArrayView1::from(&[] as &[f64]),
            scarp_height: ArrayView1::from(&[] as &[f64]),
            wall: None,
        };
        // Bearing of grid point 0 is 22.5 degrees; halfway between the last (337.5) and first (22.5) bearings is 0 degrees
        assert!((d.field_at(22.5_f64.to_radians(), 1.0) - 0.0).abs() < 1e-12);
        assert!((d.field_at(0.0, 1.0) - 3.5).abs() < 1e-12);
        assert_eq!(d.field_at(0.3, 2.5), 0.0);
    }

    #[test]
    fn zero_field_and_round_rim_give_the_basic_profile() {
        let crater = test_crater(5e3, -1500.0, 1000.0, 1200.0, 400.0);
        let (rim, field) = data(16, 11, 0.0);
        let d = RealMoonField {
            rim: rim.view(),
            field: field.view(),
            rho_max: 2.0,
            scarp_radius: ArrayView1::from(&[] as &[f64]),
            scarp_height: ArrayView1::from(&[] as &[f64]),
            wall: None,
        };
        let r = Array1::linspace(0.0, 20e3, 101);
        let b = Array1::linspace(0.0, 6.0, 101);
        let z = Array1::zeros(101);
        let a =
            realmoon_profile(r.view(), b.view(), z.view(), &crater, &d, true, true).unwrap();
        let e = basicmoon_profile(r.view(), z.view(), &crater, true, true).unwrap();
        for (x, y) in a.iter().zip(e.iter()) {
            assert!((x - y).abs() < 1e-9);
        }
    }
}
