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
//! proportion to their parts of the basic profile's height there: the ejecta thickness is the basic one times
//! (1 + D / h), where h is the height of the basic profile, or zero where the field takes the topography below the reference
//! surface, and the crater takes the rest of the field. The ejecta is so never negative, and crater formation never clips it.
//! (Given to the ejecta alone, the variation of the rim, which is mostly uplift, would make the ejecta negative just outside the
//! rim, and clipped to zero there it would leave flat patches with a sharp edge at the rim.) Below the reference surface (the
//! crater's bowl) the field belongs to the crater.

use crate::ArrayResult;
use crate::morphology::basicmoon::{BasicMoonCrater, basicmoon_profile};
use numpy::ndarray::prelude::*;
use rayon::prelude::*;

/// The rim outline and field of variation of a crater.
pub struct RealMoonField<'a> {
    /// R_rim(theta) / R at the `n_theta` bearings of the field.
    pub rim: ArrayView1<'a, f64>,
    /// The field D (m), `n_theta` x `n_rho`.
    pub field: ArrayView2<'a, f64>,
    /// The largest warped radius of the field.
    pub rho_max: f64,
}

impl RealMoonField<'_> {
    fn n_theta(&self) -> usize {
        self.field.nrows()
    }

    /// The indices of the four bearings around `bearing` (radians), periodic, and the position between the middle two.
    #[inline]
    fn theta_stencil(&self, bearing: f64) -> ([usize; 4], f64) {
        let n = self.n_theta() as i64;
        let x = (bearing.to_degrees().rem_euclid(360.0)) * n as f64 / 360.0 - 0.5;
        let i = x.floor();
        let t = x - i;
        let i = i as i64;
        let idx = [-1, 0, 1, 2].map(|d| (i + d).rem_euclid(n) as usize);
        (idx, t)
    }

    /// The rim outline R_rim / R at a bearing (radians), interpolated with a periodic cubic.
    #[inline]
    pub fn rim_at(&self, bearing: f64) -> f64 {
        let (idx, t) = self.theta_stencil(bearing);
        let w = cubic_weights(t);
        (0..4).map(|a| w[a] * self.rim[idx[a]]).sum()
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
        let (idx, t) = self.theta_stencil(bearing);
        let wt = cubic_weights(t);
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
    // The warped radius of each point, and the distance at which the basic profile is evaluated
    let rho: Vec<f64> = radial_distances
        .iter()
        .zip(bearings.iter())
        .map(|(&r, &b)| r / (radius * data.rim_at(b)))
        .collect();
    let warped = Array1::from_iter(rho.iter().map(|&x| x * radius));
    // The central peak is built on the surface by the central peak model (see `central_peak`), not by the basic profile's
    // Gaussian peak
    let mut basic = *crater;
    basic.peak_height = 0.0;
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
                (e * (1.0 + field / t)).max(0.0) - e
            } else {
                0.0
            };
            if include_ejecta { ejecta_change } else { field - ejecta_change }
        })
        .collect();
    h.iter_mut().zip(add.iter()).for_each(|(v, a)| *v += a);
    Ok(h)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data(n_theta: usize, n_rho: usize, value: f64) -> (Array1<f64>, Array2<f64>) {
        (
            Array1::from_elem(n_theta, 1.0),
            Array2::from_elem((n_theta, n_rho), value),
        )
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
        };
        // Bearing of grid point 0 is 22.5 degrees; halfway between the last (337.5) and first (22.5) bearings is 0 degrees
        assert!((d.field_at(22.5_f64.to_radians(), 1.0) - 0.0).abs() < 1e-12);
        assert!((d.field_at(0.0, 1.0) - 3.5).abs() < 1e-12);
        assert_eq!(d.field_at(0.3, 2.5), 0.0);
    }

    #[test]
    fn zero_field_and_round_rim_give_the_basic_profile() {
        let crater = BasicMoonCrater {
            diameter: 10e3,
            radius: 5e3,
            floor_elevation: -1500.0,
            floor_radius: 1000.0,
            wall_curvature: 0.2,
            floor_blend: 0.5,
            rim_width: 1200.0,
            rim_height: 400.0,
            frac_ejrim: 0.3,
            ejprofile: -3.0,
            peak_height: 0.0,
            peak_width: 0.0,
            peak_ring_radius: 0.0,
            peak_center_distance: 0.0,
            peak_center_bearing: 0.0,
            scarp_radius: 0.0,
            scarp_height: 0.0,
            scarp_width: 0.0,
        };
        let (rim, field) = data(16, 11, 0.0);
        let d = RealMoonField {
            rim: rim.view(),
            field: field.view(),
            rho_max: 2.0,
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
