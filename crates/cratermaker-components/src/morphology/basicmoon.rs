use crate::ArrayResult;
use crate::morphology::{blend_step, smooth_clamp, smooth_positive_part, smootherstep};
use crate::rootfind::{BRENTQ_MAXITER, BRENTQ_RTOL, BRENTQ_XTOL, BrentqResult, brentq};
use itertools::Itertools;
use libm::erf;
use numpy::ndarray::prelude::*;
use pyo3::FromPyObject;
use rand::SeedableRng;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha12Rng;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::f64::consts::{PI, SQRT_2, TAU};

const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.90;

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject, Clone, Copy, Debug)]
pub struct BasicMoonCrater {
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
}

// Computes a either a crater and/or ejecta 1D profile array from input radial distances and reference elevations using the model of
// Minton et al. (2026)
//
//
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater` - BasicMoonCrater struct containing crater dimensions and morphology parameters.
/// * `include_crater` - Whether to include the crater profile in the output (true/false).
/// * `include_ejecta` - Whether to include the ejecta profile in the output (true/false).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Errors if the input arrays have mismatched lengths.
pub fn basicmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &BasicMoonCrater,
    rings: &Option<Vec<BasicMoonCrater>>,
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

    // Each point is independent, so evaluate them in parallel. The minimum chunk length keeps small regions (the vast majority
    // of craters) from paying more in thread coordination than they gain.
    let profile: Vec<f64> = (0..radial_distances.len())
        .into_par_iter()
        .with_min_len(PROFILE_PAR_MIN_LEN)
        .map(|i| {
            let r = radial_distances[i];
            let h = basicmoon_profile_one(
                r,
                reference_elevations[i],
                crater,
                rings,
                include_crater,
                include_ejecta,
            );
            if r <= crater.radius {
                h.max(min_elevation)
            } else {
                h
            }
        })
        .collect();
    Ok(Array1::from_vec(profile))
}

/// Minimum number of points per parallel work item in `basicmoon_profile`.
const PROFILE_PAR_MIN_LEN: usize = 2048;

// Computes a either a crater and/or ejecta 1D profile at a single point using the model of  Minton et al. (2026)
//
//
/// # Arguments
///
/// * `r` - Radial distance from crater center (in meters).
/// * `href` - the reference surface at the point r (in meters)
/// * `crater` - BasicMoonCrater struct containing crater dimensions and morphology parameters.
/// * `include_crater` - Whether to include the crater profile in the output (true/false).
/// * `include_ejecta` - Whether to include the ejecta profile in the output (true/false).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Errors if the input arrays have mismatched lengths.
pub fn basicmoon_profile_one(
    r: f64,
    href: f64,
    crater: &BasicMoonCrater,
    rings: &Option<Vec<BasicMoonCrater>>,
    include_crater: bool,
    include_ejecta: bool,
) -> f64 {
    let (offsets, ring_r) = match rings {
        Some(rings) => (ring_elevation_offsets(crater, rings), vec![r; rings.len()]),
        None => (Vec::new(), Vec::new()),
    };
    basicmoon_profile_one_local(r, &ring_r, href, crater, rings, &offsets, include_crater, include_ejecta)
}

/// The elevation offset of each ring: the elevation of the feature it sits in (the crater for the outermost ring, otherwise
/// the next ring out) at the ring's radius.
pub fn ring_elevation_offsets(crater: &BasicMoonCrater, rings: &[BasicMoonCrater]) -> Vec<f64> {
    rings
        .iter()
        .enumerate()
        .map(|(i, ring)| {
            if i == 0 {
                elevation_offset_func(crater, ring)
            } else {
                elevation_offset_func(&rings[i - 1], ring)
            }
        })
        .collect()
}

/// The elevation of a crater's own profile (without rings) at radial distance `r`.
#[inline]
pub fn crater_elevation(crater: &BasicMoonCrater, r: f64) -> f64 {
    crater_profile_function(
        r,
        crater.radius,
        crater.floor_elevation,
        crater.floor_radius,
        crater.wall_curvature,
        crater.floor_blend,
        crater.rim_width,
        crater.rim_height,
        crater.frac_ejrim,
        crater.ejprofile,
        crater.peak_height,
        crater.peak_width,
        crater.peak_ring_radius,
    )
}

/// `basicmoon_profile_one` for models in which the crater and each of its rings are evaluated at their own radial coordinate
/// and the elevation offsets of the rings are given (see `ring_elevation_offsets`).
///
/// # Arguments
///
/// * `r` - The radial coordinate at which the main crater is evaluated.
/// * `ring_r` - The radial coordinate at which each ring is evaluated.
/// * `ring_offsets` - The elevation offset of each ring.
/// * The others are as in `basicmoon_profile_one`.
#[allow(clippy::too_many_arguments)]
pub fn basicmoon_profile_one_local(
    r: f64,
    ring_r: &[f64],
    href: f64,
    crater: &BasicMoonCrater,
    rings: &Option<Vec<BasicMoonCrater>>,
    ring_offsets: &[f64],
    include_crater: bool,
    include_ejecta: bool,
) -> f64 {
    // Prevent invalid values of parameters
    let radius = crater.radius.max(1e-8);
    let floor_elevation = crater.floor_elevation.min(-1e-8);
    let floor_radius = crater.floor_radius.max(1e-8);
    let peak_height = crater.peak_height.max(0.0);
    let peak_width = crater.peak_width.max(0.0);
    let peak_ring_radius = crater.peak_ring_radius.max(0.0);

    // Guard against NaN in inputs by falling back to safe defaults
    let wall_curvature = if crater.wall_curvature.is_nan() {
        0.5
    } else {
        crater.wall_curvature.max(0.0).min(1.0)
    };
    let floor_blend = if crater.floor_blend.is_nan() {
        0.5
    } else {
        crater.floor_blend.max(0.0).min(1.0)
    };
    let rim_width = if crater.rim_width.is_nan() {
        radius * 0.1
    } else {
        crater.rim_width.max(radius * 1e-3)
    };
    let rim_height = if crater.rim_height.is_nan() {
        0.0
    } else {
        crater.rim_height.max(0.0)
    };
    let frac_ejrim = if crater.frac_ejrim.is_nan() {
        0.5
    } else {
        crater.frac_ejrim.max(0.0).min(1.0)
    };
    let ejprofile = if crater.ejprofile.is_nan() {
        -3.0
    } else {
        crater.ejprofile.max(-6.0).min(-0.25)
    };

    let mut hcrat = crater_profile_function(
        r,
        radius,
        floor_elevation,
        floor_radius,
        wall_curvature,
        floor_blend,
        rim_width,
        rim_height,
        frac_ejrim,
        ejprofile,
        peak_height,
        peak_width,
        peak_ring_radius,
    );

    // Craters without rings use a smooth split of the profile into its crater and ejecta parts (see smooth_crater_ejecta_split).
    // Craters with rings combine the main crater and its rings with the composition below.
    if !matches!(rings, Some(v) if !v.is_empty()) {
        let hej = smooth_crater_ejecta_split(
            r, hcrat, radius, rim_width, rim_height, frac_ejrim, ejprofile,
        );
        let crater_part = if include_crater { hcrat - hej } else { 0.0 };
        let ejecta_part = if include_ejecta { hej } else { 0.0 };
        let final_val = href + crater_part + ejecta_part;
        return if final_val.is_nan() { href } else { final_val };
    }

    match rings {
        Some(rings) => {
            for (i, ring) in rings.iter().enumerate() {
                let elevation_offset = ring_offsets[i];
                let mut ring_floor_elevation = crater.floor_elevation - elevation_offset;
                ring_floor_elevation = ring_floor_elevation.min(-1e-8);

                let ring_radius = ring.radius.max(1e-8);
                let ring_floor_radius = ring.floor_radius.max(1e-8);
                let ring_rim_width = ring.rim_width.max(ring_radius * 1e-3);
                let ring_rim_height = ring.rim_height.max(0.0);
                // Guard against NaN in inputs by falling back to safe defaults
                let ring_wall_curvature = if ring.wall_curvature.is_nan() {
                    0.5
                } else {
                    ring.wall_curvature.max(0.0).min(1.0)
                };
                let ring_frac_ejrim = if ring.frac_ejrim.is_nan() {
                    0.5
                } else {
                    ring.frac_ejrim.max(0.0).min(1.0)
                };
                let ring_ejprofile = if ring.ejprofile.is_nan() {
                    -3.0
                } else {
                    ring.ejprofile.max(-6.0).min(-0.25)
                };

                let hring = crater_profile_function(
                    ring_r[i],
                    ring_radius,
                    ring_floor_elevation,
                    ring_floor_radius,
                    ring_wall_curvature,
                    floor_blend,
                    ring_rim_width,
                    ring_rim_height,
                    ring_frac_ejrim,
                    ring_ejprofile,
                    peak_height,
                    peak_width,
                    peak_ring_radius,
                ) + elevation_offset;
                // Outside its rim a ring is flat at its offset (the elevation of the feature around it at the ring's radius).
                // That surface only belongs inside the rim of that feature: beyond it, it would be a plateau standing above the
                // ejecta wherever the ring sits high on the wall. Inside the rim, the wall of the surrounding feature rises
                // above it towards the crest, so dropping the ring there leaves the profile continuous.
                let inside_parent = if i == 0 { r <= radius } else { ring_r[i - 1] <= rings[i - 1].radius };
                if inside_parent {
                    hcrat = hcrat.max(hring);
                }
            }
        }
        None => (),
    }
    let mut hej = ejecta_profile_function(r, radius, frac_ejrim, rim_height, ejprofile, rim_width);

    match rings {
        Some(rings) => {
            for (i, ring) in rings.iter().enumerate() {
                let ring_radius = ring.radius.max(1e-8);
                let ring_rim_width = ring.rim_width.max(ring_radius * 1e-3);
                let ring_rim_height = ring.rim_height.max(0.0);
                let ring_frac_ejrim = ring.frac_ejrim.max(0.0).min(1.0);
                let ring_ejprofile = ring.ejprofile.max(-6.0).min(-0.25);

                let mut hejring = ejecta_profile_function(
                    ring_r[i],
                    ring_radius,
                    ring_frac_ejrim,
                    ring_rim_height,
                    ring_ejprofile,
                    ring_rim_width,
                );
                let elevation_offset = ring_offsets[i];
                let hrel = hcrat - elevation_offset;
                if ring_r[i] <= ring.radius && hrel >= 0.0 {
                    let hup = (1.0 - ring.frac_ejrim) * ring.rim_height;
                    let max_clamp = (ring_frac_ejrim * ring_rim_height).max(0.0);
                    if !hup.is_nan() && !max_clamp.is_nan() {
                        hejring += (hrel - hup).clamp(0.0, max_clamp);
                    }
                }
                hej += hejring;
            }
        }
        None => (),
    }

    // Add the upper portion of the rim to the ejecta
    if r <= crater.radius && hcrat >= 0.0 {
        let hup = rimfunc(
            r,
            crater.radius,
            (1.0 - crater.frac_ejrim) * crater.rim_height,
            0.0,
            0.5 * crater.rim_width,
        );
        let max_clamp = (frac_ejrim * rim_height).max(0.0);
        if !hup.is_nan() && !max_clamp.is_nan() {
            hej += (hcrat - hup).clamp(0.0, max_clamp);
        }
    }

    if include_crater {
        if r > crater.radius || (hcrat > 0.0 && hej > 0.0) {
            hcrat = (hcrat - hej).max(0.0);
        }
    } else {
        hcrat = 0.0;
    }
    if !include_ejecta {
        hej = 0.0;
    }

    let final_val = href + hcrat + hej;
    if final_val.is_nan() { href } else { final_val }
}

/// Half-width of the smoothing at each end of the clamp that assigns the upper rim to the ejecta, as a fraction of the ejecta
/// rim height.
const SPLIT_CLAMP_SMOOTHING: f64 = 0.5;

/// Half-width of the blend across the rim from the inside to the outside form of the ejecta part, as a fraction of the rim
/// width.
const SPLIT_RIM_BLEND: f64 = 0.25;

/// Returns the ejecta part of a crater profile `hcrat` at radius `r`; the crater part is `hcrat` minus this.
///
/// Inside the rim, the upper part of the rim (above a narrower Gaussian of height `(1 - frac_ejrim) * rim_height`) is ejecta, up
/// to the ejecta rim height `frac_ejrim * rim_height`. Outside the rim, the ejecta part is the ejecta blanket. The clamp is
/// smoothed and the two forms are blended across the rim, so both parts of the profile are smooth and they always sum to
/// `hcrat`.
#[inline]
fn smooth_crater_ejecta_split(
    r: f64,
    hcrat: f64,
    radius: f64,
    rim_width: f64,
    rim_height: f64,
    frac_ejrim: f64,
    ejprofile: f64,
) -> f64 {
    let he = frac_ejrim * rim_height;
    let hup = rimfunc(
        r,
        radius,
        (1.0 - frac_ejrim) * rim_height,
        0.0,
        0.5 * rim_width,
    );
    let inside = smooth_clamp(hcrat - hup, 0.0, he, SPLIT_CLAMP_SMOOTHING * he);
    let outside = ejecta_rim_function(r, radius, frac_ejrim, rim_height, ejprofile, rim_width);
    blend_step(r, inside, outside, radius, SPLIT_RIM_BLEND * rim_width)
}

#[inline]
fn elevation_offset_func(crater: &BasicMoonCrater, ring: &BasicMoonCrater) -> f64 {
    // Returns the offset in elevation for the ring, which is taken to be the elevation of the crater at the ring radius.
    crater_elevation(crater, ring.radius)
}

/// Computes the integrand of the axisymmetric volume integral, r * h(r), at an array of radial distances.
///
/// This is used to estimate crater and ejecta volumes by quadrature. The profile is evaluated with the crater's
/// parameters passed through unmodified, exactly as the scalar `crater_profile_function` and `ejecta_profile_function`
/// Python bindings do.
///
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from the crater center (in meters).
/// * `crater` - The crater whose profile is integrated.
/// * `include_crater` - Include the crater (excavation and rim) part of the profile.
/// * `include_ejecta` - Include the ejecta part of the profile.
///
/// # Returns
///
/// * `r * (h_crater - h_ejecta)` if only the crater is included, `r * h_ejecta` if only the ejecta is included,
///   `r * h_crater` (which already contains the ejecta) if both are included, and zeros if neither is.
pub fn volume_integrand(
    radial_distances: ArrayView1<'_, f64>,
    crater: &BasicMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> Array1<f64> {
    let hcrater = |r: f64| {
        crater_profile_function(
            r,
            crater.radius,
            crater.floor_elevation,
            crater.floor_radius,
            crater.wall_curvature,
            crater.floor_blend,
            crater.rim_width,
            crater.rim_height,
            crater.frac_ejrim,
            crater.ejprofile,
            crater.peak_height,
            crater.peak_width,
            crater.peak_ring_radius,
        )
    };
    let hejecta = |r: f64| {
        ejecta_profile_function(
            r,
            crater.radius,
            crater.frac_ejrim,
            crater.rim_height,
            crater.ejprofile,
            crater.rim_width,
        )
    };
    radial_distances.mapv(|r| match (include_crater, include_ejecta) {
        (true, false) => r * (hcrater(r) - hejecta(r)),
        (false, true) => r * hejecta(r),
        (true, true) => r * hcrater(r),
        (false, false) => 0.0,
    })
}

/// Calculates the elevation of a crater as a function of distance from the center.
///
/// This function applies a polynomial profile for the crater interior (r < 1.0) and a rim dropoff
/// function for the exterior (r ≥ 1.0). It is based on the crater profile model described in
/// Fassett and Thomson (2014).
///
/// Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and
/// the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271.
/// https://doi.org/10.1002/2014JE004698
///
/// This function is split off from the `profile` function for clarity, and is not intended to be
/// called directly.
///
/// # Arguments
///
/// * `r` - Radial distance from crater center
/// * `radius` - Radius of the crater rim (in meters).
/// * `hf` - Elevation of the crater floor relative to the reference plane.
/// * `rf` - Radius of the crater floor.
/// * `wc` - Wall curvature parameter: 0 is straight walls, 1 is very curvy walls.
/// * `fb` - Floor blend parameter: 0 is abrupt floor-to-wall transition, 1 is very gentle floor-to-wall transition
/// * `rw` - Width of the crater rim
/// * `hr` - Height of the crater rim above the reference plane.
/// * `fe` - Fraction of the rim that is made of ejecta.
/// * `pej` - Exponent for the ejecta dropoff function.
/// * `hc` - Height of the central peak above the floor.
/// * `rc` - Width of the central peak.
/// * `ro` - Radial offset of the central peak from the crater center .
///
/// # Returns
///
/// * Elevation values at distance `r`.
#[inline]
pub fn crater_profile_function(
    r: f64,
    radius: f64,
    hf: f64,
    rf: f64,
    wc: f64,
    fb: f64,
    rw: f64,
    hr: f64,
    fe: f64,
    pej: f64,
    hc: f64,
    rc: f64,
    ro: f64,
) -> f64 {
    // Every piece below is smooth and the pieces are joined so that elevation, slope, and curvature are continuous, so the
    // profile has no breaks. (A smootherstep blend is only as smooth as the functions it blends, so none of them may switch
    // formulas or be clipped inside a blend zone.)
    let beta = (1.0 + wc * 9.0).clamp(1.0, 10.0);
    let he = fe * hr;

    // The rim and ejecta, blended in from the interior across [radius - rw, radius + rw]. Beyond that zone only this term is
    // needed (and the wall formula's exponentials would eventually overflow far from the crater).
    let phi_rim = smootherstep((r - (radius - rw)) / (2.0 * rw));
    let hrim =
        rimfunc(r, radius, hr - he, 0.0, rw) + ejecta_rim_function(r, radius, fe, hr, pej, rw);
    if phi_rim >= 1.0 {
        return hrim;
    }

    // The interior: a flat floor, a quintic fillet across the floor blend zone [rf - rfw, rf + rfw] that matches the
    // elevation, slope, and curvature of the floor at its inner end and of the wall at its outer end, then the wall. The slope
    // rises steadily through the fillet, so the floor-to-wall transition has no shoulders. Any central peak is added on top.
    let rfw = (if rf < 0.5 * radius { rf } else { radius - rf } * fb).max(1e-8 * radius);
    let (fillet_start, fillet_end) = (rf - rfw, rf + rfw);
    let base = if r <= fillet_start {
        hf
    } else if r < fillet_end {
        let (w_end, slope_end, curv_end) =
            wallfunc_derivatives(fillet_end, radius, rf, hr, hf, beta);
        quintic_hermite(
            (r - fillet_start) / (fillet_end - fillet_start),
            fillet_end - fillet_start,
            (hf, 0.0, 0.0),
            (w_end, slope_end, curv_end),
        )
    } else {
        wallfunc(r, radius, rf, hr, hf, beta)
    };
    let interior = base + floorfunc(r, rc, hc, ro, 0.0);

    // The wall formula keeps rising past the rim, so cap the interior at the rim height with a smooth minimum whose radial
    // width is a fraction of the rim width. Without it the rim blend would carry the rising wall outward, moving the crest
    // off the rim; with a hard cap the crest would be a corner. The cap level is raised slightly so the crest height is kept.
    let (_, wall_slope_at_rim, _) = wallfunc_derivatives(radius, radius, rf, hr, hf, beta);
    let cap_width = RIM_CAP_WIDTH * rw * wall_slope_at_rim.abs();
    let cap_level = hr + 0.15625 * cap_width;
    let hwall = interior - smooth_positive_part(interior - cap_level, cap_width);

    (1.0 - phi_rim) * hwall + phi_rim * hrim
}

/// Radial half-width of the smooth cap that limits the wall to the rim height, as a fraction of the rim width. It sets how
/// gradually the wall rounds over into the rim crest.
const RIM_CAP_WIDTH: f64 = 0.05;

/// Quintic Hermite interpolation on `u` in [0, 1] over an interval of length `length`, matching the (value, slope, curvature)
/// given at each end, so that a curve joined with it is continuous in elevation, slope, and curvature.
#[inline]
fn quintic_hermite(u: f64, length: f64, start: (f64, f64, f64), end: (f64, f64, f64)) -> f64 {
    let (h0, s0, k0) = start;
    let (h1, s1, k1) = end;
    let (u2, u3) = (u * u, u * u * u);
    let (u4, u5) = (u3 * u, u3 * u2);
    let h_basis0 = 1.0 - 10.0 * u3 + 15.0 * u4 - 6.0 * u5;
    let h_basis1 = 10.0 * u3 - 15.0 * u4 + 6.0 * u5;
    let s_basis0 = u - 6.0 * u3 + 8.0 * u4 - 3.0 * u5;
    let s_basis1 = -4.0 * u3 + 7.0 * u4 - 3.0 * u5;
    let k_basis0 = 0.5 * u2 - 1.5 * u3 + 1.5 * u4 - 0.5 * u5;
    let k_basis1 = 0.5 * u3 - u4 + 0.5 * u5;
    h0 * h_basis0
        + h1 * h_basis1
        + length * (s0 * s_basis0 + s1 * s_basis1)
        + length * length * (k0 * k_basis0 + k1 * k_basis1)
}

#[inline]
fn floorfunc(r: f64, rc: f64, hc: f64, ro: f64, hf: f64) -> f64 {
    if rc > f64::EPSILON {
        hc * (-((r - ro) / rc).powi(2)).exp() + hf
    } else {
        hf
    }
}

/// The crater wall, rising from `hf` at the floor radius `rf` to `hr` at the rim `radius`.
///
/// The formula is evaluated analytically on both sides of that range (it is smooth everywhere), because it is blended with
/// the floor and rim terms across zones that extend past `rf` and `radius`.
#[inline]
fn wallfunc(r: f64, radius: f64, rf: f64, hr: f64, hf: f64, beta: f64) -> f64 {
    let r0 = (r - rf) / (radius - rf);
    let c = (hr - hf) * ((-beta / 2.0).exp() + 1.0) / (beta.exp() - 1.0);
    (c * ((beta * r0).exp() - beta.exp()) / (1.0 + (beta * (r0 - 0.5)).exp())) + hr
}

/// The wall elevation and its first and second derivatives with respect to `r`, computed analytically.
#[inline]
fn wallfunc_derivatives(
    r: f64,
    radius: f64,
    rf: f64,
    hr: f64,
    hf: f64,
    beta: f64,
) -> (f64, f64, f64) {
    let span = radius - rf;
    let x = (r - rf) / span;
    let c = (hr - hf) * ((-beta / 2.0).exp() + 1.0) / (beta.exp() - 1.0);
    let ex = (beta * x).exp();
    let n = ex - beta.exp();
    let (n1, n2) = (beta * ex, beta * beta * ex);
    let ed = (beta * (x - 0.5)).exp();
    let d = 1.0 + ed;
    let (d1, d2) = (beta * ed, beta * beta * ed);
    let q = n / d;
    let q1 = (n1 * d - n * d1) / (d * d);
    let q2 = (n2 * d - n * d2) / (d * d) - 2.0 * d1 * (n1 * d - n * d1) / (d * d * d);
    (c * q + hr, c * q1 / span, c * q2 / (span * span))
}

#[inline]
fn rimfunc(r: f64, radius: f64, hr: f64, he: f64, rw: f64) -> f64 {
    let hmax = if r > radius { hr - he } else { hr };
    hmax * (-(r - radius).powi(2) / (2.0 * rw.powi(2))).exp()
}

/// Computes the ejecta profile scaling at a given radial distance.
///
/// This function returns a value based on a power-law decay that modifies the
/// ejecta distribution intensity outside the crater rim. Values inside the crater return zero.
///
/// # Arguments
///
/// * `r` - Radial distance from the crater center (in meters).
/// * `radius` - Radius of the crater (in meters).
/// * `ejrim` - Rim elevation parameter used to scale the profile.
/// *  `ejprofile` - Exponent for the power-law decay of the ejecta profile.
///
/// # Returns
///
/// * Scaled profile value representing the ejecta contribution at distance `r_actual`.
#[inline]
pub fn ejecta_profile_function(
    r: f64,
    radius: f64,
    fe: f64,
    hr: f64,
    ejprofile: f64,
    rw: f64,
) -> f64 {
    if r > radius {
        ejecta_rim_function(r, radius, fe, hr, ejprofile, rw)
    } else {
        0.0
    }
}

/// The ejecta part of the rim, defined smoothly on both sides of the rim.
///
/// It blends a Gaussian rim of height `fe * hr` into the power-law ejecta blanket across [radius - rw, radius + rw]. Outside
/// the rim it is the ejecta thickness (`ejecta_profile_function` is this function cut off at the rim). Inside the rim it is
/// the ejecta portion of the rim's height, which is needed so that the rim and the crater/ejecta split stay smooth across r =
/// radius.
#[inline]
fn ejecta_rim_function(r: f64, radius: f64, fe: f64, hr: f64, ejprofile: f64, rw: f64) -> f64 {
    let phi = smootherstep((r - (radius - rw)) / (2.0 * rw));
    let ejrim = fe * hr;
    let hg = rimfunc(r, radius, ejrim, 0.0, rw);
    if phi <= 0.0 {
        return hg;
    }

    // Guard division by zero or exceptionally small numbers close to center
    let r_ratio = (r / radius).max(1e-12);
    let hp = ejrim * r_ratio.powf(ejprofile);

    // Safeguard against Infinity and NaN propagation
    let hp_safe = if hp.is_finite() { hp } else { 0.0 };

    (1.0 - phi) * hg + phi * hp_safe
}

/// Computes the ray intensity contribution for a single ejecta point.
///
/// # Arguments
///
/// * `r` - Radial distance from the crater center (in meters).
/// * `bearing0` - Initial bearing angle (radians).
/// * `crater_radius` - Radius of the crater (in meters).
/// * `rmin` - Minimum normalized radial distance.
/// * `rmax` - Maximum normalized radial distance.
/// * `ray_center_bearingi` - Precomputed ray azimuth angles (in radians).
/// * `random_numbers` - Set of random numbers used for pattern variation.
///
/// # Returns
///
/// * Computed ray intensity value for this point.
fn ray_intensity_point(
    r: f64,
    bearing0: f64,
    crater_radius: f64,
    rmin: f64,
    rmax: f64,
    ray_center_bearingi: &[f64],
    random_numbers: &[f64],
) -> f64 {
    if r < crater_radius {
        return 0.0;
    }
    let mut sum = 0.0;
    for j in 0..NPATT as usize {
        let rn = random_numbers[j];
        let bearing = (bearing0 + rn * TAU) % TAU;
        let r_pattern = r / crater_radius - rn;
        sum += FRAYREDUCTION.powi(j as i32 - 1)
            * ray_intensity_func(
                r_pattern,
                bearing,
                rmin,
                rmax,
                ray_center_bearingi,
                2.348 * crater_radius.powf(0.006),
            );
    }
    sum
}

/// Computes the ray-modulated ejecta intensity values for a set of input points.
///
/// Each intensity value is the sum of contributions from multiple rays, modulated by random
/// angular perturbations and decaying with distance.
///
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances (meters).
/// * `bearings` - 1D array of initial bearing angles (radians, clockwise north).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A vector of normalized ray-modulated intensity values.
pub fn ray_intensity(
    radial_distances: ArrayView1<'_, f64>,
    bearings: ArrayView1<'_, f64>,
    crater_diameter: f64,
    seed: u64,
) -> ArrayResult {
    if radial_distances.len() != bearings.len() {
        return Err("radial_distances and reference_elevations must have same length".into());
    }
    let crater_radius = crater_diameter / 2.0;
    let rmax = 100.0;
    let rmin = 1.0;

    let mut rng = ChaCha12Rng::seed_from_u64(seed);

    // Distribute ray patterns evenly around the crater
    let mut ray_center_bearingi = (0..NRAYMAX)
        .map(|i| PI * 2.0 * (i + 1) as f64 / NRAYMAX as f64)
        .collect_vec();
    ray_center_bearingi.shuffle(&mut rng);

    let random_numbers: Vec<f64> = rng.random_iter().take(NPATT as usize).collect_vec();

    let mut intensity: Vec<f64> = (0..radial_distances.len())
        .into_par_iter()
        .map(|i| {
            ray_intensity_point(
                radial_distances[i],
                bearings[i],
                crater_radius,
                rmin,
                rmax,
                &ray_center_bearingi,
                &random_numbers,
            )
        })
        .collect();
    let max_val = intensity
        .iter()
        .zip(radial_distances.iter())
        .filter_map(|(&val, &r)| if r > crater_radius { Some(val) } else { None })
        .fold(f64::MIN, |a, b| a.max(b));
    intensity = intensity
        .into_iter()
        .zip(radial_distances.iter())
        .map(|(intensity, &r)| {
            if r > crater_radius {
                intensity / max_val
            } else {
                intensity
            }
        })
        .collect_vec();
    Ok(Array1::from_vec(intensity))
}

/// Computes the intensity contribution of all rays for a single radial/angular location.
///
/// Applies a distance-dependent decay of ray count, and for each ray evaluates its contribution
/// based on length and angular width using a modified Gaussian model.
///
/// # Arguments
///
/// * `r` - Normalized radial distance (unitless).
/// * `bearing` - Angle from crater center (radians).
/// * `rmin` - Minimum normalized radius.
/// * `rmax` - Maximum normalized radius (truncation).
/// * `ray_center_bearingi` - Ray azimuth directions.
/// * `minray` - Minimum ejecta blanket radius (scaling for ray length).
///
/// # Returns
///
/// * Total ray intensity contribution at the point.
fn ray_intensity_func(
    r: f64,
    bearing: f64,
    rmin: f64,
    rmax: f64,
    ray_center_bearingi: &[f64],
    minray: f64,
) -> f64 {
    const RAYP: f64 = 4.0; // Power law exponent for ray number decay
    const RAY_WIDTH_EXPONENT: f64 = 1.25; // Exponent for angular width decay
    if !r.is_finite() || r <= 0.0 || r > rmax {
        return 0.0;
    } else if r < 1.0 {
        1.0
    } else {
        let tmp = (NRAYMAX as f64).powf(RAYP)
            - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln();
        let n = if tmp < 0.0 {
            NRAYMAX // "Nrays" in Minton et al. (2019)
        } else {
            // Exponential decay of ray number with distance
            ((((NRAYMAX as f64).powf(RAYP)
                - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln())
            .powf(1.0 / RAYP))
            .floor() as i32)
                .clamp(1, NRAYMAX)
        };
        (0..NRAYMAX)
            .into_par_iter()
            .map(|i| {
                let length = minray
                    * ((rmax / minray).ln()
                        * ((((NRAYMAX - i + 2) as f64).powf(RAYP) - 1.0)
                            / ((NRAYMAX as f64).powf(RAYP) - 1.0)))
                        .exp();
                if r > length {
                    0.0 // Skip this ray contribution
                } else {
                    let w = (rmax / length).powf(1.0);
                    let rw = PI / (w * NRAYMAX as f64)
                        * (rmin / r)
                        * (1.0
                            - (1.0 - w / rmin) * (1.0 - (r / rmin).powf(RAY_WIDTH_EXPONENT)).exp());
                    ejecta_ray_func(bearing, ray_center_bearingi[i as usize], r, n, rw)
                }
            })
            .sum()
    }
}

/// Computes the contribution of a single ray to the ejecta field using a Gaussian profile.
///
/// The ray is centered at `ray_center_bearing` with angular width `w`, evaluated at angle `bearing`.
/// The function ensures normalization across rays and returns zero for values below machine epsilon.
///
/// # Arguments
///
/// * `bearing` - Sample angle (radians).
/// * `ray_center_bearing` - Ray center angle (radians).
/// * `r` - Radial distance (meters).
/// * `n` - Number of rays contributing at this distance.
/// * `w` - Angular width of the ray.
///
/// # Returns
///
/// * Ray intensity value at the specified angle and distance.
fn ejecta_ray_func(bearing: f64, ray_center_bearing: f64, r: f64, n: i32, w: f64) -> f64 {
    let c = w / r;
    let b = ray_center_bearing;
    let dbearing = f64::min(2.0 * PI - (bearing - b).abs(), (bearing - b).abs());
    let logval = -dbearing.powi(2) / (2.0 * c.powi(2));
    if logval < crate::VSMALL.ln() {
        0.0
    } else {
        let denom = erf(PI / (2.0 * SQRT_2 * c));
        if denom.abs() <= crate::VSMALL {
            return 0.0;
        }
        let a = (2.0 * PI).sqrt() / (n as f64 * c * denom);
        a * logval.exp()
    }
}

/// Computes the maximum extent of a crater or its ejecta blanket: the radius where the rim (`ejecta == false`) or ejecta
/// (`ejecta == true`) profile falls to `minimum_thickness`.
///
/// The radius is searched for in `[lower_limit, upper_limit]`. If the profile is already below the threshold at the lower limit,
/// the lower limit is returned; if it is still above it at the upper limit, the upper limit is returned. Otherwise the crossing
/// is found with [`brentq`], a port of `scipy.optimize.root_scalar(..., method="brentq")`. If the search does not converge, the
/// crater radius is returned.
///
/// # Arguments
///
/// * `crater` - The crater.
/// * `minimum_thickness` - The threshold height (m).
/// * `ejecta` - Use the ejecta profile if true, or the rim (crater) profile if false.
/// * `lower_limit`, `upper_limit` - The search range (m).
///
/// # Returns
///
/// The maximum extent (m), or an error if the profile is NaN or does not change sign across the search range.
pub fn rmax(
    crater: &BasicMoonCrater,
    minimum_thickness: f64,
    ejecta: bool,
    lower_limit: f64,
    upper_limit: f64,
) -> Result<f64, String> {
    let profile = |r: f64| -> f64 {
        if ejecta {
            ejecta_profile_function(
                r,
                crater.radius,
                crater.frac_ejrim,
                crater.rim_height,
                crater.ejprofile,
                crater.rim_width,
            )
        } else {
            crater_profile_function(
                r,
                crater.radius,
                crater.floor_elevation,
                crater.floor_radius,
                crater.wall_curvature,
                crater.floor_blend,
                crater.rim_width,
                crater.rim_height,
                crater.frac_ejrim,
                crater.ejprofile,
                crater.peak_height,
                crater.peak_width,
                crater.peak_ring_radius,
            )
        }
    };
    let invert = |r: f64| profile(r) - minimum_thickness;

    if invert(lower_limit) < 0.0 {
        return Ok(lower_limit);
    }
    if invert(upper_limit) > 0.0 {
        return Ok(upper_limit);
    }
    match brentq(
        |r| Ok(invert(r)),
        lower_limit,
        upper_limit,
        BRENTQ_XTOL,
        BRENTQ_RTOL,
        BRENTQ_MAXITER,
    )? {
        BrentqResult::Converged(r) => Ok(r),
        BrentqResult::NotConverged(_) => Ok(crater.radius),
        BrentqResult::SignError => Err("f(a) and f(b) must have different signs".to_string()),
    }
}

/// Computes the degradation function, which defines the topographic degradation that each crater contributes to the surface.
///
/// This combines the mare-scale model of Minton et al. (2019) for small craters with the highlands-scale model of
/// Riedel et al. (2020) for large craters, joined by a smoothly broken power law at a breakpoint radius of 500 m.
///
/// The arithmetic mirrors the original Python implementation term by term (same operation order and `powf` calls) so that
/// results match it bit for bit.
///
/// # Arguments
///
/// * `diameter` - Final crater diameter (meters).
/// * `fe` - Degradation function size factor.
///
/// # Returns
///
/// * The degradation state (diffusivity times time, m²) contributed by the crater.
pub fn degradation_function(diameter: f64, fe: f64) -> f64 {
    let kd1 = degradation_breakpoint_coefficient(fe);
    degradation_smooth_broken(diameter / 2.0, kd1)
}

/// Vectorized form of [`degradation_function`] over an array of diameters.
///
/// # Arguments
///
/// * `diameter` - 1D array of final crater diameters (meters).
/// * `fe` - Degradation function size factor.
///
/// # Returns
///
/// * A 1D array of degradation states, one per diameter.
pub fn degradation_function_array(diameter: ArrayView1<'_, f64>, fe: f64) -> Array1<f64> {
    let kd1 = degradation_breakpoint_coefficient(fe);
    diameter.mapv(|d| degradation_smooth_broken(d / 2.0, kd1))
}

// Breakpoint radius (m) between the mare-scale and highlands-scale degradation power laws.
const DEG_RB: f64 = 0.5e3;
// Mare-scale and highlands-scale power-law exponents.
const DEG_PSI_1: f64 = 2.0;
const DEG_PSI_2: f64 = 1.2;
// Smoothing width of the broken power law.
const DEG_DELTA: f64 = 1.0e0;

/// The amplitude of the smoothly broken power law, derived from the mare-scale degradation function of Minton et al. (2019)
/// eq. (32) evaluated at the breakpoint radius.
fn degradation_breakpoint_coefficient(fe: f64) -> f64 {
    let kv1 = 0.30;
    let neq1 = 0.0084;
    let eta = 3.2;
    let gamma = 2.0;
    let beta = 2.0;
    let kd1 = kv1
        * (PI * fe.powf(2.0) * neq1 * (gamma * beta / ((eta - 2.0) * (beta + gamma - eta))))
            .powf(-gamma / (eta - beta));
    let psi = gamma * ((eta - 2.0) / (eta - beta));
    let kdmare_rb = kd1 * DEG_RB.powf(psi);
    kdmare_rb / (1.0 + (DEG_PSI_1 - DEG_PSI_2) / DEG_PSI_1).powf(2.0)
}

/// Smoothly broken power law in crater radius `r` with amplitude `a`.
#[inline]
fn degradation_smooth_broken(r: f64, a: f64) -> f64 {
    let xr = r / DEG_RB;
    a * xr.powf(DEG_PSI_1)
        * (0.5 * (1.0 + xr.powf(1.0 / DEG_DELTA))).powf((DEG_PSI_2 - DEG_PSI_1) / DEG_DELTA)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degradation_function_scalar_matches_array() {
        let d = Array1::from(vec![1e-8, 1.0, 999.9, 1000.0, 1000.1, 2.5e4, 3.0e5]);
        let arr = degradation_function_array(d.view(), 100.0);
        for (i, &di) in d.iter().enumerate() {
            assert_eq!(arr[i].to_bits(), degradation_function(di, 100.0).to_bits());
        }
    }

    #[test]
    fn degradation_function_is_power_law_far_from_break() {
        // Far below the break the slope in log-log space tends to psi_1, far above it tends to psi_2.
        let slope = |d1: f64, d2: f64| {
            (degradation_function(d2, 100.0) / degradation_function(d1, 100.0)).ln()
                / (d2 / d1).ln()
        };
        assert!((slope(1e-3, 2e-3) - DEG_PSI_1).abs() < 1e-5);
        assert!((slope(1e8, 2e8) - DEG_PSI_2).abs() < 1e-5);
    }

    #[test]
    fn wallfunc_derivatives_match_finite_differences() {
        let (radius, rf, hr, hf) = (1000.0, 180.0, 75.0, -350.0);
        for &beta in &[1.0, 2.8, 5.0, 10.0] {
            for &r in &[0.0, 90.0, 180.0, 260.0, 500.0, 900.0, 1000.0, 1300.0] {
                let (w, w1, w2) = wallfunc_derivatives(r, radius, rf, hr, hf, beta);
                let h = 1e-3;
                let (wp, wm) = (
                    wallfunc(r + h, radius, rf, hr, hf, beta),
                    wallfunc(r - h, radius, rf, hr, hf, beta),
                );
                assert!((w - wallfunc(r, radius, rf, hr, hf, beta)).abs() < 1e-9);
                assert!(
                    (w1 - (wp - wm) / (2.0 * h)).abs() < 1e-6 * (1.0 + w1.abs()),
                    "slope beta={beta} r={r}"
                );
                assert!(
                    (w2 - (wp - 2.0 * w + wm) / (h * h)).abs() < 1e-3 * (1.0 + w2.abs()),
                    "curvature beta={beta} r={r}"
                );
            }
        }
    }

    #[test]
    fn quintic_hermite_matches_end_conditions() {
        let (length, start, end) = (170.0, (-350.0, 0.0, 0.0), (-310.0, 0.42, 0.0013));
        let f = |u: f64| quintic_hermite(u, length, start, end);
        let du = 1e-5;
        for (u, (h, s, k)) in [(0.0, start), (1.0, end)] {
            let slope = (f(u + du) - f(u - du)) / (2.0 * du * length);
            let curv = (f(u + du) - 2.0 * f(u) + f(u - du)) / (du * du * length * length);
            assert!((f(u) - h).abs() < 1e-9);
            assert!((slope - s).abs() < 1e-6, "slope at u={u}: {slope} vs {s}");
            assert!((curv - k).abs() < 1e-4, "curvature at u={u}: {curv} vs {k}");
        }
    }

    #[test]
    fn simple_crater_profile_is_smooth() {
        // No breaks in slope anywhere: the change in slope between neighbouring samples stays at the level set by curvature
        let crater = BasicMoonCrater {
            diameter: 2000.0,
            radius: 1000.0,
            floor_elevation: -347.0,
            floor_radius: 151.0,
            wall_curvature: 0.2,
            floor_blend: 0.56,
            rim_width: 354.0,
            rim_height: 70.7,
            frac_ejrim: 0.61,
            ejprofile: -3.58,
            peak_height: 0.0,
            peak_width: 0.0,
            peak_ring_radius: 0.0,
            peak_center_distance: 0.0,
            peak_center_bearing: 0.0,
        };
        let dr = 0.05;
        for (include_crater, include_ejecta) in [(true, true), (true, false), (false, true)] {
            let h: Vec<f64> = (0..60000)
                .map(|i| {
                    basicmoon_profile_one(
                        i as f64 * dr,
                        0.0,
                        &crater,
                        &None,
                        include_crater,
                        include_ejecta,
                    )
                })
                .collect();
            let slope: Vec<f64> = h.windows(2).map(|w| (w[1] - w[0]) / dr).collect();
            let max_step = slope
                .windows(2)
                .map(|w| (w[1] - w[0]).abs())
                .fold(0.0, f64::max);
            assert!(
                max_step < 1e-3,
                "slope jumps by {max_step} between samples ({include_crater}, {include_ejecta})"
            );
        }
    }
}
