pub mod basicmoon;
pub mod central_peak;
pub mod form_crater;
pub mod realmoon;

// The second-order smootherstep function of Ken Perlin.
//
//
/// # Arguments
///
/// * `x` - The parameter to smooth.
///
/// # Returns
///
/// * 0 for x<0 and 1 for x > 1 with a smooth and symmetric transition between them.
///
#[inline]
pub fn smootherstep(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if x <= 1.0 {
        6.0 * x.powi(5) - 15.0 * x.powi(4) + 10.0 * x.powi(3)
    } else {
        1.0
    }
}

// Blends two values across a single step transition boundary.
//
//
/// # Arguments
///
/// * `x` - The evaluation distance or coordinate.
/// * `yinner` - The value returned deep inside the inner zone (where r <= xtransition - xwidth).
/// * `youter` - The value returned deep in the outer zone (where r >= xtransition + xwidth).
/// * `xtransition` - The location of the step transition boundary midpoint.
/// * `xwidth` - Half-width of the transition zone across which interpolation occurs.
///
/// # Returns
///
/// * A smoothly blended value transition from yinner to youter using the second-order smootherstep function.
///
#[inline]
pub fn blend_step(x: f64, yinner: f64, youter: f64, xtransition: f64, xwidth: f64) -> f64 {
    let t = (x - (xtransition - xwidth)) / (2.0 * xwidth);
    let phi = smootherstep(t);
    (1.0 - phi) * yinner + phi * youter
}

// A smooth (C3) version of the positive part max(x, 0).
//
// Its derivative is the smootherstep function across [-w, w], so it is exactly 0 for x <= -w, exactly x for x >= w, and never
// below max(x, 0) in between (it equals 0.15625 * w at x = 0).
//
/// # Arguments
///
/// * `x` - The value to take the positive part of.
/// * `w` - Half-width of the smoothing zone around 0. If `w <= 0`, this is exactly `max(x, 0)`.
///
/// # Returns
///
/// * The smoothed positive part of `x`.
#[inline]
pub fn smooth_positive_part(x: f64, w: f64) -> f64 {
    if w <= 0.0 {
        return x.max(0.0);
    }
    if x <= -w {
        0.0
    } else if x >= w {
        x
    } else {
        // Integral of smootherstep(t) from 0 to t, scaled to the zone width
        let t = (x + w) / (2.0 * w);
        2.0 * w * (t.powi(6) - 3.0 * t.powi(5) + 2.5 * t.powi(4))
    }
}

// A smooth (C3) version of clamp(x, lo, hi), built from two smoothed positive parts.
//
// It is exactly `lo` for x <= lo - w, exactly `hi` for x >= hi + w, and exactly x for lo + w <= x <= hi - w.
//
/// # Arguments
///
/// * `x` - The value to clamp.
/// * `lo`, `hi` - The clamp bounds (`lo <= hi`).
/// * `w` - Half-width of the smoothing zone at each bound. If `w <= 0`, this is exactly `clamp(x, lo, hi)`.
///
/// # Returns
///
/// * The smoothly clamped value.
#[inline]
pub fn smooth_clamp(x: f64, lo: f64, hi: f64, w: f64) -> f64 {
    lo + smooth_positive_part(x - lo, w) - smooth_positive_part(x - hi, w)
}
