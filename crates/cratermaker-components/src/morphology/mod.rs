pub mod basicmoon;
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


// Blends a perturbed value to a background mean value wityin an annular band near a boundary.
//
//
/// # Arguments
///
/// * `r` - The evaluation distance or coordinate.
/// * `yinband` - The perturbed value returned exactly on the boundary (where r == xtransition).
/// * `youtband` - The unperturbed background value returned far from the boundary (where |r - xtransition| >= xwidth).
/// * `xtransition` - The location of the boundary midpoint around which the band is centered.
/// * `xwidth` - The total distance from the boundary at which the perturbation decays completely to the mean value.
///
/// # Returns
///
/// * A smoothly blended value that retains the yinband perturbation on the boundary and transitions symmetrically to youtband away from it.
///
#[inline]
pub fn blend_band(x: f64, yinband: f64, youtband: f64, xtransition: f64, xwidth: f64) -> f64 {
    let dist = (x - xtransition).abs();
    let t = (dist / xwidth).clamp(0.0, 1.0);
    let phi = smootherstep(t);
    (1.0 - phi) * yinband + phi * youtband
}