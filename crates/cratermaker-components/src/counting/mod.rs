use std::f64::consts::{FRAC_PI_2, PI, TAU};
use numpy::ndarray::prelude::*; 
use ndarray_linalg::error::LinalgError;
use ndarray_linalg::{Inverse, Eig};
use rayon::prelude::*;
use crate::ArrayResult;
use crate::surface::LocalSurfaceView;
use crate::crater::Crater;

/// Fits a crater rim ellipse to a local surface region.
/// This function performs an iterative fitting procedure to refine
/// the parameters of a crater rim ellipse based on the local surface data.
/// 
/// # Arguments
/// 
/// * `region` - A view of the local surface region containing the crater.
/// * `x0`, `y0` - Initial center of the crater ellipse in crater-centered coordinates.
/// * `crater` - Initial crater parameters.
/// * `score_quantile` - Quantile of rim scores to consider.
/// * `fit_center` - Whether to fit the crater center or keep it fixed.
/// * `fit_ellipse` - Whether to fit an ellipse. If false, a circle is fitted instead. 
/// * `gradmult`, `curvmult`, `heightmult` - Multipliers for gradient, curvature, and height-based scoring.
/// 
/// # Returns
/// 
/// * On success, returns a tuple `(x0_fit, y0_fit, a_fit, b_fit, o_fit)`
///  containing the fitted crater center coordinates, semi-major axis,
/// semi-minor axis, and orientation angle (in radians).
/// 
/// # Errors
/// 
/// * Returns `Err(String)` if any error occurs during the fitting process.
/// 
pub fn fit_one_rim(
    region: &LocalSurfaceView<'_>,
    x0: f64,
    y0: f64,
    crater: &Crater,
    fit_center: bool,
    fit_ellipse: bool,
    score_quantile: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) ->Result<(f64, f64, f64, f64, f64, numpy::ndarray::Array1<f64>), String> { 
    let x = region.face_proj_x.as_ref().ok_or("face_proj_x required")?;
    let y = region.face_proj_y.as_ref().ok_or("face_proj_y required")?;
    let mut x0_fit = x0;
    let mut y0_fit = y0;
    let a_fit: f64;
    let b_fit:f64;
    let o_fit:f64;
    let _wrms: f64;

    // Score the rim using the current multipliers
    let rimscore = score_rim(
        region,
        x0, 
        y0, 
        &crater, 
        score_quantile, 
        gradmult, 
        curvmult, 
        heightmult
    ).map_err(|e| e.to_string())?;

    // Fit an ellipse to the weighted points using either the floating or fixed center fitter, depending on user input
    if fit_ellipse {
        if fit_center {
            (x0_fit, y0_fit, a_fit, b_fit, o_fit, _wrms) = fit_one_ellipse(
                x.view(),
                y.view(),
                rimscore.view()
            ).map_err(|e| e.to_string())?;
        } else {
            (a_fit, b_fit, o_fit, _wrms) = fit_one_ellipse_fixed_center(
                x.view(),
                y.view(),
                rimscore.view(),
                x0,
                y0,
            ).map_err(|e| e.to_string())?;

        }
    } else {
        if fit_center {
            (x0_fit, y0_fit, a_fit, _wrms) = fit_one_circle(
                x.view(),
                y.view(),
                rimscore.view()
            ).map_err(|e| e.to_string())?;
        } else {
            (a_fit, _wrms) = fit_one_circle_fixed_center(
                x.view(),
                y.view(),
                rimscore.view(),
                x0,
                y0,
            ).map_err(|e| e.to_string())?;
        }
        b_fit = a_fit;
        o_fit = 0.0;
    }

    Ok((x0_fit, y0_fit, a_fit, b_fit, o_fit, rimscore))
}

/// Computes the signed radial distance from points to an ellipse.
///
/// For each input point `(x[i], y[i])`, this function:
///
/// 1. Computes the polar coordinates `(r, θ)` of the point relative to the
///    ellipse center `(x0, y0)`.
/// 2. Computes the radius of the ellipse along the same direction, using the
///    semi-major axis `a`, semi-minor axis `b`, and orientation.
/// 3. Returns the signed difference `r - r_ellipse(θ)`.
///
/// By this convention:
///
/// * Values **> 0** indicate points outside the ellipse.
/// * Values **< 0** indicate points inside the ellipse.
/// * Values **≈ 0** lie close to the ellipse boundary.
///
/// The orientation is given in radians. 
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `a` - Semi-major axis of the ellipse.
/// * `b` - Semi-minor axis of the ellipse.
/// * `orientation` - Orientation of the ellipse in radians (see above).
/// * `x0` - x-coordinate of the ellipse center.
/// * `y0` - y-coordinate of the ellipse center.
///
/// # Returns
///
/// On success, returns an array of signed radial distances, one value
/// per input point.
///
/// # Errors
///
/// Currently this function always returns `Ok(...)`. The `ArrayResult`
/// return type allows future extensions to return descriptive errors
/// (for example, if the input arrays have inconsistent lengths).
pub fn radial_distance_to_crater_rim(
    x: &ArrayView1<'_, f64>,
    y: &ArrayView1<'_, f64>,
    crater: &Crater,
    x0: f64, 
    y0:f64
) -> ArrayResult {
    let phi = TAU - FRAC_PI_2 - crater.measured_orientation.to_radians();
    let a = crater.measured_semimajor_axis;
    let b = crater.measured_semiminor_axis;

    let result_vec: Vec<f64> = (0..x.len())
        .into_par_iter()
        .map(|i| {
            let xi = x[i];
            let yi = y[i];
            let dx = xi - x0;
            let dy = yi - y0;
            let r = (dx * dx + dy * dy).sqrt();
            let theta = dy.atan2(dx);
            let alpha = theta - phi;
            let ca = alpha.cos();
            let sa = alpha.sin();
            r - (a * b) / ((b * ca).powi(2) + (a * sa).powi(2)).sqrt()
        })
        .collect();

    Ok(Array1::from(result_vec))
}

/// Converts conic ellipse coefficients to geometric parameters.
///
/// Given the coefficients of a general quadratic form
///
/// `a x² + b x y + c y² + d x + f y + g = 0`
///
/// (note the asymmetric naming: the parameter `f` here corresponds to the
/// `y`-linear coefficient, and `g` is the constant term), this function
/// computes the geometric parameters of the ellipse:
///
/// * Center `(x0, y0)`
/// * Semi-major axis `aₚ`
/// * Semi-minor axis `bₚ`
/// * Eccentricity `eₚ`
/// * Orientation angle `orientation`
///
/// The formulas are based on the reference from MathWorld:
/// “Ellipse” in terms of a Cartesian conic representation, with the
/// conversion adjusted for the coefficient convention used here.
///
/// If the discriminant of the conic indicates that the coefficients do not
/// describe an ellipse (e.g., a hyperbola or parabola), the function
/// will panic.
///
/// # Arguments
///
/// * `a` - Quadratic coefficient for `x²`.
/// * `b` - Quadratic coefficient for `x y`.
/// * `c` - Quadratic coefficient for `y²`.
/// * `d` - Linear coefficient for `x`.
/// * `f` - Linear coefficient for `y`.
/// * `g` - Constant term.
///
/// # Returns
///
/// A tuple `(x0, y0, ap, bp, ep, orientation)` where:
///
/// * `x0`, `y0` - Center of the ellipse.
/// * `ap` - Semi-major axis length.
/// * `bp` - Semi-minor axis length.
/// * `orientation` - Orientation angle of the ellipse in radians,
///   wrapped into `[0, π)`.
///
/// # Panics
///
/// Panics if the provided coefficients do **not** represent an ellipse,
/// i.e., if the conic discriminant satisfies `b² - a c ≥ 0`.
#[inline]
fn ellipse_coefficients_to_parameters(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    f: f64,
    g: f64) -> (f64, f64, f64, f64, f64) {

    // We use the formulas from https://mathworld.wolfram.com/Ellipse.html which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.  
    // Therefore, B, D and F must be scaled.
    let b = b /2.0;
    let d = d / 2.0; 
    let f = f / 2.0;
    let den = b*b - a*c;
    if den > 0.0 {
        panic!("The provided coefficients do not represent an ellipse.");
    }
    let x0 = (c * d - b * f) / den;
    let y0 = (a * f - b * d) / den;
    let numerator = 2.0 * (a * f * f + c * d * d + g * b * b - 2.0 * b * d * f - a * c * g);
    let fac = ((a - c).powi(2) + 4.0 * b * b).sqrt();

    // Compute unsorted semi-major and semi-minor axes
    let _ap = (numerator / (den * (fac - (a + c)))) .sqrt();
    let _bp = (numerator / (den * (-fac - (a + c)))).sqrt();
    // Sort the semi-major and semiminor axis but keep track of orientation
    let (ap, bp, width_gt_height) = if _bp >= _ap {
        (_bp, _ap, true)
    } else {
        (_ap, _bp, false)
    };


    // compute the orientation with respect to the x-axis and deal with special cases
    let mut phi = if b.abs() < 1e-12 {
            if a < c {
                0.0
            } else {
                FRAC_PI_2
            }
        } else {
            0.5 * (2.0 * b).atan2(a - c)
        }; 
    if a.abs() > c.abs() {
        phi += FRAC_PI_2;
    }
    if !width_gt_height {
        phi += FRAC_PI_2;
    }
    let orientation = phi % PI;
    

    (x0, y0, ap, bp, orientation)
}

/// Converts ellipse geometric parameters back to conic coefficients.
///
/// Given the center `(x0, y0)`, semi-axes `a` and `b`, and an orientation
/// angle, this function returns the coefficients `(a, b, c, d, e, f)`
/// of the general quadratic form
///
/// `a x² + b x y + c y² + d x + e y + f = 0`
///
/// corresponding to that ellipse.
///
/// The orientation is given in radians. Internally, the function uses
/// `phi = orientation - π/2` to align the geometric orientation with the
/// coefficient representation used in the other ellipse routines.
///
/// # Arguments
///
/// * `x0` - x-coordinate of the ellipse center.
/// * `y0` - y-coordinate of the ellipse center.
/// * `a` - Semi-major axis.
/// * `b` - Semi-minor axis.
/// * `orientation` - Orientation of the ellipse in radians.
///
/// # Returns
///
/// A tuple `(a, b, c, d, e, f)` such that the ellipse is described by:
///
/// `a x² + b x y + c y² + d x + e y + f = 0`
// #[inline]
// fn ellipse_parameters_to_coefficients(x0: f64, y0: f64, a: f64, b: f64, orientation: f64) -> (f64, f64, f64, f64, f64, f64) {
//     let phi = orientation; 
//     let s2 = phi.sin().powi(2);
//     let c2 = phi.cos().powi(2);
//     let sc = phi.sin() * phi.cos();
//     let a2 = a * a;
//     let b2 = b * b;

//     let a = a2 * s2 + b2 * c2;
//     let b = 2.0 * (b2 - a2) * sc;
//     let c = a2 * c2 + b2 * s2;
//     let d = -2.0 * a * x0 - b * y0;
//     let e = -b * x0 - 2.0 * c * y0;
//     let f = a * x0 * x0 + b * x0 * y0 + c * y0 * y0 - a2 * b2;

//     (a, b, c, d, e, f)
// }


/// Approximates geometric distances from points to an ellipse.
///
/// This function takes the conic coefficients of an ellipse and computes
/// a first-order approximation of the **geometric distance** from each
/// point `(x[i], y[i])` to the ellipse. The approximation is based on the
/// ratio of the ellipse’s implicit function value to the gradient norm,
/// evaluated at each point.
///
/// The coefficients are those of the general quadratic form:
///
/// `a x² + b x y + c y² + d x + e y + f = 0`
///
/// Internally, the coefficients are first converted to ellipse parameters
/// (center, semi-axes, orientation). Then, for each point, the distance
/// is approximated using the local normal direction.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `a` - Quadratic coefficient for `x²`.
/// * `b` - Quadratic coefficient for `x y`.
/// * `c` - Quadratic coefficient for `y²`.
/// * `d` - Linear coefficient for `x`.
/// * `e` - Linear coefficient for `y`.
/// * `f` - Constant term.
///
/// # Returns
///
/// A 1D array of approximate geometric distances, one per point.
/// Positive values indicate points lying “outside” the ellipse along the
/// normal direction; negative values indicate points “inside”.
///
/// # Notes
///
/// This is a **first-order** approximation and is intended for use in
/// robust ellipse fitting (e.g., computing a weighted RMS distance).
#[inline]
fn geometric_distances(
    x: ArrayView1<'_,f64>, 
    y: ArrayView1<'_,f64>, 
    a: f64, 
    b: f64, 
    c: f64, 
    d: f64, 
    e: f64, 
    f: f64
) -> Array1<f64> {
    let (x0, y0, a, b, orientation) = ellipse_coefficients_to_parameters(a, b, c, d, e, f);

    let n = x.len();
    let result_vec: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let xi = x[i];
            let yi = y[i];
            let dx = xi - x0;
            let dy = yi - y0;
            let theta = dy.atan2(dx);
            let alpha = theta - orientation; // + FRAC_PI_2;
            let ca = alpha.cos();
            let sa = alpha.sin();
            let capf = (a * b) / ((b * ca).powi(2) + (a * sa).powi(2)).sqrt();
            let capfx = 2.0 * a * xi + b * yi + d;
            let capfy = b * xi + 2.0 * c * yi + e;
            let gradnorm = (capfx * capfx + capfy * capfy).sqrt();
            capf / gradnorm
        })
        .collect();

    Array1::from(result_vec)
}

/// Computes a per-point “rim score” for a candidate crater rim ellipse.
///
/// This function is intended to evaluate how well a candidate ellipse
/// `(x0, y0, ap, bp, orientation)` matches a crater rim, based on both
/// the geometry of sample points `(x, y)` and the local radial topographic
/// gradient on the mesh.
///
/// Internally, it:
///
/// * Computes the radial gradient of `face_elevation` on the surface mesh,
///   using `compute_radial_gradient`.
/// * (Planned) Uses that gradient together with the projected rim geometry
///   to assign a score to each point, e.g., emphasizing locations where
///   the gradient pattern supports a crater rim interpretation.
///
/// # Arguments
///
/// * `region` - A view of the local surface region containing the crater.`
/// * `x0`, `y0` - Candidate ellipse center in crater-centered coordinates.
/// * `crater` - Candidate crater ellipse parameters.
/// * `quantile` - Quantile threshold for selecting high-scoring points.
/// * `gradmult` - Weighting factor for gradient-based score component.
/// * `curvmult` - Weighting factor for curvature-based score component.
/// * `heightmult` - Weighting factor for height-based score component.
///
/// # Returns
///
/// Currently returns a zero-filled array of length `x.len()`. The function
/// is a placeholder for a future scoring implementation.
///
/// # Panics
///
/// Panics if the input arrays have inconsistent lengths or shapes that
/// violate internal assumptions.
#[inline]
pub fn score_rim(
    region: &LocalSurfaceView<'_>,
    x0: f64, 
    y0: f64, 
    crater: &Crater,
    quantile: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) ->ArrayResult { 
    let n = region.n_face;
    const MIN_POINTS_FOR_FIT: usize = 3; // It will try to use at least this many points per sector in the fit
    const EXTENT_RADIUS_CUTOFF: f64 = 1.5; // Max radial extent as a multiple of crater semi-major axis
    const N_SECTORS: usize = 36; // Number of bearing sectors for scoring
    let sector_width = 360.0 / N_SECTORS as f64;

    let radial_gradient = crate::surface::compute_radial_gradient(
        region.face_elevation.view(),
        region,
    )?; 
    let radial_curvature = crate::surface::compute_radial_gradient(
        radial_gradient.view(),
        region,
    )?;
    let proj_x = region.face_proj_x.as_ref().ok_or("face_proj_x required")?;
    let proj_y = region.face_proj_y.as_ref().ok_or("face_proj_y required")?;
    let ap = crater.measured_semimajor_axis;
    let bp = crater.measured_semiminor_axis;
    let distances = radial_distance_to_crater_rim(proj_x, proj_y, crater, x0, y0)?; 

    // Region mask based on radial distance from origin
    let max_distance = EXTENT_RADIUS_CUTOFF * ap.max(bp);
    let mask_region: Array1<bool> = region.face_distance.as_ref().ok_or("face_distance required")?.mapv(|ri| ri > max_distance);

    // Distance score: closer to current best fit for rim = higher score
    let scale = crater.measured_diameter;
    let mut distscore = distances.mapv(|d| {
        let nd = (d / scale).powi(2);
        1.0 / (nd + 0.05)
    });
    for (v, &mask) in distscore.iter_mut().zip(mask_region.iter()) {
        if mask {
            *v = f64::NAN;
        }
    }
    if let Some(max_d) = nanmax(&distscore) {
        if max_d > 0.0 {
            distscore.map_inplace(|v| {
                if !v.is_nan() {
                    *v /= max_d;
                }
            });
        }
    }

    // Height score: high elevations (relative to mean) score high
    let mut heightscore = region.face_elevation.to_owned();
    if let (Some(min_h), Some(max_h)) = (nanmin(&heightscore), nanmax(&heightscore)) {
        let range = max_h - min_h;
        if range > 0.0 {
            heightscore.map_inplace(|v| {
                if !v.is_nan() {
                    *v = (*v - min_h) / range;
                }
            });
        }
    }

    // Weight the height score by the distance score
    for (v, &dscore) in heightscore.iter_mut().zip(distscore.iter()) {
        if !v.is_nan() {
            *v *= dscore;
        }
    }

    // Gradient score: low absolute gradient scores high (i.e., rim is at a local gradient extremum)
    let mut gradscore = radial_gradient.mapv(|g| g.abs() + 1e-16);
    for (v, &mask) in gradscore.iter_mut().zip(mask_region.iter()) {
        if mask {
            *v = f64::NAN;
        }
    }
    gradscore.map_inplace(|v| {
        if !v.is_nan() {
            *v = 1.0 / *v;
        }
    });

    // Weight the gradient score by the distance score
    for (v, &dscore) in gradscore.iter_mut().zip(distscore.iter()) {
        if !v.is_nan() {
            *v *= dscore;
        }
    }


    if let Some(max_g) = nanmax(&gradscore) {
        if max_g > 0.0 {
            gradscore.map_inplace(|v| {
                if !v.is_nan() {
                    *v /= max_g;
                }
            });
        }
    }

    // Curvature score:
    //    high positive curvature -> low score (i.e., valleys are bad)
    //    high negative curvature -> high score (i.e., ridges are good)
    let mut curvscore = radial_curvature.to_owned();
    for (v, &mask) in curvscore.iter_mut().zip(mask_region.iter()) {
        if mask {
            *v = f64::NAN;
        }
    }
    curvscore.map_inplace(|v| {
        if !v.is_nan() {
            if *v > 0.0 {
                *v = 0.0;
            } else {
                *v = -*v;
            }
        }
    });

    // Weight the curvature score by the distance score
    for (v, &dscore) in curvscore.iter_mut().zip(distscore.iter()) {
        if !v.is_nan() {
            *v *= dscore;
        }
    }

    if let Some(max_c) = nanmax(&curvscore) {
        if max_c > 0.0 {
            curvscore.map_inplace(|v| {
                if !v.is_nan() {
                    *v /= max_c;
                }
            });
        }
    }

    // Combine into rimscore 
    let mut rimscore = Array1::<f64>::zeros(n);
    for i in 0..n {
        let g = gradscore[i];
        let c = curvscore[i];
        let h = heightscore[i];

        if g.is_nan() || c.is_nan() || h.is_nan() {
            rimscore[i] = f64::NAN;
        } else {
            rimscore[i] =
                gradmult * g + curvmult * c + heightmult * h;
        }
    }
    // Bin into bearing sectors:
    let mut rimscore_sector_index = Array1::<usize>::from_elem(n, N_SECTORS);
    let bearing = region.face_bearing.as_ref().ok_or("face_bearing required")?; 
    // First pass: compute max rimscore per sector (ignoring NaNs)
    let mut sector_max = [f64::NEG_INFINITY; N_SECTORS];

    for i in 0..n {
        let v = rimscore[i];
        if v.is_nan() {
            continue;
        }

        let b = bearing[i];
        // bearings should already be in [0, 360), but we clamp just in case
        let mut idx = (b / sector_width).floor() as isize;
        if idx < 0 {
            idx = 0;
        } else if idx >= N_SECTORS as isize {
            idx = N_SECTORS as isize - 1;
        }
        let idx = idx as usize;
        rimscore_sector_index[i] = idx;

        if v > sector_max[idx] {
            sector_max[idx] = v;
        }
    }

    // Second pass: normalize each rimscore by its sector max
    for i in 0..n {
        let v = rimscore[i];
        if v.is_nan() {
            continue;
        }

        let idx = rimscore_sector_index[i];

        let max_v = sector_max[idx];
        if max_v.is_finite() && max_v > 0.0 && rimscore[i] > 0.0 {
            rimscore[i] = v / max_v;
        } else {
            // No valid points in this sector; treat as invalid
            rimscore[i] = f64::NAN;
        }

    }

    // Third pass: Apply quantile threshold to keep only highest scores in each sector
    let mut high_scores = Array1::<bool>::from_elem(n, false);
    for sector in 0..N_SECTORS {
        let sector_indices: Vec<usize> = (0..n)
            .filter(|&i| rimscore_sector_index[i] == sector && !rimscore[i].is_nan())
            .collect();
        if sector_indices.is_empty() {
            continue;
        }
        let mut sector_scores: Vec<f64> = sector_indices
            .iter()
            .map(|&i| rimscore[i])
            .filter(|v| !v.is_nan())
            .collect();
        sector_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let threshold_idx = (quantile * sector_scores.len() as f64) as usize;
        let threshold = sector_scores[threshold_idx];
        for &i in sector_indices.iter() {
            high_scores[i] = rimscore[i] >= threshold;
        }
    }

    let num_high = high_scores.iter().filter(|&&b| b).count();
    let num_valid = rimscore.iter().filter(|v| !v.is_nan()).count();

    if num_high < MIN_POINTS_FOR_FIT * N_SECTORS {
        if num_valid < MIN_POINTS_FOR_FIT * N_SECTORS {
            // use all valid points
            for i in 0..n {
                high_scores[i] = !rimscore[i].is_nan();
            }
        } else {
            // take top MIN_POINTS_FOR_FIT scores
            let mut vals: Vec<f64> = rimscore
                .iter()
                .copied()
                .filter(|v| !v.is_nan())
                .collect();
            vals.sort_by(|a, b| b.partial_cmp(a).unwrap()); // descending order
            let threshold = vals[MIN_POINTS_FOR_FIT * N_SECTORS - 1];

            for i in 0..n {
                let v = rimscore[i];
                high_scores[i] = !v.is_nan() && v >= threshold;
            }
        }
    }

    // Keep only high-score points, others -> NaN
    for i in 0..n {
        if !high_scores[i] {
            rimscore[i] = f64::NAN;
        }
    }

    // Rescale rimscore to [0, 1]
    if let (Some(min_r), Some(max_r)) = (nanmin(&rimscore), nanmax(&rimscore)) {
        let range = max_r - min_r;
        if range > 0.0 {
            rimscore.map_inplace(|v| {
                if !v.is_nan() {
                    *v = (*v - min_r) / range;
                }
            });
        }
    }

    Ok(rimscore)
}

/// Fits a single weighted ellipse to 2D points using direct least squares.
///
/// This function implements a **numerically stable, weighted** variant of the
/// direct least-squares ellipse fitting algorithm of Halír and Flusser:
///
/// > Halír, R. and Flusser, J. (1998).
/// > “Numerically Stable Direct Least Squares Fitting of Ellipses.”
///
/// The ellipse is represented in the general quadratic form:
///
/// `a x² + b x y + c y² + d x + e y + f = 0`
///
/// subject to the constraint that the conic is an ellipse. The algorithm:
///
/// 1. Constructs design matrices `D1` and `D2` from `x` and `y`.
/// 2. Applies per-point weights `w[i]` (via `sqrt(w)` row scaling).
/// 3. Computes scatter matrices `S1`, `S2`, `S3`.
/// 4. Forms the reduced matrix `M` incorporating the ellipse constraint.
/// 5. Solves the eigenproblem for `M` and selects the eigenvector that
///    satisfies the ellipse condition `4ac - b² > 0`.
/// 6. Recovers the full coefficient vector and the geometric parameters
///    (center, semi-axes, orientation).
/// 7. Computes the weighted RMS distance using the internal `geometric_distances` approximation.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `weights` - Non-negative weights for each point. Must have the same
///   length as `x` and `y`. Points with larger weights influence the fit more.
///   Non-positive weights contribute zero influence.
///
/// # Returns
///
/// On success, returns:
///
/// `(x0, y0, ap, bp, orientation, wrms)`
///
/// where:
///
/// * `x0`, `y0` - Center of the fitted ellipse.
/// * `ap` - Semi-major axis length.
/// * `bp` - Semi-minor axis length.
/// * `orientation` - Orientation angle in radians (wrapped to `[0, 2π)`).
/// * `wrms` - Weighted root-mean-square of the (approximate) geometric distances.
///
/// # Errors
///
/// Returns `Err(LinalgError)` if:
///
/// * One of the intermediate matrices (`S3`, constraint matrices, etc.) is singular and
///   cannot be inverted.
/// * The eigen decomposition of `M` fails.
///
/// # Panics
///
/// * Panics if `x`, `y`, and `weights` do not all have the same length.
/// * Panics if no eigenvector satisfies the ellipse condition `4ac - b² > 0`, indicating
///   that the data do not admit a valid ellipse under this formulation.
#[inline] 
pub fn fit_one_ellipse(
    x: ArrayView1<'_,f64>, 
    y: ArrayView1<'_,f64>, 
    weights: ArrayView1<'_,f64>
) -> Result<(f64, f64, f64, f64, f64, f64), LinalgError> {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), weights.len());
    let n = x.len();

    // Step 1: D1, D2
    let x2 = x.mapv(|v| v * v);
    let y2 = y.mapv(|v| v * v);
    let xy = &x * &y;

    let mut d1 = Array2::<f64>::zeros((n, 3));
    let mut d2 = Array2::<f64>::zeros((n, 3));
    d1.column_mut(0).assign(&x2);
    d1.column_mut(1).assign(&xy);
    d1.column_mut(2).assign(&y2);
    d2.column_mut(0).assign(&x);
    d2.column_mut(1).assign(&y);
    d2.column_mut(2).fill(1.0);

    // Step 2: apply weights
    let w_sqrt = weights.mapv(|w| if w > 0.0 { w.sqrt() } else { 0.0 });
    let w_col = w_sqrt.view().insert_axis(Axis(1));
    let wd1 = &w_col * &d1;
    let wd2 = &w_col * &d2;

    // Step 3: S1, S2, S3
    let s1 = wd1.t().dot(&wd1);
    let s2 = wd1.t().dot(&wd2);
    let s3 = wd2.t().dot(&wd2);

    // Step 4: T and M
    let s3_inv = s3.inv()?;
    let t = -s3_inv.dot(&s2.t());
    let m_tmp = s1 + s2.dot(&t);

    let c = arr2(&[
        [0.0, 0.0, 2.0],
        [0.0, -1.0, 0.0],
        [2.0, 0.0, 0.0],
    ]);
    let c_inv = c.inv()?;
    let m = c_inv.dot(&m_tmp);

    // Step 5: eigen-decomposition
    let (_eigvals, eigvecs) = m.eig()?;
    let eigvecs_re = eigvecs.mapv(|c| c.re);

    let mut chosen_col = None;
    for k in 0..3 {
        let a = eigvecs_re[[0, k]];
        let b = eigvecs_re[[1, k]];
        let c = eigvecs_re[[2, k]];
        let con = 4.0 * a * c - b * b;
        if con > 0.0 {
            chosen_col = Some(k);
            break;
        }
    }
    let k = chosen_col.expect("No valid ellipse eigenvector found");
    let ak = eigvecs_re.column(k).to_owned();
    let transposed_ak = t.dot(&ak);

    let a = ak[0];
    let b = ak[1];
    let c = ak[2];
    let d = transposed_ak[0];
    let f = transposed_ak[1];
    let g = transposed_ak[2];

    // Step 6: weighted RMS
    let delta = geometric_distances(x, y, a, b, c, d, f, g);
    let delta2 = delta.mapv(|d| d * d);
    let num = (&delta2 * &weights).sum();
    let den = weights.sum();
    let wrms = (num / den).sqrt();

    let (x0, y0, ap, bp, phi) = ellipse_coefficients_to_parameters(a, b, c, d, f, g,);

    // Convert back to the mapping convention for bearings (0)
    let orientation = (TAU - phi) % TAU;

    Ok((x0, y0, ap, bp, orientation, wrms))
}


/// Fits a single weighted ellipse to 2D points with a fixed center.
///
/// This uses the same Halír–Flusser direct least-squares formulation as `fit_one_ellipse`,
/// but enforces the center to be exactly `(x0, y0)` by reparameterizing the conic and
/// eliminating the linear terms implied by the fixed center.
///
/// The ellipse is represented in the general quadratic form:
///
/// `a x² + b x y + c y² + d x + e y + g = 0`
///
/// with `d` and `e` determined by `a`, `b`, `c` and the fixed center.
///
/// The algorithm:
///
/// 1. Constructs centered features based on the fixed center.
/// 2. Builds design matrices and applies weights via `sqrt(w)` row scaling.
/// 3. Forms the reduced 3×3 eigenproblem for `[a, b, c]`.
/// 4. Solves for the constant term `g` and reconstructs `d` and `e`.
/// 5. Converts the conic to geometric parameters and computes WRMS using `geometric_distances`.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `weights` - Non-negative weights for each point. Must have the same
///   length as `x` and `y`. Points with larger weights influence the fit more.
///   Non-positive weights contribute zero influence.
/// * `x0` - Fixed x-coordinate of the ellipse center.
/// * `y0` - Fixed y-coordinate of the ellipse center.
///
/// # Returns
///
/// On success, returns:
///
/// `(ap, bp, orientation, wrms)`
///
/// where:
///
/// * `ap` - Semi-major axis length.
/// * `bp` - Semi-minor axis length.
/// * `orientation` - Orientation angle in radians (wrapped to `[0, 2π)`).
/// * `wrms` - Weighted root-mean-square of the (approximate) geometric distances.
///
/// # Errors
///
/// Returns `Err(LinalgError)` if:
///
/// * The weighted normal-equation scalar for the eliminated constant term is zero
///   (e.g., all weights are zero).
/// * An intermediate matrix inversion fails.
/// * The eigen decomposition fails.
///
/// # Panics
///
/// * Panics if `x`, `y`, and `weights` do not all have the same length.
/// * Panics if no eigenvector satisfies the ellipse condition `4ac - b² > 0`.
#[inline]
pub fn fit_one_ellipse_fixed_center(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    weights: ArrayView1<'_, f64>,
    x0: f64,
    y0: f64,
) -> Result<(f64, f64, f64, f64), LinalgError> {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), weights.len());
    let n = x.len();

    // 1) Build centered features fa, fb, fc
    let fa = x.mapv(|xi| xi * xi - 2.0 * x0 * xi);
    let fb = x.iter()
        .zip(y.iter())
        .map(|(&xi, &yi)| xi * yi - y0 * xi - x0 * yi)
        .collect::<Array1<f64>>();
    let fc = y.mapv(|yi| yi * yi - 2.0 * y0 * yi);

    // 2) Design matrices: D1 = [fa, fb, fc], D2 = [1]
    let mut d1 = Array2::<f64>::zeros((n, 3));
    let mut d2 = Array2::<f64>::zeros((n, 1));
    d1.column_mut(0).assign(&fa);
    d1.column_mut(1).assign(&fb);
    d1.column_mut(2).assign(&fc);
    d2.column_mut(0).fill(1.0);

    // 3) Apply weights
    let w_sqrt = weights.mapv(|w| if w > 0.0 { w.sqrt() } else { 0.0 });
    let w_col = w_sqrt.view().insert_axis(Axis(1));
    let wd1 = &w_col * &d1; // n x 3
    let wd2 = &w_col * &d2; // n x 1

    // 4) Scatter matrices
    let s1 = wd1.t().dot(&wd1); // 3x3
    let s2 = wd1.t().dot(&wd2); // 3x1
    let s3 = wd2.t().dot(&wd2); // 1x1

    // 5) Eliminate g (constant term) via Schur complement
    let s3_scalar = s3[[0, 0]];
    if s3_scalar == 0.0 {
        // all weights zero? bail out
        eprintln!(
            "fit_one_ellipse_fixed_center: all weights are zero (s3 = {})",
            s3_scalar
        );
        return Err(LinalgError::NotStandardShape {
            obj: "fit_one_ellipse_fixed_center: all weights are zero",
            rows: weights.len() as i32,
            cols: 1,
        });
    }


    let s3_inv = 1.0 / s3_scalar;
    let t = s2.t().mapv(|v| -s3_inv * v);     // t: Array2<f64>
    let m_tmp = s1 + s2.dot(&t);     // 3x3

    // 6) Ellipse constraint matrix, same as before
    let c = arr2(&[
        [0.0, 0.0, 2.0],
        [0.0, -1.0, 0.0],
        [2.0, 0.0, 0.0],
    ]);
    let c_inv = c.inv()?;
    let m = c_inv.dot(&m_tmp);

    // 7) Eigen-decomposition and ellipse selection
    let (_eigvals, eigvecs) = m.eig()?;
    let eigvecs_re = eigvecs.mapv(|c| c.re);

    let mut chosen_col = None;
    for k in 0..3 {
        let a = eigvecs_re[[0, k]];
        let b = eigvecs_re[[1, k]];
        let c = eigvecs_re[[2, k]];
        let con = 4.0 * a * c - b * b;
        if con > 0.0 {
            chosen_col = Some(k);
            break;
        }
    }
    let k = chosen_col.expect("No valid ellipse eigenvector found");
    let q = eigvecs_re.column(k).to_owned(); // [a, b, c]

    let a = q[0];
    let b = q[1];
    let c = q[2];

    // 8) Recover g (constant term) using eliminated relation:
    //    g = -S22^{-1} S21 q = - (1/s3) * (s2^T q)
    // let s21 = s2.t(); // 1x3
    let g = -s3_inv * s2.t().dot(&q)[0];

    // 9) Reconstruct d,e from fixed center
    let d = -2.0 * a * x0 - b * y0;
    let e = -b * x0 - 2.0 * c * y0;

    // 10) Geometric parameters (we ignore x0_fit,y0_fit and trust the fixed center)
    let (_x0_fit, _y0_fit, ap, bp, phi) =
        ellipse_coefficients_to_parameters(a, b, c, d, e, g);

    // 11) Compute WRMS distance (you can keep using geometric_distances;
    //     it will be very close to using the fixed center, since coefficients
    //     are consistent with that center)
    let delta = geometric_distances(x, y, a, b, c, d, e, g);
    let delta2 = delta.mapv(|d| d * d);
    let num = (&delta2 * &weights).sum();
    let den = weights.sum();
    let wrms = (num / den).sqrt();
    let orientation = (TAU - phi) % TAU;

    Ok((ap, bp, orientation, wrms))
}

/// Fits a single weighted circle to 2D points using a weighted algebraic least squares fit.
///
/// The circle is represented in the algebraic form:
///
/// `x² + y² = 2*x0*x + 2*y0*y + c`, where `c = r² - x0² - y0²`.
///
/// This function solves a 3-parameter weighted linear least squares system for
/// `p = [x0, y0, c]` using the normal equations:
///
/// `(Aᵀ W A) p = (Aᵀ W b)`
///
/// with `A = [2x, 2y, 1]` and `b = x² + y²`. Weights are applied via `sqrt(w)` row scaling.
///
/// After solving for `(x0, y0, c)`, the radius is recovered as:
///
/// `r = sqrt(c + x0² + y0²)`
///
/// The reported `wrms` is the weighted RMS of the **geometric** radial residuals
/// `sqrt((x-x0)²+(y-y0)²) - r`.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `weights` - Non-negative weights for each point. Must have the same length as `x` and `y`.
///   Points with larger weights influence the fit more. Non-positive weights contribute zero influence.
///
/// # Returns
///
/// On success, returns:
///
/// `(x0, y0, r, wrms)`
///
/// where:
///
/// * `x0`, `y0` - Center of the fitted circle.
/// * `r` - Radius of the fitted circle.
/// * `wrms` - Weighted root-mean-square geometric radial residual.
///
/// # Errors
///
/// Returns `Err(LinalgError)` if:
///
/// * All weights are zero.
/// * The 3×3 normal-equation matrix is singular and cannot be inverted.
///
/// # Panics
///
/// Panics if `x`, `y`, and `weights` do not all have the same length.
#[inline]
pub fn fit_one_circle(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    weights: ArrayView1<'_, f64>,
) -> Result<(f64, f64, f64, f64), LinalgError> {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), weights.len());
    let n = x.len();

    // Guard against the all-zero weights case (ill-posed)
    let wsum = weights.sum();
    if wsum == 0.0 {
        eprintln!("fit_one_circle: all weights are zero (wsum = {})", wsum);
        return Err(LinalgError::NotStandardShape {
            obj: "fit_one_circle: all weights are zero",
            rows: n as i32,
            cols: 3,
        });
    }

    // Weighted algebraic circle fit.
    // Model: x^2 + y^2 = 2*x0*x + 2*y0*y + c, where c = r^2 - x0^2 - y0^2.
    // Solve (A^T W A) p = (A^T W b), with A = [2x, 2y, 1], p = [x0, y0, c], b = x^2 + y^2.

    let x2 = x.mapv(|v| v * v);
    let y2 = y.mapv(|v| v * v);
    let bvec = &x2 + &y2; // n

    let mut a_mat = Array2::<f64>::zeros((n, 3));
    a_mat.column_mut(0).assign(&x.mapv(|v| 2.0 * v));
    a_mat.column_mut(1).assign(&y.mapv(|v| 2.0 * v));
    a_mat.column_mut(2).fill(1.0);

    // Apply weights via sqrt(W) scaling (same pattern as ellipse fit)
    let w_sqrt = weights.mapv(|w| if w > 0.0 { w.sqrt() } else { 0.0 });
    let w_col = w_sqrt.view().insert_axis(Axis(1));
    let wa = &w_col * &a_mat; // n x 3
    let wb = &w_sqrt * &bvec; // n

    // Normal equations
    let ata = wa.t().dot(&wa); // 3x3
    let atb = wa.t().dot(&wb); // 3

    // Solve via explicit inverse (small 3x3)
    let ata_inv = ata.inv()?;
    let p = ata_inv.dot(&atb);

    let x0 = p[0];
    let y0 = p[1];
    let c = p[2];

    let r2 = c + x0 * x0 + y0 * y0;
    let r = if r2 > 0.0 { r2.sqrt() } else { 0.0 };

    // Weighted RMS geometric distance to the circle
    let deltas: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let dx = x[i] - x0;
            let dy = y[i] - y0;
            (dx * dx + dy * dy).sqrt() - r
        })
        .collect();
    let delta = Array1::from(deltas);
    let delta2 = delta.mapv(|d| d * d);
    let num = (&delta2 * &weights).sum();
    let den = weights.sum();
    let wrms = (num / den).sqrt();

    Ok((x0, y0, r, wrms))
}


/// Fits a single weighted circle to 2D points with a fixed center.
///
/// With the center fixed at `(x0, y0)`, the only remaining parameter is the radius `r`.
/// This function estimates `r` by minimizing the weighted geometric least squares objective:
///
/// `min_r Σ w_i (d_i - r)²`, where `d_i = sqrt((x_i-x0)² + (y_i-y0)²)`.
///
/// The minimizer is the weighted mean of the distances:
///
/// `r = (Σ w_i d_i) / (Σ w_i)`.
///
/// The reported `wrms` is the weighted RMS of the geometric radial residuals `(d_i - r)`.
///
/// # Arguments
///
/// * `x` - x-coordinates of sample points.
/// * `y` - y-coordinates of sample points. Must have the same length as `x`.
/// * `weights` - Non-negative weights for each point. Must have the same length as `x` and `y`.
///   Points with larger weights influence the fit more. Non-positive weights contribute zero influence.
/// * `x0`, `y0` - Fixed center of the circle.
///
/// # Returns
///
/// On success, returns:
///
/// `(r, wrms)`
///
/// where:
///
/// * `r` - Fitted radius.
/// * `wrms` - Weighted root-mean-square geometric radial residual.
///
/// # Errors
///
/// Returns `Err(LinalgError)` if all weights are zero.
///
/// # Panics
///
/// Panics if `x`, `y`, and `weights` do not all have the same length.
#[inline]
pub fn fit_one_circle_fixed_center(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    weights: ArrayView1<'_, f64>,
    x0: f64,
    y0: f64,
) -> Result<(f64, f64), LinalgError> {
    assert_eq!(x.len(), y.len());
    assert_eq!(x.len(), weights.len());
    let n = x.len();

    // Guard against the all-zero weights case (ill-posed)
    let wreal = weights.mapv(|w| if w > 0.0 { w } else { 0.0 });
    let wsum = wreal.sum();
    if wsum == 0.0 {
        eprintln!(
            "fit_one_circle_fixed_center: all weights are zero (wsum = {})",
            wsum
        );
        return Err(LinalgError::NotStandardShape {
            obj: "fit_one_circle_fixed_center: all weights are zero",
            rows: n as i32,
            cols: 1,
        });
    }

    // Geometric distances from fixed center
    let dists_vec: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let dx = x[i] - x0;
            let dy = y[i] - y0;
            (dx * dx + dy * dy).sqrt()
        })
        .collect();
    let dists = Array1::from(dists_vec);

    // Weighted least squares estimate for r:
    // minimize sum_i w_i (d_i - r)^2  =>  r = (sum_i w_i d_i) / (sum_i w_i)
    let r = (&dists * &wreal).sum() / wsum;

    // Weighted RMS geometric residual
    let delta = dists.mapv(|d| d - r);
    let delta2 = delta.mapv(|d| d * d);
    let wrms = ((&delta2 * &wreal).sum() / wsum).sqrt();

    Ok((r, wrms))
}


// Helper functions to reproduce numpy.nanmax, nanmin, nanquantile behavior
fn nanmax(arr: &Array1<f64>) -> Option<f64> {
    arr.iter()
        .filter(|v| !v.is_nan())
        .copied()
        .fold(None, |acc, x| {
            Some(match acc {
                Some(m) => m.max(x),
                None => x,
            })
        })
}

fn nanmin(arr: &Array1<f64>) -> Option<f64> {
    arr.iter()
        .filter(|v| !v.is_nan())
        .copied()
        .fold(None, |acc, x| {
            Some(match acc {
                Some(m) => m.min(x),
                None => x,
            })
        })
}


