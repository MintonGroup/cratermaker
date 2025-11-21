use std::f64::consts::{FRAC_PI_2, PI};
use numpy::ndarray::prelude::*; 
use ndarray_linalg::error::LinalgError;
use ndarray_linalg::{Inverse, Eig};
use rayon::prelude::*;
use crate::ArrayResult;

pub fn radial_distance_to_ellipse(
    x: ArrayView1<'_, f64>,
    y: ArrayView1<'_, f64>,
    a: f64, 
    b: f64,
    orientation: f64, 
    x0: f64, 
    y0:f64
) -> ArrayResult {
    let phi = orientation - FRAC_PI_2;

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

#[inline]
fn ellipse_coefficients_to_parameters(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    f: f64,
    g: f64) -> (f64, f64, f64, f64, f64, f64) {

    // We use the formulas from https://mathworld.wolfram.com/Ellipse.html which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.  
    // Therefore, B, D and F must be scaled.
    let b = b /2.0;
    let d = d / 2.0; 
    let f = f / 2.0;
    let den = b*b - a * c;
    if den > 0.0 {
        panic!("The provided coefficients do not represent an ellipse.");
    }
    let x0 = (c * d - b * f) / den;
    let y0 = (a * f - b * d) / den;
    let numerator = 2.0 * (a * f * f + c * d * d + g * b * b - 2.0 * b * d * f - a * c * g);
    let fac = ((a - c).powi(2) + 4.0 * b * b).sqrt();

    // Compute unsorted semi-major and semi-minor axes
    let ap = (numerator / (den * (fac - (a + c)))) .sqrt();
    let bp = (numerator / (den * (-fac - (a + c)))).sqrt();
    // Sort the semi-major and semiminor axis but keep track of orientation
    let (_a, _b, width_gt_height) = if ap >= bp {
        (ap, bp, true)
    } else {
        (bp, ap, false)
    };

    // Now compute eccentricity
    let mut r = (b / a).powi(2);
    if r > 1.0 {
        r = 1.0 / r;
    }
    let ep = (1.0 - r).sqrt();


    // compute the orientation with respect to the x-axis
    // if b is very close to zero, then phi is 0 if a < c otherwise it is pi/2
    let mut phi = if b.abs() < 1e-12 {
            if a < c {
                0.0
            } else {
                FRAC_PI_2
            }
        } else {
            0.5 * (2.0 * b).atan2(a - c)
        }; 
    if a > c {
        phi += FRAC_PI_2;
    }

    if !width_gt_height {
        phi += FRAC_PI_2;
    }
    let orientation = phi % PI;
    

    (x0, y0, ap, bp, ep, orientation)
}

#[inline]
fn ellipse_parameters_to_coefficients(x0: f64, y0: f64, a: f64, b: f64, orientation: f64) -> (f64, f64, f64, f64, f64, f64) {
    let phi = orientation - FRAC_PI_2;
    let s2 = phi.sin().powi(2);
    let c2 = phi.cos().powi(2);
    let sc = phi.sin() * phi.cos();
    let a2 = a * a;
    let b2 = b * b;

    let a = a2 * s2 + b2 * c2;
    let b = 2.0 * (b2 - a2) * sc;
    let c = a2 * c2 + b2 * s2;
    let d = -2.0 * a * x0 - b * y0;
    let e = -b * x0 - 2.0 * c * y0;
    let f = a * x0 * x0 + b * x0 * y0 + c * y0 * y0 - a2 * b2;

    (a, b, c, d, e, f)
}

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
    let (x0, y0, a, b, e, orientation) = ellipse_coefficients_to_parameters(a, b, c, d, e, f);

    let n = x.len();
    let result_vec: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let xi = x[i];
            let yi = y[i];
            let dx = xi - x0;
            let dy = yi - y0;
            let theta = dy.atan2(dx);
            let alpha = theta - orientation + FRAC_PI_2;
            let ca = alpha.cos();
            let sa = alpha.sin();
            let capf = (a * b) / ((b * ca).powi(2) + (a * sa).powi(2)).sqrt();
            let capfx = 2.0 * a * xi + b * yi + d;
            let capfy = b * xi + 2.0 * c * yi + e;
            let gradnorm = (capfx * capfx + capfy * capfy).sqrt();
            capf / gradnorm
        })
        .collect();

    let result = Array1::from(result_vec);
    result
}

#[inline]
fn _score_rim(
    x: ArrayView1<'_, f64>, 
    y: ArrayView1<'_, f64>, 
    face_elevation: ArrayView1<'_,f64>,   
    face_lon: ArrayView1<'_,f64>,
    face_lat: ArrayView1<'_,f64>,
    face_bearing: ArrayView1<'_,f64>,
    face_edge_connectivity: ArrayView2<'_, i64>,
    edge_face_connectivity: ArrayView2<'_,i64>,
    edge_face_distance: ArrayView1<'_,f64>,
    edge_length: ArrayView1<'_,f64>,
    x0: f64, 
    y0: f64, 
    ap: f64,
    bp: f64,
    orientation: f64
) ->Array1<f64> { 
    let mut score = Array1::<f64>::zeros(x.len());
    let MIN_POINTS_FOR_FIT = 100;

    let radial_gradient = crate::surface::compute_radial_gradient(
        face_elevation,
        face_lon,
        face_lat,
        face_bearing,
        face_edge_connectivity,
        edge_face_connectivity,
        edge_face_distance,
        edge_length,
    ); 
    score
}

#[inline] 
pub fn fit_one_ellipse(
    x: ArrayView1<'_,f64>, 
    y: ArrayView1<'_,f64>, 
    weights: ArrayView1<'_,f64>
) -> Result<(f64, f64, f64, f64, f64, f64, f64), LinalgError> {
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
    let w_col = weights.view().insert_axis(Axis(1));
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

    let (x0, y0, ap, bp, ep, orientation) = ellipse_coefficients_to_parameters(a, b, c, d, f, g,);

    Ok((x0, y0, ap, bp, ep, orientation, wrms))
}

