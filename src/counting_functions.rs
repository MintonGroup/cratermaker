use std::f64::consts::{FRAC_PI_2, PI};
use pyo3::{prelude::*}; 
use numpy::{PyArray1, PyReadonlyArray2, PyReadonlyArray1};
use numpy::ndarray::Array1;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

#[pyfunction]
pub fn tally<'py>(
    py: Python<'py>,
    //face_elevation: PyReadonlyArray1<'py, f64>,
    id_array: PyReadonlyArray2<'py, u32>, 
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let id_array = id_array.as_array();
    let mut id_vec = Vec::with_capacity(id_array.len());
    // for (key, value) in observed.as_ref(py).iter() {
    //     let id: u32 = key.extract()?;
    //     let crater: &PyDict = value.downcast::<PyDict>()?;

    //     let final_diameter: Option<f64> = crater.get_item("final_diameter").and_then(|v| v.extract().ok());
    //     let location: Option<(f64, f64)> = crater.get_item("location").and_then(|v| v.extract().ok());

    // }


    for id in id_array.iter() {
        id_vec.push(*id);
    }
    let id_array_flat = PyArray1::from_vec(py, id_vec);
    Ok(id_array_flat)
}

#[pyfunction]
pub fn radial_distance_to_ellipse<'py>(
    py: Python<'py>,
    x: PyReadonlyArray1<'py, f64>,
    y: PyReadonlyArray1<'py, f64>,
    a: f64, 
    b: f64,
    orientation: f64, 
    x0: f64, 
    y0:f64
)-> PyResult<Bound<'py, PyArray1<f64>>> {
    // Placeholder implementation
    let phi = orientation - FRAC_PI_2;
    let x = x.as_array();
    let y = y.as_array();
    //let mut result = ndarray::Array1::<f64>::zeros(x.len());

    let n = x.len();
    let result_vec: Vec<f64> = (0..n)
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

    let result = Array1::from(result_vec);

    Ok(PyArray1::from_owned_array(py, result))
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
    let orientation = phi + FRAC_PI_2 % PI;
    

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
fn geometric_distances(x: &Array1<f64>, y: &Array1<f64>, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> Array1<f64> {
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

// #[inline] 
// fit_ellipse(x: &ndarray::Array1<f64>, y: &ndarray::Array1<f64>, weights: &ndarray::Array1<f64>) -> (f64, f64, f64, f64, f64, f64, f64) {

//     (a, b, c, d, e, f, wrms)
// }
