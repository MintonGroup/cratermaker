use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray1,PyArray1};
use crate::surface_bindings::PyReadonlyLocalSurface;
use cratermaker_components::crater::Crater;
const _FITTING_RADIUS_RATIO: f64 = 3.0;
const _MEASURING_RADIUS_RATIO: f64 = 1.2;


#[pyfunction]
pub fn measure_degradation_state<'py>(
    py: Python<'py>,
    surface: &Bound<'py, PyAny>,
    crater: Crater, 
) -> PyResult<f64> {
    // for (key, value) in observed.as_ref(py).iter() {
    //     let id: u32 = key.extract()?;
    //     let crater: &PyDict = value.downcast::<PyDict>()?;

    //     let final_diameter: Option<f64> = crater.get_item("final_diameter").and_then(|v| v.extract().ok());
    //     let location: Option<(f64, f64)> = crater.get_item("location").and_then(|v| v.extract().ok());

    // }


    // for id in id_array.iter() {
    //     id_vec.push(*id);
    // }
    // let id_array_flat = PyArray1::from_vec(py, id_vec);
    Ok(0.0)
}

#[pyfunction]
pub fn measure_crater_depth<'py>(
    _py: Python<'py>,
    surface: &Bound<'py, PyAny>,
    crater: Crater, 
) -> PyResult<f64> {

    let region = surface.call_method1("extract_region",(crater.measured_location, _MEASURING_RADIUS_RATIO * crater.measured_radius))?;
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let bowl_depth = cratermaker_components::counting::measure_bowl_depth(
            &region_v,
            &crater,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let rim_height = cratermaker_components::counting::measure_rim_height(
            &region_v,
            &crater,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;


    Ok(rim_height - bowl_depth)
}

/// Fit a single ellipse to the provided x, y coordinates and weights.
/// 
/// # Arguments
/// 
/// * `py` - The Python GIL token.
/// * `x` - A 1D array of x coordinates.
/// * `y` - A 1D array of y coordinates.
/// * `weights` - A 1D array of weights corresponding to each (x, y) point.
/// 
/// # Returns
/// 
/// * On success, returns a tuple `(x0, y0, ap, bp, orientation, wrms)` where:
///  - `x0`, `y0`: Center of the fitted ellipse.
/// - `ap`: Semi-major axis length.
/// - `bp`: Semi-minor axis length.
/// - `orientation`: Orientation angle of the ellipse in radians.
/// - `wrms`: Weighted root mean square error of the fit.
/// 
/// # Errors
/// 
/// * Returns `Err(PyValueError)` if the fitting process fails.
/// 
#[pyfunction]
pub fn fit_one_ellipse<'py>(
    _py: Python<'py>,
    x: PyReadonlyArray1<'_,f64>, 
    y: PyReadonlyArray1<'_,f64>, 
    weights: PyReadonlyArray1<'_,f64>
) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let weights_v = weights.as_array();
    let (x0, y0, ap, bp, orientation, wrms) = cratermaker_components::counting::fit_one_ellipse(
            x_v,
            y_v,
            weights_v,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok((x0, y0, ap, bp, orientation, wrms))
}


/// Fit a single ellipse to the provided x, y coordinates and weights with a fixed center.
/// 
/// # Arguments
/// 
/// * `py` - The Python GIL token.
/// * `x` - A 1D array of x coordinates.
/// * `y` - A 1D array of y coordinates.
/// * `weights` - A 1D array of weights corresponding to each (x, y) point.
/// * `x0` - Fixed x coordinate of the ellipse center.
/// * `y0` - Fixed y coordinate of the ellipse center.
/// 
/// # Returns
/// 
/// * On success, returns a tuple `(ap, bp, orientation, wrms)` where:
///  - `ap`: Semi-major axis length.
/// - `bp`: Semi-minor axis length.
/// - `orientation`: Orientation angle of the ellipse in radians.
/// - `wrms`: Weighted root mean square error of the fit.
/// 
/// # Errors
/// 
/// * Returns `Err(PyValueError)` if the fitting process fails.
/// 
#[pyfunction]
pub fn fit_one_ellipse_fixed_center<'py>(
    _py: Python<'py>,
    x: PyReadonlyArray1<'_,f64>, 
    y: PyReadonlyArray1<'_,f64>, 
    weights: PyReadonlyArray1<'_,f64>,
    x0: f64,
    y0: f64,
) -> PyResult<(f64, f64, f64, f64)> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let weights_v = weights.as_array();
    let (ap, bp, orientation, wrms) = cratermaker_components::counting::fit_one_ellipse_fixed_center(
            x_v,
            y_v,
            weights_v,
            x0,
            y0,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok((ap, bp, orientation, wrms))
}


/// Score the rim of a crater on the provided surface.
/// 
/// # Arguments
/// 
/// * `py` - The Python GIL token.
/// * `surface` - A reference to the surface object.
/// * `crater` - The crater object containing location and radius.
/// * `quantile` - The quantile value for scoring.
/// * `gradmult` - Gradient multiplier for scoring.
/// * `curvmult` - Curvature multiplier for scoring.
/// * `heightmult` - Height multiplier for scoring.
/// 
/// # Returns
/// 
/// * On success, returns a reference to the modified surface object with rim scores added.
/// 
/// # Errors
/// 
/// * Returns `Err(PyValueError)` if any error occurs during the scoring process.
///
#[pyfunction]
pub fn score_rim<'py>(
    py: Python<'py>,
    surface: &Bound<'py, PyAny>,
    crater: Crater, 
    quantile: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) -> PyResult<Bound<'py, PyAny>>  {

    let region = surface.call_method1("extract_region",(crater.measured_location, _FITTING_RADIUS_RATIO * crater.measured_radius))?;
    // Ensure face projections are set
    region.call_method0("set_face_proj")?;
    let transformer = region.getattr("from_surface").unwrap();
    let x0y0 = transformer.call_method1("transform",(crater.measured_location.0, crater.measured_location.1))?;
    let (x0, y0): (f64, f64) = x0y0.extract()?;

    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let result = cratermaker_components::counting::score_rim(
            &region_v,
            x0,
            y0,
            &crater,
            quantile,
            gradmult,
            curvmult,
            heightmult
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;

    region.call_method1("add_data",("rimscore", PyArray1::from_owned_array(py, result.clone()), "Rim Score", "dimensionless", true, true, f64::NAN))?;

    Ok(region)
}


/// Fit the rim of a crater on the provided surface.
/// 
/// # Arguments
/// 
/// * `py` - The Python GIL token.
/// * `surface` - A reference to the surface object.
/// * `crater` - The crater object containing location and radius.
/// * `tol` - Tolerance for convergence of the fitting algorithm.
/// * `nloops` - Maximum number of iterations to perform.
/// * `score_quantile` - Quantile of rim scores to consider.
/// * `fit_center` - Whether to fit the crater center or keep it fixed.
/// * `fit_ellipse` - Whether to fit an ellipse to the rim. If False, it will fit a circle.
/// 
/// # Returns
///
/// * On success, returns a tuple `(location_fit, a_fit, b_fit, o_fit)` where:
/// - `location_fit`: Fitted location of the crater center.
/// - `a_fit`: Fitted semi-major axis length.
/// - `b_fit`: Fitted semi-minor axis length. (same as a_fit if fit_ellipse is False)
/// - `o_fit`: Fitted orientation angle in radians. (0 if fit_ellipse is False)
/// 
/// # Errors
/// 
/// * Returns `Err(PyValueError)` if any error occurs during the fitting process.
#[pyfunction]
pub fn fit_rim<'py>(
    py: Python<'py>,
    surface: &Bound<'py, PyAny>,
    crater: Crater, 
    tol: f64,
    nloops: usize,
    score_quantile: f64,
    fit_center: bool,
    fit_ellipse: bool,
) -> PyResult<(Bound<'py, PyAny>, f64, f64, f64)>  {
    let mut region: Bound<'_, PyAny>;
    let mut location_fit: Bound<'_, PyAny>;
    let mut crater_fit = crater.clone();
    let mut x0_fit:f64;
    let mut y0_fit:f64;
    let mut a_fit:f64;
    let mut b_fit:f64;
    let mut o_fit:f64;
    let mut rimscore: numpy::ndarray::Array1<f64>;

    let mut i = 0;
    loop {

        // Extract a region surrounding the current best fit location and best-fit radius
        region = surface.call_method1("extract_region",(crater_fit.measured_location, _FITTING_RADIUS_RATIO * crater_fit.measured_radius))?;
        // Ensure face projections are set
        region.call_method0("set_face_proj")?;
        let transformer = region.getattr("from_surface").unwrap();
        let x0y0 = transformer.call_method1("transform",(crater_fit.measured_location.0, crater_fit.measured_location.1))?;
        let (x0, y0): (f64, f64) = x0y0.extract()?;

        let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
        let region_v = region_py.as_views();

        // Update the multipliers depending on the iteration. The gradient and curvature multipliers are reduced while the hieght multiplier is increased as the fit is refined. 
        let gradmult =  1.0 / (i as f64 + 1.0);
        let curvmult = 1.0 / (i as f64 + 0.1);
        let heightmult = (i+1) as f64;

        (x0_fit, y0_fit, a_fit, b_fit, o_fit, rimscore) = cratermaker_components::counting::fit_one_rim(
                &region_v,
                x0,
                y0,
                &crater_fit,
                fit_center,
                fit_ellipse,
                score_quantile, 
                gradmult, 
                curvmult, 
                heightmult
            )
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        // Put the original face elevation back using overwrite=true
        let transformer = region.getattr("to_surface").unwrap();
        location_fit = transformer.call_method1("transform", (x0_fit, y0_fit))?;

        let delta_a = (a_fit - crater_fit.measured_semimajor_axis).abs() / crater_fit.radius;
        let delta_b = (b_fit - crater_fit.measured_semiminor_axis).abs() / crater_fit.radius;
        if delta_a < tol && delta_b < tol {
            break;
        }
        i += 1;
        if i >= nloops {
            break;
        }

        crater_fit.measured_location = location_fit.extract()?;
        crater_fit.measured_semimajor_axis = a_fit;
        crater_fit.measured_semiminor_axis = b_fit;
        crater_fit.measured_orientation = o_fit.to_degrees();

    }

    region.call_method1("add_data",("rimscore", PyArray1::from_owned_array(py, rimscore.clone()), "Rim Score", "dimensionless", true, true, f64::NAN))?; 


    Ok((location_fit, a_fit, b_fit, o_fit))
}