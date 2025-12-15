use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray1,PyReadonlyArray2,PyArray1};
use crate::surface_bindings::PyReadonlyLocalSurface;
use cratermaker_components::crater::Crater;
const _EXTENT_RADIUS_RATIO: f64 = 2.5;


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
/// * `distmult` - Distance multiplier for scoring.
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
    distmult: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) -> PyResult<Bound<'py, PyAny>>  {
    let region = surface.call_method1("extract_region",(crater.measured_location, _EXTENT_RADIUS_RATIO * crater.measured_radius))?;
    let transformer = region.getattr("from_surface").unwrap();
    let x0y0 = transformer.call_method1("transform",(crater.measured_location.0, crater.measured_location.1))?;
    let (x0, y0): (f64, f64) = x0y0.extract()?;

    // Scoring is best if the surface is modified to remove any regional slope.
    // Save the original face elevation so we can put it back at the end
    let face_elevation: Vec<f64> = region.getattr("face_elevation")?.extract()?;
    // Get the reference surface and subract it from the local surface
    let reference_radius: f64 = region.getattr("region_radius")?.extract()?;
    let mut reference_elevation: Vec<f64> = region.
                                        call_method1("get_reference_surface",(reference_radius,))?
                                        .extract()?;
    reference_elevation = reference_elevation.iter().map(|&v| -v).collect();
    region.call_method1("update_elevation", (reference_elevation,))?;
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let result = cratermaker_components::counting::score_rim(
            &region_v,
            x0,
            y0,
            &crater,
            quantile,
            distmult,
            gradmult,
            curvmult,
            heightmult
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;

    region.call_method1("add_data",("rimscore", PyArray1::from_owned_array(py, result.clone()), "Rim Score", "dimensionless", true, true, f64::NAN))?;

    // Put the original face elevation back using overwrite=true
    region.call_method1("update_elevation",(face_elevation,true))?;

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
/// 
/// # Returns
///
/// * On success, returns a tuple `(location_fit, a_fit, b_fit, o_fit)` where:
/// - `location_fit`: Fitted location of the crater center.
/// - `a_fit`: Fitted semi-major axis length.
/// - `b_fit`: Fitted semi-minor axis length.
/// - `o_fit`: Fitted orientation angle in radians.
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
) -> PyResult<(Bound<'py, PyAny>, f64, f64, f64)>  {
    let region = surface.call_method1("extract_region",(crater.location, _EXTENT_RADIUS_RATIO * crater.radius))?;
    let transformer = region.getattr("from_surface").unwrap();
    let x0y0 = transformer.call_method1("transform",(crater.location.0, crater.location.1))?;
    let (x0, y0): (f64, f64) = x0y0.extract()?;

    // Scoring is best if the surface is modified to remove any regional slope.
    // Save the original face elevation so we can put it back at the end
    let face_elevation: Vec<f64> = region.getattr("face_elevation")?.extract()?;
    // Get the reference surface and subract it from the local surface
    let reference_radius: f64 = region.getattr("region_radius")?.extract()?;
    let mut reference_elevation: Vec<f64> = region.
                                        call_method1("get_reference_surface",(reference_radius,))?
                                        .extract()?;
    reference_elevation = reference_elevation.iter().map(|&v| -v).collect();
    region.call_method1("update_elevation", (reference_elevation,))?;
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let (x0_fit, y0_fit, a_fit, b_fit, o_fit, rimscore) = cratermaker_components::counting::fit_rim(
            &region_v,
            x0,
            y0,
            &crater,
            tol,
            nloops,
            score_quantile,
            fit_center
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let transformer = region.getattr("to_surface").unwrap();
    let location_fit = transformer.call_method1("transform", (x0_fit, y0_fit))?;
    region.call_method1("add_data",("rimscore", PyArray1::from_owned_array(py, rimscore.clone()), "Rim Score", "dimensionless", true, true, f64::NAN))?; 

    // Put the original face elevation back using overwrite=true
    region.call_method1("update_elevation",(face_elevation,true))?;

    Ok((location_fit, a_fit, b_fit, o_fit))
}