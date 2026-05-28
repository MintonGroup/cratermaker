use cratermaker_components::morphology::basicmoon::BasicMoonCrater;
use cratermaker_components::morphology::realmoon::RealMoonCrater;
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Computes a crater profile elevation array from input radial distances and reference elevations.
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_radius` - Total radius of the crater (in meters).
/// * `floor_elevation` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_radius` - Radius of the crater floor (in meters).
/// * `wall_curvature` - Parameter controlling the curvature of the crater wall (>1 for more curvature)
/// * `rim_width` - Width of the crater rim (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters).
/// * `rimdrop` - Exponent for the rim dropoff function (typically -4.0 to -6.0)
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
/// * `peak_height` - Height of the central peak above the crater floor (in meters).
/// * `peak_width` - Width of the central peak (in meters).
/// * `peak_offset` - Radial offset of the central peak from the crater center (in meters).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
#[pyfunction]
pub fn basicmoon_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: BasicMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let result = cratermaker_components::morphology::basicmoon::basicmoon_profile(
        radial_distances_v,
        reference_elevations_v,
        &crater,
        include_crater,
        include_ejecta,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}

/// Computes a ray-modulated ejecta intensity field.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distances` - 1D array of radial distances from crater center.
/// * `initial_bearing` - 1D array of bearing angles (radians).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A NumPy array of ray-modulated intensity values.
#[pyfunction]
pub fn ray_intensity<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    seed: u64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let initial_bearing_v = initial_bearing.as_array();
    let result = cratermaker_components::morphology::basicmoon::ray_intensity(
        radial_distances_v,
        initial_bearing_v,
        crater_diameter,
        seed,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


/// Computes a crater profile elevation array from input radial distances and reference elevations using the realistic moon model of Du et al. (2024a,b).
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_radius` - Total radius of the crater (in meters).
/// * `floor_elevation` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_radius` - Radius of the crater floor (in meters).
/// * `wall_curvature` - Parameter controlling the curvature of the crater wall (>1 for more curvature)
/// * `rim_width` - Width of the crater rim (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters).
/// * `rimdrop` - Exponent for the rim dropoff function (typically -4.0 to -6.0)
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
/// * `peak_height` - Height of the central peak above the crater floor (in meters).
/// * `peak_width` - Width of the central peak (in meters).
/// * `peak_offset` - Radial offset of the central peak from the crater center (in meters).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
#[pyfunction]
pub fn realmoon_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: RealMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let result = cratermaker_components::morphology::realmoon::realmoon_profile(
        radial_distances_v,
        reference_elevations_v,
        &crater,
        include_crater,
        include_ejecta,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}
