use cratermaker_components::morphology::basicmoon::BasicMoonCrater;
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
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `bearings` - 1D array of bearing angles (radians, clockwise north). This is not used for this model, but is included for compatability
/// * `reference_elevations` - 1D array of reference elevations corresponding to each radius.
/// * `crater` - A BasicMoonCrater struct containing the crater's properties.
/// * `include_crater` - Boolean indicating whether to include the crater profile.
/// * `include_ejecta` - Boolean indicating whether to include the ejecta profile.
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
    bearings: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: BasicMoonCrater,
    rings: Option<Vec<BasicMoonCrater>>,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let _bearings_v = bearings.as_array(); // Suppresses the unused argument warning.
    let result = cratermaker_components::morphology::basicmoon::basicmoon_profile(
        radial_distances_v,
        reference_elevations_v,
        &crater,
        &rings,
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
/// * `earings` - 1D array of bearing angles (radians, clockwise north).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A NumPy array of ray-modulated intensity values.
#[pyfunction]
pub fn ray_intensity<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    bearings: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    seed: u64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let bearings_v = bearings.as_array();
    let result = cratermaker_components::morphology::basicmoon::ray_intensity(
        radial_distances_v,
        bearings_v,
        crater_diameter,
        seed,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}

#[pyfunction]
pub fn crater_profile_function<'py>(
    _py: Python<'py>,
    crater: BasicMoonCrater,
    r: f64,
) -> PyResult<f64> {
    let result = cratermaker_components::morphology::basicmoon::crater_profile_function(
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
    );
    Ok(result)
}

#[pyfunction]
pub fn ejecta_profile_function<'py>(
    _py: Python<'py>,
    crater: BasicMoonCrater,
    r: f64,
) -> PyResult<f64> {
    let result = cratermaker_components::morphology::basicmoon::ejecta_profile_function(
        r,
        crater.radius,
        crater.frac_ejrim,
        crater.ejprofile,
        crater.rim_width,
    );
    Ok(result)
}
