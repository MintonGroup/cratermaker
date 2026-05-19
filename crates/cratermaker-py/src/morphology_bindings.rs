use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject)]
pub struct Crater {
    pub diameter: f64,
}

/// Morphological parameters for generating and modifying lunar surface craters.
///
/// Includes floor geometry, rim height, and whether to apply ray modulation.
#[derive(FromPyObject)]
pub struct BasicMoonMorphology {
    pub floor_elevation: f64,
    pub floor_diameter: f64,
    pub rim_elevation: f64,
    pub ejrim: f64,
    pub crater: Crater,
}

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
/// * `diameter` - Total diameter of the crater (in meters).
/// * `floor_elevation` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters).
/// * `rimdrop` - Exponent for the rim dropoff function (typically -6.0)
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
/// * `ejprofile` - Exponent for the power-law decay of the ejecta profile (typically -3.0)
/// * `fassett_yang_fraction` - Weighting factor (0.0 to 1.0) for blending between the Fassett (2020) and Yang (2021) profiles.
/// * `morphology_subtype` - Subtype of crater morphology to use for the Yang (2021) profile ("normal", "central mound", "flat-bottomed", or "concentric").
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
#[pyfunction]
pub fn crater_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    floor_elevation: f64,
    floor_diameter: f64,
    rim_elevation: f64,
    rimdrop: f64,
    ejrim: f64,
    ejprofile: f64,
    fassett_yang_fraction: f64,
    morphology_subtype: &str,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let result = cratermaker_components::morphology::basicmoon::crater_profile(
        radial_distances_v,
        reference_elevations_v,
        crater_diameter,
        floor_elevation,
        floor_diameter,
        rim_elevation,
        rimdrop,
        ejrim,
        ejprofile,
        fassett_yang_fraction,
        morphology_subtype,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}

/// Computes only the radial ejecta profile without ray modulation.
///
/// This is a simple power-law decay of ejecta intensity with radial distance.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distances` - 1D array of radial distances from crater center.
/// * `crater_diameter` - Diameter of the crater (meters).
/// * `ejrim` - Profile scaling factor.
/// * `ejprofile` - Exponent for the power-law decay of the ejecta profile (typically -3.0).
///
/// # Returns
///
/// * A NumPy array of ejecta profile values.
#[pyfunction]
pub fn ejecta_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
    ejprofile: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v = radial_distances.as_array();
    let result = cratermaker_components::morphology::basicmoon::ejecta_profile(
        radial_distances_v,
        crater_diameter,
        ejrim,
        ejprofile,
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
