use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray1,PyArray1};


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
pub struct SimpleMoonMorphology {
    pub floor_depth: f64,
    pub floor_diameter: f64,
    pub rim_height: f64,
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
/// * `floor_depth` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_height` - Height of the crater rim above mean surface level (in meters).
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
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
    diameter: f64,
    floor_depth: f64,
    floor_diameter: f64,
    rim_height: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distances_v= radial_distances.as_array();
    let reference_elevations_v= reference_elevations.as_array();
    let result =  cratermaker_components::morphology::crater_profile(
                radial_distances_v,
                reference_elevations_v,
                diameter,
                floor_depth,
                floor_diameter,
                rim_height,
                ejrim
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
/// * `radial_distance` - 1D array of radial distances from crater center.
/// * `crater_diameter` - Diameter of the crater (meters).
/// * `ejrim` - Profile scaling factor.
///
/// # Returns
///
/// * A NumPy array of ejecta profile values.
#[pyfunction]
pub fn ejecta_profile<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distance_v = radial_distance.as_array();
    let result =  cratermaker_components::morphology::ejecta_profile(
                radial_distance_v,
                crater_diameter,
                ejrim
            )
            .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))    
}


/// Computes a ray-modulated ejecta intensity field.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distance` - 1D array of radial distances from crater center.
/// * `initial_bearing` - 1D array of bearing angles (degrees).
/// * `crater_diameter` - Crater diameter (meters).
///
/// # Returns
///
/// * A NumPy array of ray-modulated intensity values.
#[pyfunction]
pub fn ray_intensity<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    seed: u64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let radial_distance_v = radial_distance.as_array();
    let initial_bearing_v = initial_bearing.as_array();
    let result = cratermaker_components::morphology::ray_intensity(
        radial_distance_v,
        initial_bearing_v,
        crater_diameter,
        seed,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}
