use cratermaker_components::morphology::realmoon::RealMoonCrater;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;




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


#[pyfunction]
pub fn get_1d_psd_from_control_points<'py>(
    py: Python<'py>,
    control_points: HashMap<String, f64>,
    npoints: usize,
    add_noise: bool,
    rng_seed: u64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let result = cratermaker_components::morphology::realmoon::get_1d_psd_from_control_points(
        &control_points,
        npoints,
        add_noise,
        rng_seed,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray2::from_owned_array(py, result))
}

///
/// Generates a surface profile based on a 1D power spectral density (PSD) and optional phase information, simulating a crater surface with specified roughness characteristics.
///
/// # Arguments
/// * `crater_radius` - The radius of the crater (in meters), which scales the amplitude
/// * `ymean` - The mean elevation of the surface (in meters), which serves as a baseline for the profile.
/// * `psd` - A 2D array where the first column contains wavelengths and the second column contains power values, defining the roughness characteristics of the surface.
/// * `theta` - A 1D array of polar angles (in radians) at which to compute the profile, typically ranging from 0 to 2π.
/// * `phases` - An optional 1D array of phase values (in radians) corresponding to each frequency in the PSD. If not provided, random phases will be generated.
/// * `rng_seed` - The random seed for reproducibility when generating random phases if `phases` is not provided.
/// # Returns
/// * A 1D array of values corresponding to the input bearings, representing the linear profile generated from the PSD and phase information.
///
#[pyfunction]
pub fn profile_from_psd<'py>(
    py: Python<'py>,
    crater_radius: f64,
    ymean: f64,
    psd: PyReadonlyArray2<'py, f64>,
    theta: PyReadonlyArray1<'py, f64>,
    phases: Option<PyReadonlyArray1<'py, f64>>,
    rng_seed: u64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let psd_v = psd.as_array();
    let theta_v = theta.as_array();
    let phases_v = phases.as_ref().map(|p| p.as_array());

    let result = cratermaker_components::morphology::realmoon::profile_from_psd(
        crater_radius,
        ymean,
        psd_v,
        theta_v,
        phases_v,
        rng_seed,
    )
    .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}   

