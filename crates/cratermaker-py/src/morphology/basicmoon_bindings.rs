use crate::morphology::form_crater_bindings::run_form_craters;
use cratermaker_components::morphology::basicmoon::BasicMoonCrater;
use cratermaker_components::morphology::form_crater::BasicMoonModel;
use numpy::ndarray::Array2;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1};
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
        crater.rim_height,
        crater.ejprofile,
        crater.rim_width,
    );
    Ok(result)
}

/// Computes the volume integrand r * h(r) of a crater profile at an array of radial distances.
///
/// This is the vectorized equivalent of evaluating `r * (crater_profile_function - ejecta_profile_function)` (crater only),
/// `r * ejecta_profile_function` (ejecta only), or `r * crater_profile_function` (both) point by point, and is used to
/// estimate crater volumes by quadrature.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `radial_distances` - 1D array of radial distances from the crater center (in meters).
/// * `crater` - A BasicMoonCrater struct containing the crater's properties.
/// * `include_crater` - Include the crater part of the profile.
/// * `include_ejecta` - Include the ejecta part of the profile.
///
/// # Returns
///
/// * A NumPy array of the integrand evaluated at each radial distance.
#[pyfunction]
pub fn volume_integrand<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    crater: BasicMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let result = cratermaker_components::morphology::basicmoon::volume_integrand(
        radial_distances.as_array(),
        &crater,
        include_crater,
        include_ejecta,
    );
    Ok(PyArray1::from_owned_array(py, result))
}

/// Computes the integrand of the profile volume (see `volume_integrand`) for each of several values of the rim thickness
/// fraction `frac_ejrim` of a crater, in one call (the volume conservation of `BasicMoonCrater.maker` tries a grid of them).
///
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from the crater center (in meters).
/// * `crater` - The crater; its `frac_ejrim` is replaced by each of `fracs`.
/// * `fracs` - The values of `frac_ejrim`.
/// * `include_crater`, `include_ejecta` - The parts of the profile included.
///
/// # Returns
///
/// * A 2D NumPy array (one row per value of `frac_ejrim`) of the integrand at each radial distance.
#[pyfunction]
pub fn volume_integrand_fracs<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    crater: BasicMoonCrater,
    fracs: PyReadonlyArray1<'py, f64>,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let r = radial_distances.as_array();
    let fracs = fracs.as_array();
    let mut out = Array2::<f64>::zeros((fracs.len(), r.len()));
    for (row, &frac) in out.rows_mut().into_iter().zip(fracs.iter()) {
        let trial = BasicMoonCrater { frac_ejrim: frac, ..crater };
        let values = cratermaker_components::morphology::basicmoon::volume_integrand(r, &trial, include_crater, include_ejecta);
        row.into_iter().zip(values.iter()).for_each(|(o, v)| *o = *v);
    }
    Ok(PyArray2::from_owned_array(py, out))
}

/// Computes the degradation function for a single crater diameter.
///
/// # Arguments
///
/// * `diameter` - Final crater diameter (meters).
/// * `fe` - Degradation function size factor.
///
/// # Returns
///
/// * The degradation state contributed by the crater.
#[pyfunction]
pub fn degradation_function(diameter: f64, fe: f64) -> f64 {
    cratermaker_components::morphology::basicmoon::degradation_function(diameter, fe)
}

/// Computes the degradation function for an array of crater diameters.
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `diameter` - 1D array of final crater diameters (meters).
/// * `fe` - Degradation function size factor.
///
/// # Returns
///
/// * A NumPy array of degradation states, one per diameter.
#[pyfunction]
pub fn degradation_function_array<'py>(
    py: Python<'py>,
    diameter: PyReadonlyArray1<'py, f64>,
    fe: f64,
) -> Bound<'py, PyArray1<f64>> {
    let result = cratermaker_components::morphology::basicmoon::degradation_function_array(
        diameter.as_array(),
        fe,
    );
    PyArray1::from_owned_array(py, result)
}

/// Forms a batch of basic lunar craters with non-overlapping regions, updating the global surface arrays in place.
///
/// Each crater goes through the whole formation sequence (crater shape, excavated volume, ejecta with volume conservation,
/// burial and degradation diffusion, and optional slope collapse); see `cratermaker_components::morphology::form_crater`. The
/// craters are formed in parallel with the GIL released, which gives the same result as forming them one after another because
/// their regions do not overlap.
///
/// # Arguments
///
/// * `craters` - The BasicMoonCraters to form.
/// * `regions` - The Python `LocalSurface` each crater is formed on.
/// * `face_indices`, `node_indices` - Global indices of each region's faces and nodes.
/// * `ray_seeds` - The ray pattern seed of each crater, or None for ejecta without rays.
/// * `face_elevation`, `node_elevation`, `ejecta_thickness` - The global surface arrays, updated in place.
/// * `conserve_volume` - Scale each crater's ejecta so that its volume matches the excavated volume.
/// * `critical_slope` - Critical slope (rise over run) for slope collapse, or None to skip slope collapse.
///
/// # Returns
///
/// * A list of `(excavated_volume, conservation_factor)` for each crater, where `conservation_factor` is None if the crater was
///   not formed because its excavated volume was not negative.
#[pyfunction]
#[pyo3(signature = (craters, regions, face_indices, node_indices, ray_seeds, face_elevation, node_elevation, ejecta_thickness, conserve_volume, critical_slope=None))]
#[allow(clippy::too_many_arguments)]
pub fn form_craters<'py>(
    py: Python<'py>,
    craters: Vec<BasicMoonCrater>,
    regions: Vec<Bound<'py, PyAny>>,
    face_indices: Vec<PyReadonlyArray1<'py, i64>>,
    node_indices: Vec<PyReadonlyArray1<'py, i64>>,
    ray_seeds: Vec<Option<u64>>,
    face_elevation: Bound<'py, PyArray1<f64>>,
    node_elevation: Bound<'py, PyArray1<f64>>,
    ejecta_thickness: Bound<'py, PyArray1<f64>>,
    conserve_volume: bool,
    critical_slope: Option<f64>,
) -> PyResult<Vec<(f64, Option<f64>)>> {
    let models: Vec<BasicMoonModel> = craters
        .iter()
        .map(|crater| BasicMoonModel { crater })
        .collect();
    run_form_craters(
        py,
        &models,
        &regions,
        &face_indices,
        &node_indices,
        &ray_seeds,
        &face_elevation,
        &node_elevation,
        &ejecta_thickness,
        conserve_volume,
        critical_slope,
    )
}

/// Computes the maximum extent of a crater or its ejecta blanket: the radius where the rim or ejecta profile falls to
/// `minimum_thickness`, searched for in `[lower_limit, upper_limit]`.
///
/// # Arguments
///
/// * `crater` - The crater.
/// * `minimum_thickness` - The threshold height (m).
/// * `feature` - Either "crater" (the rim profile) or "ejecta" (the ejecta profile).
/// * `lower_limit`, `upper_limit` - The search range (m).
///
/// # Returns
///
/// * The maximum extent (m).
#[pyfunction]
pub fn rmax(
    crater: BasicMoonCrater,
    minimum_thickness: f64,
    feature: &str,
    lower_limit: f64,
    upper_limit: f64,
) -> PyResult<f64> {
    let ejecta = match feature {
        "ejecta" => true,
        "crater" => false,
        _ => {
            return Err(PyValueError::new_err(
                "Unknown feature type. Choose either 'crater' or 'ejecta'",
            ));
        }
    };
    cratermaker_components::morphology::basicmoon::rmax(
        &crater,
        minimum_thickness,
        ejecta,
        lower_limit,
        upper_limit,
    )
    .map_err(PyValueError::new_err)
}
