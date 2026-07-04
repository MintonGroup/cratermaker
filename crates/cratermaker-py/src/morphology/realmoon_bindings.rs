use cratermaker_components::morphology::realmoon::RealMoonCrater;
use interp::{InterpMode, interp_slice};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;

// Mirrors the RealMoonCrater struct in cratermaker-components and provides read-only access to its fields from Python.
pub struct PyReadonlyRealMoonCrater<'py> {
    pub diameter: f64,
    pub radius: f64,
    pub floor_elevation: f64,
    pub floor_radius: f64,
    pub wall_curvature: f64,
    pub rim_width: f64,
    pub rim_elevation: f64,
    pub rimdrop: f64,
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_ring_radius: f64,
    pub peak_center_distance: f64,
    pub peak_center_bearing: f64,
    pub elevation_offset: f64,
    pub rim_radius_psd: PyReadonlyArray2<'py, f64>,
    pub floor_radius_psd: PyReadonlyArray2<'py, f64>,
}
impl<'py> PyReadonlyRealMoonCrater<'py> {
    /// Build from a Python PyReadonlyCrater object
    pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        Ok(Self {
            diameter: obj.getattr("diameter")?.extract()?,
            radius: obj.getattr("radius")?.extract()?,
            floor_elevation: obj.getattr("floor_elevation")?.extract()?,
            floor_radius: obj.getattr("floor_radius")?.extract()?,
            wall_curvature: obj.getattr("wall_curvature")?.extract()?,
            rim_width: obj.getattr("rim_width")?.extract()?,
            rim_elevation: obj.getattr("rim_elevation")?.extract()?,
            rimdrop: obj.getattr("rimdrop")?.extract()?,
            ejrim: obj.getattr("ejrim")?.extract()?,
            ejprofile: obj.getattr("ejprofile")?.extract()?,
            peak_height: obj.getattr("peak_height")?.extract()?,
            peak_width: obj.getattr("peak_width")?.extract()?,
            peak_ring_radius: obj.getattr("peak_ring_radius")?.extract()?,
            peak_center_distance: obj.getattr("peak_center_distance")?.extract()?,
            peak_center_bearing: obj.getattr("peak_center_bearing")?.extract()?,
            elevation_offset: obj.getattr("elevation_offset")?.extract()?,
            rim_radius_psd: obj.getattr("rim_radius_psd")?.extract()?,
            floor_radius_psd: obj.getattr("floor_radius_psd")?.extract()?,
        })
    }
    /// Convert to cratermaker-components PyReadonlyLocalSurface with array views
    pub fn as_views(&self) -> RealMoonCrater<'_> {
        RealMoonCrater {
            diameter: self.diameter,
            radius: self.radius,
            floor_elevation: self.floor_elevation,
            floor_radius: self.floor_radius,
            wall_curvature: self.wall_curvature,
            rim_width: self.rim_width,
            rim_elevation: self.rim_elevation,
            rimdrop: self.rimdrop,
            ejrim: self.ejrim,
            ejprofile: self.ejprofile,
            peak_height: self.peak_height,
            peak_width: self.peak_width,
            peak_ring_radius: self.peak_ring_radius,
            peak_center_distance: self.peak_center_distance,
            peak_center_bearing: self.peak_center_bearing,
            elevation_offset: self.elevation_offset,
            rim_radius_psd: self.rim_radius_psd.as_array(),
            floor_radius_psd: self.floor_radius_psd.as_array(),
        }
    }
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
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `bearings` - 1D array of bearing angles (radians, clockwise north).
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
///
#[pyfunction]
pub fn realmoon_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    bearings: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: Bound<'py, PyAny>,
    rings: Option<Vec<Bound<'py, PyAny>>>,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let crater_py = PyReadonlyRealMoonCrater::from_py(&crater)?;

    // Convert the optional rings vector into the views that we can pass to the profile function.
    // This has to be done in two steps in order to get the PyReaonly version first as an immutable, then to extract its views.
    let rings_py: Option<Vec<PyReadonlyRealMoonCrater>> = match rings {
        Some(rings) => {
            let mut temp_vec = Vec::new();
            for ring in rings.iter() {
                temp_vec.push(PyReadonlyRealMoonCrater::from_py(&ring)?);
            }
            Some(temp_vec)
        }
        None => None,
    };

    let rings_v: Option<Vec<RealMoonCrater>> = match &rings_py {
        Some(py_vec) => Some(py_vec.iter().map(|py_ring| py_ring.as_views()).collect()),
        None => None,
    };

    let crater_v = crater_py.as_views();
    let radial_distances_v = radial_distances.as_array();
    let bearings_v = bearings.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let result = cratermaker_components::morphology::realmoon::realmoon_profile(
        radial_distances_v,
        bearings_v,
        reference_elevations_v,
        &crater_v,
        &rings_v,
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
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let psd_v = psd.as_array();
    let theta_v = theta.as_array();

    let (profile, psd_theta) =
        cratermaker_components::morphology::realmoon::compute_profile_from_psd(
            crater_radius,
            ymean,
            psd_v,
        );

    let result = interp_slice(
        &psd_theta,
        &profile,
        &theta_v.to_vec(),
        &InterpMode::default(),
    );
    Ok(PyArray1::from_vec(py, result))
}
