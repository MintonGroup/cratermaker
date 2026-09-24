use cratermaker_components::morphology::realmoon::RealMoonCrater;
use interp::{InterpMode, interp_slice};
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;
use std::f64::consts::TAU;

// Mirrors the PSD1D struct in cratermaker-components and provides read-only access to its fields from Python.
pub struct PyPSD1D<'py> {
    pub nprofile: usize,
    pub npsd: usize,
    pub normalization_length: f64,
    pub mean: f64,
    pub pix: f64,
    pub wavelength: PyReadonlyArray1<'py, f64>,
    pub amplitude: PyReadonlyArray1<'py, f64>,
    pub phase: PyReadonlyArray1<'py, f64>,
}

impl<'py> PyPSD1D<'py> {
    // Extracts the Python attributes and matches them up with the corresponding Rust struct components.
    pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        let nprofile: usize = obj.getattr("nprofile")?.extract()?;
        let npsd: usize = obj.getattr("nprofile")?.extract()?;
        let normalization_length: f64 = obj.getattr("normalization_length")?.extract()?;
        let mean: f64 = obj.getattr("mean")?.extract()?;
        let pix: f64 = obj.getattr("pix")?.extract()?;
        Ok(Self {
            nprofile,
            npsd,
            normalization_length,
            mean,
            pix,
            wavelength: obj.getattr("wavelength")?.extract()?,
            amplitude: obj.getattr("amplitude")?.extract()?,
            phase: obj.getattr("phase")?.extract()?,
        })
    }

    // Converts all of the PyArray objects inside the struct to memory views of the arrays
    pub fn as_views(&self) -> cratermaker_components::morphology::realmoon::PSD1DView<'_> {
        cratermaker_components::morphology::realmoon::PSD1DView {
            nprofile: self.nprofile,
            npsd: self.npsd,
            normalization_length: self.normalization_length,
            mean: self.mean,
            pix: self.pix,
            wavelength: self.wavelength.as_array(),
            amplitude: self.amplitude.as_array(),
            phase: self.phase.as_array(),
        }
    }
}
// Mirrors the RealMoonCrater struct in cratermaker-components and provides read-only access to its fields from Python.
pub struct PyReadonlyRealMoonCrater<'py> {
    pub diameter: f64,
    pub radius: f64,
    pub floor_elevation: f64,
    pub floor_radius: f64,
    pub wall_curvature: f64,
    pub floor_blend: f64,
    pub rim_width: f64,
    pub rim_height: f64,
    pub frac_ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_ring_radius: f64,
    pub peak_center_distance: f64,
    pub peak_center_bearing: f64,
    pub rim_radius_psd: PyPSD1D<'py>,
    pub floor_radius_psd: PyPSD1D<'py>,
}
impl<'py> PyReadonlyRealMoonCrater<'py> {
    /// Build from a Python PyReadonlyCrater object
    pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        let rim_psd_any = obj.getattr("rim_radius_psd")?;
        let floor_psd_any = obj.getattr("floor_radius_psd")?;

        let rim_radius_psd = PyPSD1D::from_py(&rim_psd_any)?;
        let floor_radius_psd = PyPSD1D::from_py(&floor_psd_any)?;
        Ok(Self {
            diameter: obj.getattr("diameter")?.extract()?,
            radius: obj.getattr("radius")?.extract()?,
            floor_elevation: obj.getattr("floor_elevation")?.extract()?,
            floor_radius: obj.getattr("floor_radius")?.extract()?,
            wall_curvature: obj.getattr("wall_curvature")?.extract()?,
            floor_blend: obj.getattr("floor_blend")?.extract()?,
            rim_width: obj.getattr("rim_width")?.extract()?,
            rim_height: obj.getattr("rim_height")?.extract()?,
            frac_ejrim: obj.getattr("frac_ejrim")?.extract()?,
            ejprofile: obj.getattr("ejprofile")?.extract()?,
            peak_height: obj.getattr("peak_height")?.extract()?,
            peak_width: obj.getattr("peak_width")?.extract()?,
            peak_ring_radius: obj.getattr("peak_ring_radius")?.extract()?,
            peak_center_distance: obj.getattr("peak_center_distance")?.extract()?,
            peak_center_bearing: obj.getattr("peak_center_bearing")?.extract()?,
            rim_radius_psd,
            floor_radius_psd,
        })
    }
    /// Convert to cratermaker-components RealMoonCrater struct with array views instead of readonly PyAarrays
    pub fn as_views(&self) -> RealMoonCrater<'_> {
        RealMoonCrater {
            diameter: self.diameter,
            radius: self.radius,
            floor_elevation: self.floor_elevation,
            floor_radius: self.floor_radius,
            wall_curvature: self.wall_curvature,
            floor_blend: self.floor_blend,
            rim_width: self.rim_width,
            rim_height: self.rim_height,
            frac_ejrim: self.frac_ejrim,
            ejprofile: self.ejprofile,
            peak_height: self.peak_height,
            peak_width: self.peak_width,
            peak_ring_radius: self.peak_ring_radius,
            peak_center_distance: self.peak_center_distance,
            peak_center_bearing: self.peak_center_bearing,
            rim_radius_psd: self.rim_radius_psd.as_views(),
            floor_radius_psd: self.floor_radius_psd.as_views(),
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
    mean: f64,
    nprofile: usize,
    add_noise: bool,
    rng_seed: u64,
) -> PyResult<(
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let (wavelength, amplitude, phase) =
        cratermaker_components::morphology::realmoon::get_1d_psd_from_control_points(
            &control_points,
            mean,
            nprofile,
            add_noise,
            rng_seed,
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    let wavelength = PyArray1::from_owned_array(py, wavelength);
    let amplitude = PyArray1::from_owned_array(py, amplitude);
    let phase = PyArray1::from_owned_array(py, phase);

    Ok((wavelength, amplitude, phase))
}

///
/// Generates a surface profile based on a 1D amplitude spectral density (PSD) and optional phase information, simulating a crater surface with specified roughness characteristics.
///
/// # Arguments
/// * `psd` - The PSD struct containing wavelength, amplitude, phase, and mean values.
/// * `theta` - A 1D array of polar angles (in radians) at which to compute the profile, typically ranging from 0 to 2π.
///
/// # Returns
///
/// * A 1D array of values corresponding to the input theta angles, representing the linear profile generated from the PSD.
///
#[pyfunction]
pub fn profile_from_psd<'py>(
    py: Python<'py>,
    psd: Bound<'py, PyAny>,
    theta: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let psd_py = PyPSD1D::from_py(&psd)?;
    let psd_v = psd_py.as_views();
    let theta_v = theta.as_array();

    let (psd_theta, profile) =
        cratermaker_components::morphology::realmoon::profile_from_psd(&psd_v);
    let n = psd_theta.len();
    if n == 0 {
        return Ok(PyArray1::from_vec(py, Vec::new()));
    }

    let mut padded_theta = Vec::with_capacity(n + 2);
    let mut padded_profile = Vec::with_capacity(n + 2);
    padded_theta.push(psd_theta[n - 1] - TAU);
    padded_profile.push(profile[n - 1]);
    padded_theta.extend_from_slice(&psd_theta);
    padded_profile.extend_from_slice(&profile);
    padded_theta.push(psd_theta[0] + TAU);
    padded_profile.push(profile[0]);

    let result = interp_slice(
        &padded_theta,
        &padded_profile,
        &theta_v.to_vec(),
        &InterpMode::default(),
    );
    Ok(PyArray1::from_vec(py, result))
}

#[pyfunction]
pub fn psd_from_profile<'py>(
    py: Python<'py>,
    profile: PyReadonlyArray1<'py, f64>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let profile_v = profile.as_array();

    let (amplitude, phase) =
        cratermaker_components::morphology::realmoon::psd_from_profile(&profile_v);

    Ok((
        PyArray1::from_vec(py, amplitude),
        PyArray1::from_vec(py, phase),
    ))
}
