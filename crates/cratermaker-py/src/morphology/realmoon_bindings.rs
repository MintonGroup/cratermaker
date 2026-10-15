use cratermaker_components::morphology::realmoon::RealMoonCrater;
use interp::{InterpMode, interp_slice};
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;
use std::f64::consts::TAU;
use paste::paste;

// Mirrors the PSD1D struct in cratermaker-components and provides read-only access to its fields from Python.
pub struct PyPSD1D<'py> {
    pub nprofile: usize,
    pub nfreq: usize,
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
        let nfreq: usize = obj.getattr("nfreq")?.extract()?;
        let normalization_length: f64 = obj.getattr("normalization_length")?.extract()?;
        let mean: f64 = obj.getattr("mean")?.extract()?;
        let pix: f64 = obj.getattr("pix")?.extract()?;
        Ok(Self {
            nprofile,
            nfreq,
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
            nfreq: self.nfreq,
            normalization_length: self.normalization_length,
            mean: self.mean,
            pix: self.pix,
            wavelength: self.wavelength.as_array(),
            amplitude: self.amplitude.as_array(),
            phase: self.phase.as_array(),
        }
    }
}

// This will build a boilerplate struct that contains the PSD objects for all of the Python crater object properties with the pattern [name]_psd.
macro_rules! define_crater_structs {
    (
        floats: [ $($float_field:ident),* ],
        psds: [ $($psd_field:ident),* ]
    ) => {
        paste! {
            pub struct PyReadonlyRealMoonCrater<'py> {
                $( pub $float_field: f64, )*
                $( pub [<$psd_field _psd>]: PyPSD1D<'py>, )*
            } 

            impl<'py> PyReadonlyRealMoonCrater<'py> {
                pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
                    Ok(Self {
                        $( $float_field: obj.getattr(stringify!($float_field))?.extract()?, )*
                        $(
                            [<$psd_field _psd>]: PyPSD1D::from_py(
                                &obj.getattr(stringify!([<$psd_field _psd>]))?
                            )?,
                        )*
                    })
                }

                pub fn as_views(&self) -> cratermaker_components::morphology::realmoon::RealMoonCrater<'_> {
                    cratermaker_components::morphology::realmoon::RealMoonCrater {
                        $( $float_field: self.$float_field, )*
                        $(
                            [<$psd_field _psd>]: self.[<$psd_field _psd>].as_views(),
                        )*
                    }
                }
            }
        }
    };
}

define_crater_structs! {
    floats: [
        diameter, radius, floor_elevation, floor_radius, wall_curvature,
        floor_blend, rim_width, rim_height, frac_ejrim, ejprofile,
        peak_height, peak_width, peak_ring_radius, peak_center_distance,
        peak_center_bearing
    ],
    psds: [ radius, floor_radius, wall_curvature, floor_blend, rim_width, rim_height, frac_ejrim, ejprofile ]
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
/// * `bearing` - A 1D array of polar angles (in radians) at which to compute the profile, typically ranging from 0 to 2π.
///
/// # Returns
///
/// * A 1D array of values corresponding to the input bearing angles, representing the linear profile generated from the PSD.
///
#[pyfunction]
pub fn profile_from_psd<'py>(
    py: Python<'py>,
    psd: Bound<'py, PyAny>,
    bearing: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let psd_py = PyPSD1D::from_py(&psd)?;
    let psd_v = psd_py.as_views();
    let bearing_v = bearing.as_array();

    let (psd_bearing, profile) =
        cratermaker_components::morphology::realmoon::profile_from_psd(&psd_v);
    let n = psd_bearing.len();
    if n == 0 {
        return Ok(PyArray1::from_vec(py, Vec::new()));
    }

    let mut padded_bearing = Vec::with_capacity(n + 2);
    let mut padded_profile = Vec::with_capacity(n + 2);
    padded_bearing.push(psd_bearing[n - 1] - TAU);
    padded_profile.push(profile[n - 1]);
    padded_bearing.extend_from_slice(&psd_bearing);
    padded_profile.extend_from_slice(&profile);
    padded_bearing.push(psd_bearing[0] + TAU);
    padded_profile.push(profile[0]);

    let result = interp_slice(
        &padded_bearing,
        &padded_profile,
        &bearing_v.to_vec(),
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
