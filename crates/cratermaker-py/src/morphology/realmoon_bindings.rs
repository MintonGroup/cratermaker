use cratermaker_components::morphology::realmoon::RealMoonCrater;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyAttributeError, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;


// Mirrors the RealMoonCrater struct in cratermaker-components and provides read-only access to its fields from Python.
pub struct PyReadonlyRealMoonCrater<'py> {
    pub id: u32,
    pub diameter: f64,
    pub radius: f64,
    pub semimajor_axis: f64,
    pub semiminor_axis: f64,
    pub orientation: f64,
    pub transient_diameter: f64,
    pub projectile_diameter: f64,
    pub projectile_velocity: f64,
    pub projectile_angle: f64,
    pub projectile_density: f64,
    pub location: (f64, f64),
    pub morphology_type: String,
    pub measured_semimajor_axis: f64,
    pub measured_semiminor_axis: f64,
    pub measured_orientation: f64,
    pub measured_diameter: f64,
    pub measured_radius: f64,
    pub measured_location: (f64, f64),
    pub time: Option<f64>,
    pub floor_elevation: f64,
    pub floor_radius: f64,
    pub wall_curvature: f64,
    pub rim_width: f64,
    pub rim_elevation: f64,
    pub rim_flank_radius: f64,
    pub rimdrop: f64,
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_offset: f64,
    pub rim_radius_rng_seed: u64,
    pub rim_flank_radius_rng_seed: u64,
    pub rim_elevation_rng_seed: u64,
    pub floor_elevation_rng_seed: u64,
    pub wall_texture_rng_seed: u64,
    pub ejecta_texture_rng_seed: u64,
    pub floor_texture_rng_seed: u64,
    pub rim_radius_psd: PyReadonlyArray1<'py, f64>,
    pub rim_elevation_psd: PyReadonlyArray1<'py, f64>,
    pub rim_flank_radius_psd: PyReadonlyArray1<'py, f64>,
}

fn getattr_optional<'py, T>(obj: &Bound<'py, PyAny>, name: &str) -> PyResult<Option<T>>
where
    T: for<'a> pyo3::FromPyObject<'a, 'py>,
    for<'a> <T as pyo3::FromPyObject<'a, 'py>>::Error: Into<pyo3::PyErr>,
{
    let py = obj.py();
    match obj.getattr(name) {
        Ok(val) => {
            if val.is_none() {
                Ok(None)
            } else {
                Ok(Some(val.extract::<T>().map_err(Into::into)?))
            }
        }
        Err(err) => {
            // If the attribute simply doesn't exist on the Python object, treat it as optional.
            if err.is_instance_of::<PyAttributeError>(py) {
                Ok(None)
            } else {
                Err(err)
            }
        }
    }
}

impl<'py> PyReadonlyRealMoonCrater<'py> {
    /// Build from a Python PyReadonlyLocalSurface object
    pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        Ok(Self {
            id: obj.getattr("id")?.extract()?,
            diameter: obj.getattr("diameter")?.extract()?,
            radius: obj.getattr("radius")?.extract()?,
            semimajor_axis: obj.getattr("semimajor_axis")?.extract()?,
            semiminor_axis: obj.getattr("semiminor_axis")?.extract()?,  
            orientation: obj.getattr("orientation")?.extract()?,
            transient_diameter: obj.getattr("transient_diameter")?.extract()?,
            projectile_diameter: obj.getattr("projectile_diameter")?.extract()?,
            projectile_velocity: obj.getattr("projectile_velocity")?.extract()?,
            projectile_angle: obj.getattr("projectile_angle")?.extract()?,
            projectile_density: obj.getattr("projectile_density")?.extract()?,
            location: obj.getattr("location")?.extract()?,
            morphology_type: obj.getattr("morphology_type")?.extract()?,
            measured_semimajor_axis: obj.getattr("measured_semimajor_axis")?.extract()?,
            measured_semiminor_axis: obj.getattr("measured_semiminor_axis")?.extract()?,
            measured_orientation: obj.getattr("measured_orientation")?.extract()?,
            measured_diameter: obj.getattr("measured_diameter")?.extract()?,
            measured_radius: obj.getattr("measured_radius")?.extract()?,
            measured_location: obj.getattr("measured_location")?.extract()?,
            time: getattr_optional(obj, "time")?,
            floor_elevation: obj.getattr("floor_elevation")?.extract()?,
            floor_radius: obj.getattr("floor_radius")?.extract()?,
            wall_curvature: obj.getattr("wall_curvature")?.extract()?,
            rim_width: obj.getattr("rim_width")?.extract()?,
            rim_elevation: obj.getattr("rim_elevation")?.extract()?,
            rim_flank_radius: obj.getattr("rim_flank_radius")?.extract()?,
            rimdrop: obj.getattr("rimdrop")?.extract()?,
            ejrim: obj.getattr("ejrim")?.extract()?,
            ejprofile: obj.getattr("ejprofile")?.extract()?,
            peak_height: obj.getattr("peak_height")?.extract()?,
            peak_width: obj.getattr("peak_width")?.extract()?,
            peak_offset: obj.getattr("peak_offset")?.extract()?,
            rim_radius_rng_seed: obj.getattr("rim_radius_rng_seed")?.extract()?,
            rim_flank_radius_rng_seed: obj.getattr("rim_flank_radius_rng_seed")?.extract()?,
            rim_elevation_rng_seed: obj.getattr("rim_elevation_rng_seed")?.extract()?,
            floor_elevation_rng_seed: obj.getattr("floor_elevation_rng_seed")?.extract()?,
            wall_texture_rng_seed: obj.getattr("wall_texture_rng_seed")?.extract()?,
            ejecta_texture_rng_seed: obj.getattr("ejecta_texture_rng_seed")?.extract()?,
            floor_texture_rng_seed: obj.getattr("floor_texture_rng_seed")?.extract()?,
            rim_radius_psd: obj.getattr("rim_radius_psd")?.extract()?,
            rim_elevation_psd: obj.getattr("rim_elevation_psd")?.extract()?,
            rim_flank_radius_psd: obj.getattr("rim_flank_radius_psd")?.extract()?,
        })
    }
    /// Convert to cratermaker-components PyReadonlyLocalSurface with array views
    pub fn as_views(&self) -> RealMoonCrater<'_> {
        RealMoonCrater{
            id: self.id,
            diameter: self.diameter,
            radius: self.radius,
            semimajor_axis: self.semimajor_axis,
            semiminor_axis: self.semiminor_axis,
            orientation: self.orientation,
            transient_diameter: self.transient_diameter,
            projectile_diameter: self.projectile_diameter,
            projectile_velocity: self.projectile_velocity,
            projectile_angle: self.projectile_angle,
            projectile_density: self.projectile_density,
            location: self.location,
            morphology_type: self.morphology_type.clone(),
            measured_semimajor_axis: self.measured_semimajor_axis,
            measured_semiminor_axis: self.measured_semiminor_axis,
            measured_orientation: self.measured_orientation,
            measured_diameter: self.measured_diameter,
            measured_radius: self.measured_radius,
            measured_location: self.measured_location,
            time: self.time,
            floor_elevation: self.floor_elevation,
            floor_radius: self.floor_radius,
            wall_curvature: self.wall_curvature,
            rim_width: self.rim_width,
            rim_elevation: self.rim_elevation,
            rim_flank_radius: self.rim_flank_radius,
            rimdrop: self.rimdrop,
            ejrim: self.ejrim,
            ejprofile: self.ejprofile,
            peak_height: self.peak_height,
            peak_width: self.peak_width,
            peak_offset: self.peak_offset,
            rim_radius_rng_seed: self.rim_radius_rng_seed,
            rim_flank_radius_rng_seed: self.rim_flank_radius_rng_seed,
            rim_elevation_rng_seed: self.rim_elevation_rng_seed,
            floor_elevation_rng_seed: self.floor_elevation_rng_seed,
            wall_texture_rng_seed: self.wall_texture_rng_seed,
            ejecta_texture_rng_seed: self.ejecta_texture_rng_seed,
            floor_texture_rng_seed: self.floor_texture_rng_seed,
            rim_radius_psd: self.rim_radius_psd.as_array(),
            rim_elevation_psd: self.rim_elevation_psd.as_array(),
            rim_flank_radius_psd: self.rim_flank_radius_psd.as_array(),
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
    bearings: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: Bound<'py, PyAny>,
    include_crater: bool,
    include_ejecta: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let crater_py = PyReadonlyRealMoonCrater::from_py(&crater)?;
    let crater_v = crater_py.as_views();
    let radial_distances_v = radial_distances.as_array();
    let bearings_v = bearings.as_array();
    let reference_elevations_v = reference_elevations.as_array();
    let result = cratermaker_components::morphology::realmoon::realmoon_profile(
        radial_distances_v,
        bearings_v,
        reference_elevations_v,
        &crater_v,
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

