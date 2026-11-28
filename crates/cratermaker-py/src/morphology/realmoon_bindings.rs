use crate::morphology::form_crater_bindings::{PyCraterRegion, run_form_craters};
use cratermaker_components::morphology::central_peak::{CentralPeak, PeakModel};
use cratermaker_components::morphology::basicmoon::BasicMoonCrater;
use cratermaker_components::morphology::form_crater::RealMoonModel;
use cratermaker_components::morphology::realmoon::{RealMoonField, WallLevels};
use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Computes the profile of a realistic lunar crater: the basic profile warped by the rim outline, plus the field of variation.
///
/// # Arguments
///
/// * `radial_distances` - Distances from the crater center (m).
/// * `bearings` - Bearings from the crater center (radians, clockwise from north).
/// * `reference_elevations` - The surface the profile is drawn relative to (m).
/// * `crater` - The basic crater.
/// * `rim` - R_rim(theta) / R at the bearings of the field.
/// * `field` - The field of variation (m), n_theta x n_rho, on bearings (i + 1/2) 360 / n_theta degrees and warped radii from 0
///   to `rho_max`.
/// * `rho_max` - The largest warped radius of the field.
/// * `include_crater`, `include_ejecta` - Which parts of the profile to include.
/// * `scarp_radius`, `scarp_height` - The variation of a multiring basin's ring scarp around the basin, as factors of its radius
///   and height (see `RealMoonMorphology.scarp_outline_data`); None or empty for a uniform scarp.
/// * `wall_nominal`, `wall_levels`, `wall_params` - The stacked levels of the wall of a complex crater (see `WallLevels`): the
///   nominal radii of the levels, their radius, scarp strength and elevation offset (3 x n_levels x n_theta), and the wall
///   turbulence (its standard deviation in m, its wavelength in rim radii and its seed); None or empty for none.
#[pyfunction]
#[pyo3(signature = (radial_distances, bearings, reference_elevations, crater, rim, field, rho_max, include_crater, include_ejecta, scarp_radius=None, scarp_height=None, wall_nominal=None, wall_levels=None, wall_params=None))]
#[allow(clippy::too_many_arguments)]
pub fn realmoon_profile<'py>(
    py: Python<'py>,
    radial_distances: PyReadonlyArray1<'py, f64>,
    bearings: PyReadonlyArray1<'py, f64>,
    reference_elevations: PyReadonlyArray1<'py, f64>,
    crater: BasicMoonCrater,
    rim: PyReadonlyArray1<'py, f64>,
    field: PyReadonlyArray2<'py, f64>,
    rho_max: f64,
    include_crater: bool,
    include_ejecta: bool,
    scarp_radius: Option<PyReadonlyArray1<'py, f64>>,
    scarp_height: Option<PyReadonlyArray1<'py, f64>>,
    wall_nominal: Option<PyReadonlyArray1<'py, f64>>,
    wall_levels: Option<PyReadonlyArray3<'py, f64>>,
    wall_params: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let empty: [f64; 0] = [];
    let data = RealMoonField {
        rim: rim.as_array(),
        field: field.as_array(),
        rho_max,
        scarp_radius: scarp_radius.as_ref().map_or(ArrayView1::from(&empty[..]), |a| a.as_array()),
        scarp_height: scarp_height.as_ref().map_or(ArrayView1::from(&empty[..]), |a| a.as_array()),
        wall: wall_view(wall_nominal.as_ref(), wall_levels.as_ref(), wall_params.as_ref())?,
    };
    let result = cratermaker_components::morphology::realmoon::realmoon_profile(
        radial_distances.as_array(),
        bearings.as_array(),
        reference_elevations.as_array(),
        &crater,
        &data,
        include_crater,
        include_ejecta,
    )
    .map_err(PyValueError::new_err)?;
    Ok(PyArray1::from_owned_array(py, result))
}

/// Forms a batch of realistic lunar craters whose regions do not overlap, in parallel (see `basicmoon_bindings.form_craters`).
///
/// # Arguments
///
/// * `craters` - The basic craters.
/// * `rims`, `fields`, `rho_max` - The rim outline and field of each crater (see `realmoon_profile`).
/// * `wall_nominals`, `wall_levels`, `wall_params` - The wall levels of each crater (see `realmoon_profile`; empty for none).
/// * The other arguments are those of `basicmoon_bindings.form_craters`.
#[pyfunction]
#[pyo3(signature = (craters, rims, fields, scarp_radii, scarp_heights, wall_nominals, wall_levels, wall_params, rho_max, regions, face_indices, node_indices, ray_seeds, face_elevation, node_elevation, ejecta_thickness, conserve_volume, critical_slope=None))]
#[allow(clippy::too_many_arguments)]
pub fn form_craters<'py>(
    py: Python<'py>,
    craters: Vec<BasicMoonCrater>,
    rims: Vec<PyReadonlyArray1<'py, f64>>,
    fields: Vec<PyReadonlyArray2<'py, f64>>,
    scarp_radii: Vec<PyReadonlyArray1<'py, f64>>,
    scarp_heights: Vec<PyReadonlyArray1<'py, f64>>,
    wall_nominals: Vec<PyReadonlyArray1<'py, f64>>,
    wall_levels: Vec<PyReadonlyArray3<'py, f64>>,
    wall_params: Vec<PyReadonlyArray1<'py, f64>>,
    rho_max: f64,
    regions: Vec<Bound<'py, PyAny>>,
    face_indices: Vec<PyReadonlyArray1<'py, i64>>,
    node_indices: Vec<PyReadonlyArray1<'py, i64>>,
    ray_seeds: Vec<Option<u64>>,
    face_elevation: Bound<'py, numpy::PyArray1<f64>>,
    node_elevation: Bound<'py, numpy::PyArray1<f64>>,
    ejecta_thickness: Bound<'py, numpy::PyArray1<f64>>,
    conserve_volume: bool,
    critical_slope: Option<f64>,
) -> PyResult<Vec<(f64, Option<f64>)>> {
    let n = craters.len();
    if rims.len() != n || fields.len() != n || scarp_radii.len() != n || scarp_heights.len() != n || wall_nominals.len() != n
        || wall_levels.len() != n || wall_params.len() != n
    {
        return Err(PyValueError::new_err(
            "craters, rims, fields, scarp outlines and wall levels must have the same length",
        ));
    }
    let models: Vec<RealMoonModel> = (0..n)
        .map(|i| {
            Ok(RealMoonModel {
                crater: &craters[i],
                field: RealMoonField {
                    rim: rims[i].as_array(),
                    field: fields[i].as_array(),
                    rho_max,
                    scarp_radius: scarp_radii[i].as_array(),
                    scarp_height: scarp_heights[i].as_array(),
                    wall: wall_view(Some(&wall_nominals[i]), Some(&wall_levels[i]), Some(&wall_params[i]))?,
                },
            })
        })
        .collect::<PyResult<_>>()?;
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

/// Builds the central peak or peak ring of a realistic lunar crater on a region (see `central_peak.rs`): the peak, less
/// the layer of the same volume removed from the floor.
///
/// # Arguments
///
/// * `region` - The Python `LocalSurface` the peak is built on, centered on the crater.
/// * `height` - The height of the highest point of the peak (m).
/// * `width` - The standard deviation of the peak's Gaussian envelope (m).
/// * `ring_radius` - The radius of the envelope's ring (m); 0 for a central peak.
/// * `floor_radius` - The radius of the floor the peak's volume is removed from (m); 0 to leave the floor unchanged.
/// * `seed` - The seed of the peak's noise fields.
/// * `center_distance`, `center_bearing` - The offset of the peak's center from the crater center (m, and degrees clockwise
///   from north); if either is None, the offset is drawn from the seed.
/// * `elongation`, `elongation_bearing` - The elongation of the peak (long to short axis ratio) and the bearing of its long axis
///   (degrees clockwise from north); if either is None, they are drawn from the seed.
/// * The remaining arguments override the parameters of the peak model (`central_peak.rs`, `PeakModel`), for its calibration;
///   by default they are Cratermaker's.
///
/// # Returns
///
/// * The change of the elevation of each face of the region (m).
#[pyfunction]
#[pyo3(signature = (region, height, width, ring_radius, floor_radius, seed, center_distance=None, center_bearing=None, elongation=None, elongation_bearing=None, threshold=None, noise_wavelength=None, octaves=None,
    persistence=None, ramp_exponent=None, critical_angle=None, angle_spread=None, angle_noise_wavelength=None, support=None, offset_sigma=None,
    elongation_sigma=None))]
#[allow(clippy::too_many_arguments)]
pub fn central_peak<'py>(
    py: Python<'py>,
    region: Bound<'py, PyAny>,
    height: f64,
    width: f64,
    ring_radius: f64,
    floor_radius: f64,
    seed: u64,
    center_distance: Option<f64>,
    center_bearing: Option<f64>,
    elongation: Option<f64>,
    elongation_bearing: Option<f64>,
    threshold: Option<f64>,
    noise_wavelength: Option<f64>,
    octaves: Option<usize>,
    persistence: Option<f64>,
    ramp_exponent: Option<f64>,
    critical_angle: Option<f64>,
    angle_spread: Option<f64>,
    angle_noise_wavelength: Option<f64>,
    support: Option<f64>,
    offset_sigma: Option<f64>,
    elongation_sigma: Option<f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let d = PeakModel::default();
    let model = PeakModel {
        threshold: threshold.unwrap_or(d.threshold),
        noise_wavelength: noise_wavelength.unwrap_or(d.noise_wavelength),
        octaves: octaves.unwrap_or(d.octaves),
        persistence: persistence.unwrap_or(d.persistence),
        ramp_exponent: ramp_exponent.unwrap_or(d.ramp_exponent),
        critical_angle: critical_angle.unwrap_or(d.critical_angle),
        angle_spread: angle_spread.unwrap_or(d.angle_spread),
        angle_noise_wavelength: angle_noise_wavelength.unwrap_or(d.angle_noise_wavelength),
        support: support.unwrap_or(d.support),
        offset_sigma: offset_sigma.unwrap_or(d.offset_sigma),
        elongation_sigma: elongation_sigma.unwrap_or(d.elongation_sigma),
    };
    let region = PyCraterRegion::from_py(&region)?;
    let peak = CentralPeak {
        height,
        width,
        ring_radius,
        floor_radius,
        center: center_distance.zip(center_bearing),
        elongation: elongation.zip(elongation_bearing),
        seed,
    };
    let view = region.as_view();
    let change = py
        .detach(|| cratermaker_components::morphology::central_peak::central_peak(&view, &peak, &model))
        .map_err(PyValueError::new_err)?;
    Ok(PyArray1::from_owned_array(py, change))
}

/// The wall levels from their arrays, or None if they are absent or empty.
fn wall_view<'a>(
    nominal: Option<&'a PyReadonlyArray1<'_, f64>>,
    levels: Option<&'a PyReadonlyArray3<'_, f64>>,
    params: Option<&'a PyReadonlyArray1<'_, f64>>,
) -> PyResult<Option<WallLevels<'a>>> {
    let (Some(nominal), Some(levels), Some(params)) = (nominal, levels, params) else {
        return Ok(None);
    };
    let nominal = nominal.as_array();
    if nominal.is_empty() {
        return Ok(None);
    }
    let params = params.as_array();
    if params.len() != 3 {
        return Err(PyValueError::new_err("wall_params must hold the turbulence amplitude, wavelength and seed"));
    }
    Ok(Some(WallLevels {
        nominal,
        levels: levels.as_array(),
        turbulence: params[0],
        wavelength: params[1],
        seed: params[2] as u32,
    }))
}
