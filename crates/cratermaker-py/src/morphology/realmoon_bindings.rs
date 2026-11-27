use crate::morphology::form_crater_bindings::run_form_craters;
use cratermaker_components::morphology::basicmoon::BasicMoonCrater;
use cratermaker_components::morphology::form_crater::RealMoonModel;
use cratermaker_components::morphology::realmoon::RealMoonField;
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
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
#[pyfunction]
#[pyo3(signature = (radial_distances, bearings, reference_elevations, crater, rim, field, rho_max, include_crater, include_ejecta))]
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
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let data = RealMoonField {
        rim: rim.as_array(),
        field: field.as_array(),
        rho_max,
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
/// * The other arguments are those of `basicmoon_bindings.form_craters`.
#[pyfunction]
#[pyo3(signature = (craters, rims, fields, rho_max, regions, face_indices, node_indices, ray_seeds, face_elevation, node_elevation, ejecta_thickness, conserve_volume, critical_slope=None))]
#[allow(clippy::too_many_arguments)]
pub fn form_craters<'py>(
    py: Python<'py>,
    craters: Vec<BasicMoonCrater>,
    rims: Vec<PyReadonlyArray1<'py, f64>>,
    fields: Vec<PyReadonlyArray2<'py, f64>>,
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
    if rims.len() != n || fields.len() != n {
        return Err(PyValueError::new_err(
            "craters, rims and fields must have the same length",
        ));
    }
    let models: Vec<RealMoonModel> = (0..n)
        .map(|i| RealMoonModel {
            crater: &craters[i],
            field: RealMoonField {
                rim: rims[i].as_array(),
                field: fields[i].as_array(),
                rho_max,
            },
        })
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
