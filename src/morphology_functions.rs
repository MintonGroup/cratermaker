use numpy::{PyReadonlyArray1, PyReadwriteArray1};
use pyo3::prelude::*;

use crate::ejecta_functions;

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject)]
pub struct Crater {
    pub final_diameter: f64,
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
    pub dorays: bool,
}

/// View into a region of the surface mesh, consisting of node and face indices.
///
/// Used to localize crater effects to a subset of the full mesh.
#[derive(FromPyObject)]
pub struct SurfaceView<'py> {
    pub node_indices: PyReadonlyArray1<'py, i64>,
    pub face_indices: PyReadonlyArray1<'py, i64>,
}

/// Applies ejecta thickness and ray modulation to a regional surface mesh.
///
/// Given radial distance and bearing values for each node and face in a selected mesh region,
/// this function calculates ejecta thickness using a radial profile and optionally modulates
/// it with ray patterns. Elevations are updated in-place for affected surface regions.
///
/// # Arguments
///
/// * `morphology` - Parameters controlling ejecta geometry and ray usage.
/// * `region_view` - Mesh indices of affected nodes and faces.
/// * `node_crater_distance` - Radial distances for nodes from crater center.
/// * `face_crater_distance` - Radial distances for faces from crater center.
/// * `node_crater_bearing` - Angular bearings for nodes from crater center.
/// * `face_crater_bearing` - Angular bearings for faces from crater center.
/// * `ejecta_truncation` - Maximum extent of ejecta distribution.
/// * `node_elevation` - Elevation values of mesh nodes (modified in-place).
/// * `face_elevation` - Elevation values of mesh faces (modified in-place).
/// * `ejecta_thickness` - Accumulated ejecta thickness per face (modified in-place).
/// * `ray_intensity` - Ray modulation values per face (modified in-place if `dorays` is true).
///
/// # Returns
///
/// * `Ok(())` on success.
#[pyfunction]
pub fn form_ejecta<'py>(
    morphology: SimpleMoonMorphology,
    region_view: SurfaceView,
    node_crater_distance: PyReadonlyArray1<f64>,
    face_crater_distance: PyReadonlyArray1<f64>,
    node_crater_bearing: PyReadonlyArray1<f64>,
    face_crater_bearing: PyReadonlyArray1<f64>,
    ejecta_truncation: f64,
    mut node_elevation: PyReadwriteArray1<f64>,
    mut face_elevation: PyReadwriteArray1<f64>,
    mut ejecta_thickness: PyReadwriteArray1<f64>,
    mut ray_intensity: PyReadwriteArray1<f64>,
) -> PyResult<()> {
    let node_thickness = ejecta_functions::distribution_internal(
        node_crater_distance.as_array(),
        node_crater_bearing.as_array(),
        morphology.crater.final_diameter,
        morphology.ejrim,
        ejecta_truncation,
        morphology.dorays,
    )?;
    let face_thickness = ejecta_functions::distribution_internal(
        face_crater_distance.as_array(),
        face_crater_bearing.as_array(),
        morphology.crater.final_diameter,
        morphology.ejrim,
        ejecta_truncation,
        morphology.dorays,
    )?;

    let mut node_elevation = node_elevation.as_array_mut();
    let mut face_elevation = face_elevation.as_array_mut();
    let mut ejecta_thickness = ejecta_thickness.as_array_mut();
    let mut ray_intensity = ray_intensity.as_array_mut();

    for (&idx, new_value) in region_view
        .node_indices
        .as_array()
        .iter()
        .zip(node_thickness)
    {
        node_elevation[idx as usize] += new_value;
    }
    for (&idx, new_value) in region_view
        .face_indices
        .as_array()
        .iter()
        .zip(face_thickness)
    {
        face_elevation[idx as usize] += new_value;
        ejecta_thickness[idx as usize] += new_value;
    }

    if morphology.dorays {
        let face_intensity = ejecta_functions::ray_intensity_internal(
            face_crater_distance.as_array(),
            face_crater_bearing.as_array(),
            morphology.crater.final_diameter,
            ejecta_truncation,
        )?;
        for (&idx, new_value) in region_view
            .face_indices
            .as_array()
            .iter()
            .zip(face_intensity)
        {
            ray_intensity[idx as usize] += new_value;
        }
    }

    Ok(())
}
