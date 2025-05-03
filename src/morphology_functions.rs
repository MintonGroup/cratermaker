use numpy::{PyReadonlyArray1, PyReadwriteArray1};
use pyo3::prelude::*;

use crate::ejecta_functions;

#[derive(FromPyObject)]
pub struct Crater {
    pub final_diameter: f64,
}

#[derive(FromPyObject)]
pub struct SimpleMoonMorphology {
    pub floordepth: f64,
    pub floor_diameter: f64,
    pub rimheight: f64,
    pub ejrim: f64,
    pub crater: Crater,
    pub dorays: bool,
}

#[derive(FromPyObject)]
pub struct SurfaceView<'py> {
    pub node_indices: PyReadonlyArray1<'py, i64>,
    pub face_indices: PyReadonlyArray1<'py, i64>,
}

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
