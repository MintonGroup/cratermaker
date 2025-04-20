use numpy::{ndarray, PyReadonlyArray1};
use pyo3::prelude::*;

#[pyfunction]
pub fn distribution<'py>(
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
    ejecta_truncation: f64,
    dorays: bool,
) {
}

#[pyfunction]
pub fn profile<'py>(radial_distance: PyReadonlyArray1<'py, f64>, crater_diameter: f64, ejrim: f64) {
}

#[pyfunction]
pub fn ray_intensity<'py>(
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
    ejecta_truncation: f64,
) {
}
