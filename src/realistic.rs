use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{prelude::*, types::PyDict};

#[pyfunction]
#[pyo3(signature = (x_array, y_array, z_array, num_octaves, anchor, **kwargs))]
pub fn apply_noise<'py>(
    py: Python<'_>,
    x_array: PyReadonlyArray1<'py, f64>,
    y_array: PyReadonlyArray1<'py, f64>,
    z_array: PyReadonlyArray1<'py, f64>,
    num_octaves: i32,
    anchor: PyReadonlyArray2<'py, f64>,
    kwargs: Option<&Bound<'_, PyDict>>,
) {
}

#[pyfunction]
pub fn realistic_crater<'py>(
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    reference_elevation_array: PyReadonlyArray1<'py, f64>,
    rim_1d_psd: PyReadonlyArray1<'py, f64>,
    floor_1d_psd: PyReadonlyArray1<'py, f64>,
    ejecta_1d_psd: PyReadonlyArray1<'py, f64>,
    floor_2d_psd: PyReadonlyArray2<'py, f64>,
    ejecta_2d_psd: PyReadonlyArray2<'py, f64>,
    wall_2d_psd: PyReadonlyArray2<'py, f64>,
) {
}
