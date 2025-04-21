use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{prelude::*, types::PyDict};

#[pyfunction]
#[pyo3(signature = (_model, _x_array, _y_array, _z_array, _num_octaves, _anchor, **_kwargs))]
pub fn apply_noise<'py>(
    _py: Python<'_>,
    _model: &str,
    _x_array: PyReadonlyArray1<'py, f64>,
    _y_array: PyReadonlyArray1<'py, f64>,
    _z_array: PyReadonlyArray1<'py, f64>,
    _num_octaves: i32,
    _anchor: PyReadonlyArray2<'py, f64>,
    _kwargs: Option<&Bound<'_, PyDict>>,
) {
    todo!()
}

#[pyfunction]
pub fn realistic_crater<'py>(
    _radial_distance: PyReadonlyArray1<'py, f64>,
    _initial_bearing: PyReadonlyArray1<'py, f64>,
    _reference_elevation_array: PyReadonlyArray1<'py, f64>,
    _rim_1d_psd: PyReadonlyArray1<'py, f64>,
    _floor_1d_psd: PyReadonlyArray1<'py, f64>,
    _ejecta_1d_psd: PyReadonlyArray1<'py, f64>,
    _floor_2d_psd: PyReadonlyArray2<'py, f64>,
    _ejecta_2d_psd: PyReadonlyArray2<'py, f64>,
    _wall_2d_psd: PyReadonlyArray2<'py, f64>,
) {
    todo!()
}
