use crate::surface::surface_bindings::PyReadonlyLocalSurface;
use cratermaker_components::crater::Crater;
use pyo3::prelude::*;
const _FITTING_RADIUS_RATIO: f64 = 2.0;

#[pyfunction]
pub fn measure_rim_height<'py>(
    _py: Python<'py>,
    region: &Bound<'py, PyAny>,
    crater: Crater,
) -> PyResult<f64> {
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let rim_height =
        match cratermaker_components::counting::measure_rim_height(&region_v, &crater) {
            Ok(v) => v,
            Err(_) => return Ok(-f64::MAX),
        };

    Ok(rim_height)
}

#[pyfunction]
pub fn measure_floor_elevation<'py>(
    _py: Python<'py>,
    region: &Bound<'py, PyAny>,
    crater: Crater,
) -> PyResult<f64> {
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();

    let floor_elevation =
        match cratermaker_components::counting::measure_floor_elevation(&region_v, &crater) {
            Ok(v) => v,
            Err(_) => return Ok(f64::MAX),
        };

    Ok(floor_elevation)
}


