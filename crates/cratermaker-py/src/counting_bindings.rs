use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray2,PyArray1};
use crate::surface_bindings::PyReadonlyLocalSurface;
use cratermaker_components::crater::Crater;


#[pyfunction]
pub fn tally<'py>(
    py: Python<'py>,
    //face_elevation: PyReadonlyArray1<'py, f64>,
    id_array: PyReadonlyArray2<'py, u32>, 
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let id_array = id_array.as_array();
    let mut id_vec = Vec::with_capacity(id_array.len());
    // for (key, value) in observed.as_ref(py).iter() {
    //     let id: u32 = key.extract()?;
    //     let crater: &PyDict = value.downcast::<PyDict>()?;

    //     let final_diameter: Option<f64> = crater.get_item("final_diameter").and_then(|v| v.extract().ok());
    //     let location: Option<(f64, f64)> = crater.get_item("location").and_then(|v| v.extract().ok());

    // }


    for id in id_array.iter() {
        id_vec.push(*id);
    }
    let id_array_flat = PyArray1::from_vec(py, id_vec);
    Ok(id_array_flat)
}
#[pyfunction]
pub fn score_rim<'py>(
    py: Python<'py>,
    region: Bound<'py, PyAny>,
    x0: f64, 
    y0: f64,
    crater: Crater, 
    quantile: f64,
    distmult: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();
    let result = cratermaker_components::counting::score_rim(
            &region_v,
            x0,
            y0,
            &crater,
            quantile,
            distmult,
            gradmult,
            curvmult,
            heightmult
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}