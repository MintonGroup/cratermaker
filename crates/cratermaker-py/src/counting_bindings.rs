use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::ndarray::prelude::*;
use numpy::{PyReadonlyArray1,PyReadonlyArray2,PyArray1};
use crate::surface_bindings::PyReadonlyLocalSurface;

#[pyfunction]
pub fn radial_distance_to_ellipse<'py>(
    py: Python<'py>,
    x: PyReadonlyArray1<'py, f64>,
    y: PyReadonlyArray1<'py, f64>,
    a: f64, 
    b: f64,
    orientation: f64, 
    x0: f64, 
    y0:f64
)-> PyResult<Bound<'py, PyArray1<f64>>> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let result = cratermaker_core::counting::radial_distance_to_ellipse(
            &x_v,
            &y_v,
            a,
            b,
            orientation,
            x0,
            y0
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}


#[pyfunction]
pub fn fit_one_ellipse<'py>(
    py: Python<'py>,
    x: PyReadonlyArray1<'py, f64>,
    y: PyReadonlyArray1<'py, f64>,
    weights: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let weights_v = weights.as_array();

    let (x0, y0, a, b, orientation, wrms) = cratermaker_core::counting::fit_one_ellipse(
            x_v, 
            y_v, 
            weights_v
        )
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to fit ellipse: {}", e)))?;

    let result = Array1::from_vec(vec![x0, y0, a, b, orientation, wrms]);
    Ok(PyArray1::from_owned_array(py, result))
}

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
    ap: f64,
    bp: f64,
    orientation: f64,
    quantile: f64,
    distmult: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();
    let result = cratermaker_core::counting::score_rim(
            &region_v,
            x0,
            y0,
            ap,
            bp,
            orientation,
            quantile,
            distmult,
            gradmult,
            curvmult,
            heightmult
        )
        .map_err(|msg| PyErr::new::<PyValueError, _>(msg))?;
    Ok(PyArray1::from_owned_array(py, result))
}