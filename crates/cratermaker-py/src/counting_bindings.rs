use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{PyReadonlyArray1,PyReadonlyArray2,PyArray1};
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
pub fn fit_one_ellipse<'py>(
    _py: Python<'py>,
    x: PyReadonlyArray1<'_,f64>, 
    y: PyReadonlyArray1<'_,f64>, 
    weights: PyReadonlyArray1<'_,f64>
) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    let x_v = x.as_array();
    let y_v = y.as_array();
    let weights_v = weights.as_array();
    let (x0, y0, ap, bp, orientation, wrms) = cratermaker_components::counting::fit_one_ellipse(
            x_v,
            y_v,
            weights_v,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok((x0, y0, ap, bp, orientation, wrms))
}


#[pyfunction]
pub fn score_rim<'py>(
    py: Python<'py>,
    surface: &Bound<'py, PyAny>,
    crater: Crater, 
    quantile: f64,
    distmult: f64,
    gradmult: f64,
    curvmult: f64,
    heightmult: f64,
) -> PyResult<Bound<'py, PyAny>>  {
    const _EXTENT_RADIUS_RATIO: f64 = 2.0;
    let region = surface.call_method1("extract_region",(crater.location, _EXTENT_RADIUS_RATIO * crater.radius))?;
    let transformer = region.getattr("from_surface").unwrap();
    let x0y0 = transformer.call_method1("transform",(crater.location.0, crater.location.1))?;
    let (x0, y0): (f64, f64) = x0y0.extract()?;
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

    region.call_method1("add_data",("rimscore", PyArray1::from_owned_array(py, result.clone()), "Rim Score", "dimensionless", true, f64::NAN))?;

    Ok(region)
}


#[pyfunction]
pub fn fit_rim<'py>(
    _py: Python<'py>,
    region: Bound<'py, PyAny>,
    crater: Crater, 
    tol: f64,
    nloops: usize,
    score_quantile: f64,
) -> PyResult<(f64, f64, f64, f64, f64)> {
    let region_py = PyReadonlyLocalSurface::from_local_surface(&region)?;
    let region_v = region_py.as_views();
    let (a, b, orientation, lon, lat) = cratermaker_components::counting::fit_rim(
            &region_v,
            &crater,
            tol,
            nloops,
            score_quantile
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok((a, b, orientation, lon, lat))
}