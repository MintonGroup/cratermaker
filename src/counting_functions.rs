use std::f64::consts::FRAC_PI_2;
use pyo3::{prelude::*}; 
use numpy::{PyArray1, PyReadonlyArray2, PyReadonlyArray1};
use ndarray::Zip;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

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
    // Placeholder implementation
    let phi = orientation - FRAC_PI_2;
    let x = x.as_array();
    let y = y.as_array();
    let mut result = ndarray::Array1::<f64>::zeros(x.len());

    Zip:: from(&mut result)
    .and_broadcast(x).and(y)
    .into_par_iter()
    .for_each(|(res, &xi, &yi)| {
        let dx = xi - x0;
        let dy = yi - y0;
        let r = (dx * dx + dy * dy).sqrt();
        let theta = dy.atan2(dx);
        let alpha = theta - phi;
        let ca = alpha.cos();
        let sa = alpha.sin();
        *res = r - (a * b) / ((b * ca).powi(2) + (a * sa).powi(2)).sqrt();
    });
    Ok(PyArray1::from_owned_array(py, result))
}
