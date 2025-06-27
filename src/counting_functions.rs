use pyo3::prelude::*;
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};

#[pyfunction]
pub fn tally_m19<'py>(
    py: Python<'py>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    unique_ids: PyReadonlyArray1<'py, u32>,
    id_array: PyReadonlyArray2<'py, u32>, 
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let face_elevation = face_elevation.as_array();
    let unique_ids = unique_ids.as_array();
    let id_array = id_array.as_array();
    // flatten the id array, then store it into a PyArray1<u32> variable and return that. This is a placeholder.
    let mut id_vec = Vec::with_capacity(id_array.len());
    for id in id_array.iter() {
        id_vec.push(*id);
    }
    let id_array_flat = PyArray1::from_vec(py, id_vec);
    Ok(id_array_flat)
}
