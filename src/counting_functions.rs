use pyo3::{prelude::*, types::PyDict, PyRef};
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};

#[pyfunction]
pub fn tally_m19<'py>(
    py: Python<'py>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    id_array: PyReadonlyArray2<'py, u32>, 
    observed: Py<PyDict>,
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let face_elevation = face_elevation.as_array();
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
