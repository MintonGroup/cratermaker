use std::{f64::consts::PI, iter::zip};

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use rayon::prelude::*;

#[pyfunction]
pub fn calculate_initial_bearing<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let lon2 = lon2.as_array();
    let lat2 = lat2.as_array();

    let res = zip(lon2, lat2)
        .par_bridge()
        .map(|(&lon2, &lat2)| {
            // Calculate differences in coordinates
            let dlon = (lon2 - lon1 + PI) % (2.0 * PI) - PI;

            // Haversine formula calculations
            let x = dlon.sin() * lat2.cos();
            let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
            let initial_bearing = f64::atan2(x, y);

            // Normalize bearing to 0 to 2*pi
            (initial_bearing + 2.0 * PI) % (2.0 * PI)
        })
        .collect();

    Ok(PyArray1::from_vec(py, res))
}
