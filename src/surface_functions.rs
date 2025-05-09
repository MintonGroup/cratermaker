use std::f64::consts::PI;

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use ndarray::Zip;
use rayon::iter::{IntoParallelIterator,ParallelIterator};

#[inline]
fn positive_mod(x: f64, m: f64) -> f64 {
    ((x % m) + m) % m
}

/// Computes the initial bearing (forward azimuth) from a fixed point to each of a set of destination points.
///
/// The bearing is calculated on a spherical surface using great-circle paths and returned in radians,
/// normalized to the range [0, 2π).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `lon1` - Longitude of the reference point, in radians.
/// * `lat1` - Latitude of the reference point, in radians.
/// * `lon2` - Longitudes of destination points, in radians.
/// * `lat2` - Latitudes of destination points, in radians.
///
/// # Returns
///
/// * A NumPy array of initial bearing angles (radians), one for each (lon2, lat2) pair.
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
    let mut result = ndarray::Array1::<f64>::zeros(lon2.len());

    Zip::from(&mut result)
        .and(lon2)
        .and(lat2)
        .into_par_iter()
        .for_each(|(out, &lon2, &lat2)| {
            // Calculate differences in coordinates
            let dlon = positive_mod(lon2 - lon1 + PI, 2.0 * PI) - PI;

            // Haversine formula calculations
            let x = dlon.sin() * lat2.cos();
            let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
            let initial_bearing = f64::atan2(x, y);

            // Normalize bearing to 0 to 2*pi
            *out = (initial_bearing + 2.0 * PI) % (2.0 * PI);
        });

    Ok(PyArray1::from_owned_array(py, result))
}
