use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};

const A: f64 = 4.0 / 11.0;
const B: f64 = -32.0 / 187.0;
const RIMDROP: f64 = 4.20; // The exponent for the uplifted rim dropoff.

#[pyfunction]
pub fn profile<'py>(
    py: Python<'py>,
    r_array: PyReadonlyArray1<'py, f64>,
    reference_elevation_array: PyReadonlyArray1<'py, f64>,
    diameter: f64,
    floordepth: f64,
    floordiam: f64,
    rimheight: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let r_array = r_array.as_array();
    let reference_elevation_array = reference_elevation_array.as_array();
    if r_array.len() != reference_elevation_array.len() {
        return Err(PyValueError::new_err(
            "input arrays must have the same length",
        ));
    }

    // Calculate the floor radius relative to the final crater radius
    let flrad = floordiam / diameter;
    let radius = diameter / 2.0;

    // Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
    let c1 = (-floordepth - rimheight)
        / (flrad - 1.0 + A * (flrad.powi(2) - 1.0) + B * (flrad.powi(3) - 1.0));
    let c0 = rimheight - c1 * (1.0 + A + B);
    let c2 = A * c1;
    let c3 = B * c1;

    let ninc = r_array.iter().filter(|&&x| x <= radius).count();
    let meanref = if ninc == 0 {
        *r_array
            .iter()
            .zip(reference_elevation_array)
            .max_by(|(&radius_a, _), (&radius_b, _)| radius_a.partial_cmp(&radius_b).unwrap())
            .unwrap()
            .1
    } else {
        r_array.iter().filter(|&&x| x <= radius).sum::<f64>() / ninc as f64
    };
    let min_elevation = meanref - floordepth;

    Ok(PyArray1::from_iter(
        py,
        reference_elevation_array
            .iter()
            .zip(r_array)
            .map(|(&elevation, &r)| {
                (
                    if r >= 1.0 {
                        elevation + (rimheight - ejrim) * (r.powf(-RIMDROP))
                    } else {
                        elevation + c0 + c1 * r + c2 * r.powi(2) + c3 * r.powi(3)
                    },
                    r,
                )
            })
            .map(|(elevation, r)| {
                if r < radius {
                    elevation.max(min_elevation)
                } else {
                    elevation
                }
            }),
    ))
}
