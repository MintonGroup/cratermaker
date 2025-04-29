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
    let radial_distances = r_array.as_array();
    let reference_elevations = reference_elevation_array.as_array();
    if radial_distances.len() != reference_elevations.len() {
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

    let ninc = radial_distances.iter().filter(|&&x| x <= radius).count();
    let meanref = if ninc == 0 {
        *radial_distances
            .iter()
            .zip(reference_elevations)
            .min_by(|(&radius_a, _), (&radius_b, _)| radius_a.partial_cmp(&radius_b).unwrap())
            .unwrap()
            .1
    } else {
        radial_distances
            .iter()
            .zip(reference_elevations)
            .filter(|(&r, _)| r <= radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let min_elevation = meanref - floordepth;

    Ok(PyArray1::from_iter(
        py,
        reference_elevations
            .iter()
            .zip(radial_distances)
            .map(|(&elevation, &radial_distance)| {
                let r = radial_distance / radius;
                (
                    if r >= 1.0 {
                        elevation + (rimheight - ejrim) * (r.powf(-RIMDROP))
                    } else {
                        elevation + c0 + c1 * r + c2 * r.powi(2) + c3 * r.powi(3)
                    },
                    radial_distance,
                )
            })
            .map(|(elevation, radial_distance)| {
                if radial_distance <= radius {
                    elevation.max(min_elevation)
                } else {
                    elevation
                }
            }),
    ))
}
