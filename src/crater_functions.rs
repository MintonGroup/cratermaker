use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};

use crate::RIMDROP;

const A: f64 = 4.0 / 11.0;
const B: f64 = -32.0 / 187.0;


/// Calculates the elevation of a crater as a function of distance from the center.
///
/// This function applies a polynomial profile for the crater interior (r < 1.0) and a rim dropoff
/// function for the exterior (r ≥ 1.0). It is based on the crater profile model described in 
/// Fassett and Thomson (2014).
/// 
/// Fassett, C.I., Thomson, B.J., 2014. Crater degradation on the lunar maria: Topographic diffusion and 
/// the rate of erosion on the Moon. J. Geophys. Res. 119, 2014JE004698-2271. 
/// https://doi.org/10.1002/2014JE004698
/// 
/// This function is split off from the `profile` function for clarity, and is not intended to be 
/// called directly.
///
/// # Arguments
///
/// * `r` - Normalized radial distance (unitless, where 1.0 corresponds to the crater rim).
/// * `elevation` - Baseline elevation before crater modification.
/// * `c0`, `c1`, `c2`, `c3` - Polynomial coefficients for the crater profile interior.
/// * `rim_height` - Height of the crater rim.
/// * `ejrim` - Rim dropoff parameter.
///
/// # Returns
///
/// * Adjusted elevation according to crater profile at distance `r`.
#[inline]
fn profile_function(r: f64, elevation: f64, c0: f64, c1: f64, c2: f64, c3: f64, rim_height: f64, ejrim: f64) -> f64 {
    if r >= 1.0 {
        elevation + (rim_height - ejrim) * r.powf(-RIMDROP)
    } else {
        elevation + c0 + c1 * r + c2 * r.powi(2) + c3 * r.powi(3)
    }
}

/// Computes a crater profile elevation array from input radial distances and reference elevations.
///
/// This function applies `profile_function` to each radial distance in the input array.
/// The coefficients c0, c1, c2, and c3 are calulated based on the crater dimensions and are based
/// on the polynomial crater profile model described in Fassett and Thomson (2014).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `r_array` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `diameter` - Total diameter of the crater (in meters).
/// * `floor_depth` - Depth of the crater floor below mean surface level (in meters).
/// * `floor_diameter` - Diameter of the crater floor (in meters).
/// * `rim_height` - Height of the crater rim above mean surface level (in meters).
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
#[pyfunction]
pub fn profile<'py>(
    py: Python<'py>,
    r_array: PyReadonlyArray1<'py, f64>,
    reference_elevation_array: PyReadonlyArray1<'py, f64>,
    diameter: f64,
    floor_depth: f64,
    floor_diameter: f64,
    rim_height: f64,
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
    let flrad = floor_diameter / diameter;
    let radius = diameter / 2.0;

    // Use polynomial crater profile similar to that of Fassett and Thomson (2014), but the parameters are set by the crater dimensions
    let c1 = (-floor_depth - rim_height)
        / (flrad - 1.0 + A * (flrad.powi(2) - 1.0) + B * (flrad.powi(3) - 1.0));
    let c0 = rim_height - c1 * (1.0 + A + B);
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
    let min_elevation = meanref - floor_depth;

    Ok(PyArray1::from_iter(
        py,
        reference_elevations
            .iter()
            .zip(radial_distances)
            .map(|(&elevation, &radial_distance)| {
                let r = radial_distance / radius;
                (
                    profile_function(r, elevation, c0, c1, c2, c3, rim_height, ejrim),
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
