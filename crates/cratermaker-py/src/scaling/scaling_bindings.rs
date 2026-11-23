use cratermaker_components::scaling::TransientScaling;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Computes the transient crater diameter made by a projectile, using the Holsapple (1993) scaling law.
///
/// # Arguments
///
/// * `projectile_diameter` - The projectile diameter (m).
/// * `k1`, `mu`, `ybar` - The material scaling coefficient, exponent, and strength (Pa).
/// * `target_density`, `projectile_density` - Densities (kg/m³).
/// * `gravity` - Target surface gravity (m/s²).
/// * `vertical_velocity` - Vertical component of the impact velocity (m/s).
///
/// # Returns
///
/// * The transient crater diameter (m).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn projectile_to_transient(
    projectile_diameter: f64,
    k1: f64,
    mu: f64,
    ybar: f64,
    target_density: f64,
    projectile_density: f64,
    gravity: f64,
    vertical_velocity: f64,
) -> f64 {
    let s = TransientScaling {
        k1,
        mu,
        ybar,
        target_density,
        projectile_density,
        gravity,
        vertical_velocity,
    };
    cratermaker_components::scaling::projectile_to_transient(projectile_diameter, &s)
}

/// Computes the diameter of the projectile that makes a transient crater of the given diameter.
///
/// # Arguments
///
/// * `transient_diameter` - The transient crater diameter (m).
/// * `k1`, `mu`, `ybar` - The material scaling coefficient, exponent, and strength (Pa).
/// * `target_density`, `projectile_density` - Densities (kg/m³).
/// * `gravity` - Target surface gravity (m/s²).
/// * `vertical_velocity` - Vertical component of the impact velocity (m/s).
///
/// # Returns
///
/// * The projectile diameter (m).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn transient_to_projectile(
    transient_diameter: f64,
    k1: f64,
    mu: f64,
    ybar: f64,
    target_density: f64,
    projectile_density: f64,
    gravity: f64,
    vertical_velocity: f64,
) -> PyResult<f64> {
    let s = TransientScaling {
        k1,
        mu,
        ybar,
        target_density,
        projectile_density,
        gravity,
        vertical_velocity,
    };
    cratermaker_components::scaling::transient_to_projectile(transient_diameter, &s)
        .map_err(PyValueError::new_err)
}
