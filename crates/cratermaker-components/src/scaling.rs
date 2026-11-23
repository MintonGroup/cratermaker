//! Crater scaling relationships between projectiles and transient craters.

use crate::rootfind::{BRENTQ_MAXITER, BRENTQ_RTOL, BRENTQ_XTOL, BrentqResult, brentq};
use std::f64::consts::PI;

/// The target, projectile, and material properties that set the size of a transient crater (Holsapple 1993).
#[derive(Clone, Copy, Debug)]
pub struct TransientScaling {
    /// Material scaling coefficient K1.
    pub k1: f64,
    /// Material scaling exponent mu.
    pub mu: f64,
    /// Target strength (Pa).
    pub ybar: f64,
    /// Target density (kg/m³).
    pub target_density: f64,
    /// Projectile density (kg/m³).
    pub projectile_density: f64,
    /// Target surface gravity (m/s²).
    pub gravity: f64,
    /// Vertical component of the impact velocity (m/s).
    pub vertical_velocity: f64,
}

/// Computes the transient crater diameter made by a projectile, using the Holsapple (1993) scaling law that combines the
/// gravity and strength regimes. The transient crater is never smaller than the projectile.
///
/// # Arguments
///
/// * `projectile_diameter` - The projectile diameter (m).
/// * `s` - The scaling properties.
///
/// # Returns
///
/// The transient crater diameter (m).
pub fn projectile_to_transient(projectile_diameter: f64, s: &TransientScaling) -> f64 {
    let projectile_radius = projectile_diameter / 2.0;
    let projectile_mass = (4.0 / 3.0) * PI * projectile_radius.powf(3.0) * s.projectile_density;
    let c1 = 1.0 + 0.5 * s.mu;
    let c2 = (-3.0 * s.mu) / (2.0 + s.mu);

    // Dimensionless quantities
    let vv2 = s.vertical_velocity.powf(2.0);
    let pitwo = (s.gravity * projectile_radius) / vv2;
    let pithree = s.ybar / (s.target_density * vv2);
    let pifour = s.target_density / s.projectile_density;
    let pivol = s.k1 * ((pitwo * pifour.powf(-1.0 / 3.0)) + pithree.powf(c1)).powf(c2);

    // Transient crater volume and diameter (depth = 1/3 diameter)
    let cvol = pivol * (projectile_mass / s.target_density);
    let transient_diameter = 2.0 * (3.0 * cvol / PI).powf(1.0 / 3.0);
    transient_diameter.max(projectile_diameter)
}

/// Computes the diameter of the projectile that makes a transient crater of the given diameter, by inverting
/// [`projectile_to_transient`] with Brent's method over projectile diameters from 1e-5 to 1.2 times the crater diameter.
///
/// # Arguments
///
/// * `transient_diameter` - The transient crater diameter (m).
/// * `s` - The scaling properties.
///
/// # Returns
///
/// The projectile diameter (m), or an error if the crater cannot be made by a projectile in the search range.
pub fn transient_to_projectile(
    transient_diameter: f64,
    s: &TransientScaling,
) -> Result<f64, String> {
    let f = |d: f64| Ok(projectile_to_transient(d, s) - transient_diameter);
    match brentq(
        f,
        1e-5 * transient_diameter,
        1.2 * transient_diameter,
        BRENTQ_XTOL,
        BRENTQ_RTOL,
        BRENTQ_MAXITER,
    )? {
        BrentqResult::Converged(d) | BrentqResult::NotConverged(d) => Ok(d),
        BrentqResult::SignError => Err(format!(
            "no projectile between {:e} and {:e} m makes a transient crater of diameter {transient_diameter:e} m",
            1e-5 * transient_diameter,
            1.2 * transient_diameter
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn moon() -> TransientScaling {
        TransientScaling {
            k1: 0.20,
            mu: 0.55,
            ybar: 7.6e6,
            target_density: 2250.0,
            projectile_density: 2500.0,
            gravity: 1.62,
            vertical_velocity: 15e3,
        }
    }

    #[test]
    fn inverse_recovers_projectile() {
        let s = moon();
        for &pd in &[1.0, 30.0, 1e3, 3e4] {
            let td = projectile_to_transient(pd, &s);
            let back = transient_to_projectile(td, &s).unwrap();
            assert!(
                (back - pd).abs() <= 1e-9 * pd,
                "pd={pd} td={td} back={back}"
            );
        }
    }
}
