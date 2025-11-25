use pyo3::prelude::*;

#[derive(FromPyObject)]
pub struct Crater {
    pub id:                 u32,
    pub semimajor_axis:     f64,
    pub semiminor_axis:     f64,
    pub orientation:        f64,
    pub transient_diameter: f64,
    pub projectile_diameter: f64,
    pub projectile_velocity: f64,
    pub projectile_angle:   f64,
    pub projectile_density: f64,
    pub location:           (f64, f64),
    pub morphology_type:    String,
    pub age:                Option<f64>,
}

