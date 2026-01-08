use pyo3::prelude::*;

#[derive(FromPyObject, Clone, Debug)]
pub struct Crater {
    pub id:                 u32,
    pub diameter:           f64,
    pub radius:             f64,
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
    pub measured_semimajor_axis: f64,
    pub measured_semiminor_axis: f64,
    pub measured_orientation: f64,
    pub measured_diameter: f64,
    pub measured_radius:   f64,
    pub measured_location: (f64, f64),
    pub measured_rim_height:  Option<f64>,
    pub measured_floor_depth: Option<f64>,
    pub measured_degradation_state: Option<f64>,
    pub age:                Option<f64>,
}

