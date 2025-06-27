#![feature(float_erf)] // Required to use f64::erf (https://github.com/rust-lang/rust/issues/136321)

pub mod simplemoon_functions;
pub mod surface_functions;
pub mod counting_functions;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use std::f64;

use pyo3::prelude::*;

const VSMALL: f64 = 10.0 * std::f64::EPSILON;
const RIMDROP: f64 = 4.20; // The exponent for the uplifted rim dropoff.
const EJPROFILE: f64 = 3.0; // The exponent for the ejecta profile

#[pymodule]
#[pyo3(name = "_cratermaker")]
mod cratermaker {
    use super::*;

    #[pymodule]
    mod simplemoon_functions {
        #[pymodule_export]
        use crate::simplemoon_functions::{crater_profile, ejecta_profile, ray_intensity};
    }

    #[pymodule]
    mod surface_functions {
        #[pymodule_export]
        use crate::surface_functions::{
            apply_diffusion, slope_collapse, calculate_bearing, interpolate_node_elevation_from_faces,
            turbulence_noise, calculate_distance, compute_slope, compute_edge_distances
        };
    }

    #[pymodule]
    mod counting_functions {
        #[pymodule_export]
        use crate::counting_functions::{
            tally_m19
        };
    }
}
