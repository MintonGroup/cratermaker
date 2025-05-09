#![feature(float_erf)] // Required to use f64::erf (https://github.com/rust-lang/rust/issues/136321)

pub mod crater_functions;
pub mod ejecta_functions;
pub mod surface_functions;

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
#[pyo3(name = "_simplemoon")]
mod simplemoon {
    use super::*;

    #[pymodule]
    mod crater_functions {
        #[pymodule_export]
        use crate::crater_functions::profile;
    }

    #[pymodule]
    mod ejecta_functions {
        #[pymodule_export]
        use crate::ejecta_functions::{profile, ray_intensity};
    }

    #[pymodule]
    mod surface_functions {
        #[pymodule_export]
        use crate::surface_functions::calculate_initial_bearing;
    }

}
