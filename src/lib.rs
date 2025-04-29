#![feature(float_erf)] // Required to use f64::erf (https://github.com/rust-lang/rust/issues/136321)

pub mod crater_functions;
pub mod ejecta_functions;
pub mod morphology_functions;
pub mod surface_functions;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use std::f64;

use pyo3::prelude::*;

const VSMALL: f64 = 10.0 * std::f64::EPSILON;

#[pymodule]
#[pyo3(name = "_cratermaker")]
mod cratermaker {
    use super::*;

    #[pymodule]
    mod crater_functions {
        #[pymodule_export]
        use crate::crater_functions::profile;
    }

    #[pymodule]
    mod ejecta_functions {
        #[pymodule_export]
        use crate::ejecta_functions::{distribution, profile, ray_intensity};
    }

    #[pymodule]
    mod surface_functions {
        #[pymodule_export]
        use crate::surface_functions::calculate_initial_bearing;
    }

    #[pymodule]
    mod morphology_functions {
        #[pymodule_export]
        use crate::morphology_functions::form_ejecta;
    }
}
