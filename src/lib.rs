#![feature(float_erf)] // Required to use f64::erf (https://github.com/rust-lang/rust/issues/136321)

pub mod crater;
pub mod ejecta;
pub mod realistic;

use std::f64;

use pyo3::prelude::*;

const VSMALL: f64 = 10.0 * std::f64::EPSILON;

#[pymodule]
mod cratermaker {
    use super::*;

    #[pymodule]
    mod crater {
        #[pymodule_export]
        use crate::crater::profile;
    }

    #[pymodule]
    mod ejecta {
        #[pymodule_export]
        use crate::ejecta::{distribution, profile, ray_intensity};
    }

    #[pymodule]
    mod realistic {
        #[pymodule_export]
        use crate::realistic::{apply_noise, realistic_crater};
    }
}
