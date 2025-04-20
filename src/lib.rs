pub mod crater;
pub mod ejecta;
pub mod realistic;

use pyo3::prelude::*;

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
        use crate::ejecta::{distribution, profile};
    }

    #[pymodule]
    mod realistic {
        #[pymodule_export]
        use crate::realistic::{apply_noise, realistic_crater};
    }
}
