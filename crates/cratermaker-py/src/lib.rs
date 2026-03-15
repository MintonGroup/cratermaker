pub mod counting_bindings;
pub mod morphology_bindings;
pub mod surface_bindings;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "bindings")]
mod cratermaker {
    use super::*;

    #[pymodule]
    mod counting_bindings {
        #[pymodule_export]
        use crate::counting_bindings::{
            fit_one_ellipse, fit_one_ellipse_fixed_center, fit_rim, measure_floor_depth,
            measure_rim_height, score_rim,
        };
    }

    #[pymodule]
    mod morphology_bindings {
        #[pymodule_export]
        use crate::morphology_bindings::{crater_profile, ejecta_profile, ray_intensity};
    }

    #[pymodule]
    mod surface_bindings {
        #[pymodule_export]
        use crate::surface_bindings::{
            apply_diffusion, compute_bearings, compute_distances, compute_edge_distances,
            compute_location_from_distance_bearing, compute_radial_gradient, compute_slope,
            interpolate_node_elevation_from_faces, slope_collapse, turbulence_noise,
        };
    }
}
