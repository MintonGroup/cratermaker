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
            fit_one_ellipse, fit_one_ellipse_fixed_center, fit_one_circle, fit_one_circle_fixed_center, fit_rim, measure_floor_elevation,
            measure_rim_elevation, score_rim,
        };
    }

    #[pymodule]
    mod morphology_bindings {
        #[pymodule_export]
        use crate::morphology_bindings::{basicmoon_profile, ray_intensity, realmoon_profile, calculate_target_1D_PSD_from_breakpoint_slope};
    }

    #[pymodule]
    mod surface_bindings {
        #[pymodule_export]
        use crate::surface_bindings::{
            apply_diffusion, compute_bearings, compute_distances, compute_edge_distances,
            compute_location_from_distance_bearing, compute_radial_gradient, compute_slope,
            interpolate_node_elevation_from_faces, reset_radial_distances, slope_collapse,
            turbulence_noise,
        };
    }
}
