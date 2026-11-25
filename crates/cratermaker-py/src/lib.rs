pub mod counting;
pub mod morphology;
pub mod scaling;
pub mod surface;

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
        use crate::counting::counting_bindings::{measure_floor_elevation, measure_rim_height};
    }

    #[pymodule]
    mod basicmoon_bindings {
        #[pymodule_export]
        use crate::morphology::basicmoon_bindings::{
            basicmoon_profile, crater_profile_function, degradation_function,
            degradation_function_array, ejecta_profile_function, form_craters, ray_intensity, rmax,
            volume_integrand,
        };
    }

    #[pymodule]
    mod realmoon_bindings {
        #[pymodule_export]
        use crate::morphology::realmoon_bindings::{form_craters, realmoon_profile};
    }

    #[pymodule]
    mod scaling_bindings {
        #[pymodule_export]
        use crate::scaling::scaling_bindings::{projectile_to_transient, transient_to_projectile};
    }

    #[pymodule]
    mod surface_bindings {
        #[pymodule_export]
        use crate::surface::surface_bindings::{
            apply_diffusion, build_local_connectivity, build_region_indices, compute_bearings,
            compute_distances, compute_distances_and_bearings, compute_edge_distances,
            compute_location_from_distance_bearing, compute_radial_gradient, compute_slope,
            interpolate_node_elevation_from_faces, reference_surface, remap_connectivity_to_local,
            reset_radial_distances, slope_collapse, turbulence_noise,
        };
    }
}
