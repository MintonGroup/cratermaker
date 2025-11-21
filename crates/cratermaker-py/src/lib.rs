pub mod morphology_bindings;
pub mod surface_bindings;
pub mod counting_bindings;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use pyo3::prelude::*;


#[pymodule]
#[pyo3(name = "_cratermaker")]
mod cratermaker {
    use super::*;

    #[pymodule]
    mod counting_bindings {
        #[pymodule_export]
        use crate::counting_bindings::{
            tally, 
            radial_distance_to_ellipse, 
            fit_one_ellipse,
            score_rim
        };
    }

    #[pymodule]
    mod morphology_bindings {
        #[pymodule_export]
        use crate::morphology_bindings::{
            crater_profile, 
            ejecta_profile, 
            ray_intensity
        };
    }

    #[pymodule]
    mod surface_bindings {
        #[pymodule_export]
        use crate::surface_bindings::{
            apply_diffusion, 
            slope_collapse, 
            compute_bearings,
            compute_distances, 
            interpolate_node_elevation_from_faces,
            turbulence_noise, 
            compute_slope, 
            compute_edge_distances, 
            compute_radial_gradient,
        };
    }


}
