use std::f64::consts::PI;

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use ndarray::Zip;
use rayon::iter::{IntoParallelIterator,ParallelIterator};

#[inline]
fn positive_mod(x: f64, m: f64) -> f64 {
    ((x % m) + m) % m
}

/// Computes the initial bearing (forward azimuth) from a fixed point to each of a set of destination points.
///
/// The bearing is calculated on a spherical surface using great-circle paths and returned in radians,
/// normalized to the range [0, 2π).
///
/// # Arguments
///
/// * `py` - Python GIL token.
/// * `lon1` - Longitude of the reference point, in radians.
/// * `lat1` - Latitude of the reference point, in radians.
/// * `lon2` - Longitudes of destination points, in radians.
/// * `lat2` - Latitudes of destination points, in radians.
///
/// # Returns
///
/// * A NumPy array of initial bearing angles (radians), one for each (lon2, lat2) pair.
#[pyfunction]
pub fn calculate_initial_bearing<'py>(
    py: Python<'py>,
    lon1: f64,
    lat1: f64,
    lon2: PyReadonlyArray1<'py, f64>,
    lat2: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let lon2 = lon2.as_array();
    let lat2 = lat2.as_array();
    let mut result = ndarray::Array1::<f64>::zeros(lon2.len());

    Zip::from(&mut result)
        .and(lon2)
        .and(lat2)
        .into_par_iter()
        .for_each(|(out, &lon2, &lat2)| {
            // Calculate differences in coordinates
            let dlon = positive_mod(lon2 - lon1 + PI, 2.0 * PI) - PI;

            // Haversine formula calculations
            let x = dlon.sin() * lat2.cos();
            let y = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
            let initial_bearing = f64::atan2(x, y);

            // Normalize bearing to 0 to 2*pi
            *out = (initial_bearing + 2.0 * PI) % (2.0 * PI);
        });

    Ok(PyArray1::from_owned_array(py, result))
}


/// View into a region of the surface mesh, consisting of node and face indices.
///
/// Used to localize crater effects to a subset of the full mesh.
#[derive(FromPyObject)]
pub struct SurfaceView<'py> {
    pub node_indices: PyReadonlyArray1<'py, i64>,
    pub face_indices: PyReadonlyArray1<'py, i64>,
}

/// Applies one explicit diffusion update step over a surface mesh with variable diffusivity.
///
/// This function computes the change in elevation for each face on the mesh using a
/// face-centered finite-volume formulation of the operator:
///     ∂h/∂t = ∇ · (κ ∇h)
/// where κ varies per face. For each face, the flux with its neighbors is computed
/// using the expression (κ_f + κ_n) * (h_n - h_f), summing over all neighbors n.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `face_areas` - Area of each face (1D array of length n_faces).
/// * `face_kappa` - Topographic diffusivity at each face (1D array of length n_faces).
/// * `face_elevation` - Elevation value at each face (1D array of length n_faces).
/// * `face_face_connectivity` - For each face, the indices of its neighboring faces (2D array).
/// * `dt` - Time step size.
///
/// # Returns
///
/// A NumPy array of shape (n_faces,) giving the elevation change per face for this update step.
#[pyfunction]
pub fn apply_diffusion_update<'py>(
    py: Python<'py>,
    face_areas: PyReadonlyArray1<'py, f64>,
    face_kappa: PyReadonlyArray1<'py, f64>,
    face_elevation: PyReadonlyArray1<'py, f64>,
    face_face_connectivity: PyReadonlyArray2<'py, i64>,
    dt: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let face_areas = face_areas.as_array();
    let face_kappa = face_kappa.as_array();
    let face_elevation = face_elevation.as_array();
    let face_face_connectivity = face_face_connectivity.as_array();

    let n_faces = face_areas.len();
    let mut dhdt = ndarray::Array1::<f64>::zeros(n_faces);

    for f in 0..n_faces {
        let k_f = face_kappa[f];
        let h_f = face_elevation[f];

        for neighbor in face_face_connectivity.row(f) {
            let f_n = *neighbor as usize;

            if f_n >= n_faces {
                continue;
            }

            let k_n = face_kappa[f_n];
            let h_n = face_elevation[f_n];

            let flux = (k_f + k_n) * (h_n - h_f);

            dhdt[f] += flux;
        }
    }

    for f in 0..n_faces {
        dhdt[f] *= dt / face_areas[f];
    }

    Ok(PyArray1::from_owned_array(py, dhdt))
}