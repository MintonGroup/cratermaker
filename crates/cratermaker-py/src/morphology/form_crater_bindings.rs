//! Python-facing plumbing shared by the morphology models' `form_crater` bindings.
//!
//! Each model's binding converts its crater object into a [`CraterProfile`] and calls [`run_form_crater`], which reads the
//! region geometry from the Python `LocalSurface` once, gathers the region's values from the global surface arrays, runs the
//! formation pipeline, and writes the results back into the global arrays in place.

use cratermaker_components::morphology::form_crater::{
    CraterProfile, CraterRegionView, CraterSurfaceState, FormCraterOptions, form_crater,
};
use numpy::ndarray::Array1;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// The geometry of a Python `LocalSurface` that crater formation needs.
pub struct PyCraterRegion<'py> {
    face_area: PyReadonlyArray1<'py, f64>,
    face_lon: PyReadonlyArray1<'py, f64>,
    face_lat: PyReadonlyArray1<'py, f64>,
    face_distance: PyReadonlyArray1<'py, f64>,
    face_bearing: PyReadonlyArray1<'py, f64>,
    face_proj_x: PyReadonlyArray1<'py, f64>,
    face_proj_y: PyReadonlyArray1<'py, f64>,
    node_distance: PyReadonlyArray1<'py, f64>,
    node_bearing: PyReadonlyArray1<'py, f64>,
    node_proj_x: PyReadonlyArray1<'py, f64>,
    node_proj_y: PyReadonlyArray1<'py, f64>,
    edge_length: PyReadonlyArray1<'py, f64>,
    edge_face_distance: PyReadonlyArray1<'py, f64>,
    face_edge_connectivity: PyReadonlyArray2<'py, i64>,
    face_node_connectivity: PyReadonlyArray2<'py, i64>,
    face_face_connectivity: PyReadonlyArray2<'py, i64>,
    node_face_connectivity: PyReadonlyArray2<'py, i64>,
    edge_face_connectivity: PyReadonlyArray2<'py, i64>,
}

impl<'py> PyCraterRegion<'py> {
    /// Reads the attributes of a Python `LocalSurface`. Accessing them computes any lazily evaluated values (distances, bearings,
    /// projected coordinates, and local connectivity) that have not been computed yet.
    pub fn from_py(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        let get1 = |name: &str| -> PyResult<PyReadonlyArray1<'py, f64>> { Ok(obj.getattr(name)?.extract()?) };
        let get2 = |name: &str| -> PyResult<PyReadonlyArray2<'py, i64>> { Ok(obj.getattr(name)?.extract()?) };
        Ok(Self {
            face_area: get1("face_area")?,
            face_lon: get1("face_lon")?,
            face_lat: get1("face_lat")?,
            face_distance: get1("face_distance")?,
            face_bearing: get1("face_bearing")?,
            face_proj_x: get1("face_proj_x")?,
            face_proj_y: get1("face_proj_y")?,
            node_distance: get1("node_distance")?,
            node_bearing: get1("node_bearing")?,
            node_proj_x: get1("node_proj_x")?,
            node_proj_y: get1("node_proj_y")?,
            edge_length: get1("edge_length")?,
            edge_face_distance: get1("edge_face_distance")?,
            face_edge_connectivity: get2("face_edge_connectivity")?,
            face_node_connectivity: get2("face_node_connectivity")?,
            face_face_connectivity: get2("face_face_connectivity")?,
            node_face_connectivity: get2("node_face_connectivity")?,
            edge_face_connectivity: get2("edge_face_connectivity")?,
        })
    }

    pub fn as_view(&self) -> CraterRegionView<'_> {
        CraterRegionView {
            face_area: self.face_area.as_array(),
            face_lon: self.face_lon.as_array(),
            face_lat: self.face_lat.as_array(),
            face_distance: self.face_distance.as_array(),
            face_bearing: self.face_bearing.as_array(),
            face_proj_x: self.face_proj_x.as_array(),
            face_proj_y: self.face_proj_y.as_array(),
            node_distance: self.node_distance.as_array(),
            node_bearing: self.node_bearing.as_array(),
            node_proj_x: self.node_proj_x.as_array(),
            node_proj_y: self.node_proj_y.as_array(),
            edge_length: self.edge_length.as_array(),
            edge_face_distance: self.edge_face_distance.as_array(),
            face_edge_connectivity: self.face_edge_connectivity.as_array(),
            face_node_connectivity: self.face_node_connectivity.as_array(),
            face_face_connectivity: self.face_face_connectivity.as_array(),
            node_face_connectivity: self.node_face_connectivity.as_array(),
            edge_face_connectivity: self.edge_face_connectivity.as_array(),
        }
    }
}

/// Converts global indices to `usize`, checking that each is within `[0, n_total)`.
fn checked_indices(indices: &PyReadonlyArray1<'_, i64>, n_total: usize, what: &str) -> PyResult<Vec<usize>> {
    indices
        .as_array()
        .iter()
        .map(|&i| {
            if i >= 0 && (i as usize) < n_total {
                Ok(i as usize)
            } else {
                Err(PyValueError::new_err(format!("{what} index {i} is out of range [0, {n_total})")))
            }
        })
        .collect()
}

/// Forms one crater on a region of the surface, updating the global surface arrays in place.
///
/// # Arguments
///
/// * `model` - The morphology model's crater profile.
/// * `region` - The Python `LocalSurface` the crater is formed on.
/// * `face_indices`, `node_indices` - Global indices of the region's faces and nodes.
/// * `face_elevation`, `node_elevation`, `ejecta_thickness` - The global surface arrays, updated in place.
/// * `options` - Volume conservation, ray, and slope collapse options.
///
/// # Returns
///
/// `(excavated_volume, conservation_factor)`, where `conservation_factor` is `None` if no ejecta was formed.
#[allow(clippy::too_many_arguments)]
pub fn run_form_crater<'py, P: CraterProfile>(
    model: &P,
    region: &Bound<'py, PyAny>,
    face_indices: PyReadonlyArray1<'py, i64>,
    node_indices: PyReadonlyArray1<'py, i64>,
    face_elevation: &Bound<'py, PyArray1<f64>>,
    node_elevation: &Bound<'py, PyArray1<f64>>,
    ejecta_thickness: &Bound<'py, PyArray1<f64>>,
    options: FormCraterOptions,
) -> PyResult<(f64, Option<f64>)> {
    let region_py = PyCraterRegion::from_py(region)?;
    let region_v = region_py.as_view();

    let mut face_elevation = face_elevation.readwrite();
    let mut node_elevation = node_elevation.readwrite();
    let mut ejecta_thickness = ejecta_thickness.readwrite();
    let face_global = face_elevation.as_slice_mut()?;
    let node_global = node_elevation.as_slice_mut()?;
    let ejecta_global = ejecta_thickness.as_slice_mut()?;
    if ejecta_global.len() != face_global.len() {
        return Err(PyValueError::new_err("ejecta_thickness must have one value per face"));
    }

    let faces = checked_indices(&face_indices, face_global.len(), "face")?;
    let nodes = checked_indices(&node_indices, node_global.len(), "node")?;

    let mut state = CraterSurfaceState {
        face_elevation: Array1::from_iter(faces.iter().map(|&f| face_global[f])),
        node_elevation: Array1::from_iter(nodes.iter().map(|&n| node_global[n])),
        ejecta_thickness: Array1::from_iter(faces.iter().map(|&f| ejecta_global[f])),
    };

    let outcome = form_crater(model, &region_v, &mut state, &options).map_err(PyValueError::new_err)?;

    for (i, &f) in faces.iter().enumerate() {
        face_global[f] = state.face_elevation[i];
        ejecta_global[f] = state.ejecta_thickness[i];
    }
    for (i, &n) in nodes.iter().enumerate() {
        node_global[n] = state.node_elevation[i];
    }
    Ok((outcome.excavated_volume, outcome.conservation_factor))
}
