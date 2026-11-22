//! Python-facing plumbing shared by the morphology models' `form_craters` bindings.
//!
//! Each model's binding converts its crater objects into [`CraterProfile`]s and calls [`run_form_craters`], which reads each
//! crater's region geometry from its Python `LocalSurface`, gathers the region's values from the global surface arrays, forms
//! the craters in parallel with the GIL released, and writes the results back into the global arrays in place.

use cratermaker_components::morphology::form_crater::{
    CraterProfile, CraterRegionView, CraterSurfaceState, FormCraterOptions, form_crater,
};
use numpy::ndarray::Array1;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::cell::RefCell;

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
        let get1 = |name: &str| -> PyResult<PyReadonlyArray1<'py, f64>> {
            Ok(obj.getattr(name)?.extract()?)
        };
        let get2 = |name: &str| -> PyResult<PyReadonlyArray2<'py, i64>> {
            Ok(obj.getattr(name)?.extract()?)
        };
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
fn checked_indices(
    indices: &PyReadonlyArray1<'_, i64>,
    n_total: usize,
    what: &str,
) -> PyResult<Vec<usize>> {
    indices
        .as_array()
        .iter()
        .map(|&i| {
            if i >= 0 && (i as usize) < n_total {
                Ok(i as usize)
            } else {
                Err(PyValueError::new_err(format!(
                    "{what} index {i} is out of range [0, {n_total})"
                )))
            }
        })
        .collect()
}

thread_local! {
    /// Reusable marks over global face or node indices for the overlap check. Always left all-false between uses.
    static OVERLAP_MARKS: RefCell<Vec<bool>> = const { RefCell::new(Vec::new()) };
}

/// Checks that no index appears in more than one of the lists.
fn check_disjoint(lists: &[Vec<usize>], n_total: usize, what: &str) -> PyResult<()> {
    if lists.len() < 2 {
        return Ok(());
    }
    OVERLAP_MARKS.with(|cell| {
        let mut marks = cell.borrow_mut();
        if marks.len() < n_total {
            marks.resize(n_total, false);
        }
        let mut overlap = None;
        'outer: for (k, list) in lists.iter().enumerate() {
            for &i in list {
                if marks[i] {
                    overlap = Some((k, i));
                    break 'outer;
                }
                marks[i] = true;
            }
        }
        for list in lists {
            for &i in list {
                marks[i] = false;
            }
        }
        match overlap {
            Some((k, i)) => Err(PyValueError::new_err(format!(
                "the regions of the craters in a batch must not overlap, but {what} {i} of crater {k} belongs to an earlier crater"
            ))),
            None => Ok(()),
        }
    })
}

/// Forms a batch of craters with non-overlapping regions, updating the global surface arrays in place.
///
/// The craters are formed in parallel with the GIL released. Each crater only reads and writes the faces and nodes of its own
/// region, and the regions are checked to be disjoint, so the result is the same as forming them one after another.
///
/// # Arguments
///
/// * `py` - Python interpreter token.
/// * `models` - The morphology model's crater profile for each crater.
/// * `regions` - The Python `LocalSurface` each crater is formed on.
/// * `face_indices`, `node_indices` - Global indices of each region's faces and nodes.
/// * `ray_seeds` - The ray pattern seed of each crater, or `None` for ejecta without rays.
/// * `face_elevation`, `node_elevation`, `ejecta_thickness` - The global surface arrays, updated in place.
/// * `conserve_volume` - Scale each crater's ejecta so that its volume matches the excavated volume.
/// * `critical_slope` - Critical slope (rise over run) for slope collapse, or `None` to skip slope collapse.
///
/// # Returns
///
/// `(excavated_volume, conservation_factor)` for each crater, where `conservation_factor` is `None` if no ejecta was formed.
#[allow(clippy::too_many_arguments)]
pub fn run_form_craters<'py, P: CraterProfile>(
    py: Python<'py>,
    models: &[P],
    regions: &[Bound<'py, PyAny>],
    face_indices: &[PyReadonlyArray1<'py, i64>],
    node_indices: &[PyReadonlyArray1<'py, i64>],
    ray_seeds: &[Option<u64>],
    face_elevation: &Bound<'py, PyArray1<f64>>,
    node_elevation: &Bound<'py, PyArray1<f64>>,
    ejecta_thickness: &Bound<'py, PyArray1<f64>>,
    conserve_volume: bool,
    critical_slope: Option<f64>,
) -> PyResult<Vec<(f64, Option<f64>)>> {
    let n = models.len();
    if regions.len() != n
        || face_indices.len() != n
        || node_indices.len() != n
        || ray_seeds.len() != n
    {
        return Err(PyValueError::new_err(
            "craters, regions, face_indices, node_indices, and ray_seeds must all have the same length",
        ));
    }

    let regions_py: Vec<PyCraterRegion> = regions
        .iter()
        .map(PyCraterRegion::from_py)
        .collect::<PyResult<_>>()?;
    let views: Vec<CraterRegionView> = regions_py.iter().map(|r| r.as_view()).collect();

    let mut face_elevation = face_elevation.readwrite();
    let mut node_elevation = node_elevation.readwrite();
    let mut ejecta_thickness = ejecta_thickness.readwrite();
    let face_global = face_elevation.as_slice_mut()?;
    let node_global = node_elevation.as_slice_mut()?;
    let ejecta_global = ejecta_thickness.as_slice_mut()?;
    if ejecta_global.len() != face_global.len() {
        return Err(PyValueError::new_err(
            "ejecta_thickness must have one value per face",
        ));
    }

    let faces: Vec<Vec<usize>> = face_indices
        .iter()
        .map(|idx| checked_indices(idx, face_global.len(), "face"))
        .collect::<PyResult<_>>()?;
    let nodes: Vec<Vec<usize>> = node_indices
        .iter()
        .map(|idx| checked_indices(idx, node_global.len(), "node"))
        .collect::<PyResult<_>>()?;
    check_disjoint(&faces, face_global.len(), "face")?;
    check_disjoint(&nodes, node_global.len(), "node")?;

    let mut states: Vec<CraterSurfaceState> = faces
        .iter()
        .zip(nodes.iter())
        .map(|(f, nd)| CraterSurfaceState {
            face_elevation: Array1::from_iter(f.iter().map(|&i| face_global[i])),
            node_elevation: Array1::from_iter(nd.iter().map(|&i| node_global[i])),
            ejecta_thickness: Array1::from_iter(f.iter().map(|&i| ejecta_global[i])),
        })
        .collect();
    let options: Vec<FormCraterOptions> = ray_seeds
        .iter()
        .map(|&ray_seed| FormCraterOptions {
            conserve_volume,
            ray_seed,
            critical_slope,
        })
        .collect();

    let outcomes = py
        .detach(|| {
            states
                .par_iter_mut()
                .enumerate()
                .map(|(k, state)| form_crater(&models[k], &views[k], state, &options[k]))
                .collect::<Result<Vec<_>, String>>()
        })
        .map_err(PyValueError::new_err)?;

    for ((state, f), nd) in states.iter().zip(faces.iter()).zip(nodes.iter()) {
        for (i, &g) in f.iter().enumerate() {
            face_global[g] = state.face_elevation[i];
            ejecta_global[g] = state.ejecta_thickness[i];
        }
        for (i, &g) in nd.iter().enumerate() {
            node_global[g] = state.node_elevation[i];
        }
    }
    Ok(outcomes
        .iter()
        .map(|o| (o.excavated_volume, o.conservation_factor))
        .collect())
}
