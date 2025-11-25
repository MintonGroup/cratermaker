use numpy::ndarray::prelude::*;

/// Represents a local region of a surface mesh with various attributes accessible as array views.
pub struct LocalSurfaceView<'a> {
    pub n_face:         usize,
    pub pix:            f64,
    pub face_area:      ArrayView1<'a, f64>,
    pub face_elevation: ArrayView1<'a, f64>,
    pub face_indices:   ArrayView1<'a, i64>,
    pub face_lon:       ArrayView1<'a, f64>,
    pub face_lat:       ArrayView1<'a, f64>,
    pub face_x:         ArrayView1<'a, f64>,
    pub face_y:         ArrayView1<'a, f64>,
    pub face_z:         ArrayView1<'a, f64>,

    pub n_node:         usize,
    pub node_elevation: ArrayView1<'a, f64>,
    pub node_indices:   ArrayView1<'a, i64>,
    pub node_lon:       ArrayView1<'a, f64>,
    pub node_lat:       ArrayView1<'a, f64>,
    pub node_x:         ArrayView1<'a, f64>,
    pub node_y:         ArrayView1<'a, f64>,
    pub node_z:         ArrayView1<'a, f64>,

    pub n_edge:                 usize,
    pub edge_indices:           ArrayView1<'a, i64>,
    pub edge_length:            ArrayView1<'a, f64>,

    pub face_edge_connectivity: ArrayView2<'a, i64>,
    pub face_node_connectivity: ArrayView2<'a, i64>,
    pub face_face_connectivity: ArrayView2<'a, i64>,
    pub node_face_connectivity: ArrayView2<'a, i64>,
    pub edge_face_connectivity: ArrayView2<'a, i64>,
    pub edge_node_connectivity: ArrayView2<'a, i64>,
    pub edge_face_distance:     ArrayView1<'a, f64>,

    // The following are optional, as they will not be present if the LocalSurfaceView represents a global Surface object
    pub face_proj_x:    Option<ArrayView1<'a, f64>>,
    pub face_proj_y:    Option<ArrayView1<'a, f64>>,
    pub face_distance:  Option<ArrayView1<'a, f64>>,
    pub face_bearing:   Option<ArrayView1<'a, f64>>,
}