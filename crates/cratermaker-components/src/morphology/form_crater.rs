//! The crater formation pipeline shared by all morphology models.
//!
//! A morphology model supplies its crater and ejecta profiles through the [`CraterProfile`] trait, and [`form_crater`] runs the
//! full sequence of surface changes for one crater on a local region of the mesh:
//!
//! 1. Fit a reference plane to the pre-impact surface and compute the crater shape relative to it.
//! 2. Compute the excavated volume.
//! 3. Clear any existing ejecta from the crater interior.
//! 4. Compute the ejecta blanket (optionally modulated by rays), scale it to conserve volume, and record its thickness.
//! 5. Degrade the pre-existing surface by ejecta burial and by the crater's degradation function (two diffusion steps).
//! 6. Add the crater and ejecta to the surface.
//! 7. Optionally collapse slopes steeper than the critical slope.
//!
//! Each step follows the original Python implementation (`Morphology.form_crater` together with
//! `BasicMoonMorphology.compute_ejecta` and the `LocalSurface` methods it calls), including the NaN handling of
//! `LocalSurface.add_data`.

use crate::ArrayResult;
use crate::morphology::basicmoon::{
    BasicMoonCrater, basicmoon_profile, degradation_function, ray_intensity,
};
use crate::morphology::realmoon::{RealMoonCrater, realmoon_profile};
use crate::surface::{
    LocalSurfaceView, apply_diffusion, interpolate_node_elevation_from_faces, reference_surface,
    slope_collapse,
};
use numpy::ndarray::concatenate;
use numpy::ndarray::prelude::*;

/// Python's `_VSMALL` (10 times the smallest positive normal f64). Diffusion steps whose largest diffusivity is below this are
/// skipped.
const PY_VSMALL: f64 = 10.0 * f64::MIN_POSITIVE;

/// Scaling factor of the ejecta burial degradation, `k = EJECTA_SOFTEN_FACTOR * thickness^2`.
const EJECTA_SOFTEN_FACTOR: f64 = 1.50;

/// Size factor passed to the degradation function for the degradation caused by each crater's ejecta.
const DEGRADATION_FE: f64 = 100.0;

/// The model-specific parts of crater formation.
///
/// A morphology model implements `profile`, and may override the ray, degradation, and burial functions, which default to the
/// basic lunar model.
pub trait CraterProfile: Sync {
    /// Final crater diameter (m).
    fn diameter(&self) -> f64;

    /// Final crater radius (m).
    fn radius(&self) -> f64;

    /// Elevations of the crater and/or ejecta profile at each point.
    ///
    /// # Arguments
    ///
    /// * `radial_distances` - Distances from the crater center (m).
    /// * `bearings` - Bearings from the crater center (radians, clockwise from north).
    /// * `reference_elevations` - The surface the profile is drawn relative to (m).
    /// * `include_crater` - Whether to include the crater profile.
    /// * `include_ejecta` - Whether to include the ejecta profile.
    fn profile(
        &self,
        radial_distances: ArrayView1<'_, f64>,
        bearings: ArrayView1<'_, f64>,
        reference_elevations: ArrayView1<'_, f64>,
        include_crater: bool,
        include_ejecta: bool,
    ) -> ArrayResult;

    /// Ray pattern intensity at each point (bearings in radians).
    fn ray_intensity(
        &self,
        radial_distances: ArrayView1<'_, f64>,
        bearings: ArrayView1<'_, f64>,
        seed: u64,
    ) -> ArrayResult {
        ray_intensity(radial_distances, bearings, self.diameter(), seed)
    }

    /// Degradation state contributed by the crater's ejecta, before it is multiplied by the ejecta intensity.
    fn degradation(&self) -> f64 {
        degradation_function(self.diameter(), DEGRADATION_FE)
    }

    /// Degradation state caused by burial under an ejecta layer of the given thickness.
    fn burial_degradation(&self, ejecta_thickness: f64) -> f64 {
        EJECTA_SOFTEN_FACTOR * (ejecta_thickness * ejecta_thickness)
    }
}

/// The basic lunar crater model (Minton et al. 2026).
pub struct BasicMoonModel<'c> {
    pub crater: &'c BasicMoonCrater,
    pub rings: &'c Option<Vec<BasicMoonCrater>>,
}

impl CraterProfile for BasicMoonModel<'_> {
    fn diameter(&self) -> f64 {
        self.crater.diameter
    }

    fn radius(&self) -> f64 {
        self.crater.radius
    }

    fn profile(
        &self,
        radial_distances: ArrayView1<'_, f64>,
        _bearings: ArrayView1<'_, f64>,
        reference_elevations: ArrayView1<'_, f64>,
        include_crater: bool,
        include_ejecta: bool,
    ) -> ArrayResult {
        basicmoon_profile(
            radial_distances,
            reference_elevations,
            self.crater,
            self.rings,
            include_crater,
            include_ejecta,
        )
    }
}

/// The lunar crater model with bearing-dependent shape parameters (Du et al. 2025).
pub struct RealMoonModel<'c, 'a> {
    pub crater: &'c RealMoonCrater<'a>,
    pub rings: &'c Option<Vec<RealMoonCrater<'a>>>,
}

impl CraterProfile for RealMoonModel<'_, '_> {
    fn diameter(&self) -> f64 {
        self.crater.diameter
    }

    fn radius(&self) -> f64 {
        self.crater.radius
    }

    fn profile(
        &self,
        radial_distances: ArrayView1<'_, f64>,
        bearings: ArrayView1<'_, f64>,
        reference_elevations: ArrayView1<'_, f64>,
        include_crater: bool,
        include_ejecta: bool,
    ) -> ArrayResult {
        realmoon_profile(
            radial_distances,
            bearings,
            reference_elevations,
            self.crater,
            self.rings,
            include_crater,
            include_ejecta,
        )
    }
}

/// The geometry of the local region a crater is formed on. Connectivity arrays use local indices, with negative values for
/// entries that fall outside the region.
pub struct CraterRegionView<'a> {
    pub face_area: ArrayView1<'a, f64>,
    pub face_lon: ArrayView1<'a, f64>,
    pub face_lat: ArrayView1<'a, f64>,
    /// Distance from the crater center (m)
    pub face_distance: ArrayView1<'a, f64>,
    /// Bearing from the crater center (degrees)
    pub face_bearing: ArrayView1<'a, f64>,
    pub face_proj_x: ArrayView1<'a, f64>,
    pub face_proj_y: ArrayView1<'a, f64>,
    pub node_distance: ArrayView1<'a, f64>,
    pub node_bearing: ArrayView1<'a, f64>,
    pub node_proj_x: ArrayView1<'a, f64>,
    pub node_proj_y: ArrayView1<'a, f64>,
    pub edge_length: ArrayView1<'a, f64>,
    pub edge_face_distance: ArrayView1<'a, f64>,
    pub face_edge_connectivity: ArrayView2<'a, i64>,
    pub face_node_connectivity: ArrayView2<'a, i64>,
    pub face_face_connectivity: ArrayView2<'a, i64>,
    pub node_face_connectivity: ArrayView2<'a, i64>,
    pub edge_face_connectivity: ArrayView2<'a, i64>,
}

impl CraterRegionView<'_> {
    pub fn n_face(&self) -> usize {
        self.face_area.len()
    }

    pub fn n_node(&self) -> usize {
        self.node_distance.len()
    }

    /// Checks that every array has the length implied by the number of faces, nodes, and edges.
    fn validate(&self) -> Result<(), String> {
        let (nf, nn, ne) = (self.n_face(), self.n_node(), self.edge_length.len());
        let checks = [
            ("face_lon", self.face_lon.len(), nf),
            ("face_lat", self.face_lat.len(), nf),
            ("face_distance", self.face_distance.len(), nf),
            ("face_bearing", self.face_bearing.len(), nf),
            ("face_proj_x", self.face_proj_x.len(), nf),
            ("face_proj_y", self.face_proj_y.len(), nf),
            ("node_bearing", self.node_bearing.len(), nn),
            ("node_proj_x", self.node_proj_x.len(), nn),
            ("node_proj_y", self.node_proj_y.len(), nn),
            ("edge_face_distance", self.edge_face_distance.len(), ne),
            (
                "face_edge_connectivity rows",
                self.face_edge_connectivity.nrows(),
                nf,
            ),
            (
                "face_node_connectivity rows",
                self.face_node_connectivity.nrows(),
                nf,
            ),
            (
                "face_face_connectivity rows",
                self.face_face_connectivity.nrows(),
                nf,
            ),
            (
                "node_face_connectivity rows",
                self.node_face_connectivity.nrows(),
                nn,
            ),
            (
                "edge_face_connectivity rows",
                self.edge_face_connectivity.nrows(),
                ne,
            ),
        ];
        for (name, len, expected) in checks {
            if len != expected {
                return Err(format!("{name} has length {len}, expected {expected}"));
            }
        }
        Ok(())
    }

    /// A view of the region for the surface kernels (diffusion, slope collapse, node interpolation), with the given current
    /// elevations. Fields those kernels do not use are left empty.
    fn kernel_view<'v>(
        &'v self,
        face_elevation: ArrayView1<'v, f64>,
        node_elevation: ArrayView1<'v, f64>,
    ) -> LocalSurfaceView<'v> {
        let empty_f: ArrayView1<'static, f64> = ArrayView1::from(&[] as &[f64]);
        let empty_i: ArrayView1<'static, i64> = ArrayView1::from(&[] as &[i64]);
        let empty_i2: ArrayView2<'static, i64> =
            ArrayView2::from_shape((0, 0), &[] as &[i64]).expect("empty shape");
        LocalSurfaceView {
            n_face: self.n_face(),
            pix: 0.0,
            face_area: self.face_area.view(),
            face_elevation,
            face_indices: empty_i,
            face_lon: self.face_lon.view(),
            face_lat: self.face_lat.view(),
            face_x: empty_f,
            face_y: empty_f,
            face_z: empty_f,
            n_node: self.n_node(),
            node_elevation,
            node_indices: empty_i,
            node_lon: empty_f,
            node_lat: empty_f,
            node_x: empty_f,
            node_y: empty_f,
            node_z: empty_f,
            n_edge: self.edge_length.len(),
            edge_indices: empty_i,
            edge_length: self.edge_length.view(),
            face_edge_connectivity: self.face_edge_connectivity.view(),
            face_node_connectivity: self.face_node_connectivity.view(),
            face_face_connectivity: self.face_face_connectivity.view(),
            node_face_connectivity: self.node_face_connectivity.view(),
            edge_face_connectivity: self.edge_face_connectivity.view(),
            edge_node_connectivity: empty_i2,
            edge_face_distance: self.edge_face_distance.view(),
            face_proj_x: None,
            face_proj_y: None,
            node_proj_x: None,
            node_proj_y: None,
            face_distance: None,
            face_bearing: None,
            region_radius: None,
            location: None,
            crater_id: None,
            desloped_face_elevation: None,
        }
    }
}

/// The surface values of a region that crater formation changes, in local indices.
pub struct CraterSurfaceState {
    pub face_elevation: Array1<f64>,
    pub node_elevation: Array1<f64>,
    pub ejecta_thickness: Array1<f64>,
}

/// Options controlling crater formation.
pub struct FormCraterOptions {
    /// Scale the ejecta so that its volume matches the excavated volume.
    pub conserve_volume: bool,
    /// Seed of the ray pattern. If `None`, the ejecta has no rays and a uniform intensity of 1.
    pub ray_seed: Option<u64>,
    /// Critical slope (rise over run) for slope collapse. If `None`, slopes are not collapsed.
    pub critical_slope: Option<f64>,
}

/// What happened when a crater was formed.
pub struct FormCraterOutcome {
    /// Volume of the crater relative to the pre-impact surface (m³, negative for a net excavation).
    pub excavated_volume: f64,
    /// Factor the ejecta was scaled by, or `None` if no ejecta was formed (zero excavated volume).
    pub conservation_factor: Option<f64>,
}

/// Forms a crater on a region of the surface, updating `state` in place.
///
/// # Arguments
///
/// * `model` - The morphology model's crater profile.
/// * `region` - The geometry of the region the crater is formed on.
/// * `state` - The current face and node elevations and face ejecta thickness of the region. Updated in place.
/// * `options` - Volume conservation, ray, and slope collapse options.
///
/// # Returns
///
/// The excavated volume and ejecta conservation factor, or an error if the inputs are inconsistent.
pub fn form_crater<P: CraterProfile + ?Sized>(
    model: &P,
    region: &CraterRegionView<'_>,
    state: &mut CraterSurfaceState,
    options: &FormCraterOptions,
) -> Result<FormCraterOutcome, String> {
    region.validate()?;
    let nf = region.n_face();
    let nn = region.n_node();
    if state.face_elevation.len() != nf
        || state.ejecta_thickness.len() != nf
        || state.node_elevation.len() != nn
    {
        return Err("the surface state arrays do not match the size of the region".to_string());
    }

    // Faces followed by nodes, as the crater and ejecta shapes are evaluated at both
    let distance = concatenate![Axis(0), region.face_distance, region.node_distance];
    let bearing =
        concatenate![Axis(0), region.face_bearing, region.node_bearing].mapv(f64::to_radians);
    let original_elevation = concatenate![
        Axis(0),
        state.face_elevation.view(),
        state.node_elevation.view()
    ];

    // 1. Crater shape relative to a reference plane fitted to the surface within one crater diameter
    let reference_radius = model.diameter();
    let mut reference_elevation = if region.node_distance.iter().any(|&d| d <= reference_radius) {
        let x = concatenate![Axis(0), region.face_proj_x, region.node_proj_x];
        let y = concatenate![Axis(0), region.face_proj_y, region.node_proj_y];
        reference_surface(
            x.view(),
            y.view(),
            original_elevation.view(),
            distance.view(),
            reference_radius,
        )?
    } else {
        original_elevation.clone()
    };
    let radius = model.radius();
    for i in 0..distance.len() {
        if distance[i] > radius {
            reference_elevation[i] = original_elevation[i];
        }
    }
    let new_elevation = model.profile(
        distance.view(),
        bearing.view(),
        reference_elevation.view(),
        true,
        false,
    )?;
    let mut elevation_change = &new_elevation - &original_elevation;

    // 2. Excavated volume
    let excavated_volume = face_volume(elevation_change.slice(s![..nf]), region.face_area);

    // 3. Remove any existing ejecta from the crater interior
    for f in inner_crater_faces(region, radius) {
        state.ejecta_thickness[f] = 0.0;
    }

    // 4-5. Ejecta blanket and the degradation it causes. NaN is truthy in Python, so only an exactly zero volume is skipped.
    let mut conservation_factor = None;
    if excavated_volume != 0.0 {
        let zeros = Array1::<f64>::zeros(distance.len());
        let mut thickness =
            model.profile(distance.view(), bearing.view(), zeros.view(), false, true)?;
        let intensity = match options.ray_seed {
            Some(seed) => {
                let intensity = model.ray_intensity(distance.view(), bearing.view(), seed)?;
                thickness *= &intensity;
                intensity
            }
            None => Array1::<f64>::ones(distance.len()),
        };
        // np.maximum(thickness, 0.0), which keeps NaN
        thickness.mapv_inplace(|t| if t >= 0.0 || t.is_nan() { t } else { 0.0 });

        let ejecta_volume = face_volume(thickness.slice(s![..nf]), region.face_area);
        let factor = if options.conserve_volume {
            python_max(-excavated_volume / ejecta_volume, 0.0)
        } else {
            1.0
        };
        thickness.mapv_inplace(|t| t * factor);
        add_data(
            state.ejecta_thickness.view_mut(),
            thickness.slice(s![..nf]),
            false,
            false,
        );

        let burial: Array1<f64> = thickness
            .slice(s![..nf])
            .mapv(|t| model.burial_degradation(t));
        diffuse(region, state, &burial)?;

        let k_degradation = model.degradation();
        let degradation: Array1<f64> = intensity.slice(s![..nf]).mapv(|i| k_degradation * i);
        diffuse(region, state, &degradation)?;

        elevation_change += &thickness;
        conservation_factor = Some(factor);
    }

    // 6. Add the crater and ejecta to the surface
    add_data(
        state.face_elevation.view_mut(),
        elevation_change.slice(s![..nf]),
        false,
        false,
    );
    add_data(
        state.node_elevation.view_mut(),
        elevation_change.slice(s![nf..]),
        false,
        false,
    );

    // 7. Slope collapse
    if let Some(critical_slope) = options.critical_slope {
        let delta = {
            let view = region.kernel_view(state.face_elevation.view(), state.node_elevation.view());
            slope_collapse(critical_slope, &view)?
        };
        apply_face_change(region, state, &delta);
    }

    Ok(FormCraterOutcome {
        excavated_volume,
        conservation_factor,
    })
}

/// Applies one diffusion step to the region, as `LocalSurface.apply_diffusion` does.
fn diffuse(
    region: &CraterRegionView<'_>,
    state: &mut CraterSurfaceState,
    kdiff: &Array1<f64>,
) -> Result<(), String> {
    if kdiff.iter().any(|&k| k < 0.0) {
        return Err("kdiff must be greater than 0.0".to_string());
    }
    // np.max, which returns NaN if any value is NaN
    let kdiff_max = kdiff.iter().fold(f64::NEG_INFINITY, |m, &k| {
        if m.is_nan() || k.is_nan() {
            f64::NAN
        } else if k > m {
            k
        } else {
            m
        }
    });
    if kdiff_max.abs() < PY_VSMALL {
        return Ok(());
    }
    let delta = {
        let view = region.kernel_view(state.face_elevation.view(), state.node_elevation.view());
        apply_diffusion(kdiff.view(), state.face_elevation.view(), &view)?
    };
    apply_face_change(region, state, &delta);
    Ok(())
}

/// Adds a face elevation change, re-interpolates the node elevations, and adds the positive part of the change to the ejecta
/// thickness, as `LocalSurface.update_elevation` followed by `add_data(..., positive_only=True)` does after a diffusion or slope
/// collapse step.
fn apply_face_change(
    region: &CraterRegionView<'_>,
    state: &mut CraterSurfaceState,
    delta: &Array1<f64>,
) {
    add_data(state.face_elevation.view_mut(), delta.view(), false, false);
    let node_elevation = {
        let view = region.kernel_view(state.face_elevation.view(), state.node_elevation.view());
        interpolate_node_elevation_from_faces(&view).expect("interpolation does not fail")
    };
    add_data(
        state.node_elevation.view_mut(),
        node_elevation.view(),
        true,
        false,
    );
    add_data(state.ejecta_thickness.view_mut(), delta.view(), false, true);
}

/// Adds (or, with `overwrite`, assigns) `data` to `target`, as `LocalSurface.add_data` does.
///
/// If every value of `data` is NaN, nothing changes. Otherwise NaN values are treated as 0 when adding, and leave the target
/// unchanged when overwriting. With `positive_only`, the result is clipped at 0 (keeping NaN, as `np.maximum` does).
fn add_data(
    mut target: ArrayViewMut1<'_, f64>,
    data: ArrayView1<'_, f64>,
    overwrite: bool,
    positive_only: bool,
) {
    debug_assert_eq!(target.len(), data.len());
    let n_nan = data.iter().filter(|v| v.is_nan()).count();
    if n_nan > 0 && n_nan == data.len() {
        return;
    }
    for (t, &d) in target.iter_mut().zip(data.iter()) {
        if d.is_nan() {
            if !overwrite {
                *t += 0.0;
            }
        } else if overwrite {
            *t = d;
        } else {
            *t += d;
        }
    }
    if positive_only {
        target.mapv_inplace(|t| if 0.0 >= t { 0.0 } else { t });
    }
}

/// Python's built-in `max(a, b)`: returns `a` unless `b > a`.
#[inline]
fn python_max(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// Volume of a face elevation change: the sum of each face's elevation change times its area.
fn face_volume(elevation: ArrayView1<'_, f64>, face_area: ArrayView1<'_, f64>) -> f64 {
    elevation.dot(&face_area)
}

/// Local indices of the faces of the crater interior, as `LocalSurface.extract_subregion(radius)` selects them: the faces within
/// `radius` of the center plus every face that shares a face, node, or edge with them, in increasing order.
///
/// The region is built from a larger radius with the same kind of buffer, so all of these faces lie inside it.
fn inner_crater_faces(region: &CraterRegionView<'_>, radius: f64) -> Vec<usize> {
    let nf = region.n_face();
    let mut inside = vec![false; nf];
    let mark = |f: i64, inside: &mut Vec<bool>| {
        if f >= 0 && (f as usize) < nf {
            inside[f as usize] = true;
        }
    };
    for seed in 0..nf {
        if region.face_distance[seed] > radius {
            continue;
        }
        inside[seed] = true;
        for &f in region.face_face_connectivity.row(seed) {
            mark(f, &mut inside);
        }
        for &node in region.face_node_connectivity.row(seed) {
            if node >= 0 {
                for &f in region.node_face_connectivity.row(node as usize) {
                    mark(f, &mut inside);
                }
            }
        }
        for &edge in region.face_edge_connectivity.row(seed) {
            if edge >= 0 {
                for &f in region.edge_face_connectivity.row(edge as usize) {
                    mark(f, &mut inside);
                }
            }
        }
    }
    (0..nf).filter(|&f| inside[f]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_data_matches_python_nan_handling() {
        let mut t = Array1::from(vec![1.0, 2.0, 3.0]);
        add_data(
            t.view_mut(),
            Array1::from(vec![f64::NAN, 1.0, -5.0]).view(),
            false,
            true,
        );
        assert_eq!(t.to_vec(), vec![1.0, 3.0, 0.0]);
        let mut t = Array1::from(vec![1.0, 2.0, 3.0]);
        add_data(
            t.view_mut(),
            Array1::from(vec![f64::NAN, 9.0, 8.0]).view(),
            true,
            false,
        );
        assert_eq!(t.to_vec(), vec![1.0, 9.0, 8.0]);
        let mut t = Array1::from(vec![1.0, 2.0]);
        add_data(
            t.view_mut(),
            Array1::from(vec![f64::NAN, f64::NAN]).view(),
            false,
            true,
        );
        assert_eq!(t.to_vec(), vec![1.0, 2.0]);
    }

    #[test]
    fn python_max_semantics() {
        assert_eq!(python_max(-1.0, 0.0), 0.0);
        assert_eq!(python_max(2.0, 0.0), 2.0);
        assert!(python_max(f64::NAN, 0.0).is_nan());
    }
}
