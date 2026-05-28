use crate::ArrayResult;
use numpy::ndarray::prelude::*;
use pyo3::FromPyObject;
use crate::morphology::basicmoon::{crater_profile_function, ejecta_profile_function};

/// Defines crater dimensions for surface modification computations.
///
/// Used to parameterize the final crater size in meters.
#[derive(FromPyObject, Clone, Debug)]
pub struct RealMoonCrater {
    pub id: u32,
    pub diameter: f64,
    pub radius: f64,
    pub semimajor_axis: f64,
    pub semiminor_axis: f64,
    pub orientation: f64,
    pub transient_diameter: f64,
    pub projectile_diameter: f64,
    pub projectile_velocity: f64,
    pub projectile_angle: f64,
    pub projectile_density: f64,
    pub location: (f64, f64),
    pub morphology_type: String,
    pub measured_semimajor_axis: f64,
    pub measured_semiminor_axis: f64,
    pub measured_orientation: f64,
    pub measured_diameter: f64,
    pub measured_radius: f64,
    pub measured_location: (f64, f64),
    pub time: Option<f64>,
    pub floor_elevation: f64,
    pub floor_radius: f64,
    pub wall_curvature: f64,
    pub rim_width: f64,
    pub rim_elevation: f64,
    pub rim_flank_radius: f64,
    pub rimdrop: f64,
    pub ejrim: f64,
    pub ejprofile: f64,
    pub peak_height: f64,
    pub peak_width: f64,
    pub peak_offset: f64,
}

// Computes a either a crater and/or ejecta 1D profile array from input radial distances and reference elevations using the model of
// Minton et al. (2026)
//
//
/// # Arguments
///
/// * `radial_distances` - 1D array of radial distances from crater center (in meters).
/// * `reference_elevation_array` - 1D array of reference elevations corresponding to each radius.
/// * `crater_radius` - Radius of the crater's rim (in meters).
/// * `floor_elevation` - Elevation of the crater floor relative to the reference surface (in meters and should be negative).
/// * `floor_radius` - Radius of the flat portion of the crater floor (in meters).
/// * `wall_curvature` - Parameter controlling the curvature of the crater wall (>1 for more curvature)
/// * `rim_width` - Width of the crater rim (in meters).
/// * `rim_elevation` - Height of the crater rim above mean surface level (in meters).
/// * `rimdrop` - Exponent for the rim dropoff function (typically -4 to -6)
/// * `ejrim` - Rim elevation adjustment parameter for the exterior dropoff.
/// * `eprofile` - Exponent for the ejecta dropoff function (typically -3)
/// * `peak_height` - Height of the central peak above the crater floor (in meters).
/// * `peak_width` - Width of the central peak (in meters).
/// * `peak_offset` - Radial offset of the central peak from the crater center (in meters). "concentric").
/// * `include_crater` - Whether to include the crater profile in the output (true/false).
/// * `include_ejecta` - Whether to include the ejecta profile in the output (true/false).
///
/// # Returns
///
/// * A NumPy array of modified elevations based on the crater model.
///
/// # Errors
///
/// Returns a `PyValueError` if the input arrays have mismatched lengths.
pub fn realmoon_profile(
    radial_distances: ArrayView1<'_, f64>,
    reference_elevations: ArrayView1<'_, f64>,
    crater: &RealMoonCrater,
    include_crater: bool,
    include_ejecta: bool,
) -> ArrayResult {
    assert_eq!(radial_distances.len(), reference_elevations.len());

    // Compute the weighted elevation profile relative to the reference plane
    let ninc = radial_distances
        .iter()
        .filter(|&&x| x <= crater.radius)
        .count();
    let meanref = if ninc == 0 {
        *radial_distances
            .iter()
            .zip(reference_elevations)
            .min_by(|(&radius_a, _), (&radius_b, _)| radius_a.partial_cmp(&radius_b).unwrap())
            .unwrap()
            .1
    } else {
        radial_distances
            .iter()
            .zip(reference_elevations)
            .filter(|(&r, _)| r <= crater.radius)
            .map(|(_, &e)| e)
            .sum::<f64>()
            / ninc as f64
    };
    let min_elevation = meanref + crater.floor_elevation;

    Ok(Array1::from_iter(
        reference_elevations
            .iter()
            .zip(radial_distances.iter().copied())
            .map(|(href, r)| {
                let mut hcrat = crater_profile_function(r, crater.radius, crater.floor_elevation, crater.floor_radius, crater.wall_curvature, crater.rim_width, crater.rim_elevation, crater.rimdrop, crater.peak_height, crater.peak_width, crater.peak_offset);
                let mut hej = ejecta_profile_function(r, crater.radius, crater.ejrim, crater.ejprofile);
                if r < crater.radius && r > crater.floor_radius {
                    hej += (hcrat - (crater.rim_elevation - crater.ejrim)).max(0.0);
                }

                if include_crater {
                    if r > crater.radius || hcrat > 0.0 {
                        hcrat = (hcrat - hej).max(0.0);
                    } 
                } else {
                    hcrat = 0.0;
                }
                if !include_ejecta {
                    hej = 0.0;
                }

                let h = href + hcrat + hej;

                if r <= crater.radius {
                    h.max(min_elevation)
                } else {
                    h
                }
            }),
    ))
}

