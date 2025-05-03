use std::f64::{
    self,
    consts::{PI, SQRT_2},
};

use itertools::Itertools;
use ndarray::ArrayView1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::{exceptions::PyValueError, prelude::*};
use rand::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::VSMALL;
const NRAYMAX: i32 = 5;
const NPATT: i32 = 8;
const FRAYREDUCTION: f64 = 0.5;
const EJPROFILE: i32 = 3;

#[pyfunction]
pub fn distribution<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
    ejecta_truncation: f64,
    dorays: bool,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    distribution_internal(
        radial_distance.as_array(),
        initial_bearing.as_array(),
        crater_diameter,
        ejrim,
        ejecta_truncation,
        dorays,
    )
    .map(|v| v.into_pyarray(py))
}

pub fn distribution_internal<'py>(
    radial_distance: ArrayView1<'py, f64>,
    initial_bearing: ArrayView1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
    ejecta_truncation: f64,
    dorays: bool,
) -> PyResult<Vec<f64>> {
    if dorays {
        let intensity = ray_intensity_internal(
            radial_distance,
            initial_bearing,
            crater_diameter,
            ejecta_truncation,
        )?;
        let crater_radius = crater_diameter / 2.0;
        Ok(intensity
            .iter()
            .zip(radial_distance)
            .map(|(&intensity, &radial_distance)| {
                if radial_distance >= crater_radius {
                    intensity * profile_func(radial_distance / crater_radius, ejrim)
                } else {
                    0.0
                }
            })
            .collect())
    } else {
        Ok(profile_internal(radial_distance, crater_diameter, ejrim))
    }
}

#[pyfunction]
pub fn profile<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    Ok(PyArray1::from_vec(
        py,
        profile_internal(radial_distance.as_array(), crater_diameter, ejrim),
    ))
}

pub fn profile_internal<'py>(
    radial_distance: ArrayView1<'py, f64>,
    crater_diameter: f64,
    ejrim: f64,
) -> Vec<f64> {
    let crater_radius = crater_diameter / 2.0;
    radial_distance
        .iter()
        .map(|&r| {
            if r >= crater_radius {
                profile_func(r / crater_radius, ejrim)
            } else {
                0.0
            }
        })
        .collect()
}

pub fn profile_func(r: f64, ejrim: f64) -> f64 {
    ejrim * r.powi(-EJPROFILE)
}

pub fn ray_intensity_internal<'py>(
    radial_distance: ArrayView1<'py, f64>,
    initial_bearing: ArrayView1<'py, f64>,
    crater_diameter: f64,
    ejecta_truncation: f64,
) -> PyResult<Vec<f64>> {
    if radial_distance.len() != initial_bearing.len() {
        return Err(PyValueError::new_err(
            "initial_bearing and radial_distance arrays must be the same size.",
        ));
    }
    let crater_radius = crater_diameter / 2.0;
    let rmax = ejecta_truncation;
    let rmin = 1.0;
    let minray = 2.348 * crater_radius.powf(0.006); // The continuous ejecta blanket radius relative to the crater radius

    let mut rng = rand::rng();

    // Distribute ray patterns evenly around the crater
    let mut thetari = (0..NRAYMAX)
        .map(|i| PI * 2.0 * (i + 1) as f64 / NRAYMAX as f64)
        .collect_vec();
    thetari.shuffle(&mut rng);

    let random_numbers: Vec<f64> = rng.random_iter().take(NPATT as usize).collect_vec();

    let mut intensity: Vec<f64> = (0..radial_distance.len())
        .into_par_iter()
        .map(|i| {
            if *radial_distance.get(i).unwrap() >= crater_radius {
                (0..NPATT as usize)
                    .into_par_iter()
                    .map(|j| {
                        let rn = random_numbers[j];
                        let theta = (initial_bearing.get(i).unwrap() + rn * 2.0 * PI) % (2.0 * PI);
                        let r_pattern = *radial_distance.get(i).unwrap() / crater_radius - rn;
                        FRAYREDUCTION.powi(j as i32)
                            * ray_intensity_func(r_pattern, theta, rmin, rmax, &thetari, minray)
                    })
                    .sum()
            } else {
                0.0
            }
        })
        .collect();
    let max_val = intensity
        .iter()
        .zip(radial_distance.iter())
        .filter_map(|(&val, &r)| if r >= crater_radius { Some(val) } else { None })
        .fold(f64::MIN, |a, b| a.max(b));
    intensity = intensity
        .into_iter()
        .zip(radial_distance)
        .map(|(intensity, &radial_distance)| {
            if radial_distance >= crater_radius {
                intensity / max_val
            } else {
                intensity
            }
        })
        .collect_vec();
    Ok(intensity)
}

#[pyfunction]
pub fn ray_intensity<'py>(
    py: Python<'py>,
    radial_distance: PyReadonlyArray1<'py, f64>,
    initial_bearing: PyReadonlyArray1<'py, f64>,
    crater_diameter: f64,
    ejecta_truncation: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let intensity = ray_intensity_internal(
        radial_distance.as_array(),
        initial_bearing.as_array(),
        crater_diameter,
        ejecta_truncation,
    )?;
    Ok(intensity.into_pyarray(py))
}

fn ray_intensity_func(
    r: f64,
    theta: f64,
    rmin: f64,
    rmax: f64,
    thetari: &[f64],
    minray: f64,
) -> f64 {
    const RAYP: f64 = 4.0;
    if !r.is_finite() || r <= 0.0 || r > rmax {
        return 0.0;
    } else if r < 1.0 {
        1.0
    } else {
        let tmp = (NRAYMAX as f64).powf(RAYP)
            - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln();
        let n = if tmp < 0.0 {
            NRAYMAX // "Nrays" in Minton et al. (2019)
        } else {
            // Exponential decay of ray number with distance
            ((((NRAYMAX as f64).powf(RAYP)
                - ((NRAYMAX as f64).powf(RAYP) - 1.0) * (r / minray).ln() / (rmax / minray).ln())
            .powf(1.0 / RAYP))
            .floor() as i32)
                .clamp(1, NRAYMAX)
        };
        (0..NRAYMAX)
            .into_par_iter()
            .map(|i| {
                let length = minray
                    * ((rmax / minray).ln()
                        * ((((NRAYMAX - i + 2) as f64).powf(RAYP) - 1.0)
                            / ((NRAYMAX as f64).powf(RAYP) - 1.0)))
                        .exp();
                if r > length {
                    return 0.0; // Don't add any material beyond the length of the ray
                }
                let w = (rmax / length).powf(1.0);
                // equation 41 Minton et al. 2019
                let rw = PI / (w * NRAYMAX as f64)
                    * (rmin / r)
                    * (1.0 - (1.0 - w / rmin) * (1.0 - (r / rmin).powi(2)).exp());
                ejecta_ray_func(theta, thetari[i as usize], r, n, rw)
            })
            .sum()
    }
}

fn ejecta_ray_func(theta: f64, thetar: f64, r: f64, n: i32, w: f64) -> f64 {
    let c = w / r;
    let b = thetar;
    let dtheta = f64::min(2.0 * PI - (theta - b).abs(), (theta - b).abs());
    let logval = -dtheta.powi(2) / (2.0 * c.powi(2));
    if logval < VSMALL.ln() {
        0.0
    } else {
        let a = (2.0 * PI).sqrt() / (n as f64 * c * (PI / (2.0 * SQRT_2 * c)).erf());
        a * logval.exp()
    }
}
