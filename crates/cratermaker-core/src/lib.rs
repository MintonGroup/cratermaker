#![feature(float_erf)] // Required to use f64::erf (https://github.com/rust-lang/rust/issues/136321)
pub mod morphology;
pub mod surface;
pub mod counting;
pub type ArrayResult = Result<numpy::ndarray::Array1<f64>, String>;
const VSMALL: f64 = 10.0 * std::f64::EPSILON;
