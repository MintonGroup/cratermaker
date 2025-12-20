pub mod morphology;
pub mod surface;
pub mod counting;
pub mod crater;
pub type ArrayResult = Result<numpy::ndarray::Array1<f64>, String>;
const VSMALL: f64 = 10.0 * std::f64::EPSILON;
