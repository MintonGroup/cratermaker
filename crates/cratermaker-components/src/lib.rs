pub mod counting;
pub mod crater;
pub mod morphology;
pub mod surface;
pub type ArrayResult = Result<numpy::ndarray::Array1<f64>, String>;
pub type ArrayResult2D = Result<numpy::ndarray::Array2<f64>, String>;
const VSMALL: f64 = 10.0 * std::f64::EPSILON;
