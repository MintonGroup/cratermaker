mod loader;

use pyo3::prelude::*;

use crate::loader::Cratermaker;

fn main() {
    Python::attach(|py| -> PyResult<()> {
        let cratermaker = Cratermaker::load(py)?;
        println!("{:#?}", cratermaker);
        Ok(())
    })
    .unwrap();
}
