use pyo3::{intern, prelude::*};

pub fn get_docstring(obj: &Bound<'_, PyAny>) -> String {
    obj.getattr(intern!(obj.py(), "__doc__"))
        .unwrap()
        .to_string()
}

pub fn get_indexable_attributes<'py>(obj: &Bound<'py, PyAny>) -> Vec<String> {
    obj.dir()
        .unwrap()
        .into_iter()
        .map(|s| s.to_string())
        .filter(|s| !s.starts_with('_'))
        .collect()
}

pub fn get_indexable_items<'py>(obj: &Bound<'py, PyAny>) -> Vec<Bound<'py, PyAny>> {
    let attrs = get_indexable_attributes(obj);

    attrs
        .into_iter()
        .map(|attr| obj.getattr(attr).unwrap())
        .collect()
}

#[derive(Debug)]
pub struct Cratermaker<'py> {
    pub module: Module<'py>,
}

impl<'py> Cratermaker<'py> {
    pub fn load(py: Python<'py>) -> PyResult<Self> {
        let cratermaker = py.import("cratermaker")?;
        Ok(Cratermaker {
            module: Module::load(cratermaker),
        })
    }
}

#[derive(Debug, Clone)]
pub struct Module<'py> {
    pub name: String,
    pub docstring: String,
    pub inner: Bound<'py, PyModule>,
    pub submodules: Vec<Module<'py>>,
}

impl<'py> Module<'py> {
    pub fn load(module: Bound<'py, PyModule>) -> Self {
        let name = module.name().unwrap().to_string();
        let docstring = get_docstring(module.as_any());
        let items = get_indexable_items(module.as_any());
        let submodules = items
            .into_iter()
            .filter_map(|item| item.cast_into::<PyModule>().ok())
            .filter(|submodule| {
                submodule
                    .name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .starts_with(&name)
            })
            .map(|submodule| Module::load(submodule))
            .collect();
        Self {
            inner: module,
            docstring,
            submodules,
            name,
        }
    }
}
