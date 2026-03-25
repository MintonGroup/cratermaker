use pyo3::{intern, prelude::*, types::PyType};

#[derive(Clone)]
pub struct Inspect<'py> {
    inspect: Bound<'py, PyModule>,
}

impl<'py> Inspect<'py> {
    pub fn load(py: Python<'py>) -> Self {
        Self {
            inspect: py.import(intern!(py, "inspect")).unwrap(),
        }
    }

    pub fn getmodule(&self, obj: &Bound<'py, PyAny>) -> Bound<'py, PyModule> {
        self.inspect
            .call_method1(intern!(obj.py(), "getmodule"), (obj,))
            .unwrap()
            .cast_into()
            .unwrap()
    }
    pub fn getdoc(&self, obj: &Bound<'py, PyAny>) -> Option<String> {
        self.inspect
            .call_method1(intern!(obj.py(), "getdoc"), (obj,))
            .unwrap()
            .extract()
            .unwrap()
    }
    pub fn getmembers(&self, obj: &Bound<'py, PyAny>) -> Vec<(String, Bound<'py, PyAny>)> {
        self.inspect
            .call_method1(intern!(obj.py(), "getmembers"), (obj,))
            .unwrap()
            .extract()
            .unwrap()
    }
}

pub fn get_indexable_members<'py>(
    inspect: &Inspect<'py>,
    obj: &Bound<'py, PyAny>,
) -> Vec<(String, Bound<'py, PyAny>)> {
    inspect
        .getmembers(obj)
        .into_iter()
        .filter(|(name, _)| !name.starts_with('_'))
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
            module: Module::load(&Inspect::load(py), cratermaker),
        })
    }
}

#[derive(Debug, Clone)]
pub struct Module<'py> {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Bound<'py, PyModule>,
    pub submodules: Vec<Module<'py>>,
    pub classes: Vec<Class<'py>>,
}

#[derive(Debug, Clone)]
pub struct Class<'py> {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Bound<'py, PyType>,
}

impl<'py> Class<'py> {
    pub fn load(inspect: &Inspect, class: Bound<'py, PyType>) -> Self {
        let name = class.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect.getdoc(class.as_any());
        let members = get_indexable_members(inspect, class.as_any());
        Self {
            name,
            docstring,
            inner: class,
        }
    }
}

impl<'py> Module<'py> {
    pub fn load(inspect: &Inspect<'py>, module: Bound<'py, PyModule>) -> Self {
        let mut name = module.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect.getdoc(module.as_any());
        let members = get_indexable_members(inspect, module.as_any());
        let submodules = members
            .iter()
            .filter_map(|(_, item)| item.cast::<PyModule>().ok())
            .filter(|submodule| {
                submodule
                    .name()
                    .unwrap()
                    .extract::<&str>()
                    .unwrap()
                    .starts_with(&name)
            })
            .map(|submodule| Module::load(inspect, submodule.clone()))
            .collect();
        let classes = members
            .iter()
            .filter_map(|(_, item)| item.cast::<PyType>().ok())
            .map(|class| Class::load(inspect, class.clone()))
            .collect();
        if let Some((_, shortname)) = name.rsplit_once('.') {
            name = shortname.to_string();
        }
        Self {
            inner: module,
            docstring,
            submodules,
            name,
            classes,
        }
    }
}
