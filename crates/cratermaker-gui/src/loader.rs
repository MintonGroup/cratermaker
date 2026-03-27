use std::{collections::HashMap, sync::Arc};

use pyo3::{
    intern,
    prelude::*,
    types::{PyFunction, PyType},
};

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
pub struct Cratermaker {
    pub module: Module,
}

impl Cratermaker {
    pub fn load(py: Python<'_>) -> PyResult<Self> {
        let cratermaker = py.import("cratermaker")?;
        Ok(Cratermaker {
            module: Module::load(&Inspect::load(py), cratermaker),
        })
    }
    pub fn get_module(&self, path: &[&str]) -> Option<&Module> {
        self.module.get_submodule_tree(path)
    }
    pub fn get_class(&self, module_path: &[&str], name: &str) -> Option<Arc<Class>> {
        self.get_module(module_path)
            .and_then(|module| module.get_class(name))
    }
}

#[derive(Debug)]
pub struct Module {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Py<PyModule>,
    pub submodules: HashMap<String, Module>,
    pub classes: HashMap<String, Arc<Class>>,
}

impl Module {
    pub fn load<'py>(inspect: &Inspect<'py>, module: Bound<'py, PyModule>) -> Self {
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
                    .to_str()
                    .unwrap()
                    .starts_with(&name)
            })
            .map(|submodule| Module::load(inspect, submodule.clone()))
            .map(|submodule| (submodule.name.clone(), submodule))
            .collect();
        let classes = members
            .iter()
            .filter_map(|(_, item)| item.cast::<PyType>().ok())
            .filter(|class| {
                class
                    .fully_qualified_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .starts_with(&name)
            })
            .map(|class| Arc::new(Class::load(inspect, class.clone())))
            .map(|class| (class.name.clone(), class))
            .collect();
        if let Some((_, shortname)) = name.rsplit_once('.') {
            name = shortname.to_string();
        }
        Self {
            inner: module.unbind(),
            docstring,
            submodules,
            name,
            classes,
        }
    }
    pub fn get_submodule(&self, name: &str) -> Option<&Module> {
        self.submodules.get(name)
    }
    pub fn get_submodule_tree(&self, path: &[&str]) -> Option<&Module> {
        path.iter().fold(Some(self), |current, &name| {
            current.and_then(|current| current.get_submodule(name))
        })
    }
    pub fn get_class(&self, name: &str) -> Option<Arc<Class>> {
        self.classes.get(name).cloned()
    }
}

#[derive(Debug)]
pub struct Class {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Py<PyType>,
    pub methods: HashMap<String, Method>,
    pub create_method: Method,
}

impl Class {
    pub fn load<'py>(inspect: &Inspect<'py>, class: Bound<'py, PyType>) -> Self {
        let name = class.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect.getdoc(class.as_any());
        let members = get_indexable_members(inspect, class.as_any());
        let mut methods: HashMap<String, Method> = members
            .iter()
            .filter_map(|(name, item)| item.cast::<PyFunction>().ok().map(|item| (name, item)))
            .map(|(name, method)| {
                (
                    name.clone(),
                    Method::load(inspect, method.as_any().clone(), name.clone()),
                )
            })
            .collect();
        let create_method = methods.remove("maker").unwrap_or_else(|| {
            let mut method = Method::load(
                inspect,
                class
                    .getattr(intern!(class.py(), "__init__"))
                    .unwrap()
                    .cast_into()
                    .unwrap(),
                "__init__".to_string(),
            );
            method.inner = class.clone().unbind().into_any();
            method
        });

        Self {
            name,
            docstring,
            inner: class.unbind(),
            methods,
            create_method,
        }
    }
}

#[derive(Debug)]
pub struct Method {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Py<PyAny>,
}

impl Method {
    pub fn load<'py>(inspect: &Inspect<'py>, method: Bound<'py, PyAny>, name: String) -> Self {
        let docstring = inspect.getdoc(&method);
        Self {
            name,
            docstring,
            inner: method.unbind(),
        }
    }
}
