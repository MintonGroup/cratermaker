use std::{collections::HashMap, sync::Arc};

use pyo3::{
    intern,
    prelude::*,
    types::{PyFunction, PyType},
};

pub mod inspect {
    use pyo3::sync::PyOnceLock;
    use pyo3::{intern, prelude::*};
    static INSPECT: PyOnceLock<Py<PyModule>> = PyOnceLock::new();
    fn inner<'py, 'a>(py: Python<'py>) -> &'a Bound<'py, PyModule> {
        INSPECT
            .get_or_init(py, || PyModule::import(py, "inspect").unwrap().unbind())
            .bind(py)
    }
    pub fn getdoc<'py>(obj: &Bound<'py, PyAny>) -> Option<String> {
        inner(obj.py())
            .call_method1(intern!(obj.py(), "getdoc"), (obj,))
            .unwrap()
            .extract()
            .unwrap()
    }
    pub fn getmembers<'py>(obj: &Bound<'py, PyAny>) -> Vec<(String, Bound<'py, PyAny>)> {
        inner(obj.py())
            .call_method1(intern!(obj.py(), "getmembers"), (obj,))
            .unwrap()
            .extract()
            .unwrap()
    }
}

pub fn get_indexable_members<'py>(obj: &Bound<'py, PyAny>) -> Vec<(String, Bound<'py, PyAny>)> {
    inspect::getmembers(obj)
        .into_iter()
        .filter(|(name, _)| !name.starts_with('_'))
        .collect()
}

#[derive(Debug)]
pub struct Cratermaker {
    pub module: Arc<Module>,
}

impl Cratermaker {
    pub fn load(py: Python<'_>) -> PyResult<Self> {
        let cratermaker = py.import("cratermaker")?;
        Ok(Cratermaker {
            module: Arc::new(Module::load(cratermaker)),
        })
    }
    pub fn get_module(&self, path: &[&str]) -> Option<Arc<Module>> {
        path.iter()
            .fold(Some(self.module.clone()), |current, &name| {
                current.and_then(|current| current.get_submodule(name))
            })
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
    pub submodules: HashMap<String, Arc<Module>>,
    pub classes: HashMap<String, Arc<Class>>,
}

impl Module {
    pub fn load<'py>(module: Bound<'py, PyModule>) -> Self {
        let mut name = module.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect::getdoc(module.as_any());
        let members = get_indexable_members(module.as_any());
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
            .map(|submodule| Module::load(submodule.clone()))
            .map(|submodule| (submodule.name.clone(), Arc::new(submodule)))
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
            .map(|class| Arc::new(Class::load(class.clone())))
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
    pub fn get_submodule(&self, name: &str) -> Option<Arc<Module>> {
        self.submodules.get(name).cloned()
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
    pub methods: HashMap<String, Arc<Method>>,
    pub properties: HashMap<String, Arc<Property>>,
    pub create_method: Arc<Method>,
}

impl Class {
    pub fn load<'py>(class: Bound<'py, PyType>) -> Self {
        let name = class.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect::getdoc(class.as_any());
        let members = get_indexable_members(class.as_any());
        let mut methods: HashMap<String, Arc<Method>> = members
            .iter()
            .filter_map(|(name, item)| item.cast::<PyFunction>().ok().map(|item| (name, item)))
            .map(|(name, method)| {
                (
                    name.clone(),
                    Arc::new(Method::load(method.as_any().clone(), name.clone())),
                )
            })
            .collect();
        let property = PyModule::import(class.py(), "builtins")
            .unwrap()
            .getattr("property")
            .unwrap();
        let properties = members
            .iter()
            .filter(|(_, item)| item.is_instance(&property).unwrap())
            .map(|(name, property)| {
                (
                    name.clone(),
                    Arc::new(Property::load(property.as_any().clone(), name.clone())),
                )
            })
            .collect();
        let create_method = methods.remove("maker").unwrap_or_else(|| {
            let mut method = Method::load(
                class
                    .getattr(intern!(class.py(), "__init__"))
                    .unwrap()
                    .cast_into()
                    .unwrap(),
                "__init__".to_string(),
            );
            method.inner = class.clone().unbind().into_any();
            Arc::new(method)
        });

        Self {
            name,
            docstring,
            inner: class.unbind(),
            methods,
            create_method,
            properties,
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
    pub fn load<'py>(method: Bound<'py, PyAny>, name: String) -> Self {
        let docstring = inspect::getdoc(&method);
        Self {
            name,
            docstring,
            inner: method.unbind(),
        }
    }
}

#[derive(Debug)]
pub struct Property {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Py<PyAny>,
}

impl Property {
    pub fn load<'py>(property: Bound<'py, PyAny>, name: String) -> Self {
        let docstring = inspect::getdoc(&property);
        Self {
            name,
            docstring,
            inner: property.unbind(),
        }
    }
    pub fn get<'py>(&self, obj: &Bound<'py, PyAny>) -> Bound<'py, PyAny> {
        obj.getattr(&self.name).unwrap()
    }
}
