use std::{collections::HashMap, ffi::CStr, sync::Arc};

use itertools::Itertools;
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
    pub fn getdoc<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Option<String>> {
        inner(obj.py())
            .call_method1(intern!(obj.py(), "getdoc"), (obj,))?
            .extract()
    }
    pub fn getmembers<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Vec<(String, Bound<'py, PyAny>)>> {
        inner(obj.py())
            .call_method1(intern!(obj.py(), "getmembers"), (obj,))?
            .extract()
    }
    pub fn signature<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        inner(obj.py()).call_method1(intern!(obj.py(), "signature"), (obj,))
    }
    pub fn getmodule<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyModule>> {
        Ok(inner(obj.py())
            .call_method1(intern!(obj.py(), "getmodule"), (obj,))?
            .cast_into()?)
    }
}

pub fn get_indexable_members<'py>(
    obj: &Bound<'py, PyAny>,
) -> PyResult<Vec<(String, Bound<'py, PyAny>)>> {
    Ok(inspect::getmembers(obj)?
        .into_iter()
        .filter(|(name, _)| !name.starts_with('_'))
        .collect())
}

#[derive(Debug, Default)]
pub struct PythonManager {
    loaded_modules: HashMap<String, Arc<Module>>,
}

impl PythonManager {
    pub fn get_loaded_module(&self, name: &str) -> Option<Arc<Module>> {
        self.loaded_modules.get(name).cloned()
    }
    pub fn get_or_import_module<'py>(
        &mut self,
        py: Python<'py>,
        name: &str,
    ) -> PyResult<Arc<Module>> {
        if let Some(module) = self.get_loaded_module(name) {
            Ok(module)
        } else {
            let module = Arc::new(Module::load(self, py.import(name)?)?);
            self.loaded_modules.insert(name.to_owned(), module.clone());
            Ok(module)
        }
    }
    pub fn get_or_load_module<'py>(
        &mut self,
        module: Bound<'py, PyModule>,
    ) -> PyResult<Arc<Module>> {
        let mod_name: String = module.name()?.extract()?;
        if let Some(module) = self.get_loaded_module(&mod_name) {
            Ok(module)
        } else {
            let module = Arc::new(Module::load(self, module)?);
            self.loaded_modules.insert(mod_name, module.clone());
            Ok(module)
        }
    }
    pub fn get_or_load_class<'py>(&mut self, typ: Bound<'py, PyType>) -> PyResult<Arc<Class>> {
        let module = self.get_or_load_module(inspect::getmodule(typ.as_any())?)?;
        Ok(module
            .get_class(typ.name()?.extract()?)
            .expect("module should contain the class"))
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
    pub fn load<'py>(py_manager: &PythonManager, module: Bound<'py, PyModule>) -> PyResult<Self> {
        let mut name = module.name()?.extract::<String>()?;
        let docstring = inspect::getdoc(module.as_any())?;
        let members = get_indexable_members(module.as_any())?;
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
            .map(|submodule| Module::load(py_manager, submodule.clone()))
            .map_ok(|submodule| (submodule.name.clone(), Arc::new(submodule)))
            .try_collect()?;
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
            .map(|class| Class::load(py_manager, Some(&module), class.clone()))
            .map_ok(|class| (class.name.clone(), Arc::new(class)))
            .try_collect()?;
        if let Some((_, shortname)) = name.rsplit_once('.') {
            name = shortname.to_string();
        }
        Ok(Self {
            inner: module.unbind(),
            docstring,
            submodules,
            name,
            classes,
        })
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
    pub fn load<'py>(
        py_manager: &PythonManager,
        module: Option<&Bound<'py, PyModule>>,
        class: Bound<'py, PyType>,
    ) -> PyResult<Self> {
        let module = module
            .cloned()
            .map(Ok)
            .unwrap_or_else(|| -> PyResult<_> { Ok(class.py().import(class.module()?)?) })?;
        let name = class.name().unwrap().extract::<String>().unwrap();
        let docstring = inspect::getdoc(class.as_any())?;
        let members = get_indexable_members(class.as_any())?;
        let mut methods: HashMap<String, Arc<Method>> = members
            .iter()
            .filter_map(|(name, item)| item.cast::<PyFunction>().ok().map(|item| (name, item)))
            .map(|(name, method)| -> PyResult<_> {
                Ok((
                    name.clone(),
                    Arc::new(Method::load(
                        py_manager,
                        method.as_any().clone(),
                        &module,
                        name.clone(),
                    )?),
                ))
            })
            .try_collect()?;
        let property = PyModule::import(class.py(), "builtins")
            .unwrap()
            .getattr("property")
            .unwrap();
        let properties = members
            .iter()
            .filter(|(_, item)| item.is_instance(&property).unwrap())
            .map(|(name, property)| -> PyResult<_> {
                Ok((
                    name.clone(),
                    Arc::new(Property::load(property.as_any().clone(), name.clone())?),
                ))
            })
            .try_collect()?;
        let create_method = methods
            .remove("maker")
            .map(Ok)
            .unwrap_or_else(|| -> PyResult<_> {
                let mut method = Method::load(
                    py_manager,
                    class
                        .getattr(intern!(class.py(), "__init__"))
                        .unwrap()
                        .cast_into()
                        .unwrap(),
                    &module,
                    "__init__".to_string(),
                )?;
                method.inner = class.clone().unbind().into_any();
                Ok(Arc::new(method))
            })?;

        Ok(Self {
            name,
            docstring,
            inner: class.unbind(),
            methods,
            create_method,
            properties,
        })
    }
}

#[derive(Debug)]
pub struct Argument {
    pub name: String,
    pub typ: Option<Arc<Class>>,
}

#[derive(Debug)]
pub struct Signature {
    pub args: Vec<Argument>,
    pub return_type: Option<Arc<Class>>,
}

#[derive(Debug)]
pub struct Method {
    pub name: String,
    pub docstring: Option<String>,
    pub signature: Signature,
    pub inner: Py<PyAny>,
}

impl Method {
    fn parse_annotation<'py>(
        module: &Bound<'py, PyModule>,
        annotation: &CStr,
    ) -> PyResult<Bound<'py, PyType>> {
        Ok(module
            .py()
            .eval(annotation, Some(&module.dict()), None)?
            .cast_into()?)
    }
    pub fn load<'py>(
        py_manager: &PythonManager,
        method: Bound<'py, PyAny>,
        module: &Bound<'py, PyModule>,
        name: String,
    ) -> PyResult<Self> {
        let py = method.py();
        let docstring = inspect::getdoc(&method)?;
        let signature_pyany = inspect::signature(&method)?;

        let args = signature_pyany
            .getattr(intern!(py, "parameters"))?
            .call_method0(intern!(py, "items"))?
            .try_iter()?
            .map(|val| -> PyResult<_> {
                let (name, parameter): (String, Bound<'py, PyAny>) = val?.extract()?;
                let annotation_py = parameter.getattr(intern!(py, "annotation"))?;
                let annotation: Option<&CStr> = annotation_py.extract().ok();

                Ok(Argument {
                    name,
                    typ: annotation.and_then(|annotation| {
                        Self::parse_annotation(module, annotation)
                            .and_then(|typ| Ok(Arc::new(Class::load(py_manager, None, typ)?)))
                            .ok()
                    }),
                })
            })
            .try_collect()?;
        let return_annotation_py =
            signature_pyany.getattr(intern!(method.py(), "return_annotation"))?;
        let return_annotation: Option<&CStr> = return_annotation_py.extract().ok();

        let return_type = return_annotation.and_then(|annotation| {
            Self::parse_annotation(module, annotation)
                .and_then(|typ| Ok(Arc::new(Class::load(py_manager, None, typ)?)))
                .ok()
        });

        Ok(Self {
            name,
            docstring,
            inner: method.unbind(),
            signature: Signature { args, return_type },
        })
    }
}

#[derive(Debug)]
pub struct Property {
    pub name: String,
    pub docstring: Option<String>,
    pub inner: Py<PyAny>,
}

impl Property {
    pub fn load<'py>(property: Bound<'py, PyAny>, name: String) -> PyResult<Self> {
        let docstring = inspect::getdoc(&property)?;
        Ok(Self {
            name,
            docstring,
            inner: property.unbind(),
        })
    }
    pub fn get<'py>(&self, obj: &Bound<'py, PyAny>) -> Bound<'py, PyAny> {
        obj.getattr(&self.name).unwrap()
    }
}
