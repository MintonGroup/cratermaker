use std::{
    collections::HashMap,
    ffi::CStr,
    sync::{Arc, RwLock},
};

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

#[derive(Debug, Default)]
pub struct PythonManager {
    /// If a module is None, it is in the process of being initialized.
    /// This is because the progam can end up in an infinite loop when a class references itself, either directly or indirectly. (eg. in return type)
    loaded_modules: HashMap<String, Option<Arc<Module>>>,
}

impl PythonManager {
    pub fn get_loaded_module(&self, name: &str) -> Option<Arc<Module>> {
        self.loaded_modules.get(name).cloned().flatten()
    }
    pub fn get_or_import_module<'py>(
        &mut self,
        py: Python<'py>,
        name: &str,
    ) -> PyResult<Arc<Module>> {
        if let Some(module) = self.get_loaded_module(name) {
            Ok(module)
        } else {
            self.loaded_modules.insert(name.to_owned(), None);
            let module = Arc::new(Module::load(self, py.import(name)?)?);
            *self.loaded_modules.get_mut(name).unwrap() = Some(module.clone());
            Ok(module)
        }
    }
    fn get_or_load_module<'py>(
        &mut self,
        module: Bound<'py, PyModule>,
    ) -> PyResult<Option<Arc<Module>>> {
        let mod_name: String = module.name()?.extract()?;
        if let Some(module) = self.loaded_modules.get(&mod_name) {
            Ok(module.clone())
        } else {
            self.loaded_modules.insert(mod_name.clone(), None);
            let module = Arc::new(Module::load(self, module)?);
            self.loaded_modules.insert(mod_name, Some(module.clone()));
            module.finish_loading(self);
            Ok(Some(module))
        }
    }
    fn get_or_load_class_inner<'py>(
        &mut self,
        class: &Bound<'py, PyType>,
    ) -> PyResult<Option<MaybeUnloadedClass>> {
        let module_py = inspect::getmodule(class.as_any())?;
        let module = self.get_or_load_module(module_py.clone())?;
        let class_name: String = class.name()?.extract()?;
        if let Some(module) = module {
            Ok(module
                .get_class(&class_name)
                .map(|class| MaybeUnloadedClass::init_loaded(class)))
        } else {
            let module_name: String = module_py.name()?.extract()?;
            Ok(Some(MaybeUnloadedClass::init_unloaded(
                module_name,
                class_name,
            )))
        }
    }
    pub fn get_or_load_class<'py>(
        &mut self,
        class: &Bound<'py, PyType>,
    ) -> PyResult<Option<Arc<Class>>> {
        Ok(self.get_or_load_class_inner(class)?.map(|class| {
            class
                .loaded()
                .expect("loaded_modules should not contain any uninitialized values")
        }))
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
    pub fn load<'py>(
        py_manager: &mut PythonManager,
        module: Bound<'py, PyModule>,
    ) -> PyResult<Self> {
        let name = module.name()?.extract::<String>()?;
        let docstring = inspect::getdoc(module.as_any())?;
        let members = inspect::getmembers(module.as_any())?;
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
            .filter(|class| class.module().unwrap().extract::<&str>().unwrap() == &name)
            .map(|class| Class::load(py_manager, Some(&module), class.clone()))
            .map_ok(|class| (class.name.clone(), Arc::new(class)))
            .try_collect()?;
        Ok(Self {
            inner: module.unbind(),
            docstring,
            submodules,
            name,
            classes,
        })
    }
    fn finish_loading(&self, py_manager: &PythonManager) {
        for module in self.submodules.values() {
            module.finish_loading(py_manager);
        }
        for class in self.classes.values() {
            class.create_method.finish_loading(py_manager);
            for method in class.methods.values() {
                method.finish_loading(py_manager);
            }
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
    pub module: String,
}

impl Class {
    pub fn load<'py>(
        py_manager: &mut PythonManager,
        module: Option<&Bound<'py, PyModule>>,
        class: Bound<'py, PyType>,
    ) -> PyResult<Self> {
        let module = module
            .cloned()
            .map(Ok)
            .unwrap_or_else(|| -> PyResult<_> { Ok(class.py().import(class.module()?)?) })?;
        let module_name = module.name()?.extract::<String>()?;
        let name = class.name()?.extract::<String>()?;
        let docstring = inspect::getdoc(class.as_any())?;
        let members = inspect::getmembers(class.as_any())?;
        let mut methods: HashMap<String, Arc<Method>> = members
            .iter()
            .filter(|(name, _)| !name.starts_with('_'))
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
            .filter(|(name, _)| !name.starts_with('_'))
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
            module: module_name,
            methods,
            create_method,
            properties,
        })
    }
}

#[derive(Debug)]
enum MaybeUnloadedClassInner {
    Loaded(Arc<Class>),
    Unloaded { module: String, class: String },
}

#[derive(Debug)]
struct MaybeUnloadedClass(RwLock<MaybeUnloadedClassInner>);

impl MaybeUnloadedClass {
    fn init_loaded(class: Arc<Class>) -> Self {
        Self(RwLock::new(MaybeUnloadedClassInner::Loaded(class)))
    }
    fn init_unloaded(module: String, class: String) -> Self {
        Self(RwLock::new(MaybeUnloadedClassInner::Unloaded {
            module,
            class,
        }))
    }
    fn loaded(&self) -> Option<Arc<Class>> {
        if let MaybeUnloadedClassInner::Loaded(class) = &*self.0.read().unwrap() {
            Some(class.clone())
        } else {
            None
        }
    }
    fn finish_loading(&self, py_manager: &PythonManager) {
        let mut lock = self.0.write().unwrap();
        if let MaybeUnloadedClassInner::Unloaded { module, class } = &*lock
            && let Some(module) = py_manager.get_loaded_module(module)
        {
            *lock = MaybeUnloadedClassInner::Loaded(module.get_class(class).unwrap());
        }
    }
}

#[derive(Debug)]
pub struct Argument {
    pub name: String,
    typ: Option<MaybeUnloadedClass>,
}

impl Argument {
    pub fn get_typ(&self) -> Option<Arc<Class>> {
        self.typ.as_ref().and_then(MaybeUnloadedClass::loaded)
    }
    fn finish_loading(&self, py_manager: &PythonManager) {
        if let Some(class) = &self.typ {
            class.finish_loading(py_manager);
        }
    }
}

#[derive(Debug)]
pub struct Signature {
    pub args: Vec<Argument>,
    return_type: Option<MaybeUnloadedClass>,
}

impl Signature {
    pub fn get_return_type(&self) -> Option<Arc<Class>> {
        self.return_type
            .as_ref()
            .and_then(MaybeUnloadedClass::loaded)
    }
    fn finish_loading(&self, py_manager: &PythonManager) {
        if let Some(class) = &self.return_type {
            class.finish_loading(py_manager);
        }
        for arg in &self.args {
            arg.finish_loading(py_manager);
        }
    }
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
    fn finish_loading<'py>(&self, py_manager: &PythonManager) {
        self.signature.finish_loading(py_manager);
    }
    pub fn load<'py>(
        py_manager: &mut PythonManager,
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
                            .and_then(|typ| py_manager.get_or_load_class_inner(&typ))
                            .ok()
                            .flatten()
                    }),
                })
            })
            .try_collect()?;
        let return_annotation_py =
            signature_pyany.getattr(intern!(method.py(), "return_annotation"))?;
        let return_annotation: Option<&CStr> = return_annotation_py.extract().ok();

        let return_type = return_annotation.and_then(|annotation| {
            Self::parse_annotation(module, annotation)
                .and_then(|typ| py_manager.get_or_load_class_inner(&typ))
                .ok()
                .flatten()
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
    pub fn get<'py>(&self, obj: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        obj.getattr(&self.name)
    }
}
