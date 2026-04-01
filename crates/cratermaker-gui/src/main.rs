pub mod context_menu;
pub mod loader;
pub mod pythonio;

use std::sync::{Arc, RwLock};

use iced::{
    Element, Length, Point, Task,
    widget::{button, column, container, opaque, row, space, text},
};
use pyo3::prelude::*;

use crate::{
    context_menu::{context_area, context_menu},
    loader::{Class, Method, PythonManager},
    pythonio::PythonIO,
};

#[derive(Debug)]
struct Variable {
    name: String,
    class: Option<Arc<Class>>,
    value: Py<PyAny>,
    focused: bool,
}

fn view_variable(var: &Arc<RwLock<Variable>>) -> Element<'_, Message> {
    let variable = var.read().unwrap();
    context_area(
        button(row![
            container(text(variable.name.clone())).align_left(Length::Fill),
            container(
                text(
                    variable
                        .class
                        .as_ref()
                        .map_or("unknown".to_string(), |class| class.name.clone())
                )
                .style(text::secondary)
            )
            .align_right(Length::Fill)
        ])
        .on_press(Message::FocusVariable(var.clone()))
        .style(if variable.focused {
            button::secondary
        } else {
            button::subtle
        }),
    )
    .on_open(|point| Message::OpenContextMenu(point, ContextMenuTarget::Variable(var.clone())))
    .into()
}

#[derive(Default)]
struct VariablesPane {
    variables: Vec<Arc<RwLock<Variable>>>,
}

impl VariablesPane {
    fn view(&self) -> Element<'_, Message> {
        column!["Variables"]
            .extend(
                self.variables
                    .iter()
                    .map(|variable| view_variable(variable)),
            )
            .spacing(6)
            .padding(2)
            .into()
    }
}

#[derive(Default)]
struct VariableInfoPane {
    variable: Option<Arc<RwLock<Variable>>>,
}

impl VariableInfoPane {
    fn view(&self) -> Element<'_, Message> {
        let col = column!["Variable Info",].spacing(6).padding(2);
        if let Some(variable) = &self.variable {
            let lock = variable.read().unwrap();
            if let Some(class) = &lock.class {
                Python::attach(|py| {
                    col.extend(class.properties.iter().map(|(name, value)| {
                        text(format!("{}: {}", name, value.get(&lock.value.bind(py)))).into()
                    }))
                })
            } else {
                col.push(text(lock.value.to_string()))
            }
        } else {
            col.push(text("No variable selected.").style(text::secondary))
        }
        .into()
    }
}

struct Toolbar;

impl Toolbar {
    fn view(&self) -> Element<'_, Message> {
        row![
            button("New simulation")
                .style(button::subtle)
                .on_press(Message::CreateSimulation)
        ]
        .spacing(2)
        .padding(2)
        .into()
    }
}

#[derive(Debug, Clone)]
enum ContextMenuTarget {
    Variable(Arc<RwLock<Variable>>),
}

impl ContextMenuTarget {
    fn contents(&self) -> Element<'_, Message> {
        match self {
            ContextMenuTarget::Variable(var) => {
                let variable = var.read().unwrap();
                if let Some(class) = variable.class.as_ref() {
                    column(class.methods.values().map(|method| {
                        button(text(method.name.clone()))
                            .on_press(Message::RunMethod(var.clone(), method.clone()))
                            .style(button::subtle)
                            .width(Length::Fill)
                            .into()
                    }))
                    .spacing(2)
                    .padding(2)
                    .into()
                } else {
                    space().into()
                }
            }
        }
    }
    fn view(&self) -> Element<'_, Message> {
        let contents = self.contents();
        opaque(
            container(contents)
                .style(container::bordered_box)
                .width(Length::Fixed(320.0)),
        )
        .into()
    }
}

struct App {
    py_manager: PythonManager,
    simulation_class: Arc<Class>,
    variables_pane: VariablesPane,
    pythonio_pane: PythonIO,
    variable_info_pane: VariableInfoPane,
    toolbar: Toolbar,
    context_menu_pos: Option<Point>,
    context_menu_target: Option<ContextMenuTarget>,
}

#[derive(Debug, Clone)]
enum Message {
    CreateSimulation,
    AddVariable(Arc<RwLock<Variable>>),
    RunMethod(Arc<RwLock<Variable>>, Arc<Method>),
    PythonIO(pythonio::Message),
    OpenContextMenu(Point, ContextMenuTarget),
    CloseContextMenu,
    FocusVariable(Arc<RwLock<Variable>>),
}

impl App {
    fn boot() -> (Self, Task<Message>) {
        Python::attach(|py| {
            let (pythonio_pane, task) = PythonIO::setup(py);
            let mut py_manager: PythonManager = Default::default();
            (
                Self {
                    simulation_class: py_manager
                        .get_or_import_module(py, "cratermaker.core.simulation")
                        .unwrap()
                        .get_class("Simulation")
                        .unwrap(),
                    py_manager,
                    variables_pane: Default::default(),
                    toolbar: Toolbar,
                    pythonio_pane,
                    context_menu_pos: None,
                    context_menu_target: None,
                    variable_info_pane: Default::default(),
                },
                task.map(Message::PythonIO),
            )
        })
    }
    fn run() -> iced::Result {
        iced::application(Self::boot, Self::update, Self::view).run()
    }
    fn update(&mut self, message: Message) -> Task<Message> {
        match message {
            Message::CreateSimulation => {
                let simulation_class = self.simulation_class.clone();

                Task::perform(
                    tokio::task::spawn_blocking(|| {
                        let value = Python::attach(|py| {
                            simulation_class.create_method.inner.call0(py).unwrap()
                        });
                        let variable = Variable {
                            name: "simulation".to_string(),
                            class: Some(simulation_class),
                            value,
                            focused: false,
                        };
                        Message::AddVariable(Arc::new(RwLock::new(variable)))
                    }),
                    |res| res.unwrap(),
                )
            }
            Message::AddVariable(variable) => {
                self.variables_pane.variables.push(variable);
                Task::none()
            }
            Message::PythonIO(message) => {
                self.pythonio_pane.update(message);
                Task::none()
            }
            Message::OpenContextMenu(point, target) => {
                self.context_menu_pos = Some(point);
                self.context_menu_target = Some(target);
                Task::none()
            }
            Message::CloseContextMenu => {
                self.context_menu_pos = None;
                self.context_menu_target = None;
                Task::none()
            }
            Message::RunMethod(variable, method) => Task::perform(
                tokio::task::spawn_blocking(move || {
                    let value = Python::attach(|py| {
                        variable
                            .read()
                            .unwrap()
                            .value
                            .call_method0(py, &method.name)
                    })
                    .unwrap();
                    let variable = Variable {
                        name: "simulation".to_string(),
                        class: method.signature.return_type.as_ref().cloned(),
                        value,
                        focused: false,
                    };
                    Message::AddVariable(Arc::new(RwLock::new(variable)))
                }),
                |res| res.unwrap(),
            ),
            Message::FocusVariable(variable) => {
                for var in &self.variables_pane.variables {
                    if var.read().unwrap().focused {
                        var.write().unwrap().focused = false;
                    }
                }
                variable.write().unwrap().focused = true;
                self.variable_info_pane.variable = Some(variable);
                Task::none()
            }
        }
    }
    fn view(&self) -> Element<'_, Message> {
        let main = column![
            container(self.toolbar.view())
                .style(container::bordered_box)
                .width(Length::Fill)
                .height(Length::Shrink),
            row![
                container(self.variables_pane.view())
                    .style(container::bordered_box)
                    .width(Length::FillPortion(1))
                    .height(Length::Fill),
                container(self.variable_info_pane.view())
                    .style(container::bordered_box)
                    .width(Length::FillPortion(2))
                    .height(Length::Fill),
            ]
            .spacing(2)
            .width(Length::Fill)
            .height(Length::FillPortion(2)),
            container(self.pythonio_pane.view().map(Message::PythonIO))
                .width(Length::Fill)
                .height(Length::FillPortion(1)),
        ]
        .spacing(2);
        context_menu(
            main,
            self.context_menu_pos,
            self.context_menu_target
                .as_ref()
                .map(ContextMenuTarget::view),
            Message::CloseContextMenu,
        )
        .into()
    }
}

#[tokio::main]
async fn main() {
    App::run().unwrap()
}
