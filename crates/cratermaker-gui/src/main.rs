pub mod context_menu;
pub mod loader;
pub mod pythonio;
pub mod surface_view;

use iced::widget::scrollable;
use std::sync::{Arc, RwLock};

use iced::{
    Element, Length, Point, Task,
    widget::{button, column, container, opaque, row, scrollable::Scrollbar, space, text},
};
use pyo3::prelude::*;

use crate::{
    context_menu::{context_area, context_menu},
    loader::{Class, Method, PythonManager},
    pythonio::PythonIO,
    surface_view::{SurfaceView, view_surface},
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
struct Variables {
    list: Vec<Arc<RwLock<Variable>>>,
    selected_variable: Option<Arc<RwLock<Variable>>>,
}

fn view_variable_pane(variables: &Variables) -> Element<'_, Message> {
    scrollable(
        column!["Variables"]
            .extend(
                variables
                    .list
                    .iter()
                    .map(|variable| view_variable(variable)),
            )
            .spacing(6)
            .padding(2),
    )
    .into()
}

fn view_variable_info(variable: Option<&Arc<RwLock<Variable>>>) -> Element<'_, Message> {
    let mut col = column!["Variable Info",].spacing(6).padding(2);
    col = if let Some(variable) = variable {
        col.push(text(variable.read().unwrap().value.to_string()))
    } else {
        col.push(text("No variable selected.").style(text::secondary))
    };
    scrollable(col)
        .direction(scrollable::Direction::Both {
            vertical: Scrollbar::new(),
            horizontal: Scrollbar::new(),
        })
        .into()
}

fn view_toolbar<'a>() -> Element<'a, Message> {
    row![
        button("New simulation")
            .style(button::subtle)
            .on_press(Message::CreateSimulation)
    ]
    .spacing(2)
    .padding(2)
    .into()
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
                .width(Length::Fixed(320.0))
                .height(Length::Shrink),
        )
        .into()
    }
}

struct App {
    py_manager: Arc<RwLock<PythonManager>>,
    variables: Variables,
    pythonio: PythonIO,
    surface_view: Option<SurfaceView>,
    context_menu_pos: Option<Point>,
    context_menu_target: Option<ContextMenuTarget>,
}

#[derive(Debug, Clone)]
enum Message {
    CreateSimulation,
    AddVariable(Arc<RwLock<Variable>>),
    RunMethod(Arc<RwLock<Variable>>, Arc<Method>),
    PythonIO(pythonio::Message),
    SurfaceView(surface_view::Message),
    OpenContextMenu(Point, ContextMenuTarget),
    CloseContextMenu,
    FocusVariable(Arc<RwLock<Variable>>),
    DoNothing,
}

impl App {
    fn boot() -> (Self, Task<Message>) {
        Python::attach(|py| {
            let (pythonio_pane, task) = PythonIO::setup(py);
            (
                Self {
                    py_manager: Default::default(),
                    variables: Default::default(),
                    pythonio: pythonio_pane,
                    context_menu_pos: None,
                    context_menu_target: None,
                    surface_view: None,
                },
                task.map(Message::PythonIO),
            )
        })
    }
    fn run() -> iced::Result {
        iced::application(Self::boot, Self::update, Self::view)
            .title("Cratermaker GUI")
            .run()
    }
    fn update(&mut self, message: Message) -> Task<Message> {
        match message {
            Message::CreateSimulation => {
                let py_manager = self.py_manager.clone();
                Task::perform(
                    tokio::task::spawn_blocking(move || {
                        let simulation_class = py_manager
                            .write()
                            .unwrap()
                            .get_or_import_module_unbound("cratermaker.core.simulation")
                            .unwrap()
                            .get_class("Simulation")
                            .expect("cratermaker.core.simulation.Simulation should exist");
                        let value = Python::attach(|py| simulation_class.inner.call0(py).unwrap());
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
                self.variables.list.push(variable);
                Task::none()
            }
            Message::PythonIO(message) => {
                self.pythonio.update(message);
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
            Message::RunMethod(variable, method) => {
                let py_manager = self.py_manager.clone();
                Task::perform(
                    tokio::task::spawn_blocking(move || {
                        Python::attach(|py| {
                            let value = variable
                                .read()
                                .unwrap()
                                .value
                                .bind(py)
                                .call_method0(&method.name)
                                .unwrap();
                            if !value.is_none() {
                                let variable = Variable {
                                    name: "result".to_string(),
                                    class: py_manager
                                        .write()
                                        .unwrap()
                                        .get_or_load_class(&value.get_type())
                                        .ok()
                                        .flatten(),
                                    value: value.unbind(),
                                    focused: false,
                                };
                                Message::AddVariable(Arc::new(RwLock::new(variable)))
                            } else {
                                Message::DoNothing
                            }
                        })
                    }),
                    |res| res.unwrap(),
                )
            }
            Message::FocusVariable(variable) => {
                for var in &self.variables.list {
                    if var.read().unwrap().focused {
                        var.write().unwrap().focused = false;
                    }
                }
                variable.write().unwrap().focused = true;

                self.surface_view = Python::attach(|py| {
                    SurfaceView::new(
                        &mut *self.py_manager.write().unwrap(),
                        variable.read().unwrap().value.bind(py).clone(),
                    )
                })
                .unwrap()
                .or(self.surface_view.take());

                self.variables.selected_variable = Some(variable);
                Task::none()
            }
            Message::DoNothing => Task::none(),
        }
    }
    fn view(&self) -> Element<'_, Message> {
        let main = column![
            container(view_toolbar())
                .style(container::bordered_box)
                .width(Length::Fill)
                .height(Length::Shrink),
            row![
                container(view_variable_pane(&self.variables))
                    .style(container::bordered_box)
                    .width(Length::FillPortion(1))
                    .height(Length::Fill),
                container(view_surface(self.surface_view.as_ref()).map(Message::SurfaceView))
                    .style(container::bordered_box)
                    .width(Length::FillPortion(3))
                    .height(Length::Fill),
                container(view_variable_info(
                    self.variables.selected_variable.as_ref()
                ))
                .style(container::bordered_box)
                .width(Length::FillPortion(2))
                .height(Length::Fill),
            ]
            .spacing(2)
            .width(Length::Fill)
            .height(Length::FillPortion(2)),
            container(self.pythonio.view().map(Message::PythonIO))
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
