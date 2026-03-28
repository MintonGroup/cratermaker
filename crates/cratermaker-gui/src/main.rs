mod context_menu;
mod loader;
mod pythonio;

use std::sync::{Arc, Weak};

use iced::{
    Element, Length, Point, Task,
    widget::{button, column, container, mouse_area, row, text},
};
use pyo3::prelude::*;

use crate::{
    context_menu::{context_area, context_menu},
    loader::{Class, Cratermaker},
    pythonio::PythonIO,
};

#[derive(Debug)]
struct Variable {
    name: String,
    class: Arc<Class>,
    value: Py<PyAny>,
}

fn view_variable(var: &Arc<Variable>) -> Element<'_, Message> {
    context_area(
        container(row![
            container(text(&var.name)).align_left(Length::Fill),
            container(text(&var.class.name).style(text::secondary)).align_right(Length::Fill)
        ])
        .style(container::rounded_box),
    )
    .on_open(|point| Message::OpenContextMenu(point, ContextMenuTarget::Variable(var.clone())))
    .into()
}

#[derive(Default)]
struct VariablePane {
    variables: Vec<Arc<Variable>>,
}

impl VariablePane {
    fn view(&self) -> Element<'_, Message> {
        column![
            "Variables",
            column(
                self.variables
                    .iter()
                    .map(|variable| view_variable(variable))
            )
            .spacing(6)
        ]
        .into()
    }
}

struct Toolbar;

impl Toolbar {
    fn view(&self) -> Element<'_, Message> {
        row![button("New simulation").on_press(Message::CreateSimulation)].into()
    }
}

#[derive(Debug, Clone)]
enum ContextMenuTarget {
    Variable(Arc<Variable>),
}

impl ContextMenuTarget {
    fn view(&self) -> Element<'_, Message> {
        match self {
            ContextMenuTarget::Variable(variable) => container(column(
                variable
                    .class
                    .methods
                    .values()
                    .map(|method| button(&method.name as &str).into()),
            ))
            .style(container::bordered_box)
            .into(),
        }
    }
}

struct App {
    cratermaker: Cratermaker,
    simulation_class: Arc<Class>,
    variable_pane: VariablePane,
    pythonio_pane: PythonIO,
    toolbar: Toolbar,
    context_menu_pos: Option<Point>,
    context_menu_target: Option<ContextMenuTarget>,
}

#[derive(Debug, Clone)]
enum Message {
    CreateSimulation,
    AddVariable(Arc<Variable>),
    PythonIO(pythonio::Message),
    OpenContextMenu(Point, ContextMenuTarget),
    CloseContextMenu,
}

impl App {
    fn boot() -> (Self, Task<Message>) {
        Python::attach(|py| {
            let (pythonio_pane, task) = PythonIO::setup(py);
            let cratermaker = Cratermaker::load(py).unwrap();
            (
                Self {
                    simulation_class: cratermaker
                        .get_class(&["core", "simulation"], "Simulation")
                        .unwrap(),
                    cratermaker,
                    variable_pane: Default::default(),
                    toolbar: Toolbar,
                    pythonio_pane,
                    context_menu_pos: None,
                    context_menu_target: None,
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
                            class: simulation_class,
                            value,
                        };
                        Message::AddVariable(Arc::new(variable))
                    }),
                    |res| res.unwrap(),
                )
            }
            Message::AddVariable(variable) => {
                self.variable_pane.variables.push(variable);
                Task::none()
            }
            Message::PythonIO(message) => {
                self.pythonio_pane.update(message);
                Task::none()
            }
            Message::OpenContextMenu(point, target) => {
                println!("opened at {:?}", point);
                self.context_menu_pos = Some(point);
                self.context_menu_target = Some(target);
                Task::none()
            }
            Message::CloseContextMenu => {
                println!("closed");
                self.context_menu_pos = None;
                self.context_menu_target = None;
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
            container(self.variable_pane.view())
                .style(container::bordered_box)
                .width(Length::Fill)
                .height(Length::FillPortion(2)),
            container(self.pythonio_pane.view().map(Message::PythonIO))
                .width(Length::Fill)
                .height(Length::FillPortion(1)),
        ];
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
