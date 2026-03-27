mod loader;
mod pythonio;

use std::sync::Arc;

use iced::{
    Element, Length, Task,
    widget::{button, column, container, row, text},
};
use pyo3::prelude::*;

use crate::{
    loader::{Class, Cratermaker},
    pythonio::PythonIO,
};

#[derive(Debug)]
struct Variable {
    name: String,
    class: Arc<Class>,
    value: Py<PyAny>,
}

impl Variable {
    fn view(&self) -> Element<'_, Message> {
        container(row![
            container(text(&self.name)).align_left(Length::Fill),
            container(text(&self.class.name).style(text::secondary)).align_right(Length::Fill)
        ])
        .style(container::rounded_box)
        .into()
    }
}

#[derive(Default)]
struct VariablePane {
    variables: Vec<Arc<Variable>>,
}

impl VariablePane {
    fn view(&self) -> Element<'_, Message> {
        column![
            "Variables",
            column(self.variables.iter().map(|variable| variable.view())).spacing(6)
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

struct App {
    cratermaker: Cratermaker,
    simulation_class: Arc<Class>,
    variable_pane: VariablePane,
    pythonio_pane: PythonIO,
    toolbar: Toolbar,
}

#[derive(Debug, Clone)]
enum Message {
    CreateSimulation,
    AddVariable(Arc<Variable>),
    PythonIO(pythonio::Message),
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
        }
    }
    fn view(&self) -> Element<'_, Message> {
        column![
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
        ]
        .into()
    }
}

#[tokio::main]
async fn main() {
    App::run().unwrap()
}
