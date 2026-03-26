mod loader;

use std::sync::Arc;

use iced::{
    Element, Length,
    widget::{button, column, container, row, text},
};
use pyo3::prelude::*;

use crate::loader::{Class, Cratermaker};

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
    variables: Vec<Variable>,
}

impl VariablePane {
    fn view(&self) -> Element<'_, Message> {
        column(self.variables.iter().map(|variable| variable.view()))
            .spacing(6)
            .into()
    }
}

struct Toolbar;

impl Toolbar {
    fn view(&self) -> Element<'_, Message> {
        row![button("new").on_press(Message::CreateObject)].into()
    }
}

struct App {
    cratermaker: Cratermaker,
    variable_pane: VariablePane,
    toolbar: Toolbar,
}

#[derive(Debug, Clone, Copy)]
enum Message {
    CreateObject,
}

impl App {
    fn run() -> iced::Result {
        iced::run(Self::update, Self::view)
    }
    fn update(&mut self, message: Message) {
        match message {
            Message::CreateObject => {
                let class = self
                    .cratermaker
                    .module
                    .submodules
                    .get("components")
                    .unwrap()
                    .submodules
                    .get("crater")
                    .unwrap()
                    .classes
                    .get("Crater")
                    .unwrap()
                    .clone();

                let value = Python::attach(|py| class.create_method.inner.call0(py).unwrap());

                self.variable_pane.variables.push(Variable {
                    name: "new object".to_string(),
                    class,
                    value,
                })
            }
        }
    }
    fn view(&self) -> Element<'_, Message> {
        column![self.toolbar.view(), self.variable_pane.view()].into()
    }
}

impl Default for App {
    fn default() -> Self {
        Python::attach(|py| Self {
            cratermaker: Cratermaker::load(py).unwrap(),
            variable_pane: Default::default(),
            toolbar: Toolbar,
        })
    }
}

fn main() {
    App::run().unwrap()
}
