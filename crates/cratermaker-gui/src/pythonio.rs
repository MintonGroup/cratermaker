use iced::Length;
use iced::widget::text_editor;
use iced::{
    Element, Task,
    widget::text_editor::{Edit, Motion},
};
use std::sync::{Arc, Mutex};
use tokio::sync::mpsc::{self, Receiver, Sender};
use tokio_stream::wrappers::ReceiverStream;

use iced::widget::text_editor::Action;
use pyo3::prelude::*;

#[pyclass]
struct StdoutHandler {
    buffer: Mutex<String>,
    channel: Sender<String>,
}

impl StdoutHandler {
    fn new() -> (Self, Receiver<String>) {
        let (tx, rx) = mpsc::channel(5);
        (
            StdoutHandler {
                buffer: Mutex::new(String::new()),
                channel: tx,
            },
            rx,
        )
    }
    fn write_to_editor(&self, text: String) {
        self.channel.blocking_send(text).unwrap();
    }
}

#[pymethods]
impl StdoutHandler {
    fn write<'py>(&self, data: &str) {
        let mut buf = self.buffer.lock().unwrap();
        let mut start = 0;

        for (i, c) in data.char_indices() {
            if c == '\n' {
                let mut line = buf.clone();
                line.push_str(&data[start..(i + 1)]);
                self.write_to_editor(line);
                buf.clear();
                start = i + 1;
            }
        }

        if start < data.len() {
            buf.push_str(&data[start..]);
        }
    }

    fn flush(&self) {
        let mut buf = self.buffer.lock().unwrap();
        if !buf.is_empty() {
            self.write_to_editor(std::mem::take(&mut buf));
        }
    }
}

pub struct PythonIO {
    content: text_editor::Content,
}

#[derive(Debug, Clone)]
pub enum Message {
    EditorAction(Action),
    StdoutLog(String),
}

impl PythonIO {
    pub fn setup<'py>(py: Python<'py>) -> (Self, Task<Message>) {
        let sys = py.import("sys").unwrap();

        let (handler, rx) = StdoutHandler::new();
        let stdout = Py::new(py, handler).unwrap();
        sys.setattr("stdout", stdout).unwrap();

        (
            Self {
                content: text_editor::Content::new(),
            },
            Task::stream(ReceiverStream::new(rx)).map(Message::StdoutLog),
        )
    }
    pub fn update(&mut self, message: Message) {
        match message {
            Message::EditorAction(action) => match action {
                Action::Edit(_) => (),
                other => self.content.perform(other),
            },
            Message::StdoutLog(text) => {
                let prev_cursor = self.content.cursor();
                self.content.perform(Action::Move(Motion::DocumentEnd));
                self.content
                    .perform(Action::Edit(Edit::Paste(Arc::new(text))));
                self.content.move_to(prev_cursor);
            }
        }
    }
    pub fn view<'a>(&'a self) -> Element<'a, Message> {
        text_editor(&self.content)
            .placeholder("Python log")
            .height(Length::Fill)
            .on_action(Message::EditorAction)
            .into()
    }
}
