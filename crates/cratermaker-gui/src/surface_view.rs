use std::ptr::slice_from_raw_parts;

use iced::{
    Element, Length, Rectangle, Renderer, Theme,
    advanced::image,
    widget::{canvas, text},
};
use pyo3::{intern, prelude::*, types::IntoPyDict};

use crate::loader::PythonManager;

#[derive(Debug, Clone)]
pub enum Message {}

pub struct SurfaceView {
    plotter: Py<PyAny>,
}

impl SurfaceView {
    pub fn new<'py>(
        py_manager: &mut PythonManager,
        object: Bound<'py, PyAny>,
    ) -> PyResult<Option<Self>> {
        if !object.hasattr(intern!(object.py(), "pyvista_plotter"))? {
            return Ok(None);
        }
        let plotter_class = py_manager
            .get_or_import_module(object.py(), "cratermaker.utils.gui_plotter")?
            .get_class("GuiPlotter")
            .expect("cratermaker.utils.gui_plotter.GuiPlotter should exist")
            .inner
            .bind(object.py())
            .to_owned();
        let plotter = plotter_class.call0()?;

        object.call_method(
            intern!(object.py(), "pyvista_plotter"),
            (),
            Some(&[("plotter", &plotter)].into_py_dict(object.py())?),
        )?;
        plotter.call_method0(intern!(object.py(), "show"))?;
        Ok(Some(Self {
            plotter: plotter.unbind(),
        }))
    }
}

pub fn view_surface<'a>(surface: Option<&'a SurfaceView>) -> Element<'a, Message> {
    match surface {
        Some(viewer) => canvas(Viewer {
            plotter: &viewer.plotter,
        })
        .width(Length::Fill)
        .height(Length::Fill)
        .into(),
        None => text("No surface selected").into(),
    }
}

struct Viewer<'a> {
    plotter: &'a Py<PyAny>,
}

#[derive(Default)]
struct ViewerState {
    scroll_diff_x: f32,
    scroll_diff_y: f32,
}

impl<'a, Message> canvas::Program<Message> for Viewer<'a> {
    type State = ViewerState;

    fn draw(
        &self,
        state: &Self::State,
        renderer: &Renderer,
        theme: &Theme,
        bounds: iced::Rectangle,
        cursor: iced::advanced::mouse::Cursor,
    ) -> Vec<canvas::Geometry<Renderer>> {
        let image = Python::attach(|py| {
            let plotter = self.plotter.bind(py);
            let (width, height) = (bounds.width.ceil() as usize, bounds.height.ceil() as usize);
            plotter
                .call_method1(intern!(py, "set_size"), &(width, height))
                .unwrap();
            plotter.call_method0(intern!(py, "render")).unwrap();
            let ptr_string = plotter
                .call_method0(intern!(py, "frame_data"))
                .unwrap()
                .extract::<String>()
                .unwrap();
            let ptr_int = usize::from_str_radix(ptr_string.split('_').nth(1).unwrap(), 16).unwrap();
            let pixels = unsafe {
                slice_from_raw_parts(ptr_int as *const u8, width * height * 4)
                    .as_ref()
                    .unwrap()
            };
            canvas::Image::new(image::Handle::from_rgba(
                width as u32,
                height as u32,
                pixels,
            ))
        });
        let mut frame = canvas::Frame::new(renderer, bounds.size());
        frame.draw_image(Rectangle::with_size(bounds.size()), image);
        vec![frame.into_geometry()]
    }

    fn update(
        &self,
        state: &mut Self::State,
        event: &iced::Event,
        bounds: Rectangle,
        cursor: iced::advanced::mouse::Cursor,
    ) -> Option<canvas::Action<Message>> {
        Python::attach(|py| {
            let interactor = self
                .plotter
                .bind(py)
                .call_method0(intern!(py, "get_interactor"))
                .unwrap();
            match event {
                iced::Event::Keyboard(event) => match event {
                    iced::keyboard::Event::KeyPressed {
                        key,
                        modified_key,
                        physical_key,
                        location,
                        modifiers,
                        text,
                        repeat,
                    } => {
                        let ctrl = modifiers.control() as i32;
                        let alt = modifiers.alt() as i32;
                        let repeat_count = *repeat as i32;
                        let keycode = key.to_latin(physical_key.clone()).unwrap_or_default();
                        interactor
                            .call_method1(
                                intern!(py, "SetKeyEventInformation"),
                                (ctrl, alt, keycode, repeat_count),
                            )
                            .unwrap();
                        interactor
                            .call_method0(intern!(py, "KeyPressEvent"))
                            .unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::keyboard::Event::KeyReleased {
                        key,
                        modified_key,
                        physical_key,
                        location,
                        modifiers,
                    } => {
                        let ctrl = modifiers.control() as i32;
                        let alt = modifiers.alt() as i32;
                        let keycode = modified_key
                            .to_latin(physical_key.clone())
                            .unwrap_or_default();
                        interactor
                            .call_method1(
                                intern!(py, "SetKeyEventInformation"),
                                (ctrl, alt, keycode),
                            )
                            .unwrap();
                        interactor
                            .call_method0(intern!(py, "KeyReleaseEvent"))
                            .unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::keyboard::Event::ModifiersChanged(modifiers) => {
                        let ctrl = modifiers.control() as i32;
                        let shift = modifiers.shift() as i32;
                        let alt = modifiers.alt() as i32;
                        interactor
                            .call_method1(intern!(py, "SetControlKey"), (ctrl,))
                            .unwrap();
                        interactor
                            .call_method1(intern!(py, "SetShiftKey"), (shift,))
                            .unwrap();
                        interactor
                            .call_method1(intern!(py, "SetAltKey"), (alt,))
                            .unwrap();
                        Some(canvas::Action::capture())
                    }
                },
                iced::Event::Mouse(event) => match event {
                    iced::mouse::Event::CursorEntered => {
                        interactor.call_method0(intern!(py, "EnterEvent")).unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::mouse::Event::CursorLeft => {
                        interactor.call_method0(intern!(py, "LeaveEvent")).unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::mouse::Event::CursorMoved { position } => {
                        interactor
                            .call_method1(
                                intern!(py, "SetEventPosition"),
                                (position.x as i32, position.y as i32),
                            )
                            .unwrap();
                        interactor
                            .call_method0(intern!(py, "MouseMoveEvent"))
                            .unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::mouse::Event::ButtonPressed(button) => {
                        interactor
                            .call_method0(match button {
                                iced::mouse::Button::Left => intern!(py, "LeftButtonPressEvent"),
                                iced::mouse::Button::Right => intern!(py, "RightButtonPressEvent"),
                                iced::mouse::Button::Middle => {
                                    intern!(py, "MiddleButtonPressEvent")
                                }
                                iced::mouse::Button::Back => intern!(py, "FourthButtonPressEvent"),
                                iced::mouse::Button::Forward => {
                                    intern!(py, "FifthButtonPressEvent")
                                }
                                _ => return None,
                            })
                            .unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::mouse::Event::ButtonReleased(button) => {
                        interactor
                            .call_method0(match button {
                                iced::mouse::Button::Left => intern!(py, "LeftButtonReleaseEvent"),
                                iced::mouse::Button::Right => {
                                    intern!(py, "RightButtonReleaseEvent")
                                }
                                iced::mouse::Button::Middle => {
                                    intern!(py, "MiddleButtonReleaseEvent")
                                }
                                iced::mouse::Button::Back => {
                                    intern!(py, "FourthButtonReleaseEvent")
                                }
                                iced::mouse::Button::Forward => {
                                    intern!(py, "FifthButtonReleaseEvent")
                                }
                                _ => return None,
                            })
                            .unwrap();
                        Some(canvas::Action::request_redraw().and_capture())
                    }
                    iced::mouse::Event::WheelScrolled { delta } => {
                        let (x, y) = match delta {
                            iced::mouse::ScrollDelta::Lines { x, y } => (*x, *y),
                            iced::mouse::ScrollDelta::Pixels { x, y } => (*x / 20.0, *y / 20.0),
                        };
                        state.scroll_diff_x += x;
                        state.scroll_diff_y += y;
                        let mut redraw = false;
                        while state.scroll_diff_x < -1.0 {
                            interactor
                                .call_method0(intern!(py, "MouseWheelLeftEvent"))
                                .unwrap();
                            state.scroll_diff_x += 1.0;
                            redraw = true;
                        }
                        while state.scroll_diff_x > 1.0 {
                            interactor
                                .call_method0(intern!(py, "MouseWheelRightEvent"))
                                .unwrap();
                            state.scroll_diff_x -= 1.0;
                            redraw = true;
                        }
                        while state.scroll_diff_y < -1.0 {
                            interactor
                                .call_method0(intern!(py, "MouseWheelBackwardEvent"))
                                .unwrap();
                            state.scroll_diff_y += 1.0;
                            redraw = true;
                        }
                        while state.scroll_diff_y > 1.0 {
                            interactor
                                .call_method0(intern!(py, "MouseWheelForwardEvent"))
                                .unwrap();
                            state.scroll_diff_y -= 1.0;
                            redraw = true;
                        }
                        Some(if redraw {
                            canvas::Action::request_redraw().and_capture()
                        } else {
                            canvas::Action::capture()
                        })
                    }
                },
                _ => None,
            }
        })
    }
}
