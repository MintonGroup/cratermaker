from collections.abc import Sequence
from typing import Literal

import pyvista
from pyvista import (
    BasePlotter,
    CameraPositionOptions,
    ColorLike,
    RenderWindowInteractor,
)
from pyvista.plotting.themes import Theme

import cratermaker.utils._vtk as vtk


class GuiPlotter(BasePlotter):
    """Plotting object to display vtk meshes in the cratermaker gui.

    Parameters
    ----------
    notebook : bool, optional
        When ``True``, the resulting plot is placed inline a jupyter
        notebook.  Assumes a jupyter console is active.  Automatically
        enables ``off_screen``.

    shape : sequence[int], optional
        Number of sub-render windows inside of the main window.
        Specify two across with ``shape=(2, 1)`` and a two by two grid
        with ``shape=(2, 2)``.  By default there is only one render
        window.  Can also accept a string descriptor as shape. E.g.:

        * ``shape="3|1"`` means 3 plots on the left and 1 on the right,
        * ``shape="4/2"`` means 4 plots on top and 2 at the bottom.

    border : bool, optional
        Draw a border around each render window.

    border_color : ColorLike, default: "k"
        Either a string, rgb list, or hex color string.  For example:

            * ``color='white'``
            * ``color='w'``
            * ``color=[1.0, 1.0, 1.0]``
            * ``color='#FFFFFF'``

    window_size : sequence[int], optional
        Window size in pixels.  Defaults to ``[1024, 768]``, unless
        set differently in the relevant theme's ``window_size``
        property.

    line_smoothing : bool, default: False
        If ``True``, enable line smoothing.

    polygon_smoothing : bool, default: False
        If ``True``, enable polygon smoothing.

    lighting : str, default: 'light kit"
        Lighting to set up for the plotter. Accepted options:

        * ``'light kit'``: a vtk Light Kit composed of 5 lights.
        * ``'three lights'``: illumination using 3 lights.
        * ``'none'``: no light sources at instantiation.

        The default is a ``'light kit'`` (to be precise, 5 separate
        lights that act like a Light Kit).

    theme : pyvista.plotting.themes.Theme, optional
        Plot-specific theme.

    image_scale : int, optional
        Scale factor when saving screenshots. Image sizes will be
        the ``window_size`` multiplied by this scale factor.

    Examples
    --------
    >>> import pyvista as pv
    >>> mesh = pv.Cube()
    >>> another_mesh = pv.Sphere()
    >>> pl = pv.Plotter()
    >>> actor = pl.add_mesh(mesh, color="red", style="wireframe", line_width=4)
    >>> actor = pl.add_mesh(another_mesh, color="blue")
    >>> pl.show()

    """

    last_update_time = 0.0

    def __init__(
        self,
        notebook: bool | None = None,
        shape: Sequence[int] | str = (1, 1),
        groups: Sequence[int] | None = None,
        row_weights: Sequence[int] | None = None,
        col_weights: Sequence[int] | None = None,
        border: bool | None = None,
        border_color: ColorLike = "k",
        border_width: float = 2.0,
        window_size: list[int] | None = None,
        line_smoothing: bool = False,  # noqa: FBT001, FBT002
        point_smoothing: bool = False,  # noqa: FBT001, FBT002
        polygon_smoothing: bool = False,  # noqa: FBT001, FBT002
        splitting_position: float | None = None,
        title: str | None = None,
        lighting: Literal["light kit", "three lights", "none"] | None = "light kit",
        theme: Theme | None = None,
        image_scale: int | None = None,
    ) -> None:
        """Initialize a vtk plotting object."""
        super().__init__(
            shape=shape,
            border=border,
            border_color=border_color,
            border_width=border_width,
            groups=groups,
            row_weights=row_weights,
            col_weights=col_weights,
            splitting_position=splitting_position,
            title=title,
            lighting=lighting,
            theme=theme,
            image_scale=image_scale,
        )
        # reset partial initialization flag
        self._initialized = False

        self.ren_win = vtk.vtkRenderWindow()
        # initialize render window
        self.render_window.SetMultiSamples(0)  # type: ignore[union-attr]
        self.render_window.SetBorders(True)  # type: ignore[union-attr]
        if line_smoothing:
            self.render_window.LineSmoothingOn()  # type: ignore[union-attr]
        if point_smoothing:
            self.render_window.PointSmoothingOn()  # type: ignore[union-attr]
        if polygon_smoothing:
            self.render_window.PolygonSmoothingOn()  # type: ignore[union-attr]

        for renderer in self.renderers:
            self.render_window.AddRenderer(renderer)  # type: ignore[union-attr]

        # Add the shadow renderer to allow us to capture interactions within
        # a given viewport
        # https://vtk.org/pipermail/vtkusers/2018-June/102030.html
        number_or_layers = self.render_window.GetNumberOfLayers()  # type: ignore[union-attr]
        current_layer = self.renderer.GetLayer()
        self.render_window.SetNumberOfLayers(number_or_layers + 1)  # type: ignore[union-attr]
        self.render_window.AddRenderer(self.renderers.shadow_renderer)  # type: ignore[union-attr]
        self.renderers.shadow_renderer.SetLayer(current_layer + 1)  # type: ignore[union-attr]
        self.renderers.shadow_renderer.SetInteractive(False)  # type: ignore[union-attr]

        self.render_window.SetUseOffScreenBuffers(True)  # type: ignore[union-attr]
        self.render_window.ShowWindowOff()  # type: ignore[union-attr]

        # Add ren win and interactor
        interactor = vtk.vtkGenericRenderWindowInteractor()
        interactor.EnableRenderOff()

        self.iren = RenderWindowInteractor(self, light_follow_camera=False, interactor=interactor)
        self.iren.set_render_window(self.render_window)
        self.enable_trackball_style()  # type: ignore[call-arg] # internally calls update_style()
        self.reset_key_events()
        self.iren.initialize()
        self.iren.add_observer("KeyPressEvent", self.key_press_event)  # type: ignore[union-attr]

        # Set camera widget based on theme. This requires that an
        # interactor be present.
        if self.theme._enable_camera_orientation_widget:
            self.add_camera_orientation_widget()

        # Set background
        self.set_background(self._theme.background)  # type: ignore[arg-type]

        # Set window size
        self._window_size_unset = False
        if window_size is None:
            self.window_size = self._theme.window_size
            if self.window_size == pyvista.plotting.themes.Theme().window_size:
                self._window_size_unset = True
        else:
            self.window_size = window_size

        if self._theme.depth_peeling.enabled and self.enable_depth_peeling():  # type: ignore[call-arg]
            for renderer in self.renderers:
                renderer.enable_depth_peeling()

        # set anti_aliasing based on theme
        if self.theme.anti_aliasing:
            self.enable_anti_aliasing(self.theme.anti_aliasing)  # type: ignore[arg-type]

        if self.theme.camera.parallel_projection:
            self.enable_parallel_projection()  # type: ignore[call-arg]

        self.parallel_scale = self.theme.camera.parallel_scale

        # some cleanup only necessary for fully initialized plotters
        self._initialized = True

    def show(  # noqa: PLR0917
        self,
        cpos: CameraPositionOptions | None = None,
    ):
        if self.render_window is None:
            msg = "This plotter has been closed and cannot be shown."
            raise RuntimeError(msg)

        # reset unless camera for the first render unless camera is set
        self.camera_position = cpos  # type: ignore[assignment]
        self._on_first_render_request()

        self.render()
        if "vtkDepthOfFieldPass" in self.renderer._render_passes._passes:
            self.render()

        self.iren.update_style()  # type: ignore[union-attr]

    def set_size(self, width: int, height: int):
        self.render_window.SetSize(width, height)  # type: ignore[union-attr]

    def get_size(self) -> tuple[int, int]:
        return self.render_window.GetSize()  # type: ignore[union-attr]

    def frame_data(self):
        self.render_window.MakeCurrent()  # type: ignore[union-attr]
        width, height = self.get_size()
        self.get_interactor().EnableRenderOn()  # type: ignore[union-attr]
        self.get_interactor().Render()  # type: ignore[union-attr]
        self.get_interactor().EnableRenderOff()  # type: ignore[union-attr]
        data = self.render_window.GetPixelData(0, 0, width - 1, height - 1, 0, 0)  # type: ignore[union-attr]
        return data

    def get_interactor(self):
        return self.iren.interactor  # type: ignore[union-attr]
