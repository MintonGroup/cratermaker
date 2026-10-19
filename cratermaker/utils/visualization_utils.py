from typing import Any, ClassVar

import pyvista as pv
from pyvista.plotting.themes import DarkTheme


class CratermakerTheme(DarkTheme):
    _default_name: ClassVar[str] = "cratermaker"

    def __init__(self):
        """Initialize the theme."""
        super().__init__()
        self.lighting_params.roughness = 1.0


def toggle_pyvista_actor(plotter, actor):
    """
    Toggle the visibility of a given actor in a PyVista plotter and updates the plotter.
    """
    actor.SetVisibility(not actor.GetVisibility())
    plotter.update()
    return


def update_pyvista_help_message(plotter, new_message: str | None = None) -> pv.Plotter:
    """
    Add a help message to a PyVista plotter with instructions for user interactions.

    Parameters
    ----------
    plotter : pyvista.Plotter
        The PyVista plotter to which the help message will be added.
    new_message : str, optional
        An additional message to prepend to the default help instructions. If None, only the default instructions will be shown.

    Returns
    -------
    pyvista.Plotter
        The updated PyVista plotter with the help message added and key event for toggling the message set up.
    """
    old_actor = plotter.actors.get("help", None)
    if old_actor is None:
        old_message = "v: Isometric view"
        old_message += "\nUp/Down: Zoom in/out"
        old_message += "\n+/-: Increase/decrease point size"
        old_message += "\nw: Wireframe view"
        old_message += "\ns: Shaded view"
        old_message += "\nv: Isometric camera view"
        old_message += "\nf: Focus and zoom in on a point"
        old_message += "\nr: Reset the camera view"
        old_message += "\nshift+s: Save a screenshot (only on BackgroundPlotter)"
        old_message += "\nshift+c: Enable interactive cell selection/picking"
        old_message += "\nshift+click or middle-click: Pan the rendering scene"
        old_message += "\nleft+click or cmd+click (Mac): Rotate the rendering scene in 3D"
        old_message += "\nctl+click: Rotate the rendering scene in 2D (view-plane)"
        old_message += "\nmouse-wheel or right-click or ctrl+click (Mac): Continuously zoom the rendering scene"
        old_message += "\nh: Toggle this help message"
        old_message += "\nq: Quit"
    else:
        plotter.remove_actor("help")
        old_message = old_actor.GetText(0)
    if new_message is None:
        help_message = old_message
    else:
        help_message = new_message + "\n" + old_message
    help_actor = pv.CornerAnnotation(0, help_message, name="help")
    help_actor.SetVisibility(False)
    plotter.add_actor(help_actor)
    plotter.add_key_event("h", lambda plotter=plotter, help_actor=help_actor: toggle_pyvista_actor(plotter, help_actor))

    return plotter
