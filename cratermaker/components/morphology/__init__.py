from __future__ import annotations

from math import pi
from typing import Any

from cratermaker._simplemoon import ejecta_functions
from cratermaker.components.surface import Surface
from cratermaker.core.crater import Crater
from cratermaker.utils.component_utils import ComponentBase, import_components
from cratermaker.utils.general_utils import format_large_units


class Morphology(ComponentBase):
    def __init__(self, crater: Crater | None = None, **kwargs: Any) -> None:
        """
        Initialize the Morphology class.

        Parameters
        ----------
        crater : Crater, optional
            The crater currently attached to the morphology model.
        **kwargs : Any
            Additional keyword arguments.

        Raises
        -------

        TypeError
            If the crater is not an instance of Crater.

        """
        super().__init__(**kwargs)
        self.crater = crater

    def __repr__(self) -> str:
        base = super().__repr__()
        if self.crater is None:
            return base
        return (
            f"{base}\n"
            f"Crater final diameter: {format_large_units(self.crater.final_diameter, quantity='length')}\n"
            f"Morphology type: {self.crater.morphology_type}"
        )

    @classmethod
    def maker(
        cls,
        morphology: str | type[Morphology] | Morphology | None = None,
        crater: Crater | None = None,
        **kwargs: Any,
    ) -> Morphology:
        """
        Initialize a component model with the given name or instance.

        Parameters
        ----------
        morphology : str or Morphology or None
            The name of the morphology model to use, or an instance of Morphology. If None, the default "simplemoon" is used.
        crater : Crater, optional
            The crater currently attached to the morphology model.
        **kwargs : Any
            Additional keyword arguments that are required for the specific morphology model being created.

        Returns
        -------
        component
            An instance of the specified component model.

        Raises
        ------
        KeyError
            If the specified morphology model name is not found in the registry.
        TypeError
            If the specified morphology model is not a string or a subclass of Morphology.
        """

        # Call the base class version of make and pass the morphology argument as the component argument
        if morphology is None:
            morphology = "simplemoon"
        morphology = super().maker(component=morphology, crater=crater, **kwargs)
        return morphology

    def form_crater(
        self, surface: Surface, crater: Crater | None = None, **kwargs
    ) -> None:
        """
        This method forms the interior of the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        surface : Surface
            The surface to be altered.
        crater : Crater
            The crater object to be formed. This is optional if it has already been added
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here).
        """
        self.crater = Crater.maker(crater, **kwargs)

        if not isinstance(surface, Surface):
            raise TypeError("surface must be an instance of Surface")
        self.node_index, self.face_index = surface.find_nearest_index(
            self.crater.location
        )

        # Test if the crater is big enough to modify the surface
        rmax = self.rmax(minimum_thickness=surface.smallest_length)
        region_view = surface.extract_region(self.crater.location, rmax)
        if region_view is None:  # The crater is too small to change the surface
            return
        crater_area = pi * rmax**2

        # Check to make sure that the face at the crater location is not smaller than the crater area
        if surface.face_areas[self.face_index] > crater_area:
            return

        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, self.crater.location
        )
        reference_face_elevation, reference_node_elevation = (
            surface.get_reference_surface(
                region_view,
                face_crater_distance,
                node_crater_distance,
                self.crater.location,
                self.crater.final_radius,
            )
        )

        try:
            node_elevation = self.crater_shape(
                node_crater_distance, reference_node_elevation
            )
            surface.node_elevation[region_view.node_indices] = node_elevation

            face_elevation = self.crater_shape(
                face_crater_distance, reference_face_elevation
            )
            surface.face_elevation[region_view.face_indices] = face_elevation
        except Exception as e:
            raise RuntimeError(f"Something went wrong with this crater!\n{self}") from e

        self.form_ejecta(surface, crater=self.crater, **kwargs)
        return

    def form_ejecta(
        self, surface: Surface, crater: Crater | None = None, **kwargs
    ) -> None:
        """
        This method forms the ejecta blanket around the crater by altering the elevation variable of the surface mesh.

        Parameters
        ----------
        surface : Surface
            The surface to be altered.
        **kwargs : dict
            Additional keyword arguments to be passed to internal functions (not used here).
        """
        if crater:
            self.crater = crater

        if not isinstance(surface, Surface):
            raise TypeError("surface must be an instance of Surface")
        self.node_index, self.face_index = surface.find_nearest_index(
            self.crater.location
        )

        # Test if the ejecta is big enough to modify the surface
        rmax = self.rmax(minimum_thickness=surface.smallest_length)
        if not self.ejecta_truncation:
            ejecta_truncation = rmax / self.crater.final_radius
        else:
            ejecta_truncation = self.ejecta_truncation
        region_view = surface.extract_region(self.crater.location, rmax)
        if region_view is None:  # The crater is too small to change the surface
            return
        ejecta_area = pi * rmax**2

        # Check to make sure that the face at the crater location is not smaller than the ejecta blanket area
        if surface.face_areas[self.face_index] > ejecta_area:
            return

        node_crater_distance, face_crater_distance = surface.get_distance(
            region_view, self.crater.location
        )
        node_crater_bearing, face_crater_bearing = surface.get_initial_bearing(
            region_view, self.crater.location
        )
        ejecta_functions.form_ejecta(
            self,
            region_view,
            node_crater_distance,
            face_crater_distance,
            node_crater_bearing,
            face_crater_bearing,
            ejecta_truncation,
            surface.node_elevation,
            surface.face_elevation,
            surface.ejecta_thickness,
            surface.ray_intensity,
        )
        return

    @property
    def crater(self):
        """
        The crater to be created.

        Returns
        -------
        Crater
        """
        return self._crater

    @crater.setter
    def crater(self, value):
        if value is not None and not isinstance(value, Crater):
            raise TypeError("crater must be an instance of Crater")
        self._crater = value


import_components(__name__, __path__, ignore_private=True)
