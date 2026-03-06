from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from cratermaker._cratermaker import counting_bindings

from cratermaker.components.crater import Crater, CraterFixed, CraterVariable
from cratermaker.constants import PairOfFloats
from cratermaker.utils.general_utils import format_large_units

if TYPE_CHECKING:
    from cratermaker.components.morphology import Morphology
    from cratermaker.components.surface import LocalSurface

# The factor by which the crater tagging region is extended beyond the final rim.
_RIM_BUFFER_FACTOR = 1.5


class MorphologyCraterVariable(CraterVariable):
    def __init__(self, morphology: Morphology | None = None, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        object.__setattr__(self, "_morphology", morphology)
        object.__setattr__(self, "_ejecta_region", None)
        object.__setattr__(self, "_crater_region", None)
        object.__setattr__(self, "_face_index", None)
        object.__setattr__(self, "_affected_face_indices", None)
        object.__setattr__(self, "_affected_node_indices", None)
        object.__setattr__(self, "_ejecta_rmax", None)
        object.__setattr__(self, "_emplaceable", None)
        object.__setattr__(self, "_Crater", None)
        return

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the crater variable properties.
        """
        dict_repr = super().as_dict()
        keys = (
            "face_index",
            "ejecta_rmax",
            "emplaceable",
            "affected_face_indices",
            "affected_node_indices",
            "morphology",
            "ejecta_region",
            "crater_region",
            "emplaceable",
        )
        for key in keys:
            dict_repr[key] = getattr(self, key)

        return dict_repr

    @property
    def _has_initialized_surface_data(self) -> bool:
        """
        A helper property to check if the morphology's surface data has been initialized, which is necessary for many of the variable properties to be computed.

        This is used to avoid unnecessary computations and potential errors when accessing properties that depend on surface data before it is available.
        """
        return self.morphology is not None and self.morphology.surface.uxds is not None

    @property
    def morphology(self) -> Morphology:
        """The morphology model associated with this crater variable object."""
        return self._morphology

    @property
    def emplaceable(self) -> bool | None:
        """Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted."""
        return self._emplaceable

    @property
    def face_index(self) -> int | None:
        """The index of the face on the surface mesh where the crater is centered."""
        return self._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        """The maximum radius of the ejecta region for this crater, as determined by the morphology model and is based on the distance that the ejecta thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object."""
        return self._ejecta_rmax

    @property
    def affected_face_indices(self) -> set[int] | None:
        """The set of face indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        """The set of node indices on the surface mesh that are affected by this crater, which is determined based on the crater region and is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously."""
        return self._affected_node_indices

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model."""
        return self._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        """The LocalSurface view extracted around the crater center that encompasses a buffered region around the cratered rim, as determined by the morphology model."""
        return self._crater_region


class MorphologyCrater(Crater):
    def __init__(self, crater: Crater | None = None, fixed_cls=CraterFixed, var_cls=MorphologyCraterVariable, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, var_cls=var_cls, **kwargs)
        return

    def __str__(self) -> str:
        output = super().__str__()
        output += f"Ejecta region maximum radius: {format_large_units(self.ejecta_rmax, quantity='length')}\n"
        output += f"\nLarge enough to be emplaced on the grid: {self.emplaceable}\n"
        if self.emplaceable:
            output += f"Face index of crater center: {self.face_index}\n"
            output += f"Crater region: {self.crater_region}\n"
            output += f"Ejecta region: {self.ejecta_region}\n"
        return output

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        location: PairOfFloats | None = None,
        relative_location: dict | None = None,
        check_redundant_inputs: bool = True,
        **kwargs: Any,
    ) -> MorphologyCrater:
        """
        Initialize a MorphologyCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a SimpleMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        location : pair of floats, optional
            The (longitude, latitude) location of the crater.
        relative_location : dict
            Set the crater's location based on the distance and bearing from a reference location. The dictionary is required to have the keys: "distance" and "bearing" with values of distance (in meters) and bearing (in degrees) from the reference point. A third key, "reference_location" should be in the form of a tuple of floats representing the longitude and latitude pair that defines the reference location. The "reference_location" is optional if the morphology.surface contains a `local_location` attribute (such as in the case of a HiResLocal surface type and its derivatives) which will be used if the key is missing or set to None, otherwise it is required. Both location and relative_location cannot be set unless check_redundant_inputs is set to False.
        check_redundant_inputs : bool, optional
            If True, check for redundant inputs such as providing both diameter and radius. Default is True.
        kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`MorphologyCrater.maker() <cratermaker.components.crater.morphologycrater.MorphologyCrater.maker>`.  Refer to its documentation for a detailed description of valid keyword arguments.
        """
        from cratermaker.components.morphology import Morphology

        morphology = Morphology.maker(morphology, **kwargs)
        if relative_location is not None:
            if location is not None and check_redundant_inputs:
                raise ValueError("Cannot provide both location and relative_location")
            if not isinstance(relative_location, dict) or "distance" not in relative_location or "bearing" not in relative_location:
                raise TypeError("relative_location must be a dict with keys 'distance' and 'bearing'.")
            location = morphology.surface.compute_location_from_distance_bearing(**relative_location)

        if crater is None:
            crater = super().maker(location=location, check_redundant_inputs=check_redundant_inputs, **kwargs)
        args = {}

        return cls(
            crater=crater,
            morphology=morphology,
            **args,
        )

    def as_dict(self, ignore_keys: list[str] | tuple[str] = (), skip_complex_data: bool = False, **kwargs) -> dict:
        """
        Return a dictionary representation of the crater properties, including morphology-specific variable properties.

        Parameters
        ----------
        ignore_keys : list[str] or tuple[str], optional
            A list or tuple of property names to ignore when creating the dictionary representation. Default is an empty tuple unless `skip_complex_data` is True, in which case it will be extended to include ("morphology", "affected_face_indices", "affected_node_indices", "ejecta_region", "ejecta_region", "crater_region").
        skip_complex_data : bool, optional
            If True, skip complex data types when creating the dictionary representation. This is useful when serializing the object for saving to a file, as it removes complex data types that may not be serializable. Default is False.
        """
        if skip_complex_data:
            ignore_keys += (
                "morphology",
                "affected_face_indices",
                "affected_node_indices",
                "ejecta_region",
                "crater_region",
            )
        dict_repr = super().as_dict(ignore_keys=ignore_keys, skip_complex_data=skip_complex_data, **kwargs)
        return dict_repr

    def remove_complex_data(self) -> None:
        """
        Remove complex data types from the crater variable properties.

        This is useful for storing a lightweight representation, as it removes complex data types that can be recomputed from the fixed properties and the morphology model when needed. The morphology model is not removed, so that the complex data can be recomputed if needed.
        """
        self._var._affected_face_indices = None
        self._var._affected_node_indices = None
        self._var._ejecta_region = None
        self._var._crater_region = None
        return

    @property
    def face_index(self) -> int | None:
        """The index of the face on the surface mesh where the crater is centered."""
        if self._var._face_index is None and self._has_initialized_surface_data:
            self._var._face_index = self.morphology.surface.find_nearest_face(self.location)
        return self._var._face_index

    @property
    def ejecta_rmax(self) -> float | None:
        """The maximum radius of the ejecta region for this crater, as determined by the morphology model and is based on the distance that the ejecta thickness falls below the value of the `smallest_length` attribute of the morphology's associated Surface object."""
        if self._var._ejecta_rmax is None and self._has_initialized_surface_data:
            self._var._ejecta_rmax = self.morphology.rmax(
                self, minimum_thickness=self.morphology.surface.smallest_length, feature="ejecta"
            )
        return self._var._ejecta_rmax

    @property
    def ejecta_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses the ejecta blanket, as determined by the morphology model.

        This is extracted from the morphology's associated Surface object based on the location of the crater and the ejecta_rmax distance. If the ejecta region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if self._var._ejecta_region is None and self.crater_region is not None:
            self._var._ejecta_region = self.morphology.surface.extract_region(
                location=self.location,
                region_radius=self.ejecta_rmax,
            )
        return self._var._ejecta_region

    @property
    def crater_region(self) -> LocalSurface | None:
        """
        The LocalSurface view extracted around the crater center that encompasses a buffered region around the cratered rim.

        This is extracted from the morphology's associated Surface object based on the location of the crater and a radius that extends by a `_RIM_BUFFER_FACTOR` constant times the crater radius. If the crater region cannot be extracted (e.g. if it is smaller than a single face of the mesh) this property will be set to None and the crater will be marked as not emplaceable.
        """
        if self._var._crater_region is None and self._has_initialized_surface_data:
            self._var._crater_region = self.morphology.surface.extract_region(
                location=self.location,
                region_radius=_RIM_BUFFER_FACTOR * self.radius,
            )
        return self._var._crater_region

    @property
    def affected_face_indices(self) -> set[int] | None:
        """
        The set of face indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_face_indices is None and self._has_initialized_surface_data and self.emplaceable:
            if self.ejecta_region is not None:
                if isinstance(self.ejecta_region.face_indices, slice):
                    self._var._affected_face_indices = set(
                        np.arange(self.morphology.surface.n_face)[self.ejecta_region.face_indices]
                    )
                else:
                    self._var._affected_face_indices = set(self.ejecta_region.face_indices)
            else:
                self._var._affected_face_indices = set()

        return self._var._affected_face_indices

    @property
    def affected_node_indices(self) -> set[int] | None:
        """
        The set of node indices on the surface mesh that are affected by this crater, which is determined based on the crater region.

        This is used by the morphology's queue manager to determine if craters overlap and can be emplaced simultaneously.
        """
        if self._var._affected_node_indices is None and self._has_initialized_surface_data and self.emplaceable:
            if self.ejecta_region is not None:
                if isinstance(self.ejecta_region.node_indices, slice):
                    self._var._affected_node_indices = set(
                        np.arange(self.morphology.surface.n_node)[self.ejecta_region.node_indices]
                    )
                else:
                    self._var._affected_node_indices = set(self.ejecta_region.node_indices)
            else:
                self._var._affected_node_indices = set()

        return self._var._affected_node_indices

    @property
    def emplaceable(self) -> bool | None:
        """Whether this crater is large enough to be emplaced on the surface mesh, which is determined based on whether the crater region could be successfully extracted."""
        if self._var._emplaceable is None and self._has_initialized_surface_data:
            self._var._emplaceable = self.crater_region is not None
        return self._var._emplaceable

    @property
    def measured_rim_height(self) -> float | None:
        """The measured rim height of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_rim_height = counting_bindings.measure_rim_height(self.crater_region, self)
        return self._var._measured_rim_height

    @property
    def measured_floor_depth(self) -> float | None:
        """The measured floor depth of the crater, which is determined based on the morphology model's crater shape and the surface elevation data in the crater region."""
        if self.crater_region is not None:
            self.crater_region.compute_desloped_face_elevation()
            self._var._measured_floor_depth = counting_bindings.measure_floor_depth(self.crater_region, self)
        return self._var._measured_floor_depth

    @property
    def measured_depth_to_diameter(self) -> float | None:
        """
        The measured depth to diameter ratio of the crater.

        This is computed from `measured_rim_height`-`measured_floor_depth`
        """
        if self.crater_region is not None:
            self.crater_region._desloped_face_elevation = None
            floor_depth = self.measured_floor_depth
            rim_height = self.measured_rim_height
            return (rim_height - floor_depth) / self.measured_diameter
        else:
            return None
