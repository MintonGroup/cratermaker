from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from cratermaker.components.crater import Crater, CraterFixed
from cratermaker.components.crater.morphologycrater import MorphologyCrater
from cratermaker.utils.general_utils import format_large_units

if TYPE_CHECKING:
    from cratermaker.components.morphology import Morphology


@dataclass(frozen=True, slots=True)
class SimpleMoonCraterFixed(CraterFixed):
    rim_height: float | None = None
    """Original rim height of the crater in meters relative to the reference surface."""
    rim_width: float | None = None
    """Original rim width of the crater in meters."""
    floor_depth: float | None = None
    """Original floor depth of the crater in meters relative to the reference surface."""
    floor_diameter: float | None = None
    """Original floor diameter of the crater in meters."""
    peak_height: float | None = None
    """Original central peak height of the crater in meters relative to the reference surface. None for simple craters."""
    ejrim: float | None = None
    """Original ejecta rim thickness of the crater in meters."""

    @property
    def depth_to_diameter(self) -> float | None:
        """
        The depth to diameter ratio of the crater.

        This is computed from `rim_height`-`floor_depth`
        """
        floor_depth = self.floor_depth
        rim_height = self.rim_height
        if floor_depth is not None and rim_height is not None:
            return (rim_height - floor_depth) / self.diameter
        else:
            return None


class SimpleMoonCrater(MorphologyCrater):
    def __init__(self, crater: Crater | None = None, fixed_cls=SimpleMoonCraterFixed, **kwargs):
        super().__init__(crater=crater, fixed_cls=fixed_cls, **kwargs)
        return

    def __str__(self) -> str:
        base = super().__str__()
        return (
            f"{base}\n"
            f"Rim height: {format_large_units(self.rim_height, quantity='length')}\n"
            f"Rim width: {format_large_units(self.rim_width, quantity='length')}\n"
            f"Floor depth: {format_large_units(self.floor_depth, quantity='length')}\n"
            f"Floor diameter: {format_large_units(self.floor_diameter, quantity='length')}\n"
            f"Central peak height: {format_large_units(self.peak_height, quantity='length') if self.peak_height else 'None'}\n"
            f"Ejecta rim thickness: {format_large_units(self.ejrim, quantity='length')}\n"
        )

    @classmethod
    def maker(
        cls,
        crater: Crater | None = None,
        morphology: Morphology | None = None,
        **kwargs: Any,
    ) -> SimpleMoonCrater:
        """
        Initialize a SimpleMoonCrater object either from an existing Crater object or from parameters.

        This generates a specialized Crater object with morphology parameters.

        Parameters
        ----------
        crater : Crater, optional
            The crater object to be converted into a SimpleMoonCrater. If None, then a new crater is created using the provided parameters.
        morphology : Morphology, optional
            The morphology model to use for generating morphology parameters.
        kwargs : Any
            The keyword arguments provided are passed down to :py:meth:`cratermaker.morphology.MorphologyCrater.maker`.  Refer to its documentation for a detailed description of valid keyword arguments.
        """
        from cratermaker.components.morphology import Morphology

        morphology = Morphology.maker(morphology, **kwargs)
        if crater is None:
            crater = super().maker(morphology=morphology, **kwargs)
        args = {}
        diameter_m = crater.diameter
        diameter_km = diameter_m * 1e-3

        if crater.morphology_type in ["simple", "transitional"]:
            args["rim_height"] = 0.043 * diameter_km**1.014 * 1e3
            args["rim_width"] = 0.257 * diameter_km**1.011 * 1e3
            args["floor_depth"] = -0.224 * diameter_km**1.010 * 1e3
            args["floor_diameter"] = 0.200 * diameter_km**1.143 * 1e3
            args["peak_height"] = None
        elif crater.morphology_type in ["complex", "peakring", "multiring"]:
            args["rim_height"] = 0.236 * diameter_km**0.399 * 1e3
            args["rim_width"] = 0.467 * diameter_km**0.836 * 1e3
            args["floor_depth"] = -1.044 * diameter_km**0.301 * 1e3
            args["floor_diameter"] = min(0.187 * diameter_km**1.249 * 1e3, 0.9 * diameter_m)
            args["peak_height"] = 0.032 * diameter_km**0.900 * 1e3
        else:
            raise ValueError(f"Unknown morphology type: {crater.morphology_type}")

        args["ejrim"] = 0.14 * (diameter_m * 0.5) ** 0.74

        return cls(
            crater=crater,
            morphology=morphology,
            **args,
        )
