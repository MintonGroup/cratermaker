from __future__ import annotations

import csv
from abc import abstractmethod
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import uxarray as uxr
from numpy.random import Generator

from cratermaker.core.crater import Crater
from cratermaker.utils import export
from cratermaker.utils.component_utils import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_NAME = "crater_id"
_TALLY_LONG_NAME = "Unique crater identification number"

_N_LAYER = (
    8  # The number of layers used for tagging faces with crater ids. This allows a single face to contain multiple crater ids
)
_MIN_FACE_FOR_COUNTING = 5
_RIM_BUFFER_FACTOR = 1.2  # The factor by which the crater taggin region is extended beyond the final rim.


class Counting(ComponentBase):
    _registry: dict[str, Counting] = {}

    _CRATER_DIR = "crater_data"

    """
    Base class for all crater counting models. It defines the interface for tallying the observable craters on a surface.

    Parameters
    ----------
    surface : Surface | LocalSurface
        The surface or local surface view to be counted. 
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        surface: Surface | LocalSurface,
        **kwargs: Any,
    ):
        self.surface = surface
        rng = kwargs.pop("rng", surface.rng)
        simdir = kwargs.pop("simdir", surface.simdir)
        super().__init__(rng=rng, simdir=simdir, **kwargs)
        object.__setattr__(self, "_true_crater_list", [])

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
        **kwargs: Any,
    ) -> Counting:
        """
        Initialize a crater counting model based on the provided name or class.

        Parameters
        ----------
        counting : str, Counting or None, default=None
            The name of the counting model to initialize. If None, the default model is used.
        surface : Surface | LocalSurface
            The surface or local surface view to be counted.
        **kwargs : Any
            Additional keyword arguments.

        Returns
        -------
        Counting
            An instance of the specified counting model.

        Raises
        ------
        KeyError
            If the specified counting model name is not found in the registry.
        TypeError
            If the specified counting model is not a string or a subclass of Scaling.
        """
        if counting is None:
            counting = "minton2019"

        if surface is None:
            raise ValueError("surface must be provided")

        counting = super().maker(
            component=counting,
            surface=surface,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nSurface: {self.surface}\n"

    def reset(self):
        """
        Remove all craters count records from the surface.
        """
        dims = ("n_face", "layer")
        data = np.zeros((self.surface.n_face, self.n_layer), dtype=np.uint32)
        uxda = uxr.UxDataArray(
            data=data,
            dims=dims,
            attrs={"long_name": _TALLY_LONG_NAME},
            name=_TALLY_NAME,
            uxgrid=self.surface.uxgrid,
        )

        self.surface._uxds[_TALLY_NAME] = uxda
        self._observed = {}
        return

    def add(self, crater: Crater, region: LocalSurface | None = None):
        """
        Add a crater to the surface.

        Parameters
        ----------
        crater : Crater
            The crater to be added to the surface.
        region : LocalSurface, optional
            A LocalSurface region that contains the crater inside. If not supplied, then the associated surface property is used.
        """
        from cratermaker.components.surface import LocalSurface

        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if _TALLY_NAME not in self.surface.uxds:
            self.reset()

        self.true_crater_list.append(crater)
        self.observed[crater.id] = crater
        region_radius = _RIM_BUFFER_FACTOR * crater.final_radius
        # Tag a region just outside crater rim with the id
        if region is None:
            count_region = self.surface.extract_region(location=crater.location, region_radius=region_radius)
        elif isinstance(region, LocalSurface):
            count_region = region.extract_subregion(subregion_radius=region_radius)
        else:
            raise TypeError("region must be a LocalSurface or None")

        if count_region and count_region.n_face >= _MIN_FACE_FOR_COUNTING:
            insert_layer = -1
            for i in reversed(range(self.n_layer)):
                if np.any(self.surface.uxds[_TALLY_NAME].isel(layer=i).data[count_region.face_indices] > 0):
                    # Gather the unique id values for the current layer
                    unique_ids = np.unique(self.surface.uxds[_TALLY_NAME].data[count_region.face_indices, i])
                    removes = [
                        id for id, v in self.observed.items() if v.id in unique_ids and v.final_diameter < crater.final_diameter
                    ]

                    # For every id that appears in the removes list, set it to 0 in the data array
                    if removes:
                        data = self.surface.uxds[_TALLY_NAME].data[count_region.face_indices, :]
                        for remove in removes:
                            data[data == remove] = 0
                        self.surface.uxds[_TALLY_NAME].data[count_region.face_indices, :] = data
                if insert_layer == -1 and np.all(self.surface.uxds[_TALLY_NAME].isel(layer=i).data[count_region.face_indices] == 0):
                    insert_layer = i
            if insert_layer == -1:
                raise ValueError("Crater counting layers are full")
            data = self.surface.uxds[_TALLY_NAME].data[count_region.face_indices, :]
            data[:, insert_layer] = crater.id
            self.surface.uxds[_TALLY_NAME].data[count_region.face_indices, :] = data

        return

    @abstractmethod
    def tally(self, region: LocalSurface | None = None) -> None: ...

    @property
    def surface(self):
        """
        Surface mesh data for the simulation. Set during initialization.
        """
        return self._surface

    @surface.setter
    def surface(self, value):
        from cratermaker.components.surface import LocalSurface, Surface

        if not isinstance(value, (Surface | LocalSurface)):
            raise TypeError("surface must be an instance of Surface or LocalSurface")
        self._surface = value

    @property
    def n_layer(self) -> int:
        """
        Number of layers in the counting model.
        """
        return _N_LAYER

    @property
    def observed(self) -> dict[int, Crater]:
        """
        List of observed craters on the surface.

        Returns
        -------
            A dict with the crater id as the key and the Crater object as the values
        """
        return self._observed

    @property
    def true_crater_list(self):
        """
        The list of craters that have been emplaced in the simulation.
        """
        return self._true_crater_list

    @property
    def output_dir(self) -> Path | None:
        """
        The output directory for the surface. If None, the surface does not have an output directory set.
        """
        if self._output_dir is None:
            self._output_dir = self.simdir / self.__class__._CRATER_DIR
        return self._output_dir

    def dump_crater_lists(self, interval_number: int = 0) -> None:
        """
        Dump the crater lists to a file and reset the true crater list.

        Parameters
        ----------
        interval_number : int, default=0
            The interval number for the output file naming.
        """
        crater_dir = self.output_dir
        crater_dir.mkdir(parents=True, exist_ok=True)
        truefilename = crater_dir / f"true_crater_list{interval_number:06d}.csv"

        # Convert current crater list to dicts, splitting location into longitude/latitude
        new_data = []
        for c in self.true_crater_list:
            d = asdict(c)
            if "location" in d:
                lon, lat = d.pop("location")
                d["longitude"] = lon
                d["latitude"] = lat
            new_data.append(d)

        # If the file already exists, read it and merge
        if truefilename.exists():
            with truefilename.open("r", newline="") as f:
                reader = csv.DictReader(f)
                existing_data = list(reader)
            combined_data = existing_data + new_data
            # Sort by final_diameter descending
            combined_data = sorted(combined_data, key=lambda d: -float(d["final_diameter"]))
        else:
            combined_data = new_data

        # Write merged data back to file
        if combined_data:
            with truefilename.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=combined_data[0].keys())
                writer.writeheader()
                writer.writerows(combined_data)
            combined_data = {k: np.array([d[k] for d in combined_data]) for k in combined_data[0]}
            export.crater_layer(combined_data, self.surface, interval_number, layer_name="True Craters")

        self._true_crater_list = []
        return


import_components(__name__, __path__, ignore_private=True)
