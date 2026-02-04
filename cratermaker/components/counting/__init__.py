from __future__ import annotations

from abc import abstractmethod
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from cratermaker._cratermaker import counting_bindings
from shapely.geometry import GeometryCollection
from shapely.ops import transform
from tqdm import tqdm
from vtk import vtkPolyData

from cratermaker import __version__ as cratermaker_version
from cratermaker.components.crater import Crater
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.utils.general_utils import get_saved_interval_numbers

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_LONG_NAME = "Unique crater identification number"

# The number of layers used for tagging faces with crater ids. This allows a single face to contain multiple crater ids
_N_LAYER = 8

# The minimum number of faces required in a region to perform crater counting
_MIN_FACE_FOR_COUNTING = 20

# The factor by which the crater tagging region is extended beyond the final rim.
_RIM_BUFFER_FACTOR = 1.2

# The factor by radius over which the local region that is extracted to evaluate the crater rim
_EXTENT_RADIUS_RATIO = 2.0

# The factor by radius over which the local region that is extracted to fit the crater rim for scoring
_FITTING_RADIUS_RATIO = 3.0

# The factor by radius over which the local region that is extracted to measure rim height and floor depth
_MEASURING_RADIUS_RATIO = 1.2


class Counting(ComponentBase):
    _registry: dict[str, Counting] = {}

    """
    Base class for all crater counting models. It defines the interface for tallying the observable craters on a surface.

    Parameters
    ----------
    surface : Surface | LocalSurface
        The surface or local surface view to be counted.
    reset : bool, optional
        Flag to indicate whether to reset the count and delete any old output files. Default is True.
    ask_overwrite : bool, optional
        If True, prompt the user for confirmation before deleting files. Default is False.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        surface: Surface | LocalSurface,
        reset: bool = True,
        ask_overwrite: bool = False,
        **kwargs: Any,
    ):
        from cratermaker.components.surface import Surface

        super().__init__(reset=reset, ask_overwrite=ask_overwrite, **kwargs)

        object.__setattr__(self, "_emplaced", [])
        object.__setattr__(self, "_observed", {})
        object.__setattr__(self, "_output_dir_name", "craters")
        object.__setattr__(self, "_output_file_prefix", "craters")
        object.__setattr__(self, "_output_file_extension", "nc")
        self._surface = Surface.maker(surface, reset=reset, ask_overwrite=ask_overwrite, **kwargs)
        self._output_file_pattern += [f"*{self._output_file_prefix}*.{self._output_file_extension}"]

        if not reset:
            observed_interval_numbers, observed_data_file_list = get_saved_interval_numbers(
                output_dir=self.output_dir,
                output_file_prefix=f"observed_{self._output_file_prefix}",
                output_file_extension=self._output_file_extension,
            )
            if observed_data_file_list:
                observed = xr.open_dataset(observed_data_file_list[-1])
                observed = self.from_xarray(observed)
                for crater in observed:
                    self._observed[crater.id] = crater
        return

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
        reset: bool = True,
        ask_overwrite: bool = False,
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
        reset : bool, optional
            Flag to indicate whether to reset the count and delete any old output files. Default is True
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before deleting files. Default is False.
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

        counting = super().maker(
            component=counting,
            surface=surface,
            reset=reset,
            ask_overwrite=ask_overwrite,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        base = super().__str__()
        return f"{base}\nSurface: {self.surface.name}"

    def reset(self, ask_overwrite: bool = False, **kwargs: Any) -> None:
        """
        Remove all craters count records from the surface and delete any output files.

        Parameters
        ----------
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before deleting files. Default is False.
        **kwargs : Any
            Additional keyword arguments for subclasses.
        """
        self.surface.add_tag(name="crater_id", long_name=_TALLY_LONG_NAME, tag=None, n_layer=_N_LAYER)
        self._emplaced = []
        self._observed = {}

        super().reset(ask_overwrite=ask_overwrite, **kwargs)
        return

    def add(self, crater: Crater, count_region: LocalSurface | None = None, **kwargs: Any):
        """
        Add a crater to the surface.

        Parameters
        ----------
        crater : Crater
            The crater to be added to the surface.
        count_region : LocalSurface, optional
            A LocalSurface region that contains the crater inside. If not supplied, then the associated surface property is used.
        **kwargs : Any
            Additional keyword arguments to pass to the surface add_tag method.
        """
        from cratermaker.components.surface import LocalSurface

        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        count_radius = _RIM_BUFFER_FACTOR * crater.final_radius
        # Tag a region just outside crater rim with the id
        if count_region is None:
            count_region = self.surface.extract_region(location=crater.location, region_radius=count_radius)
        elif isinstance(count_region, LocalSurface):
            if count_radius <= count_region.region_radius:
                count_region = count_region.extract_subregion(subregion_radius=count_radius)
        else:
            raise TypeError("region must be a LocalSurface or None")

        if count_region and count_region.n_face >= _MIN_FACE_FOR_COUNTING:
            self.surface.add_tag(
                name="crater_id",
                long_name=_TALLY_LONG_NAME,
                tag=crater.id,
                n_layer=_N_LAYER,
                region=count_region,
                **kwargs,
            )

            # Check to make sure the crater was recorded to the surface. HiResLocal surfaces may not record craters in the superdomain
            if crater.id in self.surface.uxds.crater_id:
                self.emplaced.append(crater)
                self.observed[crater.id] = crater
            else:
                return

            # Cookie cutting: remove any smaller craters that are overlapped by this new crater
            for i in range(_N_LAYER):
                unique_ids = np.unique(self.surface.uxds.crater_id.data[count_region.face_indices, i])
                if len(unique_ids) > 0:
                    # Compute cookie cutting removes list
                    removes = [
                        id for id, v in self.observed.items() if v.id in unique_ids and v.final_diameter < crater.final_diameter
                    ]

                    # For every id that appears in the removes list, set it to 0 in the data array
                    if removes:
                        for remove_id in removes:
                            self.surface.remove_tag(name="crater_id", tag=remove_id, region=count_region)
                            if not np.any(self.surface.uxds.crater_id.data == remove_id):
                                self.observed.pop(remove_id, None)

        return

    def remove(self, crater_id: int) -> None:
        """
        Remove all instances of a crater id from the surface.

        Parameters
        ----------
        crater_id : int
            The ID of the crater to be removed.
        """
        self.surface.remove_tag(name="crater_id", tag=crater_id)
        self.observed.pop(crater_id, None)
        return

    def fit_rim(self, crater: Crater, tol=0.01, nloops=10, score_quantile=0.95, fit_center=False, fit_ellipse=False) -> Crater:
        """
        Find the rim region of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to find the rim region.
        tol : float, optional
            The tolerance for the rim fitting algorithm. Default is 0.001.
        nloops : int, optional
            The number of iterations for the rim fitting algorithm. Default is 10.
        score_quantile : float, optional
            The quantile of rim scores to consider. Default is 0.95.
        fit_center : bool, optional
            If True, fit the crater center as well. Default is False.
        fit_ellipse : bool, optional
            If True, fit an ellipse to the rim, otherwise fit a circle. Default is False.

        Returns
        -------
        Crater
            A new Crater object with updated measured rim parameters.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        location, ap, bp, orientation = counting_bindings.fit_rim(
            self.surface,
            crater,
            tol,
            nloops,
            score_quantile,
            fit_center,
            fit_ellipse,
        )

        crater_fit = Crater.maker(
            crater,
            measured_semimajor_axis=ap,
            measured_semiminor_axis=bp,
            measured_orientation=np.degrees(orientation),
            measured_location=location,
        )

        return crater_fit

    def score_rim(self, crater: Crater, quantile=0.95, gradmult=1.0, curvmult=1.0, heightmult=1.0) -> None:
        """
        Score the rim region of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to score the rim region.
        quantile : float, optional
            The quantile of rim scores to consider. Default is 0.95.
        gradmult : float, optional
            Gradient multiplier for scoring. Default is 1.0.
        curvmult : float, optional
            Curvature multiplier for scoring. Default is 1.0.
        heightmult : float, optional
            Height multiplier for scoring. Default is 1.0.

        Returns
        -------
        Updates the attached surface object with the rim score for this crater.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        region = self.surface.extract_region(
            location=crater.measured_location, region_radius=_FITTING_RADIUS_RATIO * crater.measured_radius
        )
        orig_elevation = region.face_elevation.copy()
        reference_elevation = region.get_reference_surface(only_faces=True)
        region.update_elevation(-reference_elevation, overwrite=False)
        region = counting_bindings.score_rim(region, crater, quantile, gradmult, curvmult, heightmult)
        region.update_elevation(orig_elevation, overwrite=True)

        return region

    def measure_crater_depth(self, crater: Crater) -> Crater:
        """
        Measure the rim height and floor depth of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to measure the depth.

        Returns
        -------
        crater
            The updated Crater object with the new measured depth properties.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        region = self.surface.extract_region(
            location=crater.measured_location, region_radius=_MEASURING_RADIUS_RATIO * crater.measured_radius
        )
        orig_elevation = region.face_elevation.copy()
        reference_elevation = region.get_reference_surface(only_faces=True)
        region.update_elevation(-reference_elevation, overwrite=False)
        rim_height = counting_bindings.measure_rim_height(region, crater)
        floor_depth = counting_bindings.measure_floor_depth(region, crater)
        region.update_elevation(orig_elevation.data, overwrite=True)
        crater = Crater.maker(crater, measured_rim_height=rim_height, measured_floor_depth=floor_depth)
        return crater

    def tally(self, region: LocalSurface | None = None, quiet: bool = False, **kwargs: Any) -> dict[int:Crater]:
        """
        Tally the craters on the surface using.

        Parameters
        ----------
        region : LocalSurface, optional
            A LocalSurface region to count. If not supplied, then the associated surface property is used.
        quiet : bool, optional
            If True, suppress progress output. Default is False.

        Returns
        -------
        dict[int:Crater]
            A dictionary of observed craters indexed by their ID.

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        if region is None:
            region = self.surface
            id_array = self.surface.crater_id
        elif isinstance(region, LocalSurface):
            id_array = region.crater_id
        else:
            raise TypeError(f"Expected a LocalSurface, but got {type(region).__name__}.")

        unique_ids = np.unique(id_array[id_array > 0])
        remove_ids = []

        if quiet:
            iterable = unique_ids
        else:
            iterable = tqdm(
                unique_ids,
                total=len(unique_ids),
                desc="Counting craters",
                unit="craters",
                position=1,
                leave=False,
            )

        for id in iterable:
            # Check if we have orphaned crater ids for some reason and remove them
            if id not in self.observed:
                remove_ids.append(id)
                continue

            # Update the crater size measurement before computing the degradation and visibility functions
            crater = self.observed[id]
            # TODO: Make the fit_rim function more reliable before turning it on permanently
            # crater = self.fit_rim(crater=crater, fit_center=False, fit_ellipse=False, **kwargs)
            crater = self.measure_degradation_state(crater, **kwargs)
            Kd = crater.measured_degradation_state
            Kv = self.visibility_function(crater, **kwargs)
            if Kd >= Kv:
                remove_ids.append(id)
            else:
                self.observed[id] = crater  # Save the updated measurements to the observed tally

        if len(remove_ids) > 0:
            if quiet:
                iterable = remove_ids
            else:
                iterable = tqdm(
                    remove_ids,
                    total=len(remove_ids),
                    desc="Removing craters",
                    unit="craters",
                    position=2,
                    leave=False,
                )
            for id in iterable:
                self.remove(id)

        return

    @abstractmethod
    def measure_degradation_state(self, crater: Crater, **kwargs: Any) -> Crater: ...

    @abstractmethod
    def visibility_function(self, crater: Crater, Kv1: float = 0.17, gamma: float = 2.0, **kwargs: Any) -> float: ...

    @staticmethod
    def to_xarray(craters: dict[int, Crater] | list[Crater]) -> xr.Dataset:
        """
        Convert a list or dictionary of Crater objects to an xarray Dataset.

        Parameters
        ----------
        craters : dict[int, Crater] | list[Crater]
            A dictionary or list of Crater objects to convert.

        Returns
        -------
        xr.Dataset
            An xarray Dataset containing the crater data.
        """
        if len(craters) == 0:
            return xr.Dataset()
        new_data = []
        if isinstance(craters, dict):
            craters = craters.values()
        for c in craters:
            d = asdict(c)
            if "location" in d:
                lon, lat = d.pop("location")
                d["longitude"] = lon
                d["latitude"] = lat
            if "measured_location" in d:
                mlon, mlat = d.pop("measured_location")
                d["measured_longitude"] = mlon
                d["measured_latitude"] = mlat
            d = xr.Dataset(data_vars=d).set_coords("id").expand_dims(dim="id")
            d["id"].attrs["long_name"] = _TALLY_LONG_NAME
            new_data.append(d)

        return xr.concat(new_data, dim="id")

    def merge_with_file(
        self, craters: dict[int, Crater] | list[Crater], filename: Path | str, save_merged: bool = True
    ) -> xr.Dataset | None:
        """
        Merge a list or dictionary of Crater objects with an existing file.

        Parameters
        ----------
        craters : dict[int, Crater] | list[Crater]
            A dictionary or list of Crater objects to merge.
        filename : Path | str
            The path to the file to merge with.
        save_merged : bool, optional
            If True, save the merged data back to the file. Default is True.

        Returns
        -------
        xr.Dataset | None
            An xarray Dataset containing the merged crater data, or None if no data.
        """
        # Convert into an xarray dataset
        combined_data = self.to_xarray(craters)

        # If the file already exists, read it and merge
        if filename.exists():
            with xr.open_dataset(filename) as ds:
                combined_data = xr.merge([combined_data, ds])

        # Write merged data back to file
        if save_merged and combined_data:
            combined_data.to_netcdf(filename)
        return combined_data

    def save(self, interval_number: int = 0, **kwargs: Any) -> None:
        """
        Dump the crater lists to a file and reset the emplaced crater list.

        Parameters
        ----------
        interval_number : int, default=0
            The interval number for the output file naming.
        **kwargs : Any
            Additional keyword arguments (ignored)
        """
        emplaced_filename = (
            self.output_dir / f"emplaced_{self._output_file_prefix}{interval_number:06d}.{self._output_file_extension}"
        )
        observed_filename = (
            self.output_dir / f"observed_{self._output_file_prefix}{interval_number:06d}.{self._output_file_extension}"
        )

        if self.emplaced:
            self.merge_with_file(self.emplaced, emplaced_filename)
        if self.observed:
            self.merge_with_file(self.observed, observed_filename)
        return

    def export(
        self,
        craters: Crater | list[Crater] | dict[int, Crater] | str,
        name: str | None = None,
        interval_number: int = 0,
        driver: str = "GPKG",
        ask_overwrite: bool = True,
        **kwargs: Any,
    ) -> None:
        """
        Exports crater lists to a file.

        Parameters
        ----------
        craters : Crater |list[Crater] | dict[int, Crater] | str
            A crater or list or dictionary of Crater objects to export, or 'emplaced' or 'observed' to export those lists.
        name : str, default=None
            The name used for the file name or layer name. If None, uses default naming convention.
        interval_number : int, optional
            The interval number to append to the output file. Default is 0.
        driver : str, default='GPKG'
            The file format to save. Supported formats are 'VTK', 'GPKG', 'ESRI Shapefile', 'CSV', 'SCC'.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments to pass to the make_vector_file function.
        """
        if isinstance(craters, str):
            if name is None:
                name = craters.lower()
            if name == "emplaced":
                craters = self.emplaced
            elif name == "observed":
                craters = self.observed
            else:
                raise ValueError("craters string must be 'emplaced' or 'observed'")

        if isinstance(craters, Crater):
            craters = [craters]

        elif isinstance(craters, dict):
            craters = craters.values()
        elif not isinstance(craters, list):
            raise TypeError("craters must be a Crater, list of Crater, dict of Crater, or 'emplaced'/'observed' string")
        if name is not None:
            kwargs["name"] = name

        if driver.upper() in ["VTK", "VTP"]:
            self.to_vtk_file(
                craters=craters,
                interval_number=interval_number,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )
        elif driver.upper() == "CSV":
            self.to_csv_file(
                craters=craters,
                interval_number=interval_number,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )
        elif driver.upper() == "SCC":
            self.to_scc_file(
                craters=craters,
                interval_number=interval_number,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )
        else:
            self.to_vector_file(
                craters=craters,
                interval_number=interval_number,
                driver=driver,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )

        return

    def show_pyvista(
        self,
        surface: Surface | LocalSurface | None = None,
        observed_color: str = "white",
        emplaced_color: str = "red",
        interval_number: int | None = None,
        **kwargs: Any,
    ):
        """
        Passes through to the surface show_pyvista method and adds crater counts to it.

        Parameters
        ----------
        surface : Surface | LocalSurface, optional
            The surface or local surface view to be displayed. If None, uses the associated surface property
        observed_color : str, optional
            The color to use for observed craters. Default is "white".
        emplaced_color : str, optional
            The color to use for emplaced craters. Default is "red".
        interval_number : int, optional
            The interval number to load the emplaced crater data from. if None, then all emplaced data currently saved to file is used. Default is None.
        **kwargs : Any
            Additional keyword arguments to pass to the surface show_pyvista method.

        Returns
        -------
        plotter : pyvista.Plotter
            The pyvista plotter with the crater counts added.

        """
        if surface is None:
            surface = self.surface
        plotter = surface.show_pyvista(**kwargs)

        from cratermaker.utils.general_utils import toggle_pyvista_actor, update_pyvista_help_message

        emplaced_interval_numbers, emplaced_data_file_list = get_saved_interval_numbers(
            output_dir=self.output_dir,
            output_file_prefix=f"emplaced_{self._output_file_prefix}",
            output_file_extension=self._output_file_extension,
        )
        observed_interval_numbers, observed_data_file_list = get_saved_interval_numbers(
            output_dir=self.output_dir,
            output_file_prefix=f"observed_{self._output_file_prefix}",
            output_file_extension=self._output_file_extension,
        )
        if interval_number is not None:
            if interval_number in emplaced_interval_numbers:
                file_index = emplaced_interval_numbers.index(interval_number)
                emplaced_data_file_list = [emplaced_data_file_list[file_index]]
            else:
                interval_number = None
            if interval_number in observed_interval_numbers:
                file_index = observed_interval_numbers.index(interval_number)
                observed_data_file_list = [observed_data_file_list[file_index]]
            else:
                interval_number = None

        if emplaced_data_file_list:
            emplaced = xr.open_mfdataset(emplaced_data_file_list, combine="nested", parallel=True, engine="h5netcdf")
            emplaced = xr.merge([self.to_xarray(self.emplaced), emplaced])
        else:
            emplaced = self.to_xarray(self.emplaced)

        if observed_data_file_list:
            observed = xr.open_dataset(observed_data_file_list[-1])
        else:
            observed = self.to_xarray(self.observed)
        if emplaced:
            emplaced = self.from_xarray(emplaced)
            emplaced_count_actor = plotter.add_mesh(
                self.to_vtk_mesh(emplaced, use_measured_properties=False),
                line_width=2,
                color=emplaced_color,
                name="emplaced",
            )
            emplaced_count_actor.SetVisibility(False)
            plotter.add_key_event("t", lambda: toggle_pyvista_actor(plotter, emplaced_count_actor))
            plotter = update_pyvista_help_message(plotter, new_message="t: Toggle emplaced craters")
        if observed:
            observed = self.from_xarray(observed)
            observed_count_actor = plotter.add_mesh(
                self.to_vtk_mesh(observed, use_measured_properties=True), line_width=2, color=observed_color, name="observed"
            )
            observed_count_actor.SetVisibility(False)
            plotter.add_key_event("c", lambda: toggle_pyvista_actor(plotter, observed_count_actor))
            plotter = update_pyvista_help_message(plotter, new_message="c: Toggle counted craters")

        return plotter

    def show(self, engine: str = "pyvista", observed_color: str = "white", emplaced_color: str = "red", **kwargs: Any) -> None:
        """
        Passes through to the surface show method and adds crater counts to it.

        Parameters
        ----------
        engine : str, optional
            The engine to use for plotting. Currently, only "pyvista" is supported. Default is "pyvista".
        observed_color : str, optional
            The color to use for observed craters. Default is "white".
        emplaced_color : str, optional
            The color to use for emplaced craters. Default is "red".
        **kwargs : Any
            Additional keyword arguments to pass to the surface show_pyvista method.

        """
        if engine.lower() == "pyvista":
            plotter = self.show_pyvista(observed_color=observed_color, emplaced_color=emplaced_color, **kwargs)
            plotter.show()
        else:
            raise ValueError(f"Engine '{engine}' is not supported for crater counting visualization.")
        return

    @staticmethod
    def _overwrite_check(output_file: Path) -> bool:
        if output_file.exists():
            response = input(
                f"File '{output_file}' already exists. To disable this message, pass `ask_overwrite=False` to this function. Overwrite? (y/n): "
            )
            if response.lower() != "y":
                print("Operation cancelled by user.")
                return False
        output_file.unlink(missing_ok=True)
        return True

    def to_vector_file(
        self,
        craters: list[Crater],
        driver: str = "GPKG",
        interval_number: int = 0,
        name: str | None = None,
        use_measured_properties: bool = True,
        ask_overwrite: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a vector file and stores it in the default export directory.

        Parameters
        ----------
        craters : list[Crater]
            A list of Crater objects.
        driver : str, optional
            The file format to save. Supported formats are 'GPKG', 'ESRI Shapefile', etc. Default is 'GPKG'.
        interval_number : int, optional
            The interval number to append to the file name. Default is 0.
        name : str, optional
            The name of the layer in the GeoPackage file or the file name if the format does not support layers, by default "craters".
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        from cratermaker.constants import EXPORT_DRIVER_TO_EXTENSION_MAP

        def shp_key_fix(key: str) -> str:
            """
            ESRI Shapefile format limits field names to 10 characters, so this function substitues longer names with shorter alternatives, truncates the results, and sets them to upper case.
            """
            alt_names = {
                "projectile_": "proj",
                "morphology_": "morph",
                "diameter": "diam",
                "longitude": "lon",
                "latitude": "lat",
                "density": "dens",
                "velocity": "vel",
                "direction": "dir",
                "location": "loc",
                "angle": "ang",
                "transient_": "tr",
                "semimajor_axis": "a",
                "semiminor_axis": "b",
                "orientation": "orient",
                "measured_": "meas",
                "degradation_state": "kdeg",
            }
            for long, short in alt_names.items():
                if long in key:
                    key = key.replace(long, short)
            return key[:10].upper()

        if name is None:
            if hasattr(self.surface, "local"):
                name = "local_surface"
            else:
                name = "surface"
        # Common alias for Shapefile
        if driver.upper() == "SHP":
            driver = "ESRI Shapefile"
        if driver in EXPORT_DRIVER_TO_EXTENSION_MAP:
            file_extension = EXPORT_DRIVER_TO_EXTENSION_MAP[driver]
        else:
            raise ValueError("Cannot infer file extension from driver {driver}.")

        if file_extension == "shp":
            format_has_layers = False
        else:
            format_has_layers = True

        surface = self.surface
        split_antimeridian = True

        geoms = []
        attrs = []
        crater_ds = self.to_xarray(craters)
        for crater in craters:
            poly = crater.to_geoseries(
                surface=surface, split_antimeridian=split_antimeridian, use_measured_properties=use_measured_properties
            ).item()
            df = crater_ds.sel(id=[crater.id]).to_dataframe()
            if isinstance(poly, GeometryCollection):
                for p in poly.geoms:
                    geoms.append(p)
                    attrs.append(df)
            else:
                geoms.append(poly)
                attrs.append(df)

        if len(geoms) > 0:
            attrs_df = pd.concat(attrs, ignore_index=True)
        else:
            attrs_df = pd.DataFrame()
        if driver.upper() == "ESRI SHAPEFILE":
            attrs_df.rename(mapper=shp_key_fix, axis=1, inplace=True)

        gdf = gpd.GeoDataFrame(data=attrs_df, geometry=geoms, crs=surface.crs)
        if format_has_layers:
            output_file = self.export_dir / f"craters{interval_number:06d}.{file_extension}"
            print(f"Saving {name} layer to vector file: '{output_file}'...")
        else:
            output_file = self.export_dir / f"{name}{interval_number:06d}.{file_extension}"
            if ask_overwrite and not self._overwrite_check(output_file):
                return
        if driver.upper() == "ESRI SHAPEFILE":
            # Append _CRATER so that it is recognized by Craterstats
            output_file = Path(str(output_file).replace(".shp", "_CRATER.shp"))
            if hasattr(self.surface, "local"):
                # Create the _AREA file
                self.surface.local.export_region_polygon(driver=driver)

        try:
            if format_has_layers:
                gdf.to_file(output_file, layer=name)
            else:
                gdf.to_file(output_file)
        except Exception as e:
            raise RuntimeError(f"Error saving {output_file}: {e}") from e

        return

    def to_vtk_mesh(self, craters: list[Crater], use_measured_properties: bool = True, **kwargs: Any) -> vtkPolyData:
        """
        Convert the crater data to a VTK PolyData mesh.

        Parameters
        ----------
        craters : list[Crater]
            A list of Crater objects to convert.
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        **kwargs : Any
            Additional keyword arguments that are passed to the crater to_geoseries method.

        Returns
        -------
        vtkPolyData
            A VTK PolyData object representing the crater geometries.
        """
        from vtk import (
            vtkCellArray,
            vtkPoints,
            vtkPolyLine,
            vtkXMLPolyDataWriter,
        )

        def lonlat_to_xyz(R):
            def _f(lon_deg, lat_deg, z=0.0):
                lon = np.deg2rad(lon_deg)
                lat = np.deg2rad(lat_deg)
                X = (R + z) * np.cos(lat) * np.cos(lon)
                Y = (R + z) * np.cos(lat) * np.sin(lon)
                Z = (R + z) * np.sin(lat)
                return X, Y, Z

            return _f

        def polygon_xyz_coords(geom, R):
            """Yield Nx3 arrays for each exterior ring (drop closing vertex)."""
            g3d = transform(lonlat_to_xyz(R), geom)

            def _rings(g):
                if g.geom_type == "Polygon":
                    yield np.asarray(g.exterior.coords)[:-1]  # (N, 3)
                    for i in g.interiors:
                        yield np.asarray(i.coords)[:-1]
                elif g.geom_type in ("MultiPolygon", "GeometryCollection"):
                    for sub in g.geoms:
                        yield from _rings(sub)
                else:
                    yield np.asarray(g.coords)

            yield from _rings(g3d)

        surface = self.surface

        geoms = []
        for crater in craters:
            geoms.append(
                crater.to_geoseries(
                    surface=surface, split_antimeridian=False, use_measured_properties=use_measured_properties, **kwargs
                )
            )

        points = vtkPoints()
        lines = vtkCellArray()
        point_id = 0  # Keep track of the point ID across all circles
        for g in geoms:
            for ring_xyz in polygon_xyz_coords(g.item(), surface.radius):
                x, y, z = ring_xyz.T  # each is (N,)
                for i in range(len(x)):
                    points.InsertNextPoint(float(x[i]), float(y[i]), float(z[i]))
                polyline = vtkPolyLine()
                polyline.GetPointIds().SetNumberOfIds(len(x))
                for i in range(len(x)):
                    polyline.GetPointIds().SetId(i, point_id + i)
                point_id += len(x)
                lines.InsertNextCell(polyline)

        # Create a poly_data object and add points and lines to it
        poly_data = vtkPolyData()
        poly_data.SetPoints(points)
        poly_data.SetLines(lines)
        return poly_data

    def to_vtk_file(
        self,
        craters: list[Crater],
        interval_number: int = 0,
        name: str = "craters",
        ask_overwrite: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a VTK file and stores it in the default export directory.

        Notes: In order for the crater and surface to be synced up when saving to VTK/VTP format, the initial conditions (no craters) must be saved. Otherwise, saving to file only occurs if there are craters to save.

        Parameters
        ----------
        craters : list[Crater]
            A list of Crater objects.
        interval_number : int, optional
            The interval number to export. Default is 0.
        name : str, optional
            The name used for the file name, by default "craters".
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        from vtk import vtkXMLPolyDataWriter

        output_file = self.export_dir / f"{name}{interval_number:06d}.vtp"
        if ask_overwrite and not self._overwrite_check(output_file):
            return
        print(f"Saving crater data to VTK file: '{output_file}'...")
        poly_data = self.to_vtk_mesh(craters=craters)

        # Write the poly_data to a VTK file
        writer = vtkXMLPolyDataWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(poly_data)

        # Optional: set the data mode to binary to save disk space
        writer.SetDataModeToBinary()
        writer.Write()
        return

    def to_csv_file(
        self,
        craters: list[Crater],
        name: str = "craters",
        interval_number: int = 0,
        ask_overwrite: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a CSV file and stores it in the default export directory.

        Parameters
        ----------
        craters : list[Crater]
            List of Crater objects to export.
        name : str, optional
            The name used for the file name. Default is "craters".
        interval_number : int, optional
            The interval number to export. Default is 0.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        import csv

        if craters:
            output_file = self.export_dir / f"{name}{interval_number:06d}.csv"
            if ask_overwrite and not self._overwrite_check(output_file):
                return
            print(f"Saving crater data to CSV file: '{output_file}'...")
            with output_file.open(mode="w", newline="") as csvfile:
                writer = csv.writer(csvfile)
                header_written = False
                for crater in craters:
                    crater_dict = asdict(crater)

                    # Convert location fields from tuples into lon/lat
                    location = crater_dict.pop("location")
                    crater_dict["longitude"] = location[0]
                    crater_dict["latitude"] = location[1]
                    measured_location = crater_dict.pop("measured_location")
                    crater_dict["measured_longitude"] = measured_location[0]
                    crater_dict["measured_latitude"] = measured_location[1]

                    # Check if the crater is circular, and if so convert semimajor/semiminor to diameter
                    if crater_dict["measured_semimajor_axis"] == crater_dict["measured_semiminor_axis"]:
                        measured_diameter = 2.0 * crater_dict.pop("measured_semimajor_axis")
                        crater_dict.pop("measured_semiminor_axis")
                        items = list(crater_dict.items())
                        crater_dict = {"id": crater_dict["id"], "measured_diameter": measured_diameter}
                        crater_dict.update(items[1:])

                    if crater_dict["semimajor_axis"] == crater_dict["semiminor_axis"]:
                        diameter = 2.0 * crater_dict.pop("semimajor_axis")
                        crater_dict.pop("semiminor_axis")
                        items = list(crater_dict.items())
                        crater_dict = {"id": crater_dict["id"], "diameter": diameter}
                        crater_dict.update(items[1:])

                    # Pop out any None values
                    crater_dict = {k: v for k, v in crater_dict.items() if v is not None}

                    if not header_written:
                        header = list(crater_dict.keys())
                        writer.writerow(header)
                        header_written = True
                    row = [crater_dict[key] for key in header]
                    writer.writerow(row)

        return

    def from_csv_file(self, input_file: Path | str) -> list[Crater]:
        """
        Import crater data from a CSV file.

        Parameters
        ----------
        input_file : Path | str
            The path to the CSV file containing crater data.

        Returns
        -------
        list[Crater]
            A list of Crater objects imported from the CSV file.
        """
        import csv

        craters = []
        input_file = Path(input_file)
        if not input_file.exists():
            raise FileNotFoundError(f"Input file '{input_file}' does not exist.")
        with input_file.open(mode="r", newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                crater_data = {}
                crater_data["location"] = [None, None]
                crater_data["measured_location"] = [None, None]
                for key, value in row.items():
                    if key in ["id"]:
                        crater_data[key] = int(value)
                    elif key == "longitude":
                        crater_data["location"][0] = float(value)
                    elif key == "latitude":
                        crater_data["location"][1] = float(value)
                    elif key == "measured_longitude":
                        crater_data["measured_location"][0] = float(value)
                    elif key == "measured_latitude":
                        crater_data["measured_location"][1] = float(value)
                    else:
                        try:
                            crater_data[key] = float(value)
                        except ValueError:
                            crater_data[key] = value
                if None in crater_data["location"]:
                    crater_data.pop("location")
                if None in crater_data["measured_location"]:
                    crater_data.pop("measured_location")
                crater = Crater.maker(**crater_data, check_redundant_inputs=False)
                craters.append(crater)

        return craters

    def to_scc_file(
        self,
        craters: list[Crater],
        interval_number: int = 0,
        name: str = "craters",
        ask_overwrite: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to Craterstats and OpenCraterTools-compatible SCC file and stores it in the default export directory.

        Parameters
        ----------
        craters : list[Crater]
            A list of Crater objects.
        interval_number : int, optional
            The interval number to append to the file name. Default is 0.
        name : str, optional
            The name used for the file name, by default "craters".
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        import datetime

        output_file = self.export_dir / f"{name}{interval_number:06d}.scc"
        if ask_overwrite and not self._overwrite_check(output_file):
            return
        print(f"\nSaving crater data to {output_file}")
        with output_file.open(mode="w") as f:
            f.write(f"# Spatial crater count Cratermaker version {cratermaker_version}\n")
            f.write("#\n")
            f.write(f"# Exported on {datetime.datetime.now().isoformat()}\n")
            f.write("#\n")
            f.write("# Ellipsoid axes\n")
            f.write(f"a-axis radius = {self.surface.radius * 1e-3:.3f} <km>\n")
            f.write(f"b-axis radius = {self.surface.radius * 1e-3:.3f} <km>\n")
            f.write(f"c-axis radius = {self.surface.radius * 1e-3:.3f} <km>\n")

            # Start with regional area
            boundary_points = []

            if hasattr(self.surface, "local_radius") and hasattr(self.surface, "local_location"):
                f.write(f"coordinate_system_name = {self.surface.local.crs.name}\n")
                region_circle = Crater.maker(radius=self.surface.local_radius, location=self.surface.local_location)
                region_poly = region_circle.to_geoseries(
                    surface=self.surface, split_antimeridian=False, use_measured_properties=False
                ).item()
                boundary_points = list(region_poly.exterior.coords)
                area = self.surface.local.area
            else:
                f.write(f"coordinate_system_name = {self.surface.crs.name}\n")
                boundary_points = [(-180.0, -90.0), (180.0, -90.0), (180.0, 90.0), (-180.0, 90.0), (-180.0, -90.0)]
                area = self.surface.area
            f.write("# area_shapes:\n")
            f.write("unit_boundary = {vertex, sub_area, tag, lon, lat\n")
            for i, p in enumerate(boundary_points):
                f.write(f"{i}\t1\text\t{p[0]}\t{p[1]}\n")
            f.write("}\n")
            f.write("#\n")
            f.write("# area_info:\n")
            f.write(f"Total_area = {area * 1e-6} <km^2>\n")
            f.write("#\n")
            f.write("# crater_diameters\n")
            f.write("crater = {diam, fraction, lon, lat, topo_scale_factor\n")
            for crater in craters:
                f.write(f"{crater.measured_diameter * 1e-3}\t1\t{crater.measured_location[0]}\t{crater.measured_location[1]}\t 1\n")
            f.write("}\n")

        return

    def from_scc_file(self, input_file: Path | str) -> list[Crater]:
        """
        Import crater data from a Spatial Crater Count file.

        Parameters
        ----------
        input_file : Path | str
            The path to the SCC file containing crater data.

        Returns
        -------
        list[Crater]
            A list of Crater objects imported from the SCC file.
        """
        from craterstats import Spatialcount

        craters = []
        input_file = Path(input_file)
        if not input_file.exists():
            raise FileNotFoundError(f"Input file '{input_file}' does not exist.")
        if input_file.suffix != ".scc":
            raise ValueError(f"Input file '{input_file}' is not a .scc file.")
        scc = Spatialcount(filename=str(input_file))
        for diam, lon, lat in zip(scc.diam, scc.lon, scc.lat, strict=True):
            crater = Crater.maker(diameter=diam * 1e3, location=(lon, lat))
            craters.append(crater)

        return craters

    @staticmethod
    def from_xarray(dataset: xr.Dataset) -> list[Crater]:
        """
        Import crater data from an xarray Dataset.

        Parameters
        ----------
        dataset : xr.Dataset
            The xarray Dataset containing crater data.

        Returns
        -------
        list[Crater]
            A list of Crater objects imported from the xarray Dataset.
        """
        craters = []
        for id in dataset.id.data:
            crater_data = dataset.sel(id=id).to_dict()["data_vars"]
            crater_data = {k: v["data"] for k, v in crater_data.items()}
            if "longitude" in crater_data and "latitude" in crater_data:
                crater_data["location"] = (crater_data.pop("longitude"), crater_data.pop("latitude"))
            if "measured_longitude" in crater_data and "measured_latitude" in crater_data:
                crater_data["measured_location"] = (
                    crater_data.pop("measured_longitude"),
                    crater_data.pop("measured_latitude"),
                )
            for k, v in crater_data.items():
                if v is not None and np.any(np.isreal(v)) and np.any(np.isnan(v)):
                    crater_data[k] = None
            crater = Crater.maker(**crater_data, check_redundant_inputs=False)
            craters.append(crater)

        return craters

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
        Dictionary of observed craters on the surface keyed to the crater id.
        """
        return self._observed

    @property
    def emplaced(self) -> list[Crater]:
        """
        List of craters that have been emplaced in the simulation in the current interval in chronological order.
        """
        return self._emplaced


import_components(__name__, __path__)
