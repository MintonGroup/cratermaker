from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
import xarray as xr
from cratermaker._cratermaker import counting_bindings
from geopandas import GeoSeries
from matplotlib.axes import Axes
from numpy.typing import ArrayLike
from shapely.ops import transform
from tqdm import tqdm
from vtk import vtkPolyData

from cratermaker import __version__ as cratermaker_version
from cratermaker.components.crater import Crater
from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP, FloatLike
from cratermaker.core.base import ComponentBase, import_components

DRIVER_TO_EXTENSION_MAP = {"NETCDF": "nc", "SCC": "scc", "VTK": "vtp", "VTP": "vtp", "CSV": "csv", **VECTOR_DRIVER_TO_EXTENSION_MAP}


if TYPE_CHECKING:
    from cratermaker.components.morphology import Morphology, MorphologyCrater
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_LONG_NAME = "Unique crater identification number"

# The number of layers used for tagging faces with crater ids. This allows a single face to contain multiple crater ids
_N_LAYER = 8

# The minimum number of faces required in a region to perform crater counting, which corresponds to a roughly 6-8 pix diameter crater
_MIN_FACE_FOR_COUNTING = 100


class Counting(ComponentBase):
    """
    Base class for all crater counting models. It defines the interface for tallying the observable craters on a surface.
    """

    _registry: dict[str, Counting] = {}

    def __init__(
        self,
        surface: Surface | LocalSurface,
        Crater: type[Crater] | None = None,
        reset: bool = True,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        surface : Surface | LocalSurface
            The surface or local surface view to be counted.
        Crater : type[Crater], optional
            The Crater class associated with this counting model. This is used to ensure that the correct variable properties for from a specialized Crater class are available (such as one associated with a Morphology class) when importing craters from file. If not supplied, then the base Crater class is used.
        reset : bool, optional
            Flag to indicate whether to reset the count and delete any old output files. Default is True.
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.components.surface import Surface

        super().__init__(reset=reset, **kwargs)

        object.__setattr__(self, "_emplaced", [])
        object.__setattr__(self, "_observed", {})
        object.__setattr__(self, "_output_dir_name", "counting")
        object.__setattr__(self, "_output_file_prefix", "craters")
        object.__setattr__(self, "_output_file_extension", "nc")
        object.__setattr__(self, "_Crater", None)
        object.__setattr__(self, "_surface", None)
        object.__setattr__(self, "_morphology", None)
        object.__setattr__(
            self,
            "_driver_to_extension_map",
            {"SCC": "scc", "VTK": "vtk", "VTP": "vtp", "CSV": "csv", **VECTOR_DRIVER_TO_EXTENSION_MAP},
        )
        self._surface = Surface.maker(surface, reset=reset, **kwargs)
        self.Crater = Crater
        self._output_file_pattern += [
            f"observed_{self._output_file_prefix}*.{self._output_file_extension}",
            f"emplaced_{self._output_file_prefix}*.{self._output_file_extension}",
        ]

        if not reset:
            observed, emplaced = self.read_saved_output(interval=-1)
            if observed:
                interval = observed.interval.values[-1]
                observed = self.from_xarray(observed, interval=interval)
                for crater in observed:
                    self._observed[crater.id] = crater
                if emplaced:
                    self._emplaced = self.from_xarray(emplaced, interval=interval)
        return

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
        Crater: type[Crater] | None = None,
        reset: bool = True,
        **kwargs: Any,
    ) -> Counting:
        """
        Initialize a Counting model with the given name or instance.

        Parameters
        ----------
        counting : str, Counting or None, default=None
            The name of the counting model to initialize. If None, the default model is used.
        surface : Surface | LocalSurface
            The surface or local surface view to be counted.
        Crater : type[Crater], optional
            The Crater class associated with this counting model. This is used to ensure that the correct variable properties for from a specialized Crater class are available (such as one associated with a Morphology class) when importing craters from file. If not supplied, then the base Crater class is used.
        reset : bool, optional
            Flag to indicate whether to reset the count and delete any old output files. Default is True
        **kwargs : Any
            |kwargs|

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
            counting = "simplecount"

        counting = super().maker(
            component=counting,
            surface=surface,
            reset=reset,
            Crater=Crater,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        base = super().__str__()
        str_repr = f"{base}\n"
        str_repr += f"Number of observed craters: {self.n_observed}\n"
        str_repr += f"Number of emplaced craters: {self.n_emplaced}\n"
        str_repr += f"\n{self.surface}\n"
        return str_repr

    def reset(self, **kwargs: Any) -> None:
        """
        Remove all craters count records from the surface and delete any output files.

        Parameters
        ----------
        **kwargs : Any
            |kwargs|
        """
        self.surface.add_tag(name="crater_id", long_name=_TALLY_LONG_NAME, tag=None, n_layer=self.n_layer)
        self._emplaced = []
        self._observed = {}

        super().reset(**kwargs)
        return

    def add(self, crater: MorphologyCrater, **kwargs: Any):
        """
        Add a crater to the surface.

        Parameters
        ----------
        crater : Crater
            The crater to be added to the surface.
        **kwargs : Any
            |kwargs|
        """
        from cratermaker.components.morphology import MorphologyCrater

        if not isinstance(crater, MorphologyCrater):
            raise TypeError("crater must be an instance of MorphologyCrater")

        if self.surface.uxds is None:
            raise ValueError(
                "Surface must have an associated Uxarray dataset to use for counting. This is commonly caused by using a HiResLocal surface type without setting the superdomain_scale_factor."
            )

        # Tag a region just outside crater rim with the id
        crater_region = crater.crater_region

        if crater_region and crater_region.n_face >= _MIN_FACE_FOR_COUNTING:
            crater_region.add_tag(
                name="crater_id",
                long_name=_TALLY_LONG_NAME,
                tag=crater.id,
                **kwargs,
            )

            # Check to make sure the crater was recorded to the surface. HiResLocal surfaces may not record craters in the superdomain
            if crater.id in self.surface.uxds.crater_id:
                self.emplaced.append(crater)
                self.observed[crater.id] = crater
            else:
                return

            # Cookie cutting: remove any smaller craters that are overlapped by this new crater
            unique_ids = np.unique(crater_region.crater_id)
            unique_ids = unique_ids[unique_ids > 0]  # Remove the 0 id which corresponds to no crater
            if len(unique_ids) > 0:
                # Compute cookie cutting removes list
                observed = self.observed.copy()
                removes = [id for id, v in observed.items() if v.id in unique_ids and v.diameter < crater.diameter]
                # For every id that appears in the removes list, set it to 0 in the data array
                for remove_id in removes:
                    crater_region.remove_tag(name="crater_id", tag=remove_id)
                    if not np.any(
                        self.surface.uxds.crater_id.data == remove_id
                    ):  # Check to see if this crater id still appears, and if not, it's gone man.
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

    def fit_rim(
        self, crater: Crater, tol=0.01, nloops=4, score_quantile=0.95, fit_center=False, fit_ellipse=False, **kwargs
    ) -> Crater:
        """
        Find the rim region of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to find the rim region.
        tol : float, optional
            The tolerance for the rim fitting algorithm. Default is 0.01.
        nloops : int, optional
            The number of iterations for the rim fitting algorithm. Default is 4.
        score_quantile : float, optional
            The quantile of rim scores to consider. Default is 0.95.
        fit_center : bool, optional
            If True, fit the crater center as well. Default is False.
        fit_ellipse : bool, optional
            If True, fit an ellipse to the rim, otherwise fit a circle. Default is False.
        **kwargs : Any
            |kwargs|

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

        crater.measured_semimajor_axis = ap
        crater.measured_semiminor_axis = bp
        crater.measured_orientation = np.degrees(orientation)
        crater.measured_location = location

        return crater

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

        region = crater.crater_region
        region.compute_desloped_face_elevation()
        region = counting_bindings.score_rim(region, crater, quantile, gradmult, curvmult, heightmult)
        return region

    def tally(self, region: LocalSurface | None = None, measure_rim: bool = False, quiet: bool = False, **kwargs: Any) -> None:
        """
        Tally the craters on the surface using the method of Minton et al. (2019) [#]_.

        Parameters
        ----------
        region : LocalSurface, optional
            A LocalSurface region to count. If not supplied, then the associated surface property is used.
        measure_rim : bool, optional
            If True, measure the rim of each crater. Default is False.
        quiet : bool, optional
            If True, suppress progress output. Default is False.
        **kwargs : Any
            |kwargs|

        References
        ----------
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

        """
        if region is None:
            region = self.surface
            if hasattr(self.surface, "crater_id"):
                id_array = self.surface.crater_id
            else:
                return
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
            # TODO: Make the fit_rim function more reliable before turning it on by default. It is currently very slow and not reliable enough.
            if measure_rim:
                fit_center = kwargs.pop("fit_center", False)
                fit_ellipse = kwargs.pop("fit_ellipse", False)
                crater = self.fit_rim(crater=crater, fit_center=fit_center, fit_ellipse=fit_ellipse, **kwargs)
            Kd = self.measure_degradation_state(crater, **kwargs)
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

        self.remove_complex_data()

        return

    @abstractmethod
    def measure_degradation_state(self, crater: Crater, **kwargs: Any) -> float: ...

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
            d = c.as_dict(skip_complex_data=True)
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

    def _to_file(
        self,
        craters: dict[int, Crater] | list[Crater],
        filename: Path | str,
        interval: int = 0,
    ) -> xr.Dataset | None:
        """
        Merge a list or dictionary of Crater objects with an existing file.

        Parameters
        ----------
        craters : dict[int, Crater] | list[Crater]
            A dictionary or list of Crater objects to merge.
        filename : Path | str
            The path to the file to merge with.
        interval : int, optional
            The interval number. This is added to the coordinates of the dataset created from the craters before merging with the file. Default is 0.

        Returns
        -------
        xr.Dataset | None
            An xarray Dataset containing the merged crater data, or None if no data.
        """
        # Convert into an xarray dataset
        combined_data = self.to_xarray(craters)
        combined_data = combined_data.expand_dims(dim="interval").assign_coords({"interval": [interval]})

        # If the file already exists, read it and merge
        if filename.exists():
            filename.unlink()

        # Write merged data back to file
        if combined_data:
            combined_data.to_netcdf(filename)
        return combined_data

    def save(
        self,
        crater_type: Literal["observed", "emplaced", "both"] = "both",
        interval: int = 0,
        craters: list[Crater] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Dump the crater lists to a file and reset the emplaced crater list.

        Parameters
        ----------
        crater_type : str, optional
            The type of craters to export. Options are "observed", "emplaced", and "both". Default is "observed".
        interval : int, default=0
            The interval number for the output file naming.
        craters: list[Crater] | None, optional
            An arbitrary list of craters to save. If None, then the current observed and/or emplaced tallies are used. Default is None.
        **kwargs : Any
            |kwargs|
        """
        if craters is None:
            self.remove_complex_data()
            observed = self.observed
            emplaced = self.emplaced
        else:
            if crater_type == "both":
                raise ValueError(
                    "If supplying a custom list of craters to save, then crater_type must be either 'observed' or 'emplaced', not 'both'."
                )
            elif crater_type == "observed":
                observed = {crater.id: crater for crater in craters}
                emplaced = None
            elif crater_type == "emplaced":
                observed = None
                emplaced = craters

        if crater_type == "observed":
            iter = zip([observed], ["observed"], strict=True)
        elif crater_type == "emplaced":
            iter = zip([emplaced], ["emplaced"], strict=True)
        elif crater_type == "both":
            iter = zip([observed, emplaced], ["observed", "emplaced"], strict=True)
        else:
            raise ValueError(f"Invalid crater_type: {crater_type}. Must be one of 'observed', 'emplaced', or 'both'.")
        for craters, name in iter:
            if craters:
                filename = self.output_dir / f"{name}_{self.output_filename(interval)}"
                self._to_file(craters, filename, interval)
        save_args = {"interval": interval, **kwargs}
        super().save(**save_args)
        return

    def from_file(self, filename: str | Path, **kwargs: Any) -> None:
        pass

    def export(
        self,
        crater_type: Literal["observed", "emplaced", "both"] = "observed",
        interval: int | None = None,
        driver: str = "SCC",
        ask_overwrite: bool | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Exports crater lists to a file.

        Parameters
        ----------
        crater_type : str, optional
            The type of craters to export. Options are "observed", "emplaced", and "both". Default is "observed".
        interval : int | None, optional
            |interval_export|
        driver : str, default='SCC'
            The export driver used to save. The valid drivers are given by the `valid_drivers` attribute of this object.
        ask_overwrite : bool, optional
            |ask_overwrite_methods|
        **kwargs : Any
            |kwargs|
        """
        # Check if the requested driver is a supported one for crater counts, and return silently if not.
        if driver.upper() not in self.valid_drivers:
            return

        # Temporarily set the ask_overwrite attribute for the duration of the export, but reset it to its original value afterwards.
        ask_overwrite_orig = self.ask_overwrite
        if ask_overwrite is not None:
            self.ask_overwrite = ask_overwrite

        crater_names = ["observed", "emplaced"]
        output_ds = self.read_saved_output(interval=interval)
        for name, crater_ds in zip(crater_names, output_ds, strict=True):
            if crater_type == "observed" and name != "observed":
                continue
            if crater_type == "emplaced" and name != "emplaced":
                continue

            if type(crater_ds) is dict:
                interval_numbers = list(crater_ds.keys())
            elif crater_ds is None or "interval" not in crater_ds:
                if interval is None:
                    interval_numbers = [0]
                else:
                    interval_numbers = [interval]
            else:
                interval_numbers = crater_ds.interval.values
            for interval in interval_numbers:
                if driver.upper() == "NETCDF" and crater_ds is not None:
                    self.save(crater_type=name, interval=interval, craters=self.from_xarray(crater_ds, interval=interval), **kwargs)
                elif driver.upper() in ["VTK", "VTP"] and crater_ds is not None:
                    self.to_vtk_file(
                        crater_ds=crater_ds,
                        interval=interval,
                        name=name,
                        **kwargs,
                    )
                elif driver.upper() == "CSV" and crater_ds is not None:
                    self.to_csv_file(
                        crater_ds=crater_ds,
                        interval=interval,
                        name=name,
                        **kwargs,
                    )
                elif driver.upper() == "SCC":
                    if crater_ds is None:
                        crater_ds = []  # Allow for empty SCC files because they can still have the region polygon
                    self.to_scc_file(
                        crater_ds=crater_ds,
                        interval=interval,
                        name=name,
                        **kwargs,
                    )
                elif driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP and crater_ds is not None:
                    self.to_vector_file(
                        crater_ds=crater_ds,
                        interval=interval,
                        name=name,
                        driver=driver,
                        **kwargs,
                    )

        self.ask_overwrite = ask_overwrite_orig
        return

    def plot(
        self,
        interval: int | None = None,
        observed_color: str | None = "white",
        observed_original_color: str | None = None,
        emplaced_color: str | None = None,
        plot_style: Literal["map", "hillshade"] = "map",
        variable_name: str | None = None,
        cmap: str | None = None,
        show=True,
        save=True,
        ax: Axes | None = None,
        **kwargs: Any,
    ) -> Axes:
        """
        Plot an image of the local region.

        Parameters
        ----------
        interval : int | None, optional
            The interval number to load the emplaced crater data from. if None, then all emplaced data currently saved to file is used. Default is None.
        observed_color : str | None, optional
            The color to use for observed craters using their measured properties. If None, observed craters will not be plotted. Default is "white".
        observed_original_color : str | None, optional
            The color to use for observed craters using their original properties. If None, observed craters will not be plotted. Default is None.
        emplaced_color : str | None, optional
            The color to use for emplaced craters. If None, emplaced craters will not be plotted. Default is None.
        plot_style : str, optional
            The style of the plot. Options are "map" and "hillshade". In "map" mode, the variable is displayed as a colored map. In "hillshade" mode, a hillshade image is generated using "face_elevation" data. If a different variable is passed to `variable`, then the hillshade will be overlayed with that variable's data. Default is "map".
        variable_name : str | None, optional
            The variable to plot. If None is provided then "face_elevation" is used in "map" mode.
        cmap : str, optional
            The colormap to use for the plot. If None, a default colormap will be used ("cividis" by default and "grey" when plot_style=="hillshade" and variable=="face_elevation").
        show : bool, optional
            If True, the plot will be displayed. Default is True.
        save : bool, optional
            If True, the plot will be saved to the default plot directory. Default is True
        ax : matplotlib.axes.Axes, optional
            An existing Axes object to plot on. If None, a new figure and axes will be created.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        matplotlib.image.Axes object created by the surface plot method, with crater counts plotted on top if specified. If show is True, the plot will also be displayed.
            The Axes object
        """
        import matplotlib.pyplot as plt

        from cratermaker.components.surface.hireslocal import HiResLocalSurface

        crs = self.surface.crs
        split_antimeridian = True
        file_prefix = f"{self.surface.output_file_prefix}"
        # Handle the HiResLocal surface case where we may or may not be plotting the global surface
        if isinstance(self.surface, HiResLocalSurface):
            superdomain = kwargs.pop("superdomain", False)
            if not superdomain and self.surface.local is not None:
                crs = self.surface.local.crs
                split_antimeridian = False
                file_prefix = f"{self.surface.local.output_file_prefix}"

        file_prefix += f"_{self.output_file_prefix}"
        if interval is not None:
            observed, emplaced = self.read_saved_output(interval=interval)
            if observed:
                interval = observed.interval.values[-1]
                observed = self.from_xarray(observed, interval=interval)
            if emplaced:
                emplaced_interval = emplaced.interval.values[-1]
                if emplaced_interval == interval:
                    emplaced = self.from_xarray(emplaced, interval=interval)
            filename = self.plot_dir / f"{file_prefix}{interval:06d}.{self.surface.output_image_file_extension}"
        else:
            observed = [c for _, c in self.observed.items()]
            emplaced = self.emplaced
            filename = self.plot_dir / f"{file_prefix}.{self.surface.output_image_file_extension}"
        if ax is None:
            W, H = self.surface.get_raster_dims()
            _, ax = plt.subplots(figsize=(1, 1), dpi=W, frameon=False)
        ax = self.surface.plot(
            interval=interval,
            show=False,
            save=False,
            ax=ax,
            plot_style=plot_style,
            variable_name=variable_name,
            cmap=cmap,
            **kwargs,
        )

        if emplaced_color is not None and emplaced is not None and len(emplaced) > 0:
            gs = self.to_geoseries(
                craters=emplaced, use_measured_properties=False, split_antimeridian=split_antimeridian, autolim=False
            ).to_crs(crs)
            facecolor = kwargs.pop("facecolor", "none")
            edgecolor = emplaced_color
            linewidth = kwargs.pop("linewidth", 0.1)
            linestyle = kwargs.pop("linestyle", "solid")
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle)

        if observed_original_color is not None and observed is not None and len(observed) > 0:
            gs = self.to_geoseries(
                craters=observed, use_measured_properties=False, split_antimeridian=split_antimeridian, autolim=False
            ).to_crs(crs)
            facecolor = kwargs.pop("facecolor", "none")
            edgecolor = observed_original_color
            linewidth = kwargs.pop("linewidth", 0.1)
            linestyle = kwargs.pop("linestyle", "solid")
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle)

        if observed_color is not None and observed is not None and len(observed) > 0:
            gs = self.to_geoseries(
                craters=observed, use_measured_properties=True, split_antimeridian=split_antimeridian, autolim=False
            ).to_crs(crs)
            facecolor = kwargs.pop("facecolor", "none")
            edgecolor = observed_color
            linewidth = kwargs.pop("linewidth", 0.1)
            linestyle = kwargs.pop("linestyle", "solid")
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle)

        if save:
            print(f"Saving crater count plot to {filename}...")
            plt.savefig(filename, bbox_inches="tight", pad_inches=0, dpi=W)
        if show:
            plt.show()
        else:
            plt.close()
        return ax

    def show_pyvista(
        self,
        surface: Surface | LocalSurface | None = None,
        observed_color: str = "white",
        emplaced_color: str = "red",
        interval: int | None = None,
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
        interval : int, optional
            The interval number to load the emplaced crater data from. if None, then all emplaced data currently saved to file is used. Default is None.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        plotter : pyvista.Plotter
            The pyvista plotter with the crater counts added.

        """
        if surface is None:
            surface = self.surface
        plotter = surface.show_pyvista(**kwargs)

        from cratermaker.utils.general_utils import toggle_pyvista_actor, update_pyvista_help_message

        if interval is None:
            interval = -1

        observed, emplaced = self.read_saved_output(interval=interval)
        if observed:
            interval = observed.interval.values[-1]
            observed = self.from_xarray(observed, interval=interval)
            observed_count_actor = plotter.add_mesh(
                self.to_vtk_mesh(observed, use_measured_properties=True), line_width=2, color=observed_color, name="observed"
            )
            observed_count_actor.SetVisibility(False)
            plotter.add_key_event("c", lambda: toggle_pyvista_actor(plotter, observed_count_actor))
            plotter = update_pyvista_help_message(plotter, new_message="c: Toggle counted craters")
        if emplaced:
            emplaced_interval = emplaced.interval.values[-1]
            if emplaced_interval == interval:
                emplaced = self.from_xarray(emplaced, interval=interval)
                emplaced_count_actor = plotter.add_mesh(
                    self.to_vtk_mesh(emplaced, use_measured_properties=False),
                    line_width=2,
                    color=emplaced_color,
                    name="emplaced",
                )
                emplaced_count_actor.SetVisibility(False)
                plotter.add_key_event("t", lambda: toggle_pyvista_actor(plotter, emplaced_count_actor))
                plotter = update_pyvista_help_message(plotter, new_message="t: Toggle emplaced craters")

        return plotter

    def show3d(self, engine: str = "pyvista", observed_color: str = "white", emplaced_color: str = "red", **kwargs: Any) -> None:
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
            |kwargs|
        """
        if engine.lower() == "pyvista":
            plotter = self.show_pyvista(observed_color=observed_color, emplaced_color=emplaced_color, **kwargs)
            plotter.show()
        else:
            raise ValueError(f"Engine '{engine}' is not supported for crater counting visualization.")
        return

    def to_geoseries(
        self,
        craters: list[Crater] | None = None,
        use_measured_properties: bool = True,
        split_antimeridian: bool = False,
        **kwargs: Any,
    ) -> GeoSeries:
        if craters is None:
            craters = [c for _, c in self.observed.items()]
        surface = self.surface
        crater_gs = []
        for crater in tqdm(
            craters,
            total=len(craters),
            desc="Converting craters to GeoSeries polygons",
            unit="crater",
            position=0,
            leave=False,
        ):
            crater_gs.append(
                crater.to_geoseries(surface=surface, split_antimeridian=split_antimeridian, use_measured_properties=False)
            )
        return pd.concat(crater_gs)

    def _validate_export_args(
        self,
        name: Literal["observed", "emplaced"] = "observed",
        interval: int | None = None,
        crater_ds: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
    ) -> list[Crater]:
        if crater_ds is None:
            crater_list = getattr(self, name)
        if type(crater_ds) is dict:
            if type(list(crater_ds.values())[0]) is xr.Dataset:  # like multi-interval dicts keyed by interval
                crater_ds = crater_ds.get(
                    interval, []
                )  # Get the dataset for the specified interval, or an empty dataset if not found
            else:
                return list(crater_ds.values())  # like observed dicts keyed by crater id

        if type(crater_ds) is list:
            crater_list = crater_ds
        elif type(crater_ds) is dict:
            crater_list = list(crater_ds.values())
        elif type(crater_ds) is xr.Dataset:
            crater_list = self.from_xarray(crater_ds, interval=interval)
        else:
            raise ValueError(f"Unrecognized type for crater_ds: {type(crater_ds)}")

        return crater_list

    def to_vector_file(
        self,
        name: Literal["observed", "emplaced"] = "observed",
        crater_ds: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        driver: str = "GPKG",
        use_measured_properties: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a vector file and stores it in the default export directory.

        Parameters
        ----------
        name : Literal["observed", "emplaced"], optional
            The name of the crater dataset to export, either "observed" or "emplaced
        crater_ds : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the name parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
            |interval_export|
        driver : str, optional
            The file format to save. Supported formats are 'GPKG', 'ESRI Shapefile', etc. Default is 'GPKG'.
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        **kwargs : Any
            |kwargs|
        """
        return
        # TODO Fix this
        # from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP

        # crater_list = self._validate_export_args(name=name, interval=interval, crater_ds=crater_ds)

        # def shp_key_fix(key: str) -> str:
        #     """
        #     ESRI Shapefile format limits field names to 10 characters, so this function substitues longer names with shorter alternatives, truncates the results, and sets them to upper case.
        #     """
        #     alt_names = {
        #         "projectile_": "proj",
        #         "morphology_": "morph",
        #         "diameter": "diam",
        #         "longitude": "lon",
        #         "latitude": "lat",
        #         "density": "dens",
        #         "velocity": "vel",
        #         "direction": "dir",
        #         "location": "loc",
        #         "angle": "ang",
        #         "transient_": "tr",
        #         "semimajor_axis": "a",
        #         "semiminor_axis": "b",
        #         "orientation": "orient",
        #         "measured_": "meas",
        #         "degradation_state": "kdeg",
        #     }
        #     for long, short in alt_names.items():
        #         if long in key:
        #             key = key.replace(long, short)
        #     return key[:10].upper()

        # # Common alias for Shapefile
        # if driver.upper() == "SHP":
        #     driver = "ESRI Shapefile"
        # if driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP:
        #     file_extension = VECTOR_DRIVER_TO_EXTENSION_MAP[driver.upper()]
        # else:
        #     raise ValueError("Cannot infer file extension from driver {driver}.")

        # if file_extension == "shp":
        #     format_has_layers = False
        # else:
        #     format_has_layers = True

        # surface = self.surface
        # split_antimeridian = False

        # geoms = []
        # attrs = []
        # for crater in tqdm(
        #     crater_list,
        #     total=len(crater_list),
        #     desc=f"Converting {name} craters to geometries for export",
        #     unit="craters",
        #     position=0,
        #     leave=False,
        # ):
        #     poly = crater.to_geoseries(
        #         surface=surface, split_antimeridian=split_antimeridian, use_measured_properties=use_measured_properties
        #     ).item()
        #     df = crater_ds.sel(id=[crater.id]).to_dataframe()
        #     if isinstance(poly, GeometryCollection):
        #         for p in poly.geoms:
        #             geoms.append(p)
        #             attrs.append(df)
        #     else:
        #         geoms.append(poly)
        #         attrs.append(df)

        # if len(geoms) > 0:
        #     attrs_df = pd.concat(attrs, ignore_index=True)
        # else:
        #     attrs_df = pd.DataFrame()
        # if driver.upper() == "ESRI SHAPEFILE":
        #     attrs_df.rename(mapper=shp_key_fix, axis=1, inplace=True)

        # gdf = gpd.GeoDataFrame(data=attrs_df, geometry=geoms, crs=surface.crs)
        # if format_has_layers:
        #     output_file = self.export_dir / f"{self.output_file_prefix}{interval:06d}.{file_extension}"
        #     print(f"Saving {name} layer to vector file: '{output_file}'...")
        # else:
        #     output_file = self.export_dir / f"{name}_{self.output_file_prefix}{interval:06d}.{file_extension}"
        #     if not self._overwrite_check(output_file):
        #         return
        # if driver.upper() == "ESRI SHAPEFILE" and hasattr(self.surface, "local"):
        #     # Create the _AREA file
        #     self.surface.local.export_region_polygon(driver=driver)

        # try:
        #     if format_has_layers:
        #         gdf.to_file(output_file, layer=name)
        #     else:
        #         gdf.to_file(output_file)
        # except Exception as e:
        #     raise RuntimeError(f"Error saving {output_file}: {e}") from e

        # return

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
            |kwargs|

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
        name: Literal["observed", "emplaced"] = "observed",
        crater_ds: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a VTK file and stores it in the default export directory.

        Notes: In order for the crater and surface to be synced up when saving to VTK/VTP format, the initial conditions (no craters) must be saved. Otherwise, saving to file only occurs if there are craters to save.

        Parameters
        ----------
        name : Literal["observed", "emplaced"], optional
            The name of the crater dataset to export, either "observed" or "emplaced
        crater_ds : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the name parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
           |interval_export|
        **kwargs : Any
            |kwargs|
        """
        from vtk import vtkXMLPolyDataWriter

        crater_list = self._validate_export_args(name=name, interval=interval, crater_ds=crater_ds)

        filename_base = self.output_filename(interval).replace(self.output_file_extension, "vtp")
        output_file = self.export_dir / f"{name}_{filename_base}"
        if not self._overwrite_check(output_file):
            return
        print(f"Saving crater data to VTK file: '{output_file}'...")
        poly_data = self.to_vtk_mesh(craters=crater_list)

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
        name: Literal["observed", "emplaced"] = "observed",
        crater_ds: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a CSV file and stores it in the default export directory.

        Parameters
        ----------
        name : Literal["observed", "emplaced"], optional
            The name of the crater dataset to export, either "observed" or "emplaced
        crater_ds : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the name parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
            |interval_export|
        **kwargs : Any
            |kwargs|
        """
        import csv

        crater_list = self._validate_export_args(name=name, interval=interval, crater_ds=crater_ds)

        filename_base = self.output_filename(interval).replace(self.output_file_extension, "csv")
        output_file = self.export_dir / f"{name}_{filename_base}"
        if not self._overwrite_check(output_file):
            return
        print(f"Saving crater data to CSV file: '{output_file}'...")
        with output_file.open(mode="w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            header_written = False
            for crater in crater_list:
                crater_dict = crater.as_dict(skip_complex_data=True)

                # Convert location fields from tuples into lon/lat
                location = crater_dict.pop("location")
                crater_dict["longitude"] = location[0]
                crater_dict["latitude"] = location[1]
                measured_location = crater_dict.pop("measured_location")
                crater_dict["measured_longitude"] = measured_location[0]
                crater_dict["measured_latitude"] = measured_location[1]

                # Check if the crater is circular, and if so convert semimajor/semiminor to diameter
                if (
                    "measured_semimajor_axis" in crater_dict
                    and "measured_semiminor_axis" in crater_dict
                    and crater_dict["measured_semimajor_axis"] == crater_dict["measured_semiminor_axis"]
                ):
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
                crater = self.Crater.maker(**crater_data, morphology=self.morphology, check_redundant_inputs=False)
                craters.append(crater)

        return craters

    def to_scc_file(
        self,
        name: Literal["observed", "emplaced"] = "observed",
        crater_ds: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to Craterstats and OpenCraterTools-compatible SCC file and stores it in the default export directory.

        Parameters
        ----------
        name : Literal["observed", "emplaced"], optional
            The name of the crater dataset to export, either "observed" or "emplaced
        crater_ds : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the name parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
            |interval_export|
        **kwargs : Any
            |kwargs|
        """
        import datetime

        crater_list = self._validate_export_args(name=name, interval=interval, crater_ds=crater_ds)

        def overlap_fraction(crater, region_poly=None):
            if region_poly is None:
                return 1.0
            distance = self.surface.compute_distances(reference_location=self.surface.local_location, locations=[crater.location])
            if distance + crater.measured_radius > self.surface.local_radius:
                crater_poly = crater.to_geoseries(
                    surface=self.surface, split_antimeridian=False, use_measured_properties=True
                ).to_crs(self.surface.crs)
                overlap_area = crater_poly.intersection(region_poly).to_crs(self.surface.local.crs).area.item()
                return overlap_area / crater_poly.to_crs(self.surface.local.crs).area.item()
            else:
                return 1.0

        region_poly = None

        output_file = self.export_dir / f"{name}{interval:06d}.scc"
        if not self._overwrite_check(output_file):
            return
        print(f"Saving crater data to {output_file}")
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
                if region_poly is None:  # We only need to do this the first time through
                    region_circle = Crater.maker(
                        radius=self.surface.local_radius, location=self.surface.local_location
                    )  # We can get away with using just the base class for Crater here
                    region_poly = (
                        region_circle.to_geoseries(surface=self.surface, split_antimeridian=False, use_measured_properties=False)
                        .to_crs(self.surface.crs)
                        .item()
                    )
                    boundary_points = list(region_poly.exterior.coords)
                    area = self.surface.local.area
            else:
                f.write(f"coordinate_system_name = {self.surface.crs.name}\n")
                boundary_points = [(-180.0, -90.0), (180.0, -90.0), (180.0, 90.0), (-180.0, 90.0), (-180.0, -90.0)]
                area = self.surface.area
                region_poly = None
            f.write("# area_shapes:\n")
            f.write("unit_boundary = {vertex, sub_area, tag, lon, lat\n")
            for i, p in enumerate(boundary_points):
                f.write(f"{i}\t1\text\t{p[0]}\t{p[1]}\n")
            f.write("}\n")
            f.write("#\n")
            f.write("# area_info:\n")
            f.write(f"Total_area = {area * 1e-6} <km²>\n")
            f.write("#\n")
            f.write("# crater_diameters\n")
            f.write("crater = {diam, fraction, lon, lat, topo_scale_factor\n")
            for crater in crater_list:
                f.write(
                    f"{crater.measured_diameter * 1e-3}\t{overlap_fraction(crater, region_poly)}\t{crater.measured_location[0]}\t{crater.measured_location[1]}\t 1\n"
                )
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
            crater = self.Crater.maker(diameter=diam * 1e3, location=(lon, lat), morphology=self.morphology)
            craters.append(crater)

        return craters

    def from_xarray(self, dataset: xr.Dataset | dict, interval: int | None = None) -> list[Crater]:
        """
        Import crater data from an xarray Dataset.

        Parameters
        ----------
        dataset : xr.Dataset | dict
            The xarray Dataset containing crater data or a dictionary of xarray Datasets keyed by interval number.

        Returns
        -------
        list[Crater]
            A list of Crater objects imported from the xarray Dataset.
        """
        craters = []
        if type(dataset) is dict:
            if interval is None:
                dataset = dataset[-1]
            elif interval in dataset:
                dataset = dataset[interval]
            else:
                return craters
        if "interval" in dataset.coords:
            if interval is None:
                dataset = dataset.isel(interval=-1)
            elif interval in dataset.interval:
                dataset = dataset.sel(interval=interval)
            else:
                return craters
        dataset.load()
        if len(dataset) == 0:
            return craters
        for id in tqdm(dataset.id.data, desc="Converting xarray Dataset to Crater objects", unit="crater", position=0, leave=False):
            crater_data = dataset.sel(id=id).to_dict()["data_vars"]
            crater_data = {k: v["data"] for k, v in crater_data.items()}
            if np.isnan(crater_data["semimajor_axis"]):
                continue
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
            crater = self.Crater.maker(**crater_data, morphology=self.morphology, check_redundant_inputs=False)
            craters.append(crater)

        return craters

    def remove_complex_data(self):
        """
        Remove complex data from all observed and emplaced craters to free up memory. This is typically called after the tally step to clear out things like the affect node and face index sets.
        """
        for v in self.observed.values():
            v.remove_complex_data()
        for v in self.emplaced:
            v.remove_complex_data()
        return

    @property
    def surface(self):
        """The Surface object associated with this Counting object."""
        return self._surface

    @surface.setter
    def surface(self, value):
        from cratermaker.components.surface import LocalSurface, Surface

        if not isinstance(value, (Surface | LocalSurface)):
            raise TypeError("surface must be an instance of Surface or LocalSurface")
        self._surface = value

    @property
    def n_layer(self) -> int:
        """Number of layers in the counting model."""
        return _N_LAYER

    @property
    def observed(self) -> dict[int, Crater]:
        """Dictionary of observed craters on the surface keyed to the crater id."""
        return self._observed

    @property
    def emplaced(self) -> list[Crater]:
        """List of craters that have been emplaced in the simulation in the current interval in chronological order."""
        return self._emplaced

    @property
    def n_emplaced(self) -> int:
        """Number of craters that have been emplaced in the simulation in the current interval."""
        return len(self._emplaced)

    @property
    def n_observed(self) -> int:
        """Number of craters that have been observed on the surface in the current interval."""
        return len(self._observed)

    @property
    def Crater(self) -> type[Crater]:
        """The Crater class used for this counting component, which is determined by the morphology component."""
        if self._Crater is None:
            return Crater
        return self._Crater

    @Crater.setter
    def Crater(self, value):
        if value is None:
            return
        if not isinstance(value, type) or not issubclass(value, Crater):
            raise TypeError("Crater must be a subclass of the base Crater class.")
        self._Crater = value

    @property
    def morphology(self) -> Morphology:
        """The morphology component associated with this counting component, which determines the Crater class and the crater properties that are tracked in the simulation."""
        return self._morphology

    @morphology.setter
    def morphology(self, value):
        from cratermaker.components.morphology import Morphology

        if not isinstance(value, Morphology):
            raise TypeError("morphology must be an instance of the Morphology class.")
        self._morphology = value
        if not issubclass(self.Crater, value.Crater):
            self.Crater = value.Crater  # Set the Crater class to the one specified by the morphology component if they don't match
        return


def R_to_CSFD(
    R: Callable[[FloatLike | ArrayLike], FloatLike | ArrayLike],
    D: FloatLike | ArrayLike,
    Dlim: FloatLike = 1e6,
    *args: Any,
) -> FloatLike | ArrayLike:
    """
    Convert R values to cumulative N values for a given D using the R-plot function.

    Parameter
    ----------
    R : R = f(D)
        A function that computes R given D.
    D : FloatLike or ArrayLike
        diameter in units of km.
    Dlim : FloatLike
        Upper limit on the diameter over which to evaluate the integral
    args : Any
        Additional arguments to pass to the R function

    Returns
    -------
    float or ArrayLike
        The cumulative number of craters greater than D in diameter.
    """

    def _R_to_CSFD_scalar(R, D, Dlim, *args):
        # Helper function to integrate the R function
        def integrand(D):
            return R(D, *args) / D**3  # This is dN/dD

        N = 0.0
        D_i = D
        while D_i < Dlim:
            D_next = D_i * np.sqrt(2.0)
            D_mid = (D_i + D_next) / 2  # Mid-point of the bin
            bin_width = D_next - D_i
            R_value = integrand(D_mid)
            N += R_value * bin_width
            D_i = D_next  # Move to the next bin

        return N

    return _R_to_CSFD_scalar(R, D, Dlim, *args) if np.isscalar(D) else np.vectorize(_R_to_CSFD_scalar)(R, D, Dlim, *args)


def csfd_geometric_saturation(diameter: FloatLike | ArrayLike) -> FloatLike | ArrayLike:
    """
    Calculate the cumulative number of craters at geometric saturation for a given diameter.

    We use the definition of geomatric saturation from Melosh (1989) [#]_.

    Parameter
    ----------
    diameter : FloatLike or ArrayLike
        The diameter(s) for which to calculate the cumulative number of craters at geometric saturation.

    Returns
    -------
    FloatLike or ArrayLike
        The cumulative number of craters at geometric saturation for the given diameter(s).

    References
    ----------
    .. [#] Melosh, H.J., 1989. Impact cratering: A geologic process. Oxford University Press, New York, New York.

    """
    return 1.54 * diameter ** (-2)


def csfd_equilibrium(diameter: FloatLike | ArrayLike, f_geometric=0.0218) -> FloatLike | ArrayLike:
    """
    Calculate the cumulative number of craters at equilibrium for a given diameter.

    Parameter
    ----------
    diameter : FloatLike or ArrayLike
        The diameter(s) for which to calculate the cumulative number of craters at equilibrium.
    f_geometric : float, optional
        The fraction of geometric saturation at which equilibrium occurs. The default value is 0.0218, which is the value used in Minton et al. (2019) [#]_.


    Returns
    -------
    FloatLike or ArrayLike
        The cumulative number of craters at equilibrium for the given diameter(s).

    References
    ----------
    .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. https://doi.org/10.1016/j.icarus.2019.02.021

    """
    return f_geometric * csfd_geometric_saturation(diameter)


import_components(__name__, __path__)
