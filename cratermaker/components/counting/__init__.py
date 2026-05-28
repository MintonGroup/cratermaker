from __future__ import annotations

from abc import abstractmethod
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
import pyvista
import xarray as xr
from geopandas import GeoSeries
from matplotlib.axes import Axes
from numpy.typing import ArrayLike
from shapely.ops import transform
from shapely.validation import make_valid
from tqdm import tqdm
from vtk import vtkPolyData

from cratermaker import __version__ as cratermaker_version
from cratermaker.bindings import counting_bindings
from cratermaker.components.crater import _TALLY_LONG_NAME, Crater
from cratermaker.components.morphology import Morphology, MorphologyCrater
from cratermaker.constants import VECTOR_DRIVER_TO_EXTENSION_MAP, FloatLike
from cratermaker.core.base import ComponentBase, import_components

DRIVER_TO_EXTENSION_MAP = {"SCC": "scc", "VTK": "vtp", "VTP": "vtp", "CSV": "csv", **VECTOR_DRIVER_TO_EXTENSION_MAP}


if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface


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
        reset: bool = True,
        **kwargs: Any,
    ):
        """
        **Warning:** This object should not be instantiated directly. Instead, use the ``.maker()`` method.

        Parameters
        ----------
        surface : Surface | LocalSurface
            The surface or local surface view to be counted.
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
        object.__setattr__(self, "_surface", None)
        object.__setattr__(self, "_morphology", None)
        object.__setattr__(
            self,
            "_driver_to_extension_map",
            {"SCC": "scc", "VTK": "vtk", "VTP": "vtp", "CSV": "csv", **VECTOR_DRIVER_TO_EXTENSION_MAP},
        )
        self._surface = Surface.maker(surface, reset=reset, **kwargs)
        self._output_file_pattern += [
            f"observed_{self._output_file_prefix}*.{self._output_file_extension}",
            f"emplaced_{self._output_file_prefix}*.{self._output_file_extension}",
        ]

        if not reset:
            observed, emplaced = self.read_saved_output(interval=-1)
            if observed:
                interval = observed.interval.values[-1]
                observed = self.Crater.from_xarray(observed, interval=interval)
                for crater in observed:
                    self._observed[crater.id] = crater
                if emplaced:
                    self._emplaced = self.Crater.from_xarray(emplaced, interval=interval)
        return

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
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
            counting = "depthcount"

        counting = super().maker(
            component=counting,
            surface=surface,
            reset=reset,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        str_repr = super().__str__()
        str_repr += f"Number of observed craters: {self.n_observed}\n"
        str_repr += f"Number of emplaced craters: {self.n_emplaced}\n"
        str_repr += f"Surface: <{self.surface.component_name}>\n"
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

    def add(self, crater: MorphologyCrater, **kwargs: Any) -> MorphologyCrater:
        """
        Add a crater to the surface.

        Parameters
        ----------
        crater : Crater
            The crater to be added to the surface.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        MorphologyCrater:
            The emplacec crater
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")
        if not isinstance(crater, MorphologyCrater):
            crater = self.morphology.Crater.maker(crater=crater)

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
                self.observed[crater.id] = crater
            else:
                return

            rim_interior_region = crater_region.extract_subregion(subregion_radius=crater.radius)

            # Cookie cutting: remove any smaller craters that are overlapped by this new crater
            unique_ids = np.unique(rim_interior_region.crater_id)
            unique_ids = unique_ids[unique_ids > 0]  # Remove the 0 id which corresponds to no crater
            if len(unique_ids) > 0:
                # Compute cookie cutting removes list
                observed = self.observed.copy()
                removes = [id for id, v in observed.items() if v.id in unique_ids and v.diameter < crater.diameter]
                # For every id that appears in the removes list, set it to 0 in the data array
                for remove_id in removes:
                    rim_interior_region.remove_tag(name="crater_id", tag=remove_id)
                    if not np.any(
                        self.surface.uxds.crater_id.data == remove_id
                    ):  # Check to see if this crater id still appears, and if not, it's gone man.
                        self.observed.pop(remove_id, None)

        return crater

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

        if bp > ap:
            ap, bp = bp, ap
            orientation += np.pi / 2
        crater.measured_semimajor_axis = ap
        crater.measured_semiminor_axis = bp
        crater.measured_orientation = np.degrees(orientation)
        crater.measured_location = location

        return crater

    def score_rim(
        self,
        crater: Crater,
        quantile: float = 0.95,
        gradmult: float = 1.0,
        curvmult: float = 1.0,
        heightmult: float = 1.0,
        **kwargs: Any,
    ) -> None:
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
        **kwargs : Any
            |kwargs|

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
        .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. `doi: 10.1016/j.icarus.2019.02.021 <https://doi.org/10.1016/j.icarus.2019.02.021>`_

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
        from cratermaker.components.crater import _convert_tuple_vars

        if len(craters) == 0:
            return xr.Dataset()
        data = []
        if isinstance(craters, dict):
            craters = craters.values()
        for c in craters:
            d = c.as_dict(skip_complex_data=True)
            d = _convert_tuple_vars(input_dict=d, inverse=False)
            d = xr.Dataset(data_vars=d).set_coords("id").expand_dims(dim="id")
            d["id"].attrs["long_name"] = _TALLY_LONG_NAME
            data.append(d)

        return xr.concat(data, dim="id").sortby("id")

    def save(
        self,
        crater_type: Literal["observed", "emplaced", "both"] = "both",
        interval: int = 0,
        time_variables: dict | None = None,
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
        time_variables: dict | None, optional
            A dictionary of time variables to include in the output file.
        **kwargs : Any
            |kwargs|
        """
        if crater_type not in ["observed", "emplaced", "both"]:
            raise ValueError(f"Invalid crater_type: {crater_type}. Must be one of 'observed', 'emplaced', or 'both'.")
        if time_variables is not None and not isinstance(time_variables, dict):
            raise TypeError("time_variables must be a dictionary")

        self.remove_complex_data()

        if crater_type == "observed":
            iterator = zip([self.observed], ["observed"], strict=True)
        elif crater_type == "emplaced":
            iterator = zip([self.emplaced], ["emplaced"], strict=True)
        elif crater_type == "both":
            iterator = zip([self.observed, self.emplaced], ["observed", "emplaced"], strict=True)
        for craters, crater_type in iterator:
            if craters:
                filename = self.output_dir / f"{crater_type}_{self.output_filename(interval)}"
                ds = self.to_xarray(craters)
                ds = ds.expand_dims(dim="interval").assign_coords({"interval": [interval]})
                if time_variables is not None:
                    for k, v in time_variables.items():
                        ds[k] = xr.DataArray(data=[v], name=k, dims=["interval"], coords={"interval": [interval]})
                if filename.exists():
                    filename.unlink()
                ds.to_netcdf(filename)
        save_args = {"interval": interval, **kwargs}
        super().save(**save_args)
        return

    def export(
        self,
        crater_type: Literal["observed", "emplaced", "both"] = "both",
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
        if crater_type not in ["observed", "emplaced", "both"]:
            raise ValueError("crater_type must be 'observed', 'emplaced', or 'both")

        observed, emplaced = self.read_saved_output(interval=interval)
        if crater_type == "both":
            crater_types = ["observed", "emplaced"]
            crater_ds_list = [observed, emplaced]
        else:
            crater_types = [crater_type]
            if crater_type == "observed":
                crater_ds_list = [observed]
            elif crater_type == "emplaced":
                crater_ds_list = [emplaced]

        for crater_type, crater_ds in zip(crater_types, crater_ds_list, strict=True):
            if isinstance(crater_ds, dict):
                interval_numbers = list(crater_ds.keys())
            elif isinstance(crater_ds, xr.Dataset) and "interval" in crater_ds.coords:
                interval_numbers = crater_ds.interval.values
            else:
                if interval is None:
                    interval_numbers = [0]
                else:
                    interval_numbers = [interval]
            for interval in interval_numbers:
                craters = self.Crater.from_xarray(crater_ds, interval=interval) if crater_ds is not None else []
                if driver.upper() in ["VTK", "VTP"] and crater_ds is not None:
                    self.to_vtk_file(
                        craters=craters,
                        interval=interval,
                        crater_type=crater_type,
                        **kwargs,
                    )
                elif driver.upper() == "CSV" and crater_ds is not None:
                    self.to_csv_file(
                        craters=craters,
                        interval=interval,
                        crater_type=crater_type,
                        **kwargs,
                    )
                elif driver.upper() == "SCC":
                    if crater_ds is None:
                        crater_ds = []  # Allow for empty SCC files because they can still have the region polygon
                    self.to_scc_file(
                        craters=craters,
                        interval=interval,
                        crater_type=crater_type,
                        **kwargs,
                    )
                elif driver.upper() in VECTOR_DRIVER_TO_EXTENSION_MAP and crater_ds is not None:
                    self.to_vector_file(
                        craters=craters,
                        interval=interval,
                        crater_type=crater_type,
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
        label: str | None = None,
        scalebar: bool | None = None,
        colorbar: bool = True,
        show: bool = False,
        save: bool = True,
        ax: Axes | None = None,
        close_when_done: bool = True,
        minimum_plot_width: float | None = 800,
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
        label : str | None, optional
            A label for the plot. If None, no label will be added.
        scalebar : bool, optional
            If True, a scalebar will be added to the plot. Default is True.
        colorbar : bool, optional
            If True, a colorbar will be added to the plot when using "map" plot_style or "hillshade" with a variable overlay. Default is True.
        show : bool, optional
            If True, the plot will be displayed. Default is True.
        save : bool, optional
            If True, the plot will be saved to the default plot directory. Default is True
        ax : matplotlib.axes.Axes, optional
            An existing Axes object to plot on. If None, a new figure and axes will be created.
        close_when_done : bool, optional
            If True, the figure will be closed after plotting. Default is True when save is True and show is False, and False otherwise.
        minimum_plot_width : float, optional
            Because the width of the plot is determined by the number of faces, small regions will generate small plots with labels that are hard to read. This parameter sets a lower limit to the width of the image that is generated by the plot. By default it is 800. Set to None to turn it off.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        matplotlib.image.Axes object created by the surface plot method, with crater counts plotted on top if specified. If show is True, the plot will also be displayed.
            The Axes object
        """
        import matplotlib.pyplot as plt

        from cratermaker.components.surface.hireslocal import HiResLocalSurface

        if close_when_done is None:
            close_when_done = save and not show

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
        if variable_name is not None:
            file_prefix += f"_{variable_name}"
        if interval is not None:
            observed_intervals, emplaced_intervals = self.get_saved_interval_numbers()
            # Check if we are in the current interval, otherwise we have to read data from file.
            trigger_read = False
            if observed_intervals is not None and (len(observed_intervals) == 0 or observed_intervals[-1] == interval):
                observed = list(self.observed.values())
            else:
                observed = None
                trigger_read = True
            if emplaced_intervals is not None and (len(emplaced_intervals) == 0 or emplaced_intervals[-1] == interval):
                emplaced = self.emplaced
            else:
                emplaced = None
                trigger_read = True

            if trigger_read:
                observed_ds, emplaced_ds = self.read_saved_output(interval=interval)
                if observed is None and observed_ds:
                    interval = observed_ds.interval.values[-1]
                    observed = self.Crater.from_xarray(observed_ds, interval=interval)
                if emplaced is None and emplaced_ds:
                    emplaced_interval = emplaced.interval.values[-1]
                    if emplaced_interval == interval:
                        emplaced = self.Crater.from_xarray(emplaced_ds, interval=interval)
            filename = self.plot_dir / f"{file_prefix}{interval:06d}.{self.surface.output_image_file_extension}"
        else:
            observed = list(self.observed.values())
            emplaced = self.emplaced
            filename = self.plot_dir / f"{file_prefix}.{self.surface.output_image_file_extension}"

        ax = self.surface.plot(
            plot_style=plot_style,
            variable_name=variable_name,
            interval=interval,
            cmap=cmap,
            label=label,
            scalebar=scalebar,
            colorbar=colorbar,
            show=False,
            save=False,
            ax=ax,
            minimum_plot_width=minimum_plot_width,
            close_when_done=False,
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
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, zorder=100)
        if observed_original_color is not None and observed is not None and len(observed) > 0:
            gs = self.to_geoseries(
                craters=observed, use_measured_properties=False, split_antimeridian=split_antimeridian, autolim=False
            ).to_crs(crs)
            facecolor = kwargs.pop("facecolor", "none")
            edgecolor = observed_original_color
            linewidth = kwargs.pop("linewidth", 0.1)
            linestyle = kwargs.pop("linestyle", ":")
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, zorder=200)

        if observed_color is not None and observed is not None and len(observed) > 0:
            gs = self.to_geoseries(
                craters=observed, use_measured_properties=True, split_antimeridian=split_antimeridian, autolim=False
            ).to_crs(crs)
            facecolor = kwargs.pop("facecolor", "none")
            edgecolor = observed_color
            linewidth = kwargs.pop("linewidth", 0.1)
            linestyle = kwargs.pop("linestyle", "solid")
            ax = gs.plot(ax=ax, facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, linestyle=linestyle, zorder=300)

        if save:
            plt.savefig(filename, dpi=ax.figure.get_dpi())
        if show:
            plt.show()
        if close_when_done and (save or show):
            plt.close()
        return ax

    def pyvista_plotter(
        self,
        crater_type: Literal["observed", "emplaced", "both"] = "both",
        interval: int | None = None,
        enable_interactive: bool = True,
        crater_color: str | tuple[str] = ("white", "red"),
        crater_style: Literal["rings", "points", "impacts", "spheres"] = "rings",
        surface: Surface | LocalSurface | None = None,
        plotter: pyvista.Plotter | None = None,
        crater_size_scale_factor: FloatLike = 1.0,
        **kwargs: Any,
    ) -> pyvista.Plotter:
        """
        Passes through to the surface pyvista_plotter method and adds crater counts to it.

        When enable_interactive is True (the default), this will add all saved crater data (observed and emplaced) with all plot styles (rings, points, and impacts) and allow them to be activated by key presses. If False, this will only plot one set of data with one style.

        Parameters
        ----------
        crater_type: Literal["observed","emplaced","both"] | list[Crater].
            This is only used if enable_interactive is False, in which case it controls which crater data will be plotted. Default is "both".
        interval : int, optional
            The interval number to load the emplaced crater data from. if None, then all emplaced data currently saved to file is used. Default is None.
        enable_interactive : bool, optional
            If True, enables key events to toggle the visibility of all saved crater count types (observed and emplaced) with all styles (rings, points, and impacts), which are activated by keypresses ("c" for observed craters, and "t" for emplaced craters), and all are not visible on initialization. If False, only one set of data is plotted, and is visible on initialization.
        crater_color: str | tuple[str] = ("white", "yellow")
            The color to use for the plots. When used with enable_interactive or when craters is set to "both", this is a tuple where the first element is the color of observed craters and the second is the color of emplaced craters. Default is ("white,"red").
        crater_style : Literal["rings", "points"], optional
            Only used when enable_interactive is False. Sets the style of the mesh. Options are "rings", which creates polyline circles over the rim of each crater or "points" which creates a point at the center.
        surface : Surface | LocalSurface, optional
            The surface or local surface view to be displayed. If None, uses the associated surface property
        plotter : pyvista.Plotter, optional
            An existing pyvista Plotter to add the crater counts to. If None, a new Plotter will be created by the surface pyvista_plotter method. Default is None.
        crater_size_scale_factor : FloatLike, optional
            A factor to scale the size of the craters in "point", "impacts", or "spheres" styles. Default is 1.0.
        **kwargs : Any
            |kwargs|

        Returns
        -------
        plotter : pyvista.Plotter
            The pyvista plotter with the crater counts added.

        """
        from cratermaker.constants import PYVISTA_ADD_MESH_KWARGS
        from cratermaker.utils.general_utils import toggle_pyvista_actor, update_pyvista_help_message

        def _draw_craters(craters, crater_color, name, crater_style, use_measured_properties, **kwargs):
            add_mesh_kwargs = {k: v for k, v in kwargs.items() if k in PYVISTA_ADD_MESH_KWARGS}
            add_mesh_kwargs = {"color": crater_color, **add_mesh_kwargs}
            if crater_style == "rings":
                add_mesh_kwargs = {**add_mesh_kwargs, "line_width": 2}
                point_size = 1
            elif crater_style == "points":
                add_mesh_kwargs["render_points_as_spheres"] = True
                add_mesh_kwargs["style"] = "points_gaussian"
                add_mesh_kwargs["emissive"] = True
                point_size = 1
                size_scale = np.array([crater_size_scale_factor * surface.face_size[c.face_index] for c in craters])
            elif crater_style == "impacts":
                add_mesh_kwargs["render_points_as_spheres"] = False
                add_mesh_kwargs["style"] = "points_gaussian"
                add_mesh_kwargs["emissive"] = True
                add_mesh_kwargs["pbr"] = True
                size_scale = np.array([crater_size_scale_factor * c.floor_radius for c in craters])
                point_size = 1
            elif crater_style == "spheres":
                add_mesh_kwargs["render_points_as_spheres"] = True
                add_mesh_kwargs["style"] = "points_gaussian"
                add_mesh_kwargs["pbr"] = True
                point_size = 1
                size_scale = np.array([crater_size_scale_factor * c.projectile_radius for c in craters])
            mesh = self.to_vtk_mesh(
                craters=craters,
                use_measured_properties=use_measured_properties,
                crater_style=crater_style,
                crater_size_scale_factor=crater_size_scale_factor,
                **kwargs,
            )
            pdata = pyvista.PolyData(mesh)
            if crater_style in ["points", "impacts", "spheres"]:
                pdata["size_scale"] = size_scale
            actor_name = f"{crater_type}_{crater_style}"
            actor = plotter.add_mesh(pdata, name=actor_name, point_size=point_size, **add_mesh_kwargs)
            if crater_style in ["points", "impacts", "spheres"]:
                actor.mapper.scale_array = "size_scale"
            return actor

        def update_crater_style(plotter, actor_list):
            inext = 0
            for i, actor in enumerate(actor_list):
                if actor.GetVisibility():
                    actor.SetVisibility(False)
                    inext = i + 1
                    break
            if inext < len(actor_list):
                actor_list[inext].SetVisibility(True)
            plotter.update()
            return

        new_plotter = plotter is None
        surface = self.surface
        plotter = surface.pyvista_plotter(enable_interactive=enable_interactive, plotter=plotter, interval=interval, **kwargs)
        valid_crater_styles = ["rings", "points", "impacts", "spheres"]
        if not enable_interactive:
            crater_style = crater_style.lower()
            if crater_style not in valid_crater_styles:
                raise ValueError(f"Invalid crater_style: {crater_style}. Must be one of {valid_crater_styles}.")

        if interval is None:
            interval = -1

        if isinstance(crater_color, str):
            crater_color = [crater_color]
        elif not isinstance(crater_color, (tuple, list)) or (crater_type == "both" and len(crater_color) != 2):
            raise ValueError(
                "crater_color must be a tuple with the first element the color of observed craters and the second the color of emplaced craters."
            )

        if isinstance(crater_type, str):
            if crater_type not in ["observed", "emplaced", "both"]:
                raise ValueError("crater_type must be 'observed', 'emplaced', 'both', or a list of Crater objects")
            if interval == -1:
                observed = list(self.observed.values())
                emplaced = self.emplaced
            else:
                observed, emplaced = self.read_saved_output(interval=interval)
                if observed is not None and len(observed) > 0:
                    observed = self.Crater.from_xarray(observed)
                else:
                    observed = None
                if emplaced is not None and len(emplaced) > 0:
                    emplaced = self.Crater.from_xarray(emplaced)
                else:
                    emplaced = None
            if crater_type == "both":
                craterlistlist = [observed, emplaced]
                namelist = ["observed", "emplaced"]
                measured = [True, False]
                toggle = ["c", "t"]
            elif crater_type == "observed":
                craterlistlist = [observed]
                namelist = ["observed"]
                measured = [True]
                toggle = ["c"]
                crater_color = [crater_color[0]]
            elif crater_type == "emplaced":
                craterlistlist = [emplaced]
                namelist = ["emplaced"]
                measured = [False]
                toggle = ["t"]
                crater_color = [crater_color[-1]]

        if enable_interactive:
            styles = valid_crater_styles
        else:
            styles = [crater_style]

        for name, craters, use_measured_properties, col, key in zip(
            namelist, craterlistlist, measured, crater_color, toggle, strict=True
        ):
            if craters is None:
                continue
            actor_list = []
            for crater_style in styles:
                actor = _draw_craters(craters, col, name, crater_style, use_measured_properties, **kwargs)

                if enable_interactive:
                    actor.SetVisibility(False)
                    actor_list.append(actor)
            if enable_interactive and new_plotter:
                new_message = f"{key} Toggle {name} craters"
                plotter = update_pyvista_help_message(plotter, new_message=new_message)
                plotter.add_key_event(key, lambda plotter=plotter, actor_list=actor_list: update_crater_style(plotter, actor_list))

        return plotter

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
                crater.to_geoseries(
                    surface=surface, split_antimeridian=split_antimeridian, use_measured_properties=use_measured_properties
                )
            )
        return pd.concat(crater_gs)

    def _validate_export_args(
        self,
        crater_type: Literal["observed", "emplaced"] = "observed",
        interval: int | None = None,
        craters: list[Crater] | None = None,
    ) -> list[Crater]:
        if crater_type not in ["observed", "emplaced"]:
            raise ValueError("crater_type must be 'observed' or 'emplaced'")
        if craters is None:
            craters = getattr(self, crater_type)
        if type(craters) is dict:
            if type(list(craters.values())[0]) is xr.Dataset:  # like multi-interval dicts keyed by interval
                crater_ds = craters.get(
                    interval, []
                )  # Get the dataset for the specified interval, or an empty dataset if not found
            else:
                return list(crater_ds.values())  # like observed dicts keyed by crater id

        if isinstance(craters, dict):
            craters = list(craters.values())
        elif isinstance(craters, xr.Dataset):
            craters = self.Crater.from_xarray(craters, interval=interval)
        elif isinstance(craters, list):
            if not all(isinstance(c, Crater) for c in craters):
                raise TypeError("All elements of the craters list must be Crater objects.")
        else:
            raise TypeError(f"Unrecognized type for craters: {type(craters)}")

        return craters

    def to_vector_file(
        self,
        crater_type: Literal["observed", "emplaced"] = "observed",
        craters: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        driver: str = "GPKG",
        use_measured_properties: bool = True,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a vector file and stores it in the default export directory.

        Parameters
        ----------
        crater_type : Literal["observed", "emplaced"], optional
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

        # craters = self._validate_export_args(name=name, interval=interval, crater_ds=crater_ds)

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
        #     craters,
        #     total=len(craters),
        #     desc=f"Converting {crater_type} craters to geometries for export",
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
        #     print(f"Saving {crater_type} layer to vector file: '{output_file}'...")
        # else:
        #     output_file = self.export_dir / f"{crater_type}_{self.output_file_prefix}{interval:06d}.{file_extension}"
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

    def to_vtk_mesh(
        self,
        crater_type: Literal["observed", "emplaced"] | None = None,
        interval: int | None = None,
        craters: list[Crater] | None = None,
        use_measured_properties: bool = True,
        crater_style: Literal["rings", "points", "impacts", "spheres"] = "rings",
        crater_size_scale_factor: FloatLike = 1.0,
        **kwargs: Any,
    ) -> vtkPolyData:
        """
        Convert the crater data to a VTK PolyData mesh.

        Parameters
        ----------
        crater_type : Literal["observed", "emplaced"], optional
            The name of the crater dataset to convert, either "observed" or "emplaced. If None is provided, then a list of Crater objects must be provided directly to the `craters` parameter. Default is None.
        interval : int | None, optional
            The interval number to load if craters is not provided directly. If None, then the most recent interval will be used. Default is None.
        craters : list[Crater], optional
            A list of Crater objects to convert. If None is provided, then the crater dataset will be determined by the `name` parameter. Default is None.
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        crater_style : Literal["rings", "points", "impacts", "spheres"], optional
            Sets the style of the mesh. Options are "rings", which creates polyline circles over the rim of each crater, "points" which creates a small sphere at the center of each crater, and "impacts" which places a point above the floor of the center of the crater and "spheres" which creates spheres with radius equal to the projectile radius at the location of each crater. Default is "rings".
        crater_size_scale_factor : FloatLike, optional
            A factor to scale the size of the craters in "point", "impacts", or "spheres" styles, which places the center point of the actor above the surface by its radius. Default is 1.0.
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
            """Yield Nx3 arrays for each exterior ring (keep closing vertex)."""
            g3d = transform(lonlat_to_xyz(R), geom)

            def _rings(g):
                if g.geom_type == "Polygon":
                    yield np.asarray(g.exterior.coords)
                    for i in g.interiors:
                        yield np.asarray(i.coords)
                elif g.geom_type in ("MultiPolygon", "GeometryCollection"):
                    for sub in g.geoms:
                        yield from _rings(sub)
                else:
                    yield np.asarray(g.coords)

            yield from _rings(g3d)

        def draw_rings(craters, surface):
            R = surface.radius
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
                for ring_xyz in polygon_xyz_coords(g.item(), R):
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

        def draw_points(craters, surface, crater_style):
            from shapely.geometry import Point
            from vtk import VTK_FLOAT
            from vtkmodules.util.numpy_support import numpy_to_vtk

            geoms = []
            for crater in craters:
                z = surface.face_elevation[crater.face_index]
                if crater_style == "points":
                    z += crater_size_scale_factor * surface.face_size[crater.face_index] / 2
                elif crater_style == "spheres":
                    z += crater_size_scale_factor * crater.projectile_radius
                elif crater_style == "impacts":
                    z += np.sqrt(crater_size_scale_factor) * crater.floor_radius / 2
                geoms.append(Point(crater.location[0], crater.location[1], z))
            gs = GeoSeries(geoms, crs=surface.crs)

            # Convert to Cartesian coordinates using transform
            R = surface.radius
            points = vtkPoints()
            for geom in gs:
                x, y, z = lonlat_to_xyz(R)(geom.x, geom.y, geom.z)
                points.InsertNextPoint(float(x), float(y), float(z))
            poly_data = vtkPolyData()
            poly_data.SetPoints(points)
            return poly_data

        valid_crater_styles = ["rings", "points", "impacts", "spheres"]
        crater_style = crater_style.lower()
        if crater_style not in valid_crater_styles:
            vtxt = ", ".join(f'"{s}"' for s in valid_crater_styles)
            raise ValueError(f"{crater_style} is not a recognized crater_style. Valid options are {vtxt}.")

        surface = self.surface

        if craters is None:
            if crater_type is None:
                raise ValueError(
                    "If craters is not provided directly, then name must be provided to determine which crater dataset to use."
                )
            elif crater_type not in ["observed", "emplaced"]:
                raise ValueError("name must be either 'observed' or 'emplaced'")
            if interval is None:
                if crater_type == "observed":
                    craters = list(self.observed.values())
                else:
                    craters = self.emplaced
            else:
                observed, emplaced = self.read_saved_output(interval=interval)
                if crater_type == "observed":
                    if observed is None or len(observed) == 0:
                        return None
                    interval = observed.interval.values[-1]
                    craters = self.Crater.from_xarray(observed, interval=interval)
                else:
                    if emplaced is None or len(emplaced) == 0:
                        return None
                    interval = emplaced.interval.values[-1]
                    craters = self.Crater.from_xarray(emplaced, interval=interval)
            if craters is None:
                return None

        if crater_style == "rings":
            return draw_rings(craters, surface)
        else:
            return draw_points(craters, surface, crater_style)

    def to_vtk_file(
        self,
        crater_type: Literal["observed", "emplaced"] | None = None,
        interval: int | None = None,
        craters: list[Crater] | None = None,
        use_measured_properties: bool = True,
        crater_style: Literal["rings", "points", "impacts", "spheres"] = "rings",
        crater_size_scale_factor: FloatLike = 1.0,
        output_file: str | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a VTK file and stores it in the default export directory.

        Notes: In order for the crater and surface to be synced up when saving to VTK/VTP format, the initial conditions (no craters) must be saved. Otherwise, saving to file only occurs if there are craters to save.

        Parameters
        ----------
        crater_type : Literal["observed", "emplaced"], optional
            The name of the crater dataset to convert, either "observed" or "emplaced. If None is provided, then a list of Crater objects must be provided directly to the `craters` parameter. Default is None.
        interval : int | None, optional
           |interval_export|
        craters : list[Crater], optional
            A list of Crater objects to convert. If None is provided, then the crater dataset will be determined by the `crater_type` parameter. Default is None.
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        crater_style : Literal["rings", "points", "impacts", "spheres"], optional
            Sets the style of the mesh. Options are "rings", which creates polyline circles over the rim of each crater, "points" which creates a small sphere at the center of each crater, and "impacts" which places a point above the floor of the center of the crater and "spheres" which creates spheres with radius equal to the projectile radius at the location of each crater. Default is "rings".
        crater_size_scale_factor : FloatLike, optional
            A factor to scale the size of the craters in "point", "impacts", or "spheres" styles, which places the center point of the actor above the surface by its radius. Default is 1.0.
        output_file : str | None, optional
            The file path to save the VTK file to. If None, the file will be
        **kwargs : Any
            |kwargs|
        """
        from vtk import vtkXMLPolyDataWriter

        craters = self._validate_export_args(crater_type=crater_type, interval=interval, craters=craters)

        if output_file is None:
            filename_base = self.output_filename(interval).replace(self.output_file_extension, "vtp")
            output_file = self.export_dir / f"{crater_type}_{filename_base}"
        else:
            output_file = Path(output_file)
        if not self._overwrite_check(output_file):
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
        crater_type: Literal["observed", "emplaced"] = "observed",
        craters: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        output_file: str | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to a CSV file and stores it in the default export directory.

        Parameters
        ----------
        crater_type : Literal["observed", "emplaced"], optional
            The type of the crater dataset to export, either "observed" or "emplaced
        craters : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the crater_type parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
            |interval_export|
        output_file : str | None, optional
            The file path to save the CSV file to. If None, the file will be saved to the default export directory with a filename based on the crater_type and interval. Default is None.
        **kwargs : Any
            |kwargs|
        """
        import csv

        from cratermaker.components.crater import _convert_tuple_vars

        craters = self._validate_export_args(crater_type=crater_type, interval=interval, craters=craters)

        if output_file is None:
            filename_base = self.output_filename(interval).replace(self.output_file_extension, "csv")
            output_file = self.export_dir / f"{crater_type}_{filename_base}"
        else:
            output_file = Path(output_file)
        if not self._overwrite_check(output_file):
            return
        print(f"Saving crater data to CSV file: '{output_file}'...")
        with output_file.open(mode="w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            header_written = False
            for crater in craters:
                crater_dict = crater.as_dict(skip_complex_data=True)
                crater_dict = _convert_tuple_vars(input_dict=crater_dict, inverse=False)

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
                for k, v in crater_dict.items():
                    if v is None:
                        crater_dict[k] = ""

                if not header_written:
                    header = list(crater_dict.keys())
                    writer.writerow(header)
                    header_written = True
                row = [crater_dict.get(key, "") for key in header]
                writer.writerow(row)

        return

    def to_scc_file(
        self,
        crater_type: Literal["observed", "emplaced"] = "observed",
        craters: xr.Dataset | list[Crater] | dict[int, Crater] | None = None,
        interval: int | None = None,
        output_file: str | None = None,
        **kwargs,
    ) -> None:
        """
        Export the crater data to Craterstats and OpenCraterTools-compatible SCC file and stores it in the default export directory.

        Parameters
        ----------
        crater_type : Literal["observed", "emplaced"], optional
            The type of the crater dataset to export, either "observed" or "emplaced
        craters : xr.Dataset | list[Crater] | dict[int, Crater] | None, optional
            The crater data to export. Can be provided as an xarray Dataset, a list of Crater objects, or a dictionary mapping interval numbers to Crater objects. If None, the crater data will be the attribute of the class corresponding to the crater_type parameter (self.observed or self.emplaced). Default is None.
        interval : int | None, optional
            |interval_export|
        output_file : str | None, optional
            The file path to save the SCC file to. If None, the file will be saved
        **kwargs : Any
            |kwargs|
        """
        import datetime

        craters = self._validate_export_args(crater_type=crater_type, interval=interval, craters=craters)

        def overlap_fraction(crater, region_poly=None):
            if region_poly is None:
                return 1.0
            distance = self.surface.compute_distances(reference_location=self.surface.local_location, locations=[crater.location])
            if distance + crater.measured_radius > self.surface.local_radius:
                crater_poly = crater.to_geoseries(
                    surface=self.surface, split_antimeridian=False, use_measured_properties=True
                ).to_crs(self.surface.crs)
                if not crater_poly.is_valid[0]:
                    crater_poly = crater_poly.make_valid()
                overlap_area = crater_poly.intersection(region_poly).to_crs(self.surface.local.crs).area.item()
                return overlap_area / crater_poly.to_crs(self.surface.local.crs).area.item()
            else:
                return 1.0

        region_poly = None

        if output_file is None:
            output_file = self.export_dir / f"{crater_type}{interval:06d}.scc"
        else:
            output_file = Path(output_file)
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
                    region_circle = self.Crater.maker(
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
            for crater in craters:
                f.write(
                    f"{crater.measured_diameter * 1e-3}\t{overlap_fraction(crater, region_poly)}\t{crater.measured_location[0]}\t{crater.measured_location[1]}\t 1\n"
                )
            f.write("}\n")

        return

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
        if self._morphology is None:
            return Crater
        else:
            return self._morphology.Crater

    @property
    def morphology(self) -> Morphology:
        """
        The morphology component associated with this counting component, which determines the Crater class and the crater properties that are tracked in the simulation.

        Note
        ----
        The Morphology component requires a Countinng component, so this property is set by the Morphology component on initialization, rather than by the Counting component itself.
        """
        return self._morphology

    @morphology.setter
    def morphology(self, value):
        from cratermaker.components.morphology import Morphology

        if not isinstance(value, Morphology):
            raise TypeError("morphology must be an instance of the Morphology class.")
        self._morphology = value
        # Re-run any saved emplaced or observed craters through the new object's maker function so that we promote them to the new type
        for i, crater in enumerate(self.emplaced):
            self.emplaced[i] = self.Crater.maker(crater=crater)
        for id, crater in self.observed.items():
            self.observed[id] = self.Crater.maker(crater=crater)
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
    .. [#] Minton, D.A., Fassett, C.I., Hirabayashi, M., Howl, B.A., Richardson, J.E., (2019). The equilibrium size-frequency distribution of small craters reveals the effects of distal ejecta on lunar landscape morphology. Icarus 326, 63-87. `doi: 10.1016/j.icarus.2019.02.021 <https://doi.org/10.1016/j.icarus.2019.02.021>`_

    """
    return f_geometric * csfd_geometric_saturation(diameter)


import_components(__name__, __path__)
