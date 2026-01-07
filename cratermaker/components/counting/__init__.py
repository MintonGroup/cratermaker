from __future__ import annotations

from abc import abstractmethod
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
import uxarray as uxr
import xarray as xr
from cratermaker._cratermaker import counting_bindings
from shapely.geometry import GeometryCollection
from shapely.ops import transform
from vtk import vtkPolyData

from cratermaker.components.crater import Crater
from cratermaker.core.base import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_VARIABLE_NAME = "crater_id"
_TALLY_LONG_NAME = "Unique crater identification number"

_N_LAYER = (
    8  # The number of layers used for tagging faces with crater ids. This allows a single face to contain multiple crater ids
)
_MIN_FACE_FOR_COUNTING = 5
_RIM_BUFFER_FACTOR = 1.2  # The factor by which the crater tagging region is extended beyond the final rim.
_EXTENT_RADIUS_RATIO = 2.0  # The factor by radius over which the local region that is extracted to evaluate the crater rim


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
        self.surface = surface

        rng = kwargs.pop("rng", surface.rng)
        simdir = kwargs.pop("simdir", surface.simdir)
        super().__init__(rng=rng, simdir=simdir, reset=reset, ask_overwrite=ask_overwrite, **kwargs)

        object.__setattr__(self, "_emplaced", [])
        object.__setattr__(self, "_observed", {})
        object.__setattr__(self, "_output_dir_name", "craters")
        object.__setattr__(self, "_output_file_prefix", "craters")
        object.__setattr__(self, "_output_file_extension", "nc")

        self._output_file_pattern += [f"*{self._output_file_prefix}*.{self._output_file_extension}"]

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

        if surface is None:
            raise ValueError("surface must be provided")

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
        dims = ("n_face", "layer")
        data = np.zeros((self.surface.n_face, self.n_layer), dtype=np.uint32)
        uxda = uxr.UxDataArray(
            data=data,
            dims=dims,
            attrs={"long_name": _TALLY_LONG_NAME},
            name=_TALLY_VARIABLE_NAME,
            uxgrid=self.surface.uxgrid,
        )

        self.surface._uxds[_TALLY_VARIABLE_NAME] = uxda
        self._emplaced = []
        self._observed = {}

        super().reset(ask_overwrite=ask_overwrite, **kwargs)
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
        if _TALLY_VARIABLE_NAME not in self.surface.uxds:
            self.reset()

        self.emplaced.append(crater)
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
                if np.any(self.surface.uxds[_TALLY_VARIABLE_NAME].isel(layer=i).data[count_region.face_indices] > 0):
                    # Gather the unique id values for the current layer
                    unique_ids = np.unique(self.surface.uxds[_TALLY_VARIABLE_NAME].data[count_region.face_indices, i])
                    removes = [
                        id for id, v in self.observed.items() if v.id in unique_ids and v.final_diameter < crater.final_diameter
                    ]

                    # For every id that appears in the removes list, set it to 0 in the data array
                    if removes:
                        data = self.surface.uxds[_TALLY_VARIABLE_NAME].data[count_region.face_indices, :]
                        for remove in removes:
                            data[data == remove] = 0
                        self.surface.uxds[_TALLY_VARIABLE_NAME].data[count_region.face_indices, :] = data
                if insert_layer == -1 and np.all(
                    self.surface.uxds[_TALLY_VARIABLE_NAME].isel(layer=i).data[count_region.face_indices] == 0
                ):
                    insert_layer = i
            if insert_layer == -1:
                raise ValueError("Crater counting layers are full")
            data = self.surface.uxds[_TALLY_VARIABLE_NAME].data[count_region.face_indices, :]
            data[:, insert_layer] = crater.id
            self.surface.uxds[_TALLY_VARIABLE_NAME].data[count_region.face_indices, :] = data

        return

    def remove(self, crater_id: int) -> None:
        """
        Remove all instances of a crater id from the surface.

        Parameters
        ----------
        crater_id : int
            The ID of the crater to be removed.
        """
        idda = self.surface.uxds[_TALLY_VARIABLE_NAME]
        remove_mask = idda == crater_id
        self.surface.uxds[_TALLY_VARIABLE_NAME] = xr.where(remove_mask, xr.zeros_like(idda), idda, keep_attrs=True)
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

        region = counting_bindings.score_rim(self.surface, crater, quantile, gradmult, curvmult, heightmult)

        return region

    def measure_crater_depth(self, crater: Crater) -> float:
        """
        Measure the depth of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to measure the depth.

        Returns
        -------
        float
            The measured depth of the crater.
        """
        if not isinstance(crater, Crater):
            raise TypeError("crater must be an instance of Crater")

        depth = counting_bindings.measure_crater_depth(self.surface, crater)

        return depth

    @abstractmethod
    def tally(self, region: LocalSurface | None = None, **kwargs: Any) -> None: ...

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

        def _convert_and_merge(craters: dict[int, Crater] | list[Crater], filename: Path | str, layer_name: str) -> None:
            # Convert into an xarray dataset
            dsnew = self.to_xarray(craters)

            # If the file already exists, read it and merge
            if filename.exists():
                with xr.open_dataset(filename) as ds:
                    combined_data = xr.concat([ds, dsnew], dim=_TALLY_VARIABLE_NAME)
            else:
                combined_data = dsnew

            # Write merged data back to file
            if combined_data:
                combined_data.to_netcdf(filename)

        if self.emplaced:
            _convert_and_merge(self.emplaced, emplaced_filename, "emplaced_craters")
        if self.observed:
            _convert_and_merge(self.observed, observed_filename, "observed_craters")
        self._emplaced = []
        return

    def export(
        self,
        craters: Crater | list[Crater] | dict[int, Crater] | str,
        name: str | None = None,
        interval_number: int | None = None,
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
        interval_number : int, default=0
            The interval number for the output file naming.
        driver : str, default='GPKG'
            The file format to save. Supported formats are 'VTK', 'GPKG', 'ESRI Shapefile', etc.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments to pass to the make_vector_file function.
        """
        if isinstance(craters, str):
            if craters.lower() == "emplaced":
                craters = self.emplaced
                if name is None:
                    name = "emplaced_craters"
            elif craters.lower() == "observed":
                craters = self.observed
                if name is None:
                    name = "observed_craters"
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
        else:
            self.to_vector_file(
                craters=craters,
                interval_number=interval_number,
                driver=driver,
                ask_overwrite=ask_overwrite,
                **kwargs,
            )

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
        interval_number: int | None = None,
        name: str = "craters",
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
            The interval number to append to the file name. If None, then no interval number is added. Default is None.
        name : str, optional
            The name of the layer in the GeoPackage file or the file name if the format does not support layers, by default "craters".
        use_measured_properties : bool, optional
            If True, use the current measured crater properties (semimajor_axis, semiminor_axis, location, orientation) instead of the initial ones, by default True.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """

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
            }
            for long, short in alt_names.items():
                if long in key:
                    key = key.replace(long, short)
            return key[:10].upper()

        # Map of OGR drivers to file extensions
        driver_to_extension_map = {
            "PCIDSK": "pix",
            "PDS4": "xml",
            "PDF": "pdf",
            "MBTiles": "mbtiles",
            "ESRI Shapefile": "shp",
            "MapInfo File": "tab",
            "S57": "000",
            "DGN": "dgn",
            "GML": "gml",
            "GPX": "gpx",
            "KML": "kml",
            "GeoJSON": "json",
            "GeoJSONSeq": "geojsonl",
            "OGR_GMT": "gmt",
            "GPKG": "gpkg",
            "SQLite": "sqlite",
            "WAsP": "map",
            "OpenFileGDB": "gdb",
            "DXF": "dxf",
            "FlatGeobuf": "fgb",
            "PGDUMP": "sql",
            "GPSBabel": "mps",
            "JML": "jml",
            "VDV": "txt",
            "MVT": "mvt",
            "PMTiles": "pmtiles",
            "JSONFG": "json",
            "MiraMonVector": "pol",
        }

        if driver in driver_to_extension_map:
            file_extension = driver_to_extension_map[driver]
        else:
            raise ValueError("Cannot infer file extension from driver {driver}.")

        if file_extension == "shp":
            format_has_layers = False
            split_antimeridian = True
        else:
            format_has_layers = True
            split_antimeridian = True

        surface = self.surface

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
            if driver.upper() == "SHP":
                attrs_df.rename(mapper=shp_key_fix, axis=1, inplace=True)

            gdf = gpd.GeoDataFrame(data=attrs_df, geometry=geoms, crs=surface.crs)
            if format_has_layers:
                if interval_number is None:
                    output_file = self.output_dir / f"craters.{file_extension}"
                else:
                    output_file = self.output_dir / f"craters{interval_number:06d}.{file_extension}"
                print(f"Saving {name} layer to vector file: '{output_file}'...")
            else:
                if interval_number is None:
                    output_file = self.output_dir / f"{name}.{file_extension}"
                else:
                    output_file = self.output_dir / f"{name}{interval_number:06d}.{file_extension}"
                if ask_overwrite and not self._overwrite_check(output_file):
                    return
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
        interval_number: int | None = None,
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
            The interval number to export. If None, all intervals currently saved will be exported. Default is None.
        name : str, optional
            The name used for the file name, by default "craters".
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        from vtk import vtkXMLPolyDataWriter

        if interval_number is not None:
            output_file = self.output_dir / f"{name}{interval_number:06d}.vtp"
        else:
            output_file = self.output_dir / f"{name}.vtp"
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
        interval_number: int | None = None,
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
            The interval number to export. If None, all intervals currently saved will be exported. Default is None.
        ask_overwrite : bool, optional
            If True, prompt the user for confirmation before overwriting files. Default is True.
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        import csv

        if craters:
            if interval_number is None:
                output_file = self.output_dir / f"{name}.csv"
            else:
                output_file = self.output_dir / f"{name}{interval_number:06d}.csv"
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
