from __future__ import annotations

import csv
from abc import abstractmethod
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import uxarray as uxr

from cratermaker.constants import FloatLike
from cratermaker.core.base import ComponentBase, import_components
from cratermaker.core.crater import Crater
from cratermaker.utils.general_utils import parameter

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
    vector_format: str | None, optional
        The format of the output file used for the vector representation of the craters. By default, no vector output file is saved. If set, the value will be used as the file extension, and geopandas.to_file will attempt to infer the driver from this. By default, no additional vector The recommended options are "gpkg" (GeoPackage) or "shp" (ESRI shape file). Other drivers may require additional kwargs.
    **kwargs : Any
        Additional keyword arguments.
    """

    def __init__(
        self,
        surface: Surface | LocalSurface,
        vector_format: str | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_emplaced", [])
        object.__setattr__(self, "_observed", {})
        object.__setattr__(self, "_vector_format", None)
        self.surface = surface
        rng = kwargs.pop("rng", surface.rng)
        simdir = kwargs.pop("simdir", surface.simdir)
        super().__init__(rng=rng, simdir=simdir, **kwargs)
        self.vector_format = vector_format

    @classmethod
    def maker(
        cls,
        counting: str | Counting | None = None,
        surface: Surface | LocalSurface | None = None,
        vector_format: str | None = None,
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
        vector_format: str | None, optional
            The format of the output file used for the vector representation of the craters. By default, no vector output file is saved. If set, the value will be used as the file extension, and geopandas.to_file will attempt to infer the driver from this. By default, no additional vector The recommended options are "gpkg" (GeoPackage) or "shp" (ESRI shape file). Other drivers may require additional kwargs.
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
            vector_format=vector_format,
            **kwargs,
        )

        return counting

    def __str__(self) -> str:
        base = super().__str__()
        if self.vector_format is not None:
            base += f"\nOutput vector format: {self.vector_format}"
        return f"{base}\nSurface: {self.surface.name}"

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
        self._emplaced = []
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
    def emplaced(self):
        """
        The list of craters that have been emplaced in the simulation.
        """
        return self._emplaced

    @property
    def output_dir(self) -> Path | None:
        """
        The output directory for the surface. If None, the surface does not have an output directory set.
        """
        if self._output_dir is None:
            self._output_dir = self.simdir / self.__class__._CRATER_DIR
        return self._output_dir

    def save(self, interval_number: int = 0, **kwargs: Any) -> None:
        """
        Dump the crater lists to a file and reset the emplaced crater list.

        Parameters
        ----------
        interval_number : int, default=0
            The interval number for the output file naming.
        **kwargs : Any
            Additional keyword arguments to pass to the save_vector function.
        """
        crater_dir = self.output_dir
        crater_dir.mkdir(parents=True, exist_ok=True)
        emplaced_filename = crater_dir / f"emplaced_craters{interval_number:06d}.csv"
        observed_filename = crater_dir / f"observed_craters{interval_number:06d}.csv"

        def _convert_and_merge(craters: dict[int, Crater] | list[Crater], filename: Path | str, layer_name: str) -> None:
            new_data = []
            if isinstance(craters, dict):
                craters = craters.values()
            for c in craters:
                d = asdict(c)
                if "location" in d:
                    lon, lat = d.pop("location")
                    d["longitude"] = lon
                    d["latitude"] = lat
                new_data.append(d)

            # If the file already exists, read it and merge
            if filename.exists():
                with filename.open("r", newline="") as f:
                    reader = csv.DictReader(f)
                existing_data = list(reader)
                combined_data = existing_data + new_data
                # Sort by final_diameter descending
                combined_data = sorted(combined_data, key=lambda d: -float(d["final_diameter"]))
            else:
                combined_data = new_data

            # Write merged data back to file
            if combined_data:
                with filename.open("w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=combined_data[0].keys())
                    writer.writeheader()
                    writer.writerows(combined_data)
                if self.vector_format is not None:
                    self.save_vector(combined_data, interval_number, layer_name=layer_name, **kwargs)

        _convert_and_merge(self.emplaced, emplaced_filename, "emplaced_craters")
        _convert_and_merge(self.observed, observed_filename, "observed_craters")
        self._emplaced = []
        return

    def save_vector(
        self,
        crater_data: dict,
        interval_number: int = 0,
        layer_name: str = "craters",
        **kwargs,
    ) -> None:
        """
        Export the crater data to a vector file and stores it in the default export directory.

        Parameters
        ----------
        crater_list : dict
            Dictionary containing crater attributes. Must include 'final_diameter', 'longitude', and 'latitude' keys. Any additional key value pairs will be added as attributes to each crater.
        interval_number : int, optional
            The interval number to save, by default 0.
        layer_name : str, optional
            The name of the layer in the GeoPackage file, by default "craters".

        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        import geopandas as gpd
        import pandas as pd
        from pyproj import Geod
        from shapely.geometry import GeometryCollection, LineString, Polygon
        from shapely.ops import split, transform

        def _geodesic_ellipse_polygon(
            lon: float,
            lat: float,
            a: float,
            b: float,
            az_deg: float = 0.0,
            n: int = 150,
            R_pole: FloatLike = 1.0,
            R_equator: FloatLike = 1.0,
        ) -> Polygon:
            """
            Geodesic ellipse on a sphere: for each bearing theta from the center, we shoot a geodesic with distance r(theta) = (a*b)/sqrt((b*ct)^2 + (a*st)^2), then rotate all bearings by az_deg.

            Parameters
            ----------
            lon : float
                Longitude of the ellipse center in degrees.
            lat : float
                Latitude of the ellipse center in degrees.
            a : float
                Semi-major axis in meters.
            b : float
                Semi-minor axis in meters.
            az_deg : float, optional
                Azimuth rotation of the ellipse in degrees clockwise from north, by default 0.0.
            n : int, optional
                Number of points to use for the polygon, by default 150.
            R_pole : FloatLike, optional
                Planetary polar radius in units of meters, by default 1.0.
            R_equator : FloatLike, optional
                Planetary equatorial radius in units of meters, by default 1.0.

            Returns
            -------
            Returns a Shapely Polygon in lon/lat degrees.
            """
            geod = Geod(a=R_pole, b=R_equator)
            theta = np.linspace(0.0, 360.0, num=n, endpoint=False)

            # Polar radius of an axis-aligned ellipse in a Euclidean tangent plane
            ct = np.cos(np.deg2rad(theta))
            st = np.sin(np.deg2rad(theta))
            r = (a * b) / np.sqrt((b * ct) ** 2 + (a * st) ** 2)

            # Bearings (from east, CCW) rotated by azimuth
            bearings = theta + az_deg % 360.0

            # Forward geodesic for each bearing/distance
            poly_lon, poly_lat, _ = geod.fwd(lon * np.ones_like(bearings), lat * np.ones_like(bearings), bearings, r)

            # Correct for potential antimeridian crossing
            if np.ptp(poly_lon) > 180.0:
                center_sign = np.sign(lon)
                poly_lon = np.where(np.sign(poly_lon) != center_sign, poly_lon + 360.0 * center_sign, poly_lon)
                poly = Polygon(zip(poly_lon, poly_lat, strict=False))
                merdian_lon = 180.0 * center_sign
                meridian = LineString([(merdian_lon, -90.0), (merdian_lon, 90.0)])
                poly = split(poly, meridian)

                def lon_flip(lon, lat):
                    lon = np.where(np.abs(lon) >= 180.0, lon - 360.0 * np.sign(lon), lon)
                    return lon, lat

                new_geoms = []
                for p in poly.geoms:
                    if np.abs(p.centroid.x) > 180.0:
                        p = transform(lon_flip, p)
                    new_geoms.append(p)
                poly = GeometryCollection(new_geoms)

            else:
                poly = Polygon(zip(poly_lon, poly_lat, strict=False))

            return poly

        out_dir = self.output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        if "final_diameter" not in crater_data[0] or "longitude" not in crater_data[0] or "latitude" not in crater_data[0]:
            raise ValueError("crater_data must contain 'final_diameter', 'longitude', and 'latitude' keys.")

        geoms = []
        attrs = {}
        surface = self.surface

        def add_attrs(crater: dict, attrs: dict) -> dict:
            for k, v in crater.items():
                if k in attrs:
                    attrs[k].append(v)
                else:
                    attrs[k] = [v]
            return attrs

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

        for crater in crater_data:
            lon = float(crater["longitude"])
            lat = float(crater["latitude"])
            radius = float(crater["final_diameter"]) / 2.0
            poly = _geodesic_ellipse_polygon(lon, lat, a=radius, b=radius, R_pole=surface.radius, R_equator=surface.radius)
            if isinstance(poly, GeometryCollection):
                for p in poly.geoms:
                    geoms.append(p)
                    attrs = add_attrs(crater, attrs)
            else:
                geoms.append(poly)
                attrs = add_attrs(crater, attrs)

        if self.vector_format == "shp":
            attrs = {shp_key_fix(k): v for k, v in attrs.items()}
        attrs_df = pd.DataFrame({k: np.asarray(v) for k, v in attrs.items()})

        gdf = gpd.GeoDataFrame(data=attrs_df, geometry=geoms, crs=surface.crs)
        if self.vector_format == "shp":
            output_file = out_dir / f"{layer_name}{interval_number:06d}.{self.vector_format}"
        else:
            output_file = out_dir / f"craters{interval_number:06d}.{self.vector_format}"
        print(f"Saving vector file: '{output_file}'...")
        try:
            gdf.to_file(output_file, layer=layer_name)
        except Exception as e:
            raise RuntimeError(f"Error saving {output_file}: {e}") from e

        return

    @parameter
    def vector_format(self) -> str | None:
        return self._vector_format

    @vector_format.setter
    def vector_format(self, value) -> None:
        """
        The file extension for the vector representation of the craters that will be generated when the save method is called.
        """
        if not isinstance(value, (str | None)):
            raise TypeError("vector_format must be a str or None")
        self._vector_format = value
        return


import_components(__name__, __path__, ignore_private=True)
