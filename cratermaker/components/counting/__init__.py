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

from cratermaker.components.crater import Crater
from cratermaker.core.base import ComponentBase, import_components

if TYPE_CHECKING:
    from cratermaker.components.surface import LocalSurface, Surface

_TALLY_ID = "id"
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
            name=_TALLY_ID,
            uxgrid=self.surface.uxgrid,
        )

        self.surface._uxds[_TALLY_ID] = uxda
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
        if _TALLY_ID not in self.surface.uxds:
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
                if np.any(self.surface.uxds[_TALLY_ID].isel(layer=i).data[count_region.face_indices] > 0):
                    # Gather the unique id values for the current layer
                    unique_ids = np.unique(self.surface.uxds[_TALLY_ID].data[count_region.face_indices, i])
                    removes = [
                        id for id, v in self.observed.items() if v.id in unique_ids and v.final_diameter < crater.final_diameter
                    ]

                    # For every id that appears in the removes list, set it to 0 in the data array
                    if removes:
                        data = self.surface.uxds[_TALLY_ID].data[count_region.face_indices, :]
                        for remove in removes:
                            data[data == remove] = 0
                        self.surface.uxds[_TALLY_ID].data[count_region.face_indices, :] = data
                if insert_layer == -1 and np.all(self.surface.uxds[_TALLY_ID].isel(layer=i).data[count_region.face_indices] == 0):
                    insert_layer = i
            if insert_layer == -1:
                raise ValueError("Crater counting layers are full")
            data = self.surface.uxds[_TALLY_ID].data[count_region.face_indices, :]
            data[:, insert_layer] = crater.id
            self.surface.uxds[_TALLY_ID].data[count_region.face_indices, :] = data

        return

    def _fit_rim_one(self, crater, score_quantile=0.99, distmult=1.0, gradmult=1.0, curvmult=1.0, heightmult=1.0):
        pass

    def fit_rim(self, crater: Crater, tol=0.001, nloops=10, score_quantile=0.95, fit_center=False) -> Crater:
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

        Returns
        -------
        Crater
            A new Crater object with updated rim parameters.
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
        )

        crater_fit = Crater.maker(
            crater,
            measured_semimajor_axis=ap,
            measured_semiminor_axis=bp,
            measured_orientation=orientation,
            measured_location=location,
        )

        return crater_fit

    def score_rim(self, crater: Crater, quantile=0.95, distmult=1.0, gradmult=1.0, curvmult=1.0, heightmult=1.0) -> None:
        """
        Score the rim region of a crater on the surface.

        Parameters
        ----------
        crater : Crater
            The crater for which to score the rim region.
        quantile : float, optional
            The quantile of rim scores to consider. Default is 0.95.
        distmult : float, optional
            Distance multiplier for scoring. Default is 1.0.
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

        region = counting_bindings.score_rim(self.surface, crater, quantile, distmult, gradmult, curvmult, heightmult)

        return region

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
        Dictionary of observed craters on the surface keyed to the crater id.
        """
        return self._observed

    @property
    def emplaced(self) -> list[Crater]:
        """
        List of craters that have been emplaced in the simulation in the current interval in chronological order.
        """
        return self._emplaced

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
            d = xr.Dataset(data_vars=d).set_coords(_TALLY_ID).expand_dims(dim=_TALLY_ID)
            d[_TALLY_ID].attrs["long_name"] = _TALLY_LONG_NAME
            new_data.append(d)

        return xr.concat(new_data, dim=_TALLY_ID)

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
                    combined_data = xr.concat([ds, dsnew], dim=_TALLY_ID)
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

    def export(self, interval_number: int = 0, driver: str = "GPKG", **kwargs: Any) -> None:
        """
        Dump the crater lists to a file and reset the emplaced crater list.

        Parameters
        ----------
        interval_number : int, default=0
            The interval number for the output file naming.
        driver : str, default='GPKG'
            The vector file format to save. Supported formats are 'VTK', 'GPKG', 'ESRI Shapefile', etc.
        **kwargs : Any
            Additional keyword arguments to pass to the make_vector_file function.
        """
        # for crater_type in ["emplaced", "observed"]:
        #     filename = (
        #         self.output_dir / f"{crater_type}_{self._output_file_prefix}{interval_number:06d}.{self._output_file_extension}"
        #     )
        #     if filename.exists():
        #         with xr.open_dataset(filename) as ds:
        #             self.make_vector_file(ds, interval_number, layer_name=f"{crater_type}_craters", format=format, **kwargs)
        pass
        # return

    @staticmethod
    def _compute_geometric_distances_to_ellipse(x, y, coeffs):
        """
        Compute the geometric distances from points (x, y) to the ellipse defined by the coefficients.

        Parameters
        ----------
        x : array_like
            1D array of x-coordinates of data points.
        y : array_like
            1D array of y-coordinates of data points.
        coeffs : array_like
            Coefficients [a, b, c, d, e, f] defining the ellipse.

        Returns
        -------
        delta : array_like
            1D array of geometric distances from each point to the ellipse.
        """
        a, b, c, d, e, f = coeffs
        F = a * x**2 + b * x * y + c * y**2 + d * x + e * y + f

        # gradient magnitude (to convert algebraic to geometric residual)
        Fx = 2 * a * x + b * y + d
        Fy = b * x + 2 * c * y + e
        gradnorm = np.hypot(Fx, Fy)

        # geometric distances
        delta = F / gradnorm
        return delta

    @staticmethod
    def cart_to_pol(coeffs):
        """

        Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.

        The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the ellipse centre; (ap, bp) are the semi-major and semi-minor axes, respectively; e is the eccentricity; and phi is the rotation of the semi- major axis from the x-axis.

        Adapted from: https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/

        """
        # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
        # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
        # Therefore, rename and scale b, d and f appropriately.
        a = coeffs[0]
        b = coeffs[1] / 2
        c = coeffs[2]
        d = coeffs[3] / 2
        f = coeffs[4] / 2
        g = coeffs[5]

        den = b**2 - a * c
        if den > 0:
            raise ValueError("coeffs do not represent an ellipse: b^2 - 4ac must be negative!")

        # The location of the ellipse centre.
        x0, y0 = (c * d - b * f) / den, (a * f - b * d) / den

        num = 2 * (a * f**2 + c * d**2 + g * b**2 - 2 * b * d * f - a * c * g)
        fac = np.sqrt((a - c) ** 2 + 4 * b**2)
        # The semi-major and semi-minor axis lengths (these are not sorted).
        ap = np.sqrt(num / den / (fac - a - c))
        bp = np.sqrt(num / den / (-fac - a - c))

        # Sort the semi-major and semi-minor axis lengths but keep track of
        # the original relative magnitudes of width and height.
        width_gt_height = True
        if ap < bp:
            width_gt_height = False
            ap, bp = bp, ap

        # The eccentricity.
        r = (bp / ap) ** 2
        if r > 1:
            r = 1 / r
        e = np.sqrt(1 - r)

        # The angle of anticlockwise rotation of the major-axis from x-axis.
        if b == 0:
            phi = 0 if a < c else np.pi / 2
        else:
            phi = np.arctan((2.0 * b) / (a - c)) / 2
            if a > c:
                phi += np.pi / 2
        if not width_gt_height:
            # Ensure that phi is the angle to rotate to the semi-major axis.
            phi += np.pi / 2
        phi = phi % np.pi

        return x0, y0, ap, bp, e, phi

    def pol_to_coeff(x0, y0, ap, bp, e, phi):
        """
        Convert from ellipse parameters to cartesian conic coefficients.
        """
        s2 = np.sin(phi) ** 2
        c2 = np.cos(phi) ** 2
        sc = np.sin(phi) * np.cos(phi)
        a2 = ap**2
        b2 = bp**2

        a = a2 * s2 + b2 * c2
        b = 2 * (b2 - a2) * sc
        c = a2 * c2 + b2 * s2
        d = -2 * a * x0 - b * y0
        e = -b * x0 - 2 * c * y0
        f = a * x0**2 + b * x0 * y0 + c * y0**2 - a2 * b2

        return (a, b, c, d, e, f)

    def to_vector_file(
        self,
        crater_data: xr.Dataset | dict[int, Crater] | list[Crater],
        driver: str = "GPKG",
        interval_number: int | None = None,
        layer_name: str = "craters",
        **kwargs,
    ) -> None:
        """
        Export the crater data to a vector file and stores it in the default export directory.

        Notes: In order for the crater and surface to be synced up when saving to VTK/VTP format, the initial conditions (no craters) must be saved. Otherwise, saving to file only occurs if there are craters in the crater_data.

        Parameters
        ----------
        crater_data : dict
            Dictionary containing crater attributes. Must include 'final_diameter', 'longitude', and 'latitude' keys. Any additional key value pairs will be added as attributes to each crater.
        interval_number : int, optional
            The interval number to export. If None, all intervals currently saved will be exported. Default is None.
        layer_name : str, optional
            The name of the layer in the GeoPackage file, by default "craters".
        **kwargs : Any
            Additional keyword arguments that are ignored.
        """
        from vtk import (
            vtkCellArray,
            vtkPoints,
            vtkPolyData,
            vtkPolyLine,
            vtkXMLPolyDataWriter,
        )

        def lonlat_to_xyz(R):
            def _f(lon_deg, lat_deg, z=None):
                lon = np.deg2rad(lon_deg)
                lat = np.deg2rad(lat_deg)
                X = R * np.cos(lat) * np.cos(lon)
                Y = R * np.cos(lat) * np.sin(lon)
                Z = R * np.sin(lat)
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

        if isinstance(crater_data, (dict | list)):
            crater_data = self.to_xarray(crater_data)

        if format == "shp":
            format_has_layers = False
            split_antimeridian = True
        elif format == "vtp" or format == "vtk":
            format_has_layers = False
            split_antimeridian = False
        else:
            format_has_layers = True
            split_antimeridian = True

        surface = self.surface

        geoms = []
        attrs = []
        if "final_diameter" in crater_data and "longitude" in crater_data and "latitude" in crater_data:
            for id in crater_data.id:
                crater = crater_data.sel(id=[id])
                lon = float(crater["longitude"])
                lat = float(crater["latitude"])
                radius = float(crater["final_diameter"]) / 2.0
                poly = self.geodesic_ellipse_polygon(
                    lon,
                    lat,
                    a=radius,
                    b=radius,
                    R_pole=surface.radius,
                    R_equator=surface.radius,
                    split_antimeridian=split_antimeridian,
                )
                if isinstance(poly, GeometryCollection):
                    for p in poly.geoms:
                        geoms.append(p)
                        attrs.append(crater.to_dataframe())
                else:
                    geoms.append(poly)
                    attrs.append(crater.to_dataframe())

        if format == "vtp" or format == "vtk":
            output_file = self.output_dir / f"{layer_name}{interval_number:06d}.vtp"

            points = vtkPoints()
            lines = vtkCellArray()
            point_id = 0  # Keep track of the point ID across all circles
            for poly in geoms:
                for ring_xyz in polygon_xyz_coords(poly, surface.radius):
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

            # Write the poly_data to a VTK file
            writer = vtkXMLPolyDataWriter()
            writer.SetFileName(output_file)
            writer.SetInputData(poly_data)

            # Optional: set the data mode to binary to save disk space
            writer.SetDataModeToBinary()
            writer.Write()
        elif len(geoms) > 0:
            attrs_df = pd.concat(attrs, ignore_index=True)
            if format == "shp":
                attrs_df.rename(mapper=shp_key_fix, axis=1, inplace=True)

            gdf = gpd.GeoDataFrame(data=attrs_df, geometry=geoms, crs=surface.crs)
            try:
                if format_has_layers:
                    output_file = self.output_dir / f"craters{interval_number:06d}.{format}"
                    print(f"Saving {layer_name} layer to vector file: '{output_file}'...")
                    gdf.to_file(output_file, layer=layer_name)
                else:
                    output_file = self.output_dir / f"{layer_name}{interval_number:06d}.{format}"
                    print(f"Saving to vector file: '{output_file}'...")
                    gdf.to_file(output_file)
            except Exception as e:
                raise RuntimeError(f"Error saving {output_file}: {e}") from e

        return


import_components(__name__, __path__)
