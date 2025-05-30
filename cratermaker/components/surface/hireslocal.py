from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.spatial.transform import Rotation as R

from cratermaker.components.morphology import Morphology
from cratermaker.components.scaling import Scaling
from cratermaker.components.surface import Surface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import (
    format_large_units,
    parameter,
    validate_and_normalize_location,
)


@Surface.register("hireslocal")
class HiResLocalSurface(Surface):
    """
    Create a uniform grid configuration with the given pixel size.

    Parameters
    ----------
    pix : FloatLike
        The approximate face size inside the local region in meters.
    local_radius : FloatLike
        The radius of the local region in meters.
    local_location : PairOfFloats
        The longitude and latitude of the location in degrees.
    superdomain_scale_factor : FloatLike, optional
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters
        that are modeled outside the local region are those whose ejecta could just reach the boundary. If not provide, then it will
        be computed based on a provided (or default) scaling and morphology model
    target : Target, optional
        The target body or name of a known target body for the impact simulation. If none provide, it will be either set to the default,
        or extracted from the scaling model if it is provied
    reset : bool, optional
        Flag to indicate whether to reset the surface. Default is True.
    regrid : bool, optional
        Flag to indicate whether to regrid the surface. Default is False.
    simdir : str | Path
        The main project simulation directory. Defaults to the current working directory if None.

    Returns
    -------
    HiResLocalSurface
        An instance of the HiResLocalSurface clas initialized with the given point distribution
    """

    def __init__(
        self,
        pix: FloatLike,
        local_radius: FloatLike,
        local_location: PairOfFloats,
        superdomain_scale_factor: FloatLike | None = None,
        target: Target | str | None = None,
        reset: bool = False,
        regrid: bool = False,
        simdir: str | Path | None = None,
        **kwargs: Any,
    ):
        object.__setattr__(self, "_local_view", None)
        object.__setattr__(self, "_superdomain_function_slope", None)
        object.__setattr__(self, "_superdomain_function_exponent", None)
        super().__init__(target=target, simdir=simdir, **kwargs)

        self.pix = pix
        self.local_radius = local_radius
        self.local_location = local_location
        if superdomain_scale_factor is not None:
            self.superdomain_scale_factor = superdomain_scale_factor
            self._load_from_files(reset=reset, regrid=regrid, **kwargs)

        return

    def __str__(self) -> str:
        base = super().__str__()
        pix = format_large_units(self.pix, quantity="length")
        pix_min = format_large_units(self.pix_min, quantity="length")
        pix_max = format_large_units(self.pix_max, quantity="length")
        local_radius = format_large_units(self.local_radius, quantity="length")
        return (
            f"{base}\n"
            f"Local pixel size: {pix}\n"
            f"Local Radius: {local_radius}\n"
            f"Local Location: ({self.local_location[0]:.2f}°, {self.local_location[1]:.2f}°)\n"
            f"Minimum effective pixel size: {pix_min}\n"
            f"Maximum effective pixel size: {pix_max}"
        )

    def plot_hillshade(self, imagefile=None, **kwargs: Any) -> None:
        """
        Plot a hillshade image of the local region.

        Parameters
        ----------
        imagefile : str | Path, optional
            The file path to save the hillshade image. If None, the image will be displayed instead of saved.
        **kwargs : Any
            Additional keyword arguments to pass to the plotting function.
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LightSource
        from scipy.interpolate import griddata

        region_view = self.local_view
        local_radius = self.local_radius
        pix = self.pix

        # Dynamically compute image resolution and dpi
        extent_val = local_radius
        resolution = int(2 * extent_val / pix)
        dpi = resolution
        x = np.linspace(-extent_val, extent_val, resolution)
        y = np.linspace(-extent_val, extent_val, resolution)
        grid_x, grid_y = np.meshgrid(x, y)

        # Use polar coordinates from region_view
        r = region_view.face_distance
        theta = region_view.face_bearing
        x_cart = r * np.cos(theta)
        y_cart = r * np.sin(theta)

        points = np.column_stack((x_cart, y_cart))
        values = region_view.face_elevation
        grid_z = griddata(points, values, (grid_x, grid_y), method="linear")

        # Generate hillshade
        azimuth = 300.0
        solar_angle = 20.0
        ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
        hillshade = ls.hillshade(grid_z, dx=pix, dy=pix, fraction=1.0)

        # Plot hillshade with (1, 1) inch figure and dpi=resolution for exact pixel size
        fig, ax = plt.subplots(figsize=(1, 1), dpi=dpi, frameon=False)
        ax.imshow(
            hillshade,
            interpolation="nearest",
            cmap="gray",
            vmin=0.0,
            vmax=1.0,
            aspect="equal",
            extent=(-extent_val, extent_val, -extent_val, extent_val),
        )
        ax.axis("off")
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        if imagefile:
            plt.savefig(imagefile, bbox_inches="tight", pad_inches=0, dpi=dpi, **kwargs)
        else:
            plt.show(**kwargs)
        plt.close(fig)
        return

    def superdomain_function(self, r):
        """
        A function that defines the superdomain scale factor based on the distance from the local boundary.
        This is set so that smallest craters that are modeled outside the local region are those whose ejecta could just reach the boundary.
        It is a piecewise function that returns the local pixel size inside the local region and a power law function outside. The slope and exponent of the power law is linear if superdomain_scale_factor is set explicitly, otherwise it is computed by the `_set_superdomain` method.

        Parameters
        ----------
        r : FloatLike
            The distance from the local center in meters.
        Returns
        -------
        FloatLike
            The effective pixel size at the given distance from the local center.
        """

        return np.where(
            r < self.local_radius,
            self.pix,
            self.pix
            + self.superdomain_function_slope
            * (r - self.local_radius) ** self.superdomain_function_exponent,
        )

    def _set_superdomain(
        self,
        scaling: Scaling | str | None = None,
        morphology: Morphology | str | None = None,
        reset: bool = False,
        regrid: bool = False,
        **kwargs: Any,
    ):
        """
        Set the superdomain scale factor based on the scaling and morphology models.
        This sets the cell size at the antipode such that ejecta from a crater of that size just reaches the local region.

        Parameters
        ----------
        scaling : Scaling | str | None, optional
            The scaling model to use. If None, the default scaling model will be used.
        morphology : Morphology | str | None, optional
            The morphology model to use. If None, the default morphology model will be used.
        reset : bool, optional
            Flag to indicate whether to reset the surface. Default is False.
        regrid : bool, optional
            Flag to indicate whether to regrid the surface. Default is False.
        **kwargs : Any
            Additional keyword arguments to pass to the scaling and morphology models.
        """
        from scipy.optimize import curve_fit

        from cratermaker import Crater

        antipode_distance = np.pi * self.target.radius
        projectile_velocity = scaling.projectile.mean_velocity * 10

        distance = 1.0
        dvals = []
        sdvals = []
        superdomain_size = distance * 0.1
        while distance < antipode_distance:
            for final_diameter in np.logspace(
                np.log10(superdomain_size),
                np.log10(self.target.radius * 2),
                1000,
            ):
                crater = Crater.maker(
                    final_diameter=final_diameter,
                    angle=90.0,
                    scaling=scaling,
                    projectile_velocity=projectile_velocity,
                )
                rmax = morphology.rmax(crater=crater, minimum_thickness=1e-3)
                if rmax >= distance:
                    superdomain_size = crater.final_diameter
                    break
            dvals.append(distance)
            sdvals.append(superdomain_size)
            distance = distance + superdomain_size

        def plaw(x, a, b):
            return a * x**b

        try:
            popt = curve_fit(plaw, dvals, sdvals, p0=[1.0, 1.0])[0]
            self._superdomain_function_slope = popt[0]
            self._superdomain_function_exponent = popt[1]
        except Exception:
            print("Could not fit superdomain function, using default values.")
            self._superdomain_function_slope = 1.0
            self._superdomain_function_exponent = 1.0
        self._superdomain_scale_factor = self.superdomain_function(antipode_distance)

        self._load_from_files(
            reset=reset, regrid=regrid, scaling=scaling, morphology=morphology, **kwargs
        )
        return

    def _rotate_point_cloud(self, points):
        """
        Rotate a point cloud so that the point at [0,0,1] moves to (lon,lat)
        using the convention:
        - Longitude [-180, 180] degrees, increasing eastward
        - Latitude [-90, 90] degrees, increasing northward
        - (0,0) lon/lat corresponds to [0,0,1] in Cartesian space

        Parameters:
            points (np.ndarray): Nx3 array of (x,y,z) points.

        Returns:
            np.ndarray: Rotated Nx3 point cloud.
        """

        # Convert target lon, lat to radians
        lon_rad, lat_rad = (
            np.radians(self.local_location[0]),
            np.radians(self.local_location[1]),
        )

        # Compute target unit vector
        target = np.array(
            [
                np.cos(lat_rad) * np.cos(lon_rad),
                np.cos(lat_rad) * np.sin(lon_rad),
                np.sin(lat_rad),
            ]
        )

        # Original vector is the north pole
        original = np.array([0, 0, 1])

        axis = np.cross(original, target)
        axis_norm = np.linalg.norm(axis)

        if axis_norm < 1e-10:
            return points

        axis /= axis_norm
        angle = np.arccos(np.clip(np.dot(original, target), -1.0, 1.0))

        rotvec = axis * angle
        rotation = R.from_rotvec(rotvec)
        return rotation.apply(points)

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """

        def _interior_distribution(theta):
            arc_distance = self.radius * theta
            delta = self.pix / self.radius
            if arc_distance < self.pix:
                return [], theta + delta
            if np.pi * self.radius - arc_distance < self.pix:
                return [], np.pi

            max_extent = arc_distance / self.radius
            coords = np.arange(-max_extent, max_extent + delta, delta)
            lat_grid, lon_grid = np.meshgrid(coords, coords, indexing="ij")

            pix_distance_sq = (lon_grid / delta) ** 2 + (lat_grid / delta) ** 2
            inner = (arc_distance / self.pix - 1.0) ** 2
            outer = (arc_distance / self.pix) ** 2
            mask = (pix_distance_sq >= inner) & (pix_distance_sq < outer)

            lat = lat_grid[mask]
            lon = lon_grid[mask]

            cos_lat = np.cos(lat)
            sin_lat = np.sin(lat)
            cos_lon = np.cos(lon)
            sin_lon = np.sin(lon)

            z = cos_lat * cos_lon
            x = cos_lat * sin_lon
            y = sin_lat

            points = np.column_stack([x, y, z])
            return points.tolist(), theta + delta

        def _exterior_distribution(theta):
            r = self.radius
            arc_distance = r * theta
            pix_local = self.superdomain_function(arc_distance)
            if (np.pi * r - arc_distance) < pix_local:
                return [], np.pi

            theta_next = theta + pix_local / r

            sin_theta = np.sin(theta)
            cos_theta = np.cos(theta)

            n_phi = max(1, int(round(2 * np.pi * r * sin_theta / pix_local)))
            phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)

            x = sin_theta * np.cos(phi)
            y = sin_theta * np.sin(phi)
            z = np.full_like(phi, cos_theta)

            points = np.column_stack((x, y, z))
            return points.tolist(), theta_next

        print(f"Center of local region: {self.local_location}")
        print(
            f"Size of local region: {format_large_units(self.local_radius, quantity='length')}"
        )
        print(
            f"Local region pixel size: {format_large_units(self.pix, quantity='length')}"
        )

        interior_points = []
        theta = 0.0
        while theta <= self.local_radius / self.radius:
            new_points, theta = _interior_distribution(theta)
            interior_points.extend(new_points)
        print(f"Generated {len(interior_points)} points in the local region.")

        exterior_points = []
        while theta < np.pi:
            new_points, theta = _exterior_distribution(theta)
            exterior_points.extend(new_points)
        print(f"Generated {len(exterior_points)} points in the superdomain region.")
        points = interior_points + exterior_points
        decimals = max(6, -int(np.floor(np.log10(self.pix / self.radius))))

        points = np.array(points, dtype=np.float64)
        points = np.round(points, decimals=decimals)
        points = np.unique(points, axis=0)
        points = self._rotate_point_cloud(
            points
        ).T  # rotates from the north pole to local_location

        return points

    def _save_to_files(
        self,
        combine_data_files: bool = False,
        interval_number: int = 0,
        time_variables: dict | None = None,
        **kwargs,
    ) -> None:
        """
        Save the surface data to the specified directory. Each data variable is saved to a separate NetCDF file. If 'time_variables' is specified, then a one or more variables will be added to the dataset along the time dimension. If 'interval_number' is included as a key in `time_variables`, then this will be appended to the data file name.

        Parameters
        ----------
        combine_data_files : bool, optional
            If True, combine all data variables into a single NetCDF file, otherwise each variable will be saved to its own NetCDF file. Default is False.
        interval_number : int, optional
            Interval number to append to the data file name. Default is 0.
        time_variables : dict, optional
            Dictionary containing one or more variable name and value pairs. These will be added to the dataset along the time dimension. Default is None.
        """

        super()._save_to_files(
            combine_data_files=combine_data_files,
            interval_number=interval_number,
            time_variables=time_variables,
            **kwargs,
        )
        imgdir = Path(self.simdir) / "surface_images"
        imgdir.mkdir(parents=True, exist_ok=True)
        imagefile = imgdir / f"hillshade{interval_number:06d}.png"
        self.plot_hillshade(imagefile=imagefile, **kwargs)
        return

    @property
    def local_view(self):
        """
        Returns the local view of the surface.
        """
        if self._local_view is None:
            self._local_view = self.extract_region(
                location=self.local_location, region_radius=self.local_radius
            )
        return self._local_view

    @parameter
    def pix(self):
        """
        The approximate face size for a cell of the mesh.
        """
        return self._pix

    @pix.setter
    def pix(self, value: FloatLike):
        if (
            not isinstance(value, FloatLike)
            or np.isnan(value)
            or np.isinf(value)
            or value <= 0
        ):
            raise TypeError("pix must be a positive float")
        self._pix = value

    @parameter
    def local_radius(self):
        """
        The radius of the local region in meters.
        """
        return self._local_radius

    @local_radius.setter
    def local_radius(self, value: FloatLike):
        if (
            not isinstance(value, FloatLike)
            or np.isnan(value)
            or np.isinf(value)
            or value <= 0
        ):
            raise TypeError("local_radius must be a positive float")
        if value > np.pi * self.radius:
            raise ValueError(
                "local_radius must be less than pi * radius of the target body"
            )
        if value < self.pix:
            raise ValueError(
                "local_radius must be greater than or equal to pix (the approximate face size in the local region"
            )
        self._local_radius = value

    @parameter
    def local_location(self):
        """
        The longitude and latitude of the location in degrees.
        """
        return self._local_location

    @local_location.setter
    def local_location(self, value: PairOfFloats):
        if not isinstance(value, tuple) or len(value) != 2:
            raise TypeError("local_location must be a tuple of two floats")
        self._local_location = validate_and_normalize_location(value)

    @property
    def superdomain_function_slope(self):
        if self._superdomain_function_slope is None:
            d0 = self.local_radius
            d1 = np.pi * self.radius - self.local_radius
            pix_lo = self.pix
            pix_hi = self.pix * self.superdomain_scale_factor
            self._superdomain_function_slope = (pix_hi - pix_lo) / (d1 - d0)
        return self._superdomain_function_slope

    @property
    def superdomain_function_exponent(self):
        if self._superdomain_function_exponent is None:
            self._superdomain_function_exponent = 1.0
        return self._superdomain_function_exponent

    @parameter
    def superdomain_scale_factor(self):
        """
        A factor defining the ratio of cell size to the distance from the local boundary. This is set so that smallest craters that are
        modeled outside the local region are those whose ejecta could just reach the boundary.
        """
        return self._superdomain_scale_factor

    @superdomain_scale_factor.setter
    def superdomain_scale_factor(self, value: FloatLike):
        if (
            not isinstance(value, FloatLike)
            or np.isnan(value)
            or np.isinf(value)
            or value < 1.0
        ):
            raise TypeError(
                "superdomain_scale_factor must be a positive float greater than or equal to 1"
            )
        self._superdomain_scale_factor = value

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [
            self._component_name,
            self.radius,
            self.pix,
            self.local_radius,
            self.local_location,
            self.superdomain_scale_factor,
        ]
