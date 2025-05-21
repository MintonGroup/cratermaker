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
from cratermaker.core.crater import Crater
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
        super().__init__(target=target, simdir=simdir, **kwargs)
        self.pix = pix
        self.local_radius = local_radius
        self.local_location = local_location
        if superdomain_scale_factor is not None:
            self.superdomain_scale_factor = superdomain_scale_factor
            self.load_from_files(reset=reset, regrid=regrid, **kwargs)

        return

    def __str__(self) -> str:
        base = super().__str__()
        pix = format_large_units(self.pix, quantity="length")
        local_radius = format_large_units(self.local_radius, quantity="length")
        return (
            f"{base}\n"
            f"Pixel Size: {pix}\n"
            f"Local Radius: {local_radius}\n"
            f"Local Location: ({self.local_location[0]:.2f}°, {self.local_location[1]:.2f}°)\n"
            f"Superdomain Scale Factor: {self.superdomain_scale_factor:.2f}"
        )

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

    def set_superdomain(
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
        scaling = Scaling.maker(scaling, target=self.target, **kwargs)
        morphology = Morphology.maker(morphology, target=self.target, **kwargs)
        projectile_velocity = scaling.projectile.mean_velocity * 10

        antipode_distance = np.pi * self.target.radius - self.local_radius

        for final_diameter in np.logspace(
            np.log10(self.target.radius / 1e6),
            np.log10(self.target.radius * 2),
            1000,
        ):
            crater = Crater.maker(
                final_diameter=final_diameter,
                angle=90.0,
                scaling=scaling,
                projectile_velocity=projectile_velocity,
                **kwargs,
            )
            rmax = morphology.rmax(crater=crater, minimum_thickness=1e-3)
            if rmax >= antipode_distance:
                superdomain_scale_factor = crater.final_diameter / self.pix
                break

        self.superdomain_scale_factor = superdomain_scale_factor
        self.load_from_files(reset=reset, regrid=regrid, **kwargs)
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

    def generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """

        print("Generating a mesh with variable resolution faces")
        print(f"Center of local region: {self.local_location}")
        print(
            f"Size of local region: {format_large_units(self.local_radius, quantity='length')}"
        )
        print(
            f"Hires region pixel size: {format_large_units(self.pix, quantity='length')}"
        )
        print(
            f"Lores region pixel size: {format_large_units(self.pix * self.superdomain_scale_factor, quantity='length')}"
        )

        def _pix_func(distance):
            d0 = self.local_radius
            d1 = np.pi * self.radius - self.local_radius
            pix_lo = self.pix
            pix_hi = self.pix * self.superdomain_scale_factor
            slope = (pix_hi - pix_lo) / (d1 - d0)
            return np.where(
                distance <= d0,
                pix_lo,
                pix_lo + slope * (distance - d0),
            )

        points = []
        r = self.radius
        dtheta = 0.5 * self.pix / r
        theta = dtheta

        while theta < np.pi:
            arc_distance = r * theta
            pix_local = _pix_func(arc_distance)
            # Stop if the distance to the antipode is less than pix_local
            if (np.pi * r - arc_distance) < pix_local:
                break
            dphi = (
                pix_local / (r * np.sin(theta)) if np.sin(theta) > 1e-6 else 2 * np.pi
            )
            n_phi = max(1, int(round(2 * np.pi / dphi)))
            for i in range(n_phi):
                phi = i * 2 * np.pi / n_phi
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                points.append([x, y, z])
            theta += 0.5 * pix_local / r

        points = np.round(points, decimals=5)
        points = np.unique(np.array(points), axis=0)
        points = np.array(points).T
        points = self._rotate_point_cloud(points.T).T  # rotates from north pole
        return points

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
        if value > 2 * np.pi * self.radius:
            raise ValueError(
                "local_radius must be less than 2 * pi * radius of the target body"
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
