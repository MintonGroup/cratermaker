from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cratermaker.components.surface import Surface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike, PairOfFloats
from cratermaker.utils.general_utils import format_large_units, parameter


@Surface.register("arbitrary_resolution")
class ArbitraryResolutionSurface(Surface):
    """
    Create a uniform grid configuration with an arbitrary user-defined pixel size. This will not be as nice as the regular IcosphereSurface, but can be any resolution desired.

    Parameters
    ----------
    pix : float
        The approximate face size for the mesh in meters.
    target : Target, optional
        The target body or name of a known target body for the impact simulation.
    reset : bool, optional
        Flag to indicate whether to reset the surface. Default is True.
    regrid : bool, optional
        Flag to indicate whether to regrid the surface. Default is False.
    simdir : str | Path
        The main project simulation directory. Defaults to the current working directory if None.

    Returns
    -------
    ArbitraryResolutionSurface
        An instance of the ArbitraryResolutionSurface class initialized with the given pixel size.
    """

    def __init__(
        self,
        pix: FloatLike | None = None,
        target: Target | str | None = None,
        reset: bool = False,
        regrid: bool = False,
        simdir: str | Path | None = None,
        **kwargs: Any,
    ):
        super().__init__(target=target, simdir=simdir, **kwargs)
        self.pix = pix
        self._load_from_files(reset=reset, regrid=regrid, **kwargs)

    def __str__(self) -> str:
        base = super().__str__()
        pix_mean = format_large_units(self.pix_mean, quantity="length")
        pix_std = format_large_units(self.pix_std, quantity="length")
        return f"{base}\nEffective pixel size: {pix_mean} +/- {pix_std}"

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self.radius, self.pix]

    @parameter
    def pix(self):
        """
        The approximate face size for a cell of the mesh.
        """
        return self._pix

    @pix.setter
    def pix(self, value: FloatLike):
        if value is None:
            value = (
                np.sqrt(4 * np.pi * self.radius**2) * 1e-3
            )  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid
        elif not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
            raise TypeError("pix must be a positive float")
        self._pix = float(value)

    @staticmethod
    def _distribute_points(
        distance: FloatLike,
        radius: FloatLike = 1.0,
        lon_range: PairOfFloats = (-180, 180),
        lat_range: PairOfFloats = (-90, 90),
    ) -> NDArray:
        """
        Distributes points on a sphere using Deserno's algorithm (Deserno 2004).

        Parameters
        ----------
        distance : float
            Approximate distance between points, used to determine the number of points, where n = 1/distance**2 when distributed over the whole sphere.
        radius : float, optional
            Radius of the sphere. Default is 1.0
        lon_range : tuple, optional
            Range of longitudes in degrees. Default is (-180,180).
        lat_range : tuple, optional
            Range of latitudes in degrees. Default is (-90,90).

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of cartesian points on the sphere.

        References
        ----------
        - Deserno, Markus., 2004. How to generate equidistributed points on the surface of a sphere. https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf

        """

        def _sph2cart(theta, phi, r):
            """
            Converts spherical coordinates to Cartesian coordinates.

            Parameters
            ----------
            theta : float
                Inclination angle in radians.
            phi : float
                Azimuthal angle in radians.
            r : float
                Radius.
            """
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            return x, y, z

        phi_range = np.deg2rad(lon_range) + np.pi
        theta_range = np.deg2rad(lat_range) + np.pi / 2
        points = []

        n = int(4 * np.pi / distance**2)
        if n < 1:
            return

        a = 4 * np.pi / n
        d = np.sqrt(a)
        mtheta = int(np.round(np.pi / d))
        dtheta = np.pi / mtheta
        dphi = a / dtheta

        thetavals = np.pi * (np.arange(mtheta) + 0.5) / mtheta
        thetavals = thetavals[(thetavals >= theta_range[0]) & (thetavals < theta_range[1])]

        for theta in thetavals:
            mphi = max(1, int(np.round(2 * np.pi * np.sin(theta) / dphi)))
            phivals = 2 * np.pi * np.arange(mphi) / mphi
            phivals = phivals[(phivals >= phi_range[0]) & (phivals < phi_range[1])]
            for phi in phivals:
                points.append(_sph2cart(theta, phi, radius))
        if len(points) == 0:
            return

        points = np.array(points, dtype=np.float64)
        decimals = int(np.ceil(np.log10((2 * np.pi * radius) / distance)))
        points = np.unique(points.T.round(decimals=decimals), axis=0).T
        points = points.T

        return points

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """
        pix_str = format_large_units(self.pix, quantity="length")
        print(f"Generating a mesh with uniformly distributed faces of size ~{pix_str}.")
        points = self._distribute_points(distance=self.pix / self.radius)
        points[:, 0] = np.array([0, 0, 1])
        points[:, -1] = np.array([0, 0, -1])
        return points
