from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cratermaker.components.surface import Surface
from cratermaker.components.target import Target
from cratermaker.constants import FloatLike
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
        self.load_from_files(reset=reset, regrid=regrid, **kwargs)

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
        elif (
            not isinstance(value, FloatLike)
            or np.isnan(value)
            or np.isinf(value)
            or value <= 0
        ):
            raise TypeError("pix must be a positive float")
        self._pix = float(value)

    def generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.

        """

        print(
            f"Generating a mesh with uniformly distributed faces of size ~{self.pix} m."
        )
        points = self._distribute_points(distance=self.pix / self.radius)
        points[:, 0] = np.array([0, 0, 1])
        points[:, -1] = np.array([0, 0, -1])
        return points
