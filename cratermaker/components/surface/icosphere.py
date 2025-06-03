import math
from pathlib import Path
from typing import Any

from numpy.typing import NDArray

from cratermaker.components.surface import Surface
from cratermaker.components.target import Target
from cratermaker.utils.general_utils import format_large_units, parameter


@Surface.register("icosphere")
class IcosphereSurface(Surface):
    """
    Create a uniform grid configuration using an icosphere. This is the most accurate and efficient way to create a uniform grid, but is limited to a few resolutions.

    Parameters
    ----------
    gridlevel : float, default 8
        The subdivision level of the icosphere.
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
    IcosphereSurface
        An instance of the IcosphereSurface class initialized with the given pixel size.

    Notes
    -----
    The number of faces of the icosphere is given by the formula:

    .. math::
        n_face = 10 * 4^{gridlevel} + 2

    The number of nodes is given by the formula:

    .. math::
        n_node = 20 * 4^{gridlevel}
    """

    def __init__(
        self,
        gridlevel: int = 8,
        target: Target | str | None = None,
        reset: bool = False,
        regrid: bool = False,
        simdir: str | Path | None = None,
        **kwargs: Any,
    ):
        super().__init__(target=target, simdir=simdir, **kwargs)
        self.gridlevel = gridlevel
        self._load_from_files(reset=reset, regrid=regrid, **kwargs)

    def __str__(self) -> str:
        base = super().__str__()
        pix_mean = format_large_units(self.pix_mean, quantity="length")
        pix_std = format_large_units(self.pix_std, quantity="length")
        return f"{base}\nGrid Level: {self.gridlevel}\nEffective pixel size: {pix_mean} +/- {pix_std}"

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self.radius, self.gridlevel]

    def _generate_face_distribution(self, **kwargs: Any) -> NDArray:
        """
        Creates the points that define the mesh centers.

        Returns
        -------
        (3,n) ndarray of np.float64
            Array of points on a unit sphere.
        """
        from trimesh.creation import icosphere

        print(f"Generating a mesh with icosphere level {self.gridlevel}.")
        mesh = icosphere(self.gridlevel)
        points = mesh.vertices.T
        return points

    @parameter
    def gridlevel(self) -> int:
        return self._gridlevel

    @gridlevel.setter
    def gridlevel(self, value: int) -> None:
        if value < 0:
            raise ValueError("Grid level must be a non-negative integer.")
        self._gridlevel = int(value)

    @property
    def pix(self) -> float:
        """
        The approximate face size for a cell of the mesh.
        """
        if self._pix is None:
            nfaces = 10 * 4**self.gridlevel + 2
            self._pix = (4 * math.pi * self.radius**2 / nfaces) ** 0.5
        return self._pix
