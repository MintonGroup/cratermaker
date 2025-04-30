from pathlib import Path
import numpy as np
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.grid import Grid
from cratermaker.utils.general_utils import parameter

@Grid.register("arbitrary_resolution")
class ArbitraryResolutionGrid(Grid):
    """
    Create a uniform grid configuration with an arbitrary user-defined pixel size. This will not be as nice as the regular IcosphereGrid, but can be any resolution desired.
    
    Parameters
    ----------
    pix : float
        The approximate face size for the mesh in meters.
    radius: FloatLike
        The radius of the target body in meters.
    simdir : str | Path
        The main project simulation directory.

    Returns
    -------
    ArbitraryResolutionGrid
        An instance of the ArbitraryResolutionGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 pix: FloatLike | None = None, 
                 radius: FloatLike = 1.0, 
                 simdir: str | Path | None = None,
                 **kwargs: Any):
        super().__init__(radius=radius, simdir=simdir, **kwargs)
        self.pix = pix

    def __repr__(self) -> str:
        base = super().__repr__()
        return (
            f"{base}"
            f"Effective pixel size: {self.pix_mean:.2f} +/- {self.pix_std:.2f} m\n"
        )           

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self._radius, self._pix]

    @parameter
    def pix(self):
        """
        The approximate face size for a cell of the mesh.
        """
        return self._pix

    
    @pix.setter
    def pix(self, value: FloatLike):
        if value is None:
            value= np.sqrt(4 * np.pi * self.radius**2) * 1e-3  # Default mesh scale that is somewhat comparable to a 1000x1000 CTEM grid
        elif not isinstance(value, FloatLike) or np.isnan(value) or np.isinf(value) or value <= 0:
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

        print(f"Generating a mesh with uniformly distributed faces of size ~{self.pix} m.")
        points = self._distribute_points(distance=self.pix/self.radius) 
        points[:,0] = np.array([0,0,1])
        points[:,-1] = np.array([0,0,-1])
        return points
    

