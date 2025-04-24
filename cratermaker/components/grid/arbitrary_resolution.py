import os
from pathlib import Path
import numpy as np
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.grid import register_grid_type, GridMaker
from cratermaker.utils.general_utils import parameter

@register_grid_type("arbitrary_resolution")
class ArbitraryResolutionGrid(GridMaker):
    """
    Create a uniform grid configuration with an arbitrary user-defined pixel size. This will not be as nice as the regular IcosphereGrid, but can be any resolution desired.
    
    Parameters
    ----------
    pix : float
        The approximate face size for the mesh in meters.
    radius: FloatLike
        The radius of the target body in meters.
    simdir: os.PathLike
        The directory where the simulation files are stored. Default is the current working directory.  

    Returns
    -------
    ArbitraryResolutionGrid
        An instance of the ArbitraryResolutionGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 pix: FloatLike | None = None, 
                 radius: FloatLike = 1.0, 
                 simdir: os.PathLike = Path.cwd(),             
                 **kwargs: Any):
        super().__init__(radius=radius, simdir=simdir, **kwargs)
        self.pix = pix
        


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
   
    
    def generate_grid(self,
                      grid_file: os.PathLike,
                      id: str | None = None,
                      **kwargs: Any) -> tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, id=id, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        print(f"Effective pixel size: {pix_mean:.2f} +/- {pix_std:.2f} m")
        return

