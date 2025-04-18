import os
import numpy as np
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.grid import register_grid_type, GridMaker

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
        
    Returns
    -------
    ArbitraryResolutionGrid
        An instance of the ArbitraryResolutionGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 pix: FloatLike, 
                 radius: FloatLike, 
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.pix = pix
        self.radius = radius
    

    def generate_face_distribution(self) -> NDArray:
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
                      grid_hash: str | None = None,
                      **kwargs: Any) -> tuple[os.PathLike, os.PathLike]:        
        super().generate_grid(grid_file=grid_file, grid_hash=grid_hash, **kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        print(f"Effective pixel size: {pix_mean:.2f} +/- {pix_std:.2f} m")
        return

