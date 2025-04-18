import os
import numpy as np
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.grid import register_grid_type, GridMaker

@register_grid_type("icosphere")
class IcosphereGrid(GridMaker):    
    """
    Create a uniform grid configuration using an icosphere. This is the most accurate and efficient way to create a uniform grid, but is limited to a few resolutions.
    
    Parameters
    ----------
    gridlevel : float
        The subdivision level of the icosphere. The number of faces is 20 * 4**level. The default level is 8.
    radius: FloatLike
        The radius of the target body in meters.
        
    Returns
    -------
    IcosphereGrid
        An instance of the IcosphereGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 gridlevel: int = 8, 
                 radius: FloatLike = 1.0, 
                 **kwargs: Any):
        super().__init__(**kwargs)
        self.gridlevel = gridlevel
        self.radius = radius
        
        
    def generate_face_distribution(self, **kwargs: Any) -> NDArray:
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