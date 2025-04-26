from pathlib import Path
import numpy as np
from typing import Any
from numpy.typing import NDArray
from cratermaker.utils.custom_types import FloatLike
from cratermaker.components.grid import register_grid_type, Grid
from cratermaker.utils.general_utils import parameter

@register_grid_type("icosphere")
class IcosphereGrid(Grid):    
    """
    Create a uniform grid configuration using an icosphere. This is the most accurate and efficient way to create a uniform grid, but is limited to a few resolutions.
    
    Parameters
    ----------
    gridlevel : float
        The subdivision level of the icosphere. The number of faces is 20 * 4**level. The default level is 8.
    radius: FloatLike
        The radius of the target body in meters.
    simdir : str | Path
        The main project simulation directory.


    Returns
    -------
    IcosphereGrid
        An instance of the IcosphereGrid class initialized with the given pixel size. 
    """    
    
    def __init__(self, 
                 gridlevel: int = 8, 
                 radius: FloatLike = 1.0, 
                 simdir: str | Path = Path.cwd(),
                 **kwargs: Any):
        super().__init__(radius=radius, simdir=simdir, **kwargs)
        self.gridlevel = gridlevel
        
    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self._radius, self._gridlevel] 


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
   
    
    def generate_grid(self, **kwargs: Any):
        super().generate_grid(**kwargs)
        face_areas = self.grid.face_areas 
        face_sizes = np.sqrt(face_areas / (4 * np.pi))
        pix_mean = face_sizes.mean().item() * self.radius
        pix_std = face_sizes.std().item() * self.radius
        print(f"Effective pixel size: {pix_mean:.2f} +/- {pix_std:.2f} m")
        return    
    
    @parameter
    def gridlevel(self) -> int:
        return self._gridlevel

    @gridlevel.setter
    def gridlevel(self, value: int) -> None:
        if value < 0:
            raise ValueError("Grid level must be a non-negative integer.")
        self._gridlevel = int(value)

    @property
    def _hashvars(self):
        """
        The variables used to generate the hash.
        """
        return [self._component_name, self._radius, self._gridlevel]
