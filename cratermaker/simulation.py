from ._surface import _SurfaceBind
import xarray as xr
import numpy as np
import xesmf as xe

class Simulation(object):
    """
    This is a class that defines the basic Cratermaker surface object. 
    """
    def __init__(self, param_file="cratermaker.json", gridshape=(2000,1000)):
        self._surface = _SurfaceBind(gridshape)     
        lat = np.linspace(-90, 90, gridshape[0])
        lon = np.linspace(0, 360, gridshape[1])
        self.ds = xr.Dataset(
            {'elevation': (['lat', 'lon'], self.get_elevation())},
                        coords={'lat': lat, 'lon': lon}
        )
            
        
    @property
    def elevation(self):
        return self._surface.fobj.elevation
    
    @property
    def stringvar(self):
        return self._surface.fobj.stringvar

    def get_elevation(self):
        return self._surface.get_elevation()
    
    def set_elevation(self, elevation_array):
        self._surface.set_elevation(elevation_array)

