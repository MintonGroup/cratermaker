from ._surface import _SurfaceBind

class Simulation(object):
    """
    This is a class that defines the basic Cratermaker surface object. 
    """
    def __init__(self, param_file="cratermaker.json", gridshape=(1000,1000)):
        self._surface = _SurfaceBind(gridshape)     
        
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

