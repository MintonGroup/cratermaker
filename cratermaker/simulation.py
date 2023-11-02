from ._surface import _SurfaceBind
import xarray as xr
import numpy as np
import jigsawpy
import os
from numpy.random import default_rng

# Define some basic parameters
body_radius_values = {
    "Moon" : 1737.35e3
    } 

valid_body_name = ["Moon","Custom"]

class Simulation(object):
    """
    This is a class that defines the basic Cratermaker surface object. 
    """
    def __init__(self):
        self.cwd = os.getcwd()
        
        # Set default configuration options
        # TODO: Initialize with configure options as arguments or read in configuration file
        self.config = {
            "body" : "Moon", # Which body to simulation (Options are "Moon","Custom" )
            "pix"  : 6.16e3, # Approximate cell global cell size
            "cachedir" : os.path.join(self.cwd,".cache"), # Directory location of output files
            "meshname" : "body", # Name of files generated during mesh process
            "seed" : 5421109845, # RNG seed
        }
        
        # Determine the radius of the target body
        if self.config['body'] != 'Custom':
            self.config['body_radius'] = body_radius_values[self.config['body']]
           
        # Set up Jigsaw mesh files and mesh body if a file doesn't already exist
        self.meshfile = os.path.join(self.config['cachedir'],f"{self.config['meshname']}.msh")
        if not os.path.exists(self.config['cachedir']):
            os.mkdir(self.config['cachedir'])
        if not os.path.exists(self.meshfile):
            self.mesh_body()   
            
        self.rng = default_rng(seed=self.config['seed'])
        # self._surface = _SurfaceBind(gridshape)     
        # lat = np.linspace(-90, 90, gridshape[0])
        # lon = np.linspace(0, 360, gridshape[1])
        # self.ds = xr.Dataset(
        #     {'elevation': (['lat', 'lon'], self.get_elevation())},
        #                 coords={'lat': lat, 'lon': lon}
        # )
            
        
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

    def mesh_body(self):
        # This will get updated evantually after we're done testing
        meshname=self.config['meshname']
        cachedir=self.config['cachedir']
        # Set up jigsaw objects
        opts = jigsawpy.jigsaw_jig_t()
        geom = jigsawpy.jigsaw_msh_t()
        mesh = jigsawpy.jigsaw_msh_t()

        # Define file names
        opts.geom_file = os.path.join(cachedir,f"{meshname}_base.msh")
        opts.jcfg_file = os.path.join(cachedir,f"{meshname}.jig")
        opts.mesh_file = os.path.join(cachedir,f"{meshname}.msh")

        # Define basic sphere using the built-in ellipsoid mesh model
        geom.mshID = "ellipsoid-mesh"
        geom.radii = np.full(3, self.config['body_radius'], dtype=geom.REALS_t)
        jigsawpy.savemsh(opts.geom_file, geom)

        # Set mesh options
        opts.hfun_scal = "absolute"  #scaling type for mesh-size function, either "relative" or "absolute." For "relative" mesh-size values as percentages of the (mean) length of the axis-aligned bounding-box (AABB) associated with the geometry. "absolute" interprets mesh-size values as absolute measures.
        opts.hfun_hmax = self.config['pix']       # mesh-size function value (meaning depends on hfun_scal: for "relative" it is a percentage, and for "absolute" it is a dimensional value)
        opts.mesh_dims = +2          # Number of mesh dimensions (2 for surface mesh, 3 for volume)
        opts.optm_qlim = +9.5E-01    # threshold on mesh cost function above which gradient-based optimisation is attempted.
        opts.optm_iter = +32         # max. number of mesh optimisation iterations. 
        opts.optm_qtol = +1.0E-05    # tolerance on mesh cost function for convergence. Iteration on a given node is terminated if adjacent element cost-functions are improved by less than QTOL.
        opts.mesh_kern = "delfront"  # meshing kernel, choice of the standard Delaunay-refinement algorithm ('delaunay') or the Frontal-Delaunay method ('delfront').

        # Generate tesselated mesh
        jigsawpy.cmd.jigsaw(opts, mesh)

        scr2 = jigsawpy.triscr2(            # "quality" metric
            mesh.point["coord"],
            mesh.tria3["index"])

        jigsawpy.savevtk(os.path.join(cachedir,f"{meshname}.vtk"), mesh)
        self.mesh = mesh
        return