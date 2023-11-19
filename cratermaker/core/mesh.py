import trimesh
import jigsawpy
import numpy as np
import os

class Mesh(trimesh.Trimesh):
   def __init__(self, mesh_file, target, pix, *args, **kwargs):
      self.mesh_file = mesh_file
      if os.path.exists(self.mesh_file):
         vertices, faces = self.make_target_mesh(target, pix)
         super().__init__(vertices=vertices, faces=faces, *args, **kwargs)  # Initialize the base Trimesh class
      else:
         super().__init__(*args, **kwargs)  # Initialize the base Trimesh class
   
   def make_target_mesh(self,target, pix):
      """
      Generate a tessellated mesh of a sphere using the jigsawpy library and convert it to GLB format.

      This function sets up Jigsaw mesh files and a mesh body, defines a basic sphere using an ellipsoid mesh model,
      sets mesh options, generates a tessellated mesh, and converts the mesh to a trimesh object. The generated mesh 
      is then saved in GLB format.

      Returns
      -------
      None
         The function does not return a value but updates the `self.mesh` attribute with the generated trimesh object.
      """
      
      # This will get updated evantually after we're done testing
      # Set up jigsaw objects
      opts = jigsawpy.jigsaw_jig_t()
      geom = jigsawpy.jigsaw_msh_t()
      mesh = jigsawpy.jigsaw_msh_t()
      
      # Set up Jigsaw mesh files and mesh body if a file doesn't already exist
      opts.mesh_file = self.mesh_file.replace(".glb",".msh")
      opts.jcfg_file = self.mesh_file.replace(".glb",".jig")
      opts.geom_file = self.mesh_file.replace(".glb","_geom.msh")
      
      # Define basic sphere using the built-in ellipsoid mesh model
      geom.mshID = "ellipsoid-mesh"
      geom.radii = np.full(3, target.radius, dtype=geom.REALS_t)
      jigsawpy.savemsh(opts.geom_file,geom)

      # Set mesh options
      opts.hfun_scal = "absolute"  #scaling type for mesh-size function, either "relative" or "absolute." For "relative" mesh-size values as percentages of the (mean) length of the axis-aligned bounding-box (AABB) associated with the geometry. "absolute" interprets mesh-size values as absolute measures.
      opts.hfun_hmax = np.sqrt(2.0) * pix       # mesh-size function value. The "pix" value is adjusted to keep it scaled  to roughly the same as the older CTEM pix
      opts.mesh_dims = +2          # Number of mesh dimensions (2 for body mesh, 3 for volume)
      opts.optm_qlim = +9.5E-01    # threshold on mesh cost function above which gradient-based optimisation is attempted.
      opts.optm_iter = +64         # max. number of mesh optimisation iterations. 
      opts.optm_qtol = +1.0E-05    # tolerance on mesh cost function for convergence. Iteration on a given node is terminated if adjacent element cost-functions are improved by less than QTOL.
      opts.mesh_kern = "delfront"  # meshing kernel, choice of the standard Delaunay-refinement algorithm ('delaunay') or the Frontal-Delaunay method ('delfront').
      # Generate tesselated mesh
      jigsawpy.cmd.jigsaw(opts, mesh)
      
      # Convert the mesh to trimesh
      # Ensure the vertex and face data is in numpy array format
      vertices = np.array([vert[0] for vert in mesh.vert3])
      faces = np.array([face[0] for face in mesh.tria3])

      # Create a trimesh object
      
      # Save generated mesh 
      export = mesh.export(self.mesh_file)
      
      return vertices, faces
   
def load_target_mesh(mesh_file):
   """
   Load a target mesh from a file into the `self.mesh` attribute.

   This function uses the trimesh library to read a mesh file, which is expected to be in GLB format. The function
   handles the peculiarity of trimesh reading the mesh as a Scene instead of a Trimesh object and extracts the 
   actual mesh from the scene.

   Returns
   -------
   None
      The function does not return a value but updates the `self.mesh` attribute with the loaded trimesh object.
   """

   # This is not well documented, but trimesh reads the mesh in as a Scene instead of a Trimesh, so we need to extract the actual mesh
   scene = trimesh.load_mesh(mesh_file)
   
   # Assuming that the only geometry in the file is the body mesh, this should work
   mesh = next(iter(scene.geometry.values()))
   
   return mesh