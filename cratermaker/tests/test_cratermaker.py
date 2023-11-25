import numpy as np
import cratermaker

sim = cratermaker.Simulation()

noise_height = 10e4 # 10 km noise height
num_octaves=12
pers = 0.5
freq = 1.92
scale = 1.0e4
gain0 = 70.0
gain = 0.5
warp0 = 0.4
warp  = 0.35
damp0 = 1.0
damp  = 0.8
jordscale = 0.01
damp_scale  = noise_height * jordscale

anchor = sim.rng.uniform(0.0,scale,size=(num_octaves,3))

def make_noise(x,y,z):
    noise = 0.0
    noise += cratermaker.util_perlin("jordan",x,y,z,num_octaves,anchor,lacunarity=freq,gain=gain,gain0=gain0,warp0=warp0,warp=warp,damp0=damp0,damp=damp,damp_scale=damp_scale)
    return noise

# Get the vertices of the mesh
print("Stacking node arrays")
vertices = np.column_stack((sim.data.uxgrid.node_x,sim.data.uxgrid.node_y,sim.data.uxgrid.node_z))

# Normalize the position vectors
norms = np.linalg.norm(vertices, axis=1, keepdims=True)
normalized_vectors = vertices / norms

# Vectorized calculation of noise for each vertex
print("Computing noise")
noises = np.vectorize(make_noise)(normalized_vectors[:, 0], normalized_vectors[:, 1], normalized_vectors[:, 2])
print("normalized: ",np.max(normalized_vectors),np.min(normalized_vectors))
print("noises ",np.max(noises),np.min(noises))

# Update the mesh vertices
sim.data['elevation'] += noises
sim.data.to_netcdf("noise_sphere.nc")