# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "perlin.h":
    ctypedef struct PerlinArguments:
        double damp
        double damp0
        double damp_scale
        double freq
        double gain
        double gain0
        double lacunarity
        double noise_height
        double pers
        double slope
        double warp
        double warp0

    double bind_perlin_noise(const char *model, double x, double y, double z, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0)


cdef void _to_fortran_2D_double_array(cnp.ndarray[cnp.float64_t, ndim=2] src, double* dest, int rows, int cols):
    """
    Copies data from a 2D NumPy array to a 2D Fortran-style array.

    Parameters
    ----------
    src : cnp.ndarray
        The source NumPy array with 2 dimensions.
    dest : double*
        The pointer to the destination array in memory which should be already allocated.
    rows : int
        The number of rows in the destination array.
    cols : int
        The number of columns in the destination array.
    """
    # Assuming dest is already allocated with the correct size
    # NumPy arrays are row-major, while Fortran arrays are column-major,
    # so we need to transpose the indices when copying.
    for i in range(rows):
        for j in range(cols):
            dest[j * rows + i] = src[i, j]  # Transpose the index for Fortran's column-major order


def apply_noise(str model, cnp.float64_t x, cnp.float64_t y, cnp.float64_t z, int num_octaves, cnp.ndarray[cnp.float64_t, ndim=2] anchor, **kwargs):

    # Ensure memory-contiguous numpy array
    anchor = np.ascontiguousarray(anchor, dtype=np.float64)

    # Convert numpy 2D array to double**
    cdef double* f_anchor = <double*> malloc(3 * num_octaves * sizeof(double))
    _to_fortran_2D_double_array(anchor, f_anchor, num_octaves, 3)

    model = model.lower()
    kw = {
        'damp' : 0.0, 
        'damp0' : 0.0, 
        'damp_scale' : 0.0, 
        'freq' : 0.0, 
        'gain' : 0.0, 
        'gain0' : 0.0, 
        'lacunarity' : 0.0, 
        'noise_height' : 0.0, 
        'pers' : 0.0, 
        'slope' : 0.0, 
        'warp' : 0.0, 
        'warp0' : 0.0
    }
    required_kwargs = {
        'turbulence' : ['noise_height', 'freq', 'pers'],
        'billowed'   : ['noise_height', 'freq', 'pers'],
        'plaw'       : ['noise_height', 'freq', 'pers', 'slope'],
        'ridged'     : ['noise_height', 'freq', 'pers'],
        'swiss'      : ['lacunarity', 'gain', 'warp'],
        'jordan'     : ['lacunarity', 'gain0', 'gain', 'warp0', 'warp', 'damp0', 'damp', 'damp_scale'],
    }

    if model in required_kwargs:
        # If it is a valid model, make sure we've passed all of the required arguments
        if set(required_kwargs[model]).issubset(kwargs.keys()):
            # iterate over required_kwargs[model] for all keys and set the value of kw[key]to be given by kwargs[key]
            for arg in required_kwargs[model]:
                kw[arg] = kwargs[arg]
            noise = bind_perlin_noise(model, x, y, z, num_octaves, f_anchor, kw['damp'], kw['damp0'], kw['damp_scale'], kw['freq'], kw['gain'], kw['gain0'], kw['lacunarity'], kw['noise_height'], kw['pers'], kw['slope'], kw['warp'], kw['warp0'] )
        else:
            missing_args = set(required_kwargs[model]) - kwargs.keys()
            raise ValueError(f"The {model} model requires the following missing keywords: {', '.join(missing_args)}")
    else:
        raise ValueError(f"{model} is an invalid model. Valid options are {', '.join(required_kwargs.keys())}")

    # Free the allocated memory
    free(f_anchor)

    return noise  # Return the noise value
