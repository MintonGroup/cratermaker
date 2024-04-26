# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "ejecta.h":
    void bind_ejecta_profile(double *r, int num_elements, double diameter, double ejrim, double *elevation)

def profile(cnp.ndarray[cnp.float64_t, ndim=1] r_array, 
            cnp.float64_t diameter,
            cnp.float64_t ejrim):
    """
    Generate an ejecta profile.

    Parameters
    ----------
    r_array : ndarray(N,)
        radial distance from the center of the crater.
    diameter : float
        diameter of the crater.
    ejrim : float
        ejecta rim height.

    Returns
    -------
    ndarray(N,)
        computed elevation values

    Raises
    ------
    ValueError - If required parameters are missing, arrays are mismatched, or if an invalid model is specified.
    """

    cdef int num_elements = r_array.size

    # Ensure memory-contiguous numpy array
    r_array = np.ascontiguousarray(r_array, dtype=np.float64)
    elevation_array = np.empty(num_elements, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] r = r_array
    cdef cnp.float64_t[::1] elevation = elevation_array

    bind_ejecta_profile(&r[0], num_elements, diameter, ejrim, &elevation[0])

    return elevation  
