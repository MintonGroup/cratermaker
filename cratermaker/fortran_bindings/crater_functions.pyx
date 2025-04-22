# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "crater_functions.h":
    void bind_crater_profile(double *r, double *reference_elevation, int num_elements, double diameter, double floordepth, double floor_diameter, double rimheight, double ejrim, double *elevation)


def profile(cnp.ndarray[cnp.float64_t, ndim=1] r_array, 
            cnp.ndarray[cnp.float64_t, ndim=1] reference_elevation_array, 
            cnp.float64_t diameter,
            cnp.float64_t floordepth,
            cnp.float64_t floor_diameter,
            cnp.float64_t rimheight,
            cnp.float64_t ejrim):
    """
    Generate a crater profile.

    Parameters
    ----------
    r_array : ndarray(N,)
        radial distance from the center of the crater.
    reference_elevation_array : ndarray(N,)
        reference elevation at each radial distance.
    diameter : float
        diameter of the crater.
    floordepth : float
        depth of the crater floor.
    floor_diameter : float
        diameter of the crater floor.
    rimheight : float
        height of the crater rim.
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

    if not (r_array.size == reference_elevation_array.size):
        raise ValueError("input arrays must have the same length")

    cdef int num_elements = r_array.size

    # Ensure memory-contiguous numpy array
    r_array = np.ascontiguousarray(r_array, dtype=np.float64)
    reference_elevation_array = np.ascontiguousarray(reference_elevation_array, dtype=np.float64)
    elevation_array = np.empty(num_elements, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] r = r_array
    cdef cnp.float64_t[::1] reference_elevation = reference_elevation_array
    cdef cnp.float64_t[::1] elevation = elevation_array

    bind_crater_profile(&r[0], &reference_elevation[0], num_elements, diameter, floordepth, floor_diameter, rimheight, ejrim, &elevation[0])

    return elevation  
