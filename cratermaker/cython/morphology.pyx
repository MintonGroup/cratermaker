# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "morphology.h":
    void bind_morphology_profile(double *r, double *reference_elevation, int num_elements, double diameter, double floordepth, double floordiam, double rimheight, double ejrim, double RIMDROP, double *elevation)


def profile(cnp.ndarray[cnp.float64_t, ndim=1] r_array, 
            cnp.ndarray[cnp.float64_t, ndim=1] reference_elevation_array, 
            cnp.float64_t diameter,
            cnp.float64_t floordepth,
            cnp.float64_t floordiam,
            cnp.float64_t rimheight,
            cnp.float64_t ejrim,
            cnp.float64_t RIMDROP):
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
    floordiam : float
        diameter of the crater floor.
    rimheight : float
        height of the crater rim.
    ejrim : float
        ejecta riom height.
    RIMDROP : float
        rim drop exponent

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

    bind_morphology_profile(&r[0], &reference_elevation[0], num_elements, diameter, floordepth, floordiam, rimheight, ejrim, RIMDROP, &elevation[0])

    return elevation  
