# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "ejecta.h":
    void bind_ejecta_distribution(double *radial_distance, double *initial_bearing, int num_elements, double crater_diameter, double ejrim, double ejecta_truncation, int dorays, double *ejecta_thickness);
    void bind_ejecta_profile(double *radial_distance, int num_elements, double crater_diameter, double ejrim, double *elevation)
    void bind_ejecta_ray_pattern(double *radial_distance, double *initial_bearing, int num_elements, double crater_diameter, double ejrim, double ejecta_truncation, double *ejecta_thickness);


def distribution(cnp.ndarray[cnp.float64_t, ndim=1] radial_distance, 
                cnp.ndarray[cnp.float64_t, ndim=1] initial_bearing,
                cnp.float64_t crater_diameter,
                cnp.float64_t ejrim,
                bint dorays,
                cnp.float64_t ejecta_truncation):
    """
    Generate an ejecta profile.

    Parameters
    ----------
    radial_distance : ndarray(N,)
        Radial distance from the center of the crater.
    initial_bearing : ndarray(N,)
        Initial bearing from the denter of the crater in radians.
    crater_diameter : float
        Diameter of the crater.
    ejrim : float
        Ejecta rim height.
    ejecta_truncation : float
        Ejecta truncation distance relative to the crater radius.
    dorays : bool
        Set to True to compute a ray pattern. Set to False to compute a homogenous distribution. 

    Returns
    -------
    ndarray(N,)
        computed elevation values

    Raises
    ------
    ValueError - If required parameters are missing, arrays are mismatched, or if an invalid model is specified.
    """

    cdef int num_elements = radial_distance.size

    # Ensure memory-contiguous numpy array
    radial_distance = np.ascontiguousarray(radial_distance, dtype=np.float64)
    initial_bearing = np.ascontiguousarray(initial_bearing, dtype=np.float64)
    thickness_array = np.empty(num_elements, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] r = radial_distance
    cdef cnp.float64_t[::1] theta = initial_bearing
    cdef cnp.float64_t[::1] ejecta_thickness = thickness_array
    cdef int c_dorays = 1 if dorays else 0

    bind_ejecta_distribution(&r[0], &theta[0], num_elements, crater_diameter, ejrim, ejecta_truncation, c_dorays, &ejecta_thickness[0])

    return ejecta_thickness  


def profile(cnp.ndarray[cnp.float64_t, ndim=1] radial_distance, 
            cnp.float64_t crater_diameter,
            cnp.float64_t ejrim):
    """
    Generate an ejecta profile.

    Parameters
    ----------
    radial_distance : ndarray(N,)
        radial distance from the center of the crater.
    crater_diameter : float
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

    cdef int num_elements = radial_distance.size

    # Ensure memory-contiguous numpy array
    radial_distance = np.ascontiguousarray(radial_distance, dtype=np.float64)
    elevation_array = np.empty(num_elements, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] r = radial_distance
    cdef cnp.float64_t[::1] elevation = elevation_array

    bind_ejecta_profile(&r[0], num_elements, crater_diameter, ejrim, &elevation[0])

    return elevation  

def ray_pattern(cnp.ndarray[cnp.float64_t, ndim=1] radial_distance, 
                cnp.ndarray[cnp.float64_t, ndim=1] initial_bearing,
                cnp.float64_t crater_diameter,
                cnp.float64_t ejrim,
                cnp.float64_t ejecta_truncation):
    """
    Generate an ejecta profile.

    Parameters
    ----------
    radial_distance : ndarray(N,)
        radial distance from the center of the crater.
    initial_bearing : ndarray(N,)
        initial bearing from the denter of the crater in radians.
    crater_diameter : float
        diameter of the crater.
    ejrim : float
        ejecta rim height.
    ejecta_truncation : float
        ejecta truncation distance relative to the crater radius.

    Returns
    -------
    ndarray(N,)
        computed elevation values

    Raises
    ------
    ValueError - If required parameters are missing, arrays are mismatched, or if an invalid model is specified.
    """

    cdef int num_elements = radial_distance.size

    # Ensure memory-contiguous numpy array
    radial_distance = np.ascontiguousarray(radial_distance, dtype=np.float64)
    initial_bearing = np.ascontiguousarray(initial_bearing, dtype=np.float64)
    thickness_array = np.empty(num_elements, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] r = radial_distance
    cdef cnp.float64_t[::1] theta = initial_bearing
    cdef cnp.float64_t[::1] ejecta_thickness = thickness_array

    bind_ejecta_ray_pattern(&r[0], &theta[0], num_elements, crater_diameter, ejrim, ejecta_truncation, &ejecta_thickness[0])

    return ejecta_thickness  

