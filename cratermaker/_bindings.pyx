# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset 

cdef extern from "bindings.h":
    ctypedef struct surface_type:
        double *elev_data
        int elev_shape[2]

    ctypedef double c_double

    surface_type* bind_surface_init(int gridsize)
    void bind_surface_final(surface_type* obj)

cdef class Surface:
    cdef surface_type* c_surf

    def __cinit__(self, int gridsize):
        self.c_surf = bind_surface_init(gridsize)
        if self.c_surf is NULL:
            raise MemoryError("Failed to allocate surface object")

    def __dealloc__(self):
        if self.c_surf is not NULL:
            bind_surface_final(self.c_surf)
            free(self.c_surf)

    # def get_elev(self):
    #    # Create a NumPy array from the Fortran array
    #    cdef np.ndarray[np.float64_t, ndim=2] elev_array
    #    elev_array = np.PyArray_SimpleNewFromData(2, <np.npy_intp>[gridsize, gridsize], np.NPY_FLOAT64, <void*>self.c_surf.elev)
    #    return elev_array
#
#    def set_elev(self, np.ndarray[np.float64_t, ndim=2] elev_array):
#        if elev_array.shape != (gridsize, gridsize):
#            raise ValueError("Invalid shape for elev array")
#        # Copy data from elev_array to the Fortran array
#        np.copyto(elev_array, self.get_elev())
