# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset 

cdef extern from "_bindings.h":
    ctypedef packed struct c_surface_type:
        double *elev_data
        int elev_shape[2]
 
    ctypedef double c_double
 
    c_surface_type* bind_surface_init(int gridsize)
    void bind_surface_final(c_surface_type *obj)

cdef class Surface:
    cdef c_surface_type* c_surf

    def __cinit__(self, int gridsize):
        self.c_surf = bind_surface_init(gridsize)
        if self.c_surf is NULL:
            raise MemoryError("Failed to allocate surface object")

        # Manually set the shape values based on gridsize
        self.c_surf.elev_shape[0] = gridsize
        self.c_surf.elev_shape[1] = gridsize

    def __dealloc__(self):
        if self.c_surf is not NULL:
            bind_surface_final(self.c_surf)

    def get_elev(self):
        # Retrieve the shape of the elevation data
        cdef np.npy_intp shape[2]
        shape[0] = self.c_surf.elev_shape[0]
        shape[1] = self.c_surf.elev_shape[1]

        # Create a NumPy array from the Fortran array
        cdef np.ndarray[np.float64_t, ndim=2, mode="c"] elev_array
        elev_array = np.PyArray_SimpleNewFromData(2, shape, np.NPY_FLOAT64, <void*>(self.c_surf.elev_data))

        return elev_array
 
    def set_elev(self, np.ndarray[np.float64_t, ndim=2] elev_array):
        cdef np.npy_intp old_shape[2]
        cdef np.npy_intp new_shape[2]
        old_shape[0] = self.c_surf.elev_shape[0]
        old_shape[1] = self.c_surf.elev_shape[1]
        new_shape[0] = elev_array.shape[0]
        new_shape[1] = elev_array.shape[1]

        if new_shape[0] != old_shape[0] or new_shape[1] != old_shape[1]:
            raise ValueError(f"Invalid shape for elev array: {new_shape} does not match {old_shape}")

        # Get the dimensions of the elev_array
        cdef int rows = elev_array.shape[0]
        cdef int cols = elev_array.shape[1]

        # Manually copy data from the NumPy array to the Fortran array
        cdef double* c_elev_data = self.c_surf.elev_data
        for row in range(rows):
            for col in range(cols):
                c_elev_data[row * cols + col] = elev_array[row, col]
