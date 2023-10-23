# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport numpy as np
import numpy as np

cdef extern from "_bindings.h":
    ctypedef struct surface_tally:
        double *diam, *xl, *yl
        int data_size

    void initialize_type(my_type* obj, int size)
    void process_type(my_type* obj)
    void finalize_type(my_type* obj)

cdef class MyType:
    cdef my_type* c_obj

    def __cinit__(self, int size):
        self.c_obj = <my_type*> malloc(sizeof(my_type))
        initialize_type(self.c_obj, size)
        self.c_obj.data_size = size

    def __dealloc__(self):
        if self.c_obj is not NULL:
            finalize_type(self.c_obj)
            free(self.c_obj)

    def process(self):
        process_type(self.c_obj)

    property data:
        def __get__(self):
            return <np.ndarray[np.double_t, ndim=1]>np.PyArray_SimpleNewFromData(1, &self.c_obj.data_size, np.NPY_DOUBLE, <void*>self.c_obj.data)
