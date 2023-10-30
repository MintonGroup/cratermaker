# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
from cpython cimport PyUnicode_AsUTF8AndSize, PyUnicode_FromString
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memset 

cdef extern from "_surface.h":
    ctypedef struct c_surface_type:
        double *elevation
        int elevation_shape[2]
        char *stringvar
        int stringvar_len

    c_surface_type* bind_surface_init(int ny, int nx)
    void bind_surface_final(c_surface_type *obj)
    void bind_surface_set_stringvar(c_surface_type *obj, const char *c_string)
    char* bind_surface_get_stringvar(c_surface_type *obj)


cdef class Simulation_binding:
    cdef c_surface_type* fobj

    def __cinit__(self, tuple shape):
        """
        Initializes the Simulation object by calling the Fortran initializer, which will allocate the array in Fortran, set some initial values, and return a pointer that connects the Fortran derived type class variable to the Python object.

        Parameters
        ----------
            Shape of the allocatable Numpy array to initialize in Fortran 

        Returns
        -------
            Sets the fobj component variable containing the components set by the Fortran object.
        """

        # Check to make sure we are passing a correct 2D array for the shape.
        if len(shape) != 2:
            raise ValueError("Expected a tuple of length 2 for shape")

        print("Cython: calling bind_surface_init")
        self.fobj = bind_surface_init(shape[0],shape[1])  
        print("Cython: Successfully returned")

        # Do some basic checks to make sure the object variable and all its components were allocated succesfully
        if self.fobj is NULL:
            raise MemoryError("Failed to allocate Fortran object.")
        else:
            print("The Fortran object was allocated successfully ")
        print(f"self.fobj           = {<unsigned long>self.fobj          }")
        print(f"self.fobj.elevation = {<unsigned long>self.fobj.elevation}")
        print(f"self.fobj.stringvar = {<unsigned long>self.fobj.stringvar}")

        if self.fobj.elevation is NULL:
            raise MemoryError("Failed to allocate component variable 'elevation' in the Fortran object.")
        else:
            print("The component variable 'elevation' was allocated successfuly in the Fortran object")

        if self.fobj.stringvar is NULL: # <- This is where the problem lies
            raise MemoryError("Failed to allocate component variable 'stringvar' in the Fortran object.")
        else:
            print("The component variable 'stringvar' was allocated successfuly in the Fortran object")

        # Manually set the shape of the 2D component array in the Python object
        self.fobj.elevation_shape[0] = shape[0]
        self.fobj.elevation_shape[1] = shape[1]

        return


    def __dealloc__(self):
        """
        Finalizes the Fortran component variable.

        Parameters
        ----------
            None
        Returns
        -------
            Deallocates the fobj component variables from Fortran.
        """
        if self.fobj is not NULL:
            bind_surface_final(self.fobj)


    def get_elevation(self):
        """
        A getter method that retrieves the elevation allocatable array from Fortran and returns it as a Numpy array

        Parameters
        ----------
            None
        Returns
        -------
            elevation_array : [y,x] Numpy array
        """

        # Retrieve the shape of the elevation data
        cdef cnp.npy_intp shape[2]
        shape[0] = self.fobj.elevation_shape[0]
        shape[1] = self.fobj.elevation_shape[1]

        # Create a NumPy array from the Fortran array
        cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] elevation_array
        elevation_array = cnp.PyArray_SimpleNewFromData(2, shape, cnp.NPY_FLOAT64, <void*>(self.fobj.elevation))

        return elevation_array
 
    def set_elevation(self, cnp.ndarray[cnp.float64_t, ndim=2] elevation_array):
        """
        A setter method that sets the value of the elevation allocatable array in Fortran from a Numpy array

        Parameters
        ----------
            elevation_array : [y,x] Numpy array
        Returns
        -------
            None : Sets the values of self.fobj
        """
        cdef cnp.npy_intp old_shape[2]
        cdef cnp.npy_intp new_shape[2]
        old_shape[0] = self.fobj.elevation_shape[0]
        old_shape[1] = self.fobj.elevation_shape[1]
        new_shape[0] = elevation_array.shape[0]
        new_shape[1] = elevation_array.shape[1]

        if new_shape[0] != old_shape[0] or new_shape[1] != old_shape[1]:
            raise ValueError(f"Invalid shape for elevation array: {new_shape} does not match {old_shape}")

        # Get the dimensions of the elevation_array
        cdef int rows = elevation_array.shape[0]
        cdef int cols = elevation_array.shape[1]

        # Manually copy data from the NumPy array to the Fortran array
        cdef double* c_elevation = self.fobj.elevation
        for row in range(rows):
            for col in range(cols):
                c_elevation[row * cols + col] = elevation_array[row, col]

    def get_stringvar(self):
        """
        A getter method that retrieves the stringvar from Fortran and returns it as a Python string 

        Parameters
        ----------
            None
        Returns
        -------
            string : str
        """
        cdef char *c_string
        c_string = bind_surface_get_stringvar(self.fobj)
        
        if c_string == NULL:
            return None
        else:
            py_string = PyUnicode_FromString(c_string)
            # Don't forget to free the C string if allocated in Fortran
            #free(c_string)
            return py_string

    
    def set_stringvar(self, str string):
        """
        A setter method that sets the value of the stringvar in Fortran from a Python string. 

        Parameters
        ----------
            string : str 
                Input string
        Returns
        -------
            None : Sets the values of self.fobj
        """
        cdef const char *c_string
        cdef Py_ssize_t length

        c_string = PyUnicode_AsUTF8AndSize(string, &length)
        bind_surface_set_stringvar(self.fobj, c_string)
