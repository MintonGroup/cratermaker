# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
from cpython cimport PyUnicode_AsUTF8AndSize, PyUnicode_FromString
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset 
from libc.stdint cimport uintptr_t

cdef extern from "_bind.h":
    ctypedef struct surface_type:
        double *elevation
        int elevation_shape[2]
        char *stringvar
        int stringvar_len

    surface_type* bind_surface_init(int ny, int nx)
    void bind_surface_final(surface_type *obj)
    void bind_surface_set_stringvar(surface_type *obj, const char *c_string)
    char* bind_surface_get_stringvar(surface_type *obj)

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


cdef void to_fortran_2D_double_array(cnp.ndarray[cnp.float64_t, ndim=2] src, double* dest, int rows, int cols):
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



cdef class _SurfaceBind:
    cdef surface_type* fobj

    def __cinit__(self, tuple gridshape=(1000,1000)):
        """
        Initializes the Simulation object by calling the Fortran initializer, which will allocate the array in Fortran, set some initial values, and return a pointer that connects the Fortran derived type class variable to the Python object.

        Parameters
        ----------
            gridshape:
                Shape of the allocatable Numpy array to initialize in Fortran 

        Returns
        -------
            Sets the fobj component variable containing the components set by the Fortran object.
        """

        # Check to make sure we are passing a correct 2D array for the shape.
        if len(gridshape) != 2:
            raise ValueError("Expected a tuple of length 2 for gridshape")

        self.fobj = bind_surface_init(gridshape[0],gridshape[1])  

        # Do some basic checks to make sure the object variable and all its components were allocated succesfully
        if self.fobj is NULL:
            raise MemoryError("Failed to allocate Fortran object.")

        if self.fobj.elevation is NULL:
            raise MemoryError("Failed to allocate component variable 'elevation' in the Fortran object.")

        if self.fobj.stringvar is NULL: # <- This is where the problem lies
            raise MemoryError("Failed to allocate component variable 'stringvar' in the Fortran object.")

        # Manually set the shape of the 2D component array in the Python object
        self.fobj.elevation_shape[0] = gridshape[0]
        self.fobj.elevation_shape[1] = gridshape[1]

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

        to_fortran_2D_double_array(elevation_array, self.fobj.elevation, new_shape[0], new_shape[1])
        return


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


    
def util_perlin(str model, cnp.float64_t x, cnp.float64_t y, cnp.float64_t z, int num_octaves, cnp.ndarray[cnp.float64_t, ndim=2] anchor, **kwargs):

    # Ensure memory-contiguous numpy array
    anchor = np.ascontiguousarray(anchor, dtype=np.float64)

    # Convert numpy 2D array to double**
    cdef double* f_anchor = <double*> malloc(3 * num_octaves * sizeof(double))
    to_fortran_2D_double_array(anchor, f_anchor, num_octaves, 3)

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
