# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
from cpython cimport PyUnicode_AsUTF8AndSize, PyUnicode_FromString
cimport numpy as cnp
import numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset 
from libc.stdint cimport uintptr_t
from uxarray import UxDataset

cdef extern from "fortran_bind.h":
    ctypedef struct surface_type:
        double *elevation
        int elevation_shape[1]
        double *node_x
        int node_x_shape[1]
        double *node_y
        int node_y_shape[1]
        double *node_z
        int node_z_shape[1]
        double *face_x
        int face_x_shape[1]
        double *face_y
        int face_y_shape[1]
        double *face_z
        int face_z_shape[1]

    surface_type* bind_surface_init(int n_node, int n_face)
    void bind_surface_final(surface_type *obj)

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

    double bind_perlin_noise(const char *model, 
                             double *x, 
                             double *y, 
                             double *z, 
                             int num_octaves, 
                             double *anchor, 
                             double damp, 
                             double damp0, 
                             double damp_scale, 
                             double freq, 
                             double gain, 
                             double gain0, 
                             double lacunarity, 
                             double noise_height, 
                             double pers, 
                             double slope, 
                             double warp, 
                             double warp0)


cdef cnp.ndarray get_array(double* array, int length):
    """
    Generic function to create a NumPy array from a Fortran-allocated array.

    Parameters
    ----------
    array : double*
        Pointer to the Fortran-allocated array.
    length : int
        Length of the array.

    Returns
    -------
    numpy_array : cnp.ndarray
        A NumPy array view of the Fortran array.
    """
    if array is NULL or length <= 0:
        raise ValueError("Invalid array or length.")
    cdef cnp.npy_intp size = length
    return cnp.PyArray_SimpleNewFromData(1, &size, cnp.NPY_FLOAT64, <void*>array)


cdef void set_array(double* array, cnp.ndarray[cnp.float64_t, ndim=1] input_array, int length):
    """
    Generic function to set a Fortran-allocated array from a NumPy array.

    Parameters
    ----------
    array : double*
        Pointer to the Fortran-allocated array.
    input_array : cnp.ndarray
        The NumPy array to copy data from.
    length : int
        Length of the Fortran array.
    """
    if array is NULL or length <= 0:
        raise ValueError("Invalid array or length.")
    if input_array.shape[0] != length:
        raise ValueError("Length of input NumPy array does not match the length of the Fortran array.")
    for i in range(input_array.shape[0]):
        array[i] = input_array[i]


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

    def __cinit__(self, surf: UxDataset):
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

        cdef cnp.npy_intp n_node = surf.uxgrid.n_node
        cdef cnp.npy_intp n_face = surf.uxgrid.n_face

        self.fobj = bind_surface_init(n_node, n_face)

        # Do some basic checks to make sure the object variable and all its components were allocated succesfully
        if self.fobj is NULL:
            raise MemoryError("Failed to allocate Fortran object.")

        if self.fobj.elevation is NULL:
            raise MemoryError("Failed to allocate component variable 'elevation' in the Fortran object.")
        
        if self.fobj.node_x is NULL:
            raise MemoryError("Failed to allocate component variable 'node_x' in the Fortran object.")

        if self.fobj.node_y is NULL:
            raise MemoryError("Failed to allocate component variable 'node_y' in the Fortran object.")

        if self.fobj.node_z is NULL:
            raise MemoryError("Failed to allocate component variable 'node_z' in the Fortran object.")

        if self.fobj.face_x is NULL:
            raise MemoryError("Failed to allocate component variable 'face_x' in the Fortran object.")
        
        if self.fobj.face_y is NULL:
            raise MemoryError("Failed to allocate component variable 'face_y' in the Fortran object.")
        
        if self.fobj.face_z is NULL:
            raise MemoryError("Failed to allocate component variable 'face_z' in the Fortran object.")

        # Explicitly set the length of the allocatable arrays
        self.fobj.elevation_shape[0] = n_node
        self.fobj.node_x_shape[0] = n_node
        self.fobj.node_y_shape[0] = n_node
        self.fobj.node_z_shape[0] = n_node
        self.fobj.face_x_shape[0] = n_face
        self.fobj.face_y_shape[0] = n_face
        self.fobj.face_z_shape[0] = n_face

        print("Elevation length from Cython:", self.fobj.elevation_shape[0])
        print("Elevation pointer address from Cython:", <unsigned long>self.fobj.elevation)
        print("Node_x length from Cython:", self.fobj.node_x_shape[0])
        print("Node_x pointer address from Cython:", <unsigned long>self.fobj.node_x)
        print("Node_y length from Cython:", self.fobj.node_y_shape[0])
        print("Node_y pointer address from Cython:", <uintptr_t>self.fobj.node_y)
        print("Node_z length from Cython:", self.fobj.node_z_shape[0])
        print("Node_z pointer address from Cython:", <uintptr_t>self.fobj.node_z)
        print("Face_x length from Cython:", self.fobj.face_x_shape[0])
        print("Face_x pointer address from Cython:", <uintptr_t>self.fobj.face_x)
        print("Face_y length from Cython:", self.fobj.face_y_shape[0])
        print("Face_y pointer address from Cython:", <uintptr_t>self.fobj.face_y)
        print("Face_z length from Cython:", self.fobj.face_z_shape[0])
        print("Face_z pointer address from Cython:", <uintptr_t>self.fobj.face_z)

        # Copy over the data using the setter methods
        self.elevation = surf.elevation.values
        self.node_x = surf.uxgrid.node_x.values
        self.node_y = surf.uxgrid.node_y.values
        self.node_z = surf.uxgrid.node_z.values
        self.face_x = surf.uxgrid.face_x.values
        self.face_y = surf.uxgrid.face_y.values
        self.face_z = surf.uxgrid.face_z.values
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

    @property
    def elevation(self):
        """
        A getter method that retrieves the elevation allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.elevation, self.fobj.elevation_shape[0])

    @elevation.setter
    def elevation(self, cnp.ndarray[cnp.float64_t, ndim=1] elevation_array):
        """
        Property to set the elevation data.
        """
        set_array(self.fobj.elevation, elevation_array, self.fobj.elevation_shape[0])

    @property
    def node_x(self):
        """
        A getter method that retrieves the node_x allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.node_x, self.fobj.node_x_shape[0])
    
    @node_x.setter
    def node_x(self, cnp.ndarray[cnp.float64_t, ndim=1] node_x_array):
        """
        Property to set the node_x data.
        """
        set_array(self.fobj.node_x, node_x_array, self.fobj.node_x_shape[0])

    @property
    def node_y(self):
        """
        A getter method that retrieves the node_y allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.node_y, self.fobj.node_y_shape[0])
    
    @node_y.setter
    def node_y(self, cnp.ndarray[cnp.float64_t, ndim=1] node_y_array):
        """
        Property to set the node_y data.
        """
        set_array(self.fobj.node_y, node_y_array, self.fobj.node_y_shape[0])
    
    @property
    def node_z(self):
        """
        A getter method that retrieves the node_z allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.node_z, self.fobj.node_z_shape[0])
    
    @node_z.setter
    def node_z(self, cnp.ndarray[cnp.float64_t, ndim=1] node_z_array):
        """
        Property to set the node_z data.
        """
        set_array(self.fobj.node_z, node_z_array, self.fobj.node_z_shape[0])

    @property
    def face_x(self):
        """
        A getter method that retrieves the face_x allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.face_x, self.fobj.face_x_shape[0])
    
    @face_x.setter
    def face_x(self, cnp.ndarray[cnp.float64_t, ndim=1] face_x_array):
        """
        Property to set the face_x data.
        """
        set_array(self.fobj.face_x, face_x_array, self.fobj.face_x_shape[0])
    
    @property
    def face_y(self):
        """
        A getter method that retrieves the face_y allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.face_y, self.fobj.face_y_shape[0])
    
    @face_y.setter
    def face_y(self, cnp.ndarray[cnp.float64_t, ndim=1] face_y_array):
        """
        Property to set the face_y data.
        """
        set_array(self.fobj.face_y, face_y_array, self.fobj.face_y_shape[0])
    
    @property
    def face_z(self):
        """
        A getter method that retrieves the face_z allocatable array from Fortran and returns it as a Numpy array
        """
        return get_array(self.fobj.face_z, self.fobj.face_z_shape[0])
    
    @face_z.setter
    def face_z(self, cnp.ndarray[cnp.float64_t, ndim=1] face_z_array):
        """
        Property to set the face_z data.
        """
        set_array(self.fobj.face_z, face_z_array, self.fobj.face_z_shape[0])
    
    
def util_perlin(str model, 
                cnp.ndarray[cnp.float64_t, ndim=1] x, 
                cnp.ndarray[cnp.float64_t, ndim=1] y, 
                cnp.ndarray[cnp.float64_t, ndim=1] z, 
                int num_octaves, 
                cnp.ndarray[cnp.float64_t, ndim=2] anchor, 
                **kwargs):

    # Ensure memory-contiguous numpy array
    anchor = np.ascontiguousarray(anchor, dtype=np.float64)

    # Convert numpy 2D array to double**
    cdef double* f_anchor = <double*> malloc(3 * num_octaves * sizeof(double))
    cdef double* f_x = <double*> malloc(x.shape[0] * sizeof(double))
    cdef double* f_y = <double*> malloc(y.shape[0] * sizeof(double))
    cdef double* f_z = <double*> malloc(z.shape[0] * sizeof(double))

    to_fortran_2D_double_array(anchor, f_anchor, num_octaves, 3)

    set_array(f_x, x, x.shape[0])
    set_array(f_y, y, y.shape[0])
    set_array(f_z, z, z.shape[0])

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
            noise = bind_perlin_noise(model, f_x, f_y, f_z, num_octaves, f_anchor, kw['damp'], kw['damp0'], kw['damp_scale'], kw['freq'], kw['gain'], kw['gain0'], kw['lacunarity'], kw['noise_height'], kw['pers'], kw['slope'], kw['warp'], kw['warp0'] )
        else:
            missing_args = set(required_kwargs[model]) - kwargs.keys()
            raise ValueError(f"The {model} model requires the following missing keywords: {', '.join(missing_args)}")
    else:
        raise ValueError(f"{model} is an invalid model. Valid options are {', '.join(required_kwargs.keys())}")

    # Free the allocated memory
    free(f_anchor)
    free(f_x)
    free(f_y)
    free(f_z)

    return noise  # Return the noise value


