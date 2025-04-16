import yaml
import numpy as np
from numpy.typing import ArrayLike
from cratermaker.utils.custom_types import FloatLike
from typing import Callable, Union, Any
import inspect

class Parameter(property):
    """
    A property that flags whether a parameter is explicitly set by the user.
    This class extends the built-in property to allow tracking of user-set values.
    
    Parameter
    ----------
    fget : function
        The getter function for the property.
    fset : function, optional
        The setter function for the property. Default is None.
    fdel : function, optional
        The deleter function for the property. Default is None.
    doc : str, optional
        The documentation string for the property. Default is None.
    """
    def __init__(self, fget, fset=None, fdel=None, doc=None):
        super().__init__(fget, fset, fdel, doc)
        self.param_name = None  # Will be set via __set_name__

    def __set_name__(self, owner, name):
        self.param_name = name

    def __set__(self, instance, value):
        if self.fset is None:
            raise AttributeError("can't set attribute")
        # Invoke the original setter.
        self.fset(instance, value)
        # Flag that this property was explicitly set by the user.
        if not hasattr(instance, '_user_defined'):
            instance._user_defined = {}
        instance._user_defined[self.param_name] = True

    def setter(self, fset):
        return type(self)(self.fget, fset, self.fdel, self.__doc__)

    def deleter(self, fdel):
        return type(self)(self.fget, self.fset, fdel, self.__doc__)

def parameter(fget=None):
    """
    A decorator to mark a property as a user-settable parameter.
    """
    if fget is None:
        def decorator(fget):
            return Parameter(fget)
        return decorator
    else:
        return Parameter(fget)
   
def set_properties(obj,**kwargs):
    """
    Set properties of a simulation object from various sources.

    This function sets the properties of a simulation object based on the provided arguments.
    Properties can be read from a YAML file, a pre-defined catalogue, or directly passed as keyword arguments.

    Parameter
    ----------
    obj : object
        The simulation object whose properties are to be set.
    **kwargs : dict
        Keyword arguments that can include 'filename', 'catalogue', and other direct property settings.

    Notes
    -----
    The order of property precedence is: 
    1. Direct keyword arguments (kwargs).
    2. Pre-defined catalogue (specified by 'catalogue' key in kwargs).
    3. YAML file (specified by 'filename' key in kwargs).
    Properties set by kwargs override those set by 'catalogue' or 'filename'.
    """
    
    def set_properties_from_arguments(obj, **kwargs):
        for key, value in kwargs.items():
            if hasattr(obj, key) and value is not None:
                setattr(obj, key, value)   
            
    def set_properties_from_catalogue(obj, catalogue, name=None, **kwargs):
        # Check to make sure that the catalogue argument is in the valid nested dict format
        if not isinstance(catalogue, dict):
            raise ValueError("Catalogue must be a dictionary")

        for key, value in catalogue.items():
            if not isinstance(value, dict):
                raise ValueError(f"Value for key '{key}' in catalogue must be a dictionary")
        
        # Look up material in the catalogue
        if name is None:
            if len(catalogue) == 1:
                name = next(iter(catalogue)) 
            else:
                raise ValueError("A name argument must be passed if there is more than one item in the catalogue!")
        
        properties = catalogue.get(name) 
        if properties: # A match was found to the catalogue 
            set_properties_from_arguments(obj, **properties)
        else:
            set_properties_from_arguments(obj, name=name, **kwargs)
            
    def set_properties_from_file(obj, filename, name=None, **kwargs):
        with open(filename, 'r') as f:
            catalogue = yaml.safe_load(f)
            
        set_properties_from_catalogue(obj,catalogue=catalogue,name=name)
        set_properties_from_arguments(obj,name=name)
       
    filename = kwargs.get('filename') 
    if filename:
        set_properties_from_file(obj,**kwargs)
   
    catalogue = kwargs.get('catalogue') 
    if catalogue:
        set_properties_from_catalogue(obj,**kwargs)
        
    set_properties_from_arguments(obj,**kwargs)
    
    if not hasattr(obj,"name"):
        raise ValueError("The object must be given a name")
    
    return


def check_properties(obj):
     # Check for any unset properties
    missing_prop = []
    for property_name, value in obj.__dict__.items():
        if value is None:
            missing_prop.append(property_name)
       
    if len(missing_prop) == 0:
        return     
    elif len(missing_prop) == 1:
        raise ValueError(f"The required property {missing_prop[0]} has not been set")    
    else:
        raise ValueError(f"The following required properties have not been set: {missing_prop}")
    

def to_config(obj) -> dict:
    """
    Serialize the properties of this instance based on whether they have been explicitly set by the user.
    
    Only properties flagged as user-set (by the user_param decorator) are included.
    Other properties (those not flagged) are added only if they have not already been included.

    Parameter
    ---------
    obj : object
        The object whose properties are to be serialized.
    
    Returns
    -------
    dict
        A dictionary containing the selected key/value pairs from properties.
    """
    config = {}
    
    # Include properties decorated as user parameters if they were explicitly set by the user.
    for name, _ in inspect.getmembers(type(obj), lambda o: isinstance(o, Parameter)):
        if getattr(obj, "_user_defined", {}).get(name, False):
            config[name] = getattr(obj, name)
                
    return config


def create_catalogue(header,values):
    """
    Create and return a catalogue of properties or items based on the given inputs.

    This function generates a catalogue, which could be a collection of properties, configurations,
    or any other set of items, based on the provided arguments.

    Parameter
    ----------
    args : various
        The arguments that determine the contents of the catalogue. The type and number of arguments
        can vary based on the intended use of the catalogue.

    Returns
    -------
    catalogue_type
        A catalogue of items or properties. The exact type of this catalogue (e.g., dict, list, custom object)
        depends on the implementation.

    Notes
    -----
    The catalogues built by this function are the built-in catalogues for material properties and target bodie
    """   
    # Create the catalogue dictionary using the class variables
    catalogue = {
        tab[0]: dict(zip(header, tab))
        for tab in values
    }

    # Remove the first key from each dictionary in the catalogue
    for k in list(catalogue):
        del catalogue[k][header[0]]

    return catalogue 


def validate_and_convert_location(location):
    """
    Validate and convert a given location into a standard structured format.

    This function checks the input location data and converts it into a 
    consistent structured array format if it is a valid location representation.
    Valid formats for location include a tuple, a dictionary, or a structured 
    array with latitude ('lon') and longitude ('lat').

    Parameter
    ----------
    location : tuple, dict, ArrayLike
        The input location data. It can be:
        - A tuple with two elements (longitude, latitude).
        - A dictionary with keys 'lon' and 'lat'.
        - A structured numpy array with 'lon' and 'lon' fields.
        - A list or an unstructured numpy array with two elements or an array of pairs of elements.

    Returns
    -------
    ArrayLike
        A structured numpy array with the location data in the format 
        [('lon', 'f8'), ('lat', 'f8')].

    Raises
    ------
    ValueError
        If the input does not conform to one of the expected formats for location data.

    Examples
    --------
    >>> validate_and_convert_location((-120.0, -45.0))
    array([(-120., 45.)], dtype=[('lon', '<f8'), ('lat', '<f8')])

    >>> validate_and_convert_location({'lon': -120.0, 'lat': -45.0})
    array([(-120, 45.)], dtype=[('lon', '<f8'), ('lat', '<f8')])

    >>> validate_and_convert_location(np.array([(45.0, -120.0)], dtype=[('lon', 'f8'), ('lat', 'f8')]))
    array([(-120., 45.)], dtype=[('lon', '<f8'), ('lat', '<f8')])

    Notes
    -----
    The function ensures that the output is always a structured numpy array with 
    'lon' and 'lat' fields for consistent handling of location data across different
    input formats.
    """    
    # Check if it's already a tuple
    if isinstance(location, tuple) and len(location) == 2:
        return location
    
    # Check if it's already a structured array with 'lon' and 'lat'
    if isinstance(location, ArrayLike) and location.dtype.names == ('lon', 'lat'):
        return location
    
    # Check if it's a dictionary with 'lon' and 'lat' keys
    if isinstance(location, dict):
        if "lon" in location and "lat" in location:
            return np.array([(location['lon'], location['lon'])], dtype=[('lon', 'f8'), ('lat', 'f8')])
    
    # Check if it's a tuple, list, or array of the correct shape
    if isinstance(location, (tuple, list, ArrayLike)):
        if len(location) == 2:
            return np.array([(location[0], location[1])], dtype=[('lon', 'f8'), ('lat', 'f8')])
    
    raise ValueError("location must be a dict with 'lon' and 'lat', a 2-element tuple/list, or a structured array with 'lon' and 'lat'")


def normalize_coords(loc):
    """
    Normalize geographic coordinates to ensure longitude is within [-180, 180) degrees 
    and latitude within [-90, 90] degrees.

    This function takes a tuple of longitude and latitude values in degrees, normalizes 
    them to the specified ranges, and handles cases where latitude values exceed the 
    polar extremes, adjusting both latitude and longitude accordingly.

    Parameter
    ----------
    loc : tuple
        A tuple containing two elements: (longitude, latitude) in degrees. 
        Longitude and latitude can be any float values.

    Returns
    -------
    tuple
        A tuple of two elements: (normalized_longitude, normalized_latitude). 
        The normalized longitude is in the range [-180, 180) degrees, and the 
        normalized latitude is in the range [-90, 90] degrees.

    Notes
    -----
    - The longitude is normalized using a modulo operation with 360 degrees and then adjusted to the range [-180, 180).
    - Latitude values beyond 90 or below -90 degrees are adjusted by reflecting them within the range and flipping the longitude by 180 degrees, then re-normalizing it to the [-180, 180) range.

    Examples
    --------
    >>> normalize_coords((370, 95))
    (10.0, 85.0)

    >>> normalize_coords((-185, -100))
    (-5.0, 80.0)
    """
    lon, lat = loc

    # Normalize longitude to be within [-180, 180)
    normalized_lon = ((lon + 180) % 360) - 180

    # Normalize latitude
    if lat > 90:
        normalized_lat = 180 - lat 
        normalized_lon = lon - 180 # Flip the longitude 
    elif lat < -90:
        normalized_lat = -180 - lat 
        normalized_lon = lon - 180 # Flip the longitude
    else:
        normalized_lat = lat

    # Ensure latitude is within the range [-90, 90] after adjustments
    normalized_lat = np.clip(normalized_lat, -90, 90)

    return normalized_lon, normalized_lat


def R_to_CSFD(R: Callable[[Union[FloatLike, ArrayLike]], Union[FloatLike, ArrayLike]], 
              D: Union[FloatLike, ArrayLike],
              Dlim: FloatLike = 1e6,
              *args: Any,
            ) -> Union[FloatLike, ArrayLike]:
    """
    Convert R values to cumulative N values for a given D using the R-plot function.

    Parameter
    ----------
    R : R = f(D) 
        A function that computes R given D.
    D : FloatLike or ArrayLike
        diameter in units of km.
    Dlim : FloatLike
        Upper limit on the diameter over which to evaluate the integral
    *args : Any
        Additional arguments to pass to the R function

    Returns
    -------
    float or ArrayLike
        The cumulative number of craters greater than D in diameter.
    """
    
    def _R_to_CSFD_scalar(R, D, Dlim, *args):
        # Helper function to integrate the R function
        def integrand(D):
            return R(D,*args) / D**3  # This is dN/dD
        
        N = 0.0
        D_i = D
        while D_i < Dlim:
            D_next = D_i * np.sqrt(2.0) 
            D_mid = (D_i + D_next) / 2  # Mid-point of the bin
            bin_width = D_next - D_i
            R_value = integrand(D_mid)
            N += R_value * bin_width
            D_i = D_next  # Move to the next bin
    
        return N
    
    return _R_to_CSFD_scalar(R, D, Dlim, *args) if np.isscalar(D) else np.vectorize(_R_to_CSFD_scalar)(R, D, Dlim, *args)

