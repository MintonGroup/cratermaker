import json
from typing import Union
import numpy as np

float_like = Union[float, int, np.float64, np.float32, np.int64, np.int32]

def to_config(obj):
    """
    Serialize the attributes of an object into a dictionary.

    This function generates a dictionary of serializable attributes of the given object,
    excluding those specified in the object's 'config_ignore' attribute.

    Parameters
    ----------
    obj : object
        The object whose attributes are to be serialized.

    Returns
    -------
    dict
        A dictionary containing the serializable attributes of the object.

    Notes
    -----
    Only attributes that are instances of basic data types (int, float, str, list, dict, bool, None) are included.
    Attributes listed in 'config_ignore' of the object are excluded from serialization.
    """   
    # Check if the object has the attribute 'config_ignore'
    ignores = getattr(obj, 'config_ignore', [])
        
    # Generate a dictionary of serializable attributes, excluding those in 'ignores'
    return {
        k: v for k, v in obj.__dict__.items()
        if isinstance(v, (int, float, str, list, dict, bool, type(None))) and k not in ignores
    }
  
   
def set_properties(obj,**kwargs):
    """
    Set properties of a simulation object from various sources.

    This function sets the properties of a simulation object based on the provided arguments.
    Properties can be read from a JSON file, a pre-defined catalogue, or directly passed as keyword arguments.

    Parameters
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
    3. JSON file (specified by 'filename' key in kwargs).
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
            catalogue = json.load(f)
            
        set_properties_from_catalogue(obj,catalogue=catalogue,name=name)
        set_properties_from_arguments(obj,name=name)
        
    if 'filename' in kwargs:
        set_properties_from_file(obj,**kwargs)
    
    if 'catalogue' in kwargs: 
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
    
            
def create_catalogue(header,values):
    """
    Create and return a catalogue of properties or items based on the given inputs.

    This function generates a catalogue, which could be a collection of properties, configurations,
    or any other set of items, based on the provided arguments.

    Parameters
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
    # Check if it's already a tuple
    if isinstance(location, tuple) and len(location) == 2:
        return location
    
    # Check if it's already a structured array with 'lat' and 'lon'
    if isinstance(location, np.ndarray) and location.dtype.names == ('lat', 'lon'):
        return location
    
    # Check if it's a dictionary with 'lat' and 'lon' keys
    if isinstance(location, dict):
        if "lat" in location and "lon" in location:
            return np.array([(location['lat'], location['lon'])], dtype=[('lat', 'f8'), ('lon', 'f8')])
    
    # Check if it's a tuple, list, or array of the correct shape
    if isinstance(location, (tuple, list, np.ndarray)):
        if len(location) == 2:
            return np.array([(location[0], location[1])], dtype=[('lat', 'f8'), ('lon', 'f8')])
    
    raise ValueError("location must be a dict with 'lat' and 'lon', a 2-element tuple/list, or a structured array with 'lat' and 'lon'")


def normalize_coords(loc):
    """
    Normalize geographic coordinates to ensure longitude is within [-180, 180) degrees 
    and latitude within [-90, 90] degrees.

    This function takes a tuple of longitude and latitude values in degrees, normalizes 
    them to the specified ranges, and handles cases where latitude values exceed the 
    polar extremes, adjusting both latitude and longitude accordingly.

    Parameters
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

