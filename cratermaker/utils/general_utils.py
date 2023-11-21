import json
from dataclasses import fields
from typing import Union
import numpy as np

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
            if hasattr(obj, key):
                setattr(obj, key, value)   
            
    def set_properties_from_catalogue(obj, catalogue, key):
        # Look up material in the catalogue
            
        properties = catalogue.get(key) 
        if properties: # A match was found to the catalogue 
            set_properties_from_arguments(obj, **properties)
            
    def set_properties_from_file(obj, filename):
        with open(filename, 'r') as f:
            properties =  json.load(f)
            n = obj.__class__.__name__.lower()
            if n in properties:
                set_properties_from_arguments(obj,**properties[n])
        
    if 'filename' in kwargs:
        set_properties_from_file(obj,filename=kwargs['filename'])
    
    if 'catalogue' in kwargs and 'key' in kwargs:
        set_properties_from_catalogue(obj,kwargs['catalogue'],kwargs['key'])
        
    set_properties_from_arguments(obj,**kwargs)
    
    # Check for any unset properties
    for property_name, value in obj.__dict__.items():
        if value is None:
            raise ValueError(f"The property {property_name} has not been set!")    
    
    return
    
            
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
