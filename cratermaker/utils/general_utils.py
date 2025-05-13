from pathlib import Path
from typing import Any, Callable
from warnings import warn

import numpy as np
import yaml
from numpy.typing import ArrayLike

from cratermaker.constants import FloatLike


class Parameter(property):
    """
    A property descriptor that tracks user-defined properties.  This class is a subclass of the built-in property class and is used
    to create properties in a class that can be set and retrieved. It also tracks whether the property has been set by the user,
    allowing for parameters to be exported to a YAML configuration file.
    """

    def __init__(self, fget, fset=None, fdel=None, doc=None):
        super().__init__(fget, fset, fdel, doc)
        self.name = fget.__name__

    def setter(self, fset):
        def wrapped(instance, value):
            if not hasattr(instance, "_user_defined"):
                instance._user_defined = set()
            instance._user_defined.add(self.name)
            fset(instance, value)

        return Parameter(self.fget, wrapped, self.fdel, self.__doc__)


def parameter(fget=None):
    """
    A decorator to mark a property as a user-settable parameter.
    Can be used with or without parentheses.
    """
    if fget is None:

        def decorator(fget):
            return Parameter(fget)

        return decorator
    else:
        return Parameter(fget)


def _set_properties(
    obj,
    catalogue: dict | None = None,
    key: str | None = None,
    config_file: str | Path | None = None,
    **kwargs: Any,
):
    """
    Set properties of a simulation object from various sources.

    This function sets the properties of a simulation object based on the provided arguments.
    Properties can be read from a YAML file, a pre-defined catalogue, or directly passed as keyword arguments.

    Parameter
    ----------
    obj : object
        The simulation object whose properties are to be set.
    catalogue : dict, optional
        A dictionary representing a catalogue of properties. It must be in the form of a nested dict. If provided, it will be used to set properties.
    key : str, optional
        The key to look up in the catalogue. It must be provided if the catalogue is provided.
    config_file : str or Path, optional
        The path to a YAML file containing properties. If provided, it will be used to set properties.
    **kwargs : dict
        Keyword arguments that can include 'config_file', 'catalogue', and other direct property settings.

    Returns
    -------
    matched : dict
        A dictionary of properties that were successfully set on the object.
    unmatched : dict
        A dictionary of properties that were not set, either due to being None or not matching any known properties.

    Notes
    -----
    The order of property precedence is:
    1. Direct keyword arguments (kwargs).
    2. Pre-defined catalogue (specified by 'catalogue' key in kwargs).
    3. YAML file (specified by 'config_file' key in kwargs).
    Properties set by kwargs override those set by 'catalogue' or 'config_file'.
    """

    def _set_properties_from_arguments(obj, **kwargs):
        matched = {}
        unmatched = {}
        cls = type(obj)
        for key, value in kwargs.items():
            if value is None:
                continue
            param = getattr(cls, key, None)
            if (
                isinstance(param, (property, Parameter))
                and getattr(param, "fset", None) is not None
            ):
                setattr(obj, key, value)
                matched[key] = value
            else:
                unmatched[key] = value
        return matched, unmatched

    def _set_properties_from_catalogue(obj, catalogue, key, **kwargs):
        if "catalogue_key" in dir(obj):
            catalogue_key = getattr(obj, "catalogue_key")
        else:
            raise ValueError(
                "The object does not have a catalogue_key property, and therefore is not set up to receive catalogue entries."
            )
        if catalogue_key in kwargs:
            key = kwargs.pop(catalogue_key)

        if not isinstance(catalogue, dict):
            raise ValueError("Catalogue must be a dictionary")

        for k, v in catalogue.items():
            if not isinstance(v, dict):
                raise ValueError(
                    f"Value for key '{k}' in catalogue must be a dictionary"
                )

        if key not in catalogue:
            return {}, {}

        properties = catalogue.get(key)
        properties.update({catalogue_key: key})
        # Remove any items in kwargs that are already in properties
        for k in properties.keys():
            if k in kwargs:
                del kwargs[k]
        if properties:  # A match was found to the catalogue
            matched, unmatched = _set_properties_from_arguments(
                obj, **properties, **kwargs
            )
            properties.pop(
                catalogue_key
            )  # Make sure that the catlogue key doesn't stay in the properties
        return matched, unmatched

    def _set_properties_from_file(obj, config_file, key=None, **kwargs):
        try:
            with open(config_file, "r") as f:
                properties = yaml.safe_load(f)
        except Exception as e:
            warn(f"Could not read the file {config_file}.\n{e}", RuntimeWarning)
            return {}, {}
        merged = {**properties, **{k: v for k, v in kwargs.items() if v is not None}}
        if key is None:
            matched, unmatched = _set_properties_from_arguments(obj, **merged)
        else:
            if key not in properties:
                raise ValueError(f"Key '{key}' not found in the file '{config_file}'.")
            matched, unmatched = _set_properties_from_catalogue(
                obj, key=key, catalogue=properties, **kwargs
            )
        return matched, unmatched

    matched = {}
    unmatched = {}
    if config_file:
        m, u = _set_properties_from_file(
            obj, config_file=config_file, key=key, **kwargs
        )
        matched.update(m)
        unmatched.update(u)

    if catalogue:
        m, u = _set_properties_from_catalogue(
            obj, catalogue=catalogue, key=key, **kwargs
        )
        matched.update(m)
        unmatched.update(u)

    m, u = _set_properties_from_arguments(obj, **kwargs)
    matched.update(m)
    unmatched.update(u)

    # if there are any keys in unmatched that are also present in matched, remove them from unmatched
    for key in matched.keys():
        if key in unmatched:
            del unmatched[key]

    return matched, unmatched


def _create_catalogue(header, values):
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
    catalogue = {tab[0]: dict(zip(header, tab)) for tab in values}

    # Remove the first key from each dictionary in the catalogue
    for k in list(catalogue):
        del catalogue[k][header[0]]

    return catalogue


def normalize_coords(location: tuple[FloatLike, FloatLike]) -> tuple[float, float]:
    """
    Normalize geographic coordinates to ensure longitude is within [-180, 180) degrees
    and latitude within [-90, 90] degrees.

    This function takes a tuple of longitude and latitude values in degrees, normalizes
    them to the specified ranges, and handles cases where latitude values exceed the
    polar extremes, adjusting both latitude and longitude accordingly.

    Parameter
    ----------
    location : tuple
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
    lon, lat = location

    # Normalize longitude to be within [-180, 180)
    normalized_lon = ((lon + 180) % 360) - 180

    # Normalize latitude
    if lat > 90:
        normalized_lat = 180 - lat
        normalized_lon = lon - 180  # Flip the longitude
    elif lat < -90:
        normalized_lat = -180 - lat
        normalized_lon = lon - 180  # Flip the longitude
    else:
        normalized_lat = lat

    # Ensure latitude is within the range [-90, 90] after adjustments
    normalized_lat = np.clip(normalized_lat, -90, 90)

    return float(normalized_lon), float(normalized_lat)


def validate_and_normalize_location(location):
    """
    Validate and normalize a given location into a standard structured format.

    This function checks the input location data and converts it into a
    consistent structured array format if it is a valid location representation.
    Valid formats for location include a tuple, a dictionary, or a structured
    array with longitude ('lon') and latitude ('lat').

    Parameters
    ----------
    location : tuple, dict, ArrayLike
        The input location data. It can be:
        - A tuple, list, or array with two elements (longitude, latitude).
        - A dictionary with keys 'lon' and 'lat'.
        - A structured numpy array with 'lon' and 'lat' fields.

    Returns
    -------
    tuple
        longitude and latitude as a tuple of floats.

    Raises
    ------
    ValueError
        If the input does not conform to one of the expected formats for location data.

    Examples
    --------
    >>> validate_and_normalize_location((370, 95))
    (10.0, 85.0))

    >>> validate_and_normalize_location({'lat': 45.0, 'lon': 120.0})
    (-120., 45.)

    >>> validate_and_normalize_location(np.array([(-120.0, 45.0)], dtype=[('lon', 'f8'), ('lat', 'f8')]))
    (-120., 45.)

    """
    # Check if it's already a tuple

    if isinstance(location, np.ndarray) and location.dtype.names == ("lon", "lat"):
        return normalize_coords((location[0], location[1]))

    if isinstance(location, np.ndarray) and location.dtype.names == ("lat", "lon"):
        return normalize_coords((location[1], location[0]))

    if isinstance(location, (tuple, list, np.ndarray)) and len(location) == 2:
        return normalize_coords(location)

    # Check if it's a dictionary with 'lon' and 'lat' keys
    if isinstance(location, dict):
        if "lon" in location and "lat" in location:
            return normalize_coords((location["lon"], location["lat"]))

    if len(location) == 2:
        return normalize_coords((location[0], location[1]))

    raise ValueError(
        "location must a tuple, list, or ArrayLike of len==2, a dict with 'lon' and 'lat', or a structured array with 'lon' and 'lat' names"
    )


def R_to_CSFD(
    R: Callable[[FloatLike | ArrayLike], FloatLike | ArrayLike],
    D: FloatLike | ArrayLike,
    Dlim: FloatLike = 1e6,
    *args: Any,
) -> FloatLike | ArrayLike:
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
    args : Any
        Additional arguments to pass to the R function

    Returns
    -------
    float or ArrayLike
        The cumulative number of craters greater than D in diameter.
    """

    def _R_to_CSFD_scalar(R, D, Dlim, *args):
        # Helper function to integrate the R function
        def integrand(D):
            return R(D, *args) / D**3  # This is dN/dD

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

    return (
        _R_to_CSFD_scalar(R, D, Dlim, *args)
        if np.isscalar(D)
        else np.vectorize(_R_to_CSFD_scalar)(R, D, Dlim, *args)
    )


def format_large_units(
    value: float, threshold: float = 1000.0, quantity: str = "length"
) -> str:
    """
    Format a value and automatically shift units based on threshold.
    """
    if quantity == "length":
        units = ["m", "km"]
    elif quantity == "velocity":
        units = ["m/s", "km/s"]
    elif quantity == "time":
        units = ["My", "Gy"]
    elif quantity == "pressure":
        units = ["Pa", "kPa", "MPa", "GPa"]

    if value is None:
        return "N/A"

    unit_index = 0
    while unit_index + 1 < len(units) and value >= threshold:
        value /= threshold
        unit_index += 1

    if value >= 100:
        fmt = "{:.0f} {}"
    elif value >= 10:
        fmt = "{:.1f} {}"
    elif value >= 1:
        fmt = "{:.2f} {}"
    else:
        fmt = "{:.3g} {}"
    return fmt.format(value, units[unit_index])
