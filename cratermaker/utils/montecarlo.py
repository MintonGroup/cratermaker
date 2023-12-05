import numpy as np
from numpy.random import Generator
from typing import Union, Optional, Tuple
from numpy.typing import NDArray
from scipy.stats import truncnorm

def get_random_location(size: Optional[Union[int, Tuple[int, ...]]]=1, rng: Optional[Generator]=None) -> Union[np.float64, Tuple[np.float64, np.float64], np.ndarray]:
    """
    Computes random longitude and latitude values.
    
    Generates a set of latitude and longitude values that are uniformly distributed on the surface of a sphere.
    
    Parameters
    ----------
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator, optional
        An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.
    
    Returns
    -------
    (lon,lat) or ndarray[(lon,lat)] of given size
        A pair or array of pairs of longitude and latitude values in radians.
    """

    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 

    u = rng.uniform(size=size)
    v = rng.uniform(size=size)
    
    # Compute the angles theta and phi
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    
    # Convert to lon/lat
    lon = theta
    lat = phi - np.pi / 2.0
    
    if size == 1: 
        return (np.float64(lon.item()),np.float64(lat.item()))
    else:
        # Reshape lat and lon to the original size if necessary
        lon = lon.reshape(size)
        lat = lat.reshape(size)
  
        # Combine lat and lon into a structured array
        lonlat_arr = np.empty(size, dtype=[('lon', 'float'), ('lat', 'float')])
        lonlat_arr['lon'] = lon
        lonlat_arr['lat'] = lat
    
    return lonlat_arr


def get_random_impact_angle(size: Optional[Union[int, Tuple[int, ...]]]=1, rng: Optional[Generator]=None) -> Union[np.float64,NDArray[np.float64]]:
    """
    Sample impact angles from a distribution centered on 45deg.
    
    For the theory, see Shoemaker (1962) "Interpretation of lunar craters."
    
    Parameters 
    ----------
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator, optional
        An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.     
        
    Returns
    ----------
    np.float64 or ndarray of np.float 64 
        A scalar or array of impact angles (in degrees).
    """    
    
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 

    u = np.sqrt(rng.uniform(size=size))
    impact_angle = np.arcsin(u)
    if size == 1:
        return np.rad2deg(impact_angle[0])
    else:
        return np.rad2deg(impact_angle)


def get_random_size(diameters: NDArray[np.float64], cdf: NDArray[np.float64], size: Optional[Union[int, Tuple[int, ...]]]=1, rng: Optional[Generator]=None) -> Union[np.float64,NDArray[np.float64]]:
    """
    Sample diameters from a cumulative size-frequency distribution (SFD).
    
    Given an array of diameters and optionally a cumulative distribution function (CDF), this function generates new diameter values that follow the specified distribution. The SFD is treated as a continuous function that interpolates between the provided diameters, which are assumed to represent a power-law distribution.
    
    Parameters
    ----------
    diameters : array_like
        An array of diameters from which the SFD is constructed. Must be 1-dimensional.
    cdf : array_like
        The cumulative distribution function corresponding to `diameters`. Must be the same size as `diameters`. If not provided, a CDF is calculated by cumulatively summing the sorted `diameters`.
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator, optional
        An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.
    
    Returns
    -------
    np.float64 or ndarray of np.float 64 
        A scalar or array of sampled diameter values from the SFD. The size of the array is determined by the `size` parameter.
    
    Notes
    -----
    The SFD is assumed to be a continuous distribution that follows a power-law between the provided discrete diameter values. Linear interpolation in log-space is used to sample new values between the known diameters.
    A small amount of random noise (of the order 1e-3 the diameter value) is added to the final diameter value to ensure that diameter values are unlikely to be identical, even when an input CDF could lead to identical diameter values.
    
    Examples
    --------
    >>> diameters = np.array([100.,  56.,  32.,  18.,  10.])
    >>> ncumul = np.array([1.  , 0.51, 0.21, 0.06, 0.01])
    >>> sample_from_sfd(diameters, cdf=ncumul, size=4)
    array([14.80803668, 44.95292261, 29.80797715, 23.11082091])
    
    See Also
    --------
    numpy.random.Generator : The numpy random generator class used for random sampling.
    
    """

    # Check if rng has 'uniform' method which is a characteristic of numpy's random generator objects and use that to generate our values
    if rng and not isinstance(rng, Generator):
        raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
    if rng is None:
        rng = np.random.default_rng() 

    # Sort the diameters in descending order and get the cumulative distribution if it was not supplied
    sorted_indices = np.argsort(diameters)[::-1]
    sorted_diameters = diameters[sorted_indices]
    sorted_cdf = cdf[sorted_indices]

    # Normalize the cdf
    sorted_cdf /= sorted_cdf[-1]

    # Generate uniform random numbers for the entire sample size
    u = rng.uniform(low=sorted_cdf[0], high=sorted_cdf[-1], size=size)

    # Handle the situation where u is a scalar (size is None and diameters is a scalar)
    if np.isscalar(u):
        u = np.array([u])
        
    # Flatten u to work with it as a 1D array
    original_shape = u.shape
    u = u.flatten()     

    # Find the indices where the random numbers would be inserted to maintain order
    indices = np.searchsorted(sorted_cdf, u, side="left")

    # Initialize the new_diameters array
    new_diameters = np.empty(u.shape)

    # Handle edge cases
    new_diameters[indices == 0] = sorted_diameters[0]
    new_diameters[indices == len(sorted_diameters)] = sorted_diameters[-1]

    # Handle non-edge cases by linear interpolation in log space for other values
    not_edge = (indices > 0) & (indices < len(sorted_diameters))
    log_diam_low = np.log(sorted_diameters[indices[not_edge] - 1])
    log_diam_high = np.log(sorted_diameters[indices[not_edge]])
    fractions = (u[not_edge] - sorted_cdf[indices[not_edge] - 1]) / (
        sorted_cdf[indices[not_edge]] - sorted_cdf[indices[not_edge] - 1])
    log_diam_interp = log_diam_low + fractions * (log_diam_high - log_diam_low)
    new_diameters[not_edge] = np.exp(log_diam_interp)

    # Reshape new_diameters to the original shape of u
    new_diameters = new_diameters.reshape(original_shape)

    # Add a small random noise to the diameters
    noise = 1e-8 * rng.uniform(size=new_diameters.shape)
    new_diameters *= (1 + noise)
    if size == 1:
        return new_diameters[0]
    else:
        return new_diameters


def get_random_velocity(vmean: np.float64, size: Optional[Union[int, Tuple[int, ...]]]=1, rng: Optional[Generator]=None) -> Union[np.float64,NDArray[np.float64]]:
    """
    Sample impact velocities from a Maxwell-Boltzmann distribution given a mean velocity.
    
    Parameters 
    ----------
    vmean : np.float64
        The mean velocity of the distribution.
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator, optional
        An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.    
       
    Returns
    ----------
    ndarray 
        An array of impact velocities (in m/s).
    """
    
    # Check if rng has 'uniform' method which is a characteristic of numpy's random generator objects and use that to generate our values
    if rng and hasattr(rng, 'normal'):
        pass
    elif rng is None:  # Just use the basic normal random number generator
        rng = np.random.default_rng()
    else:
        raise TypeError("The 'rng' argument must be a compatible with numpy random generator or None")
    
    sigma = vmean / np.sqrt(8/np.pi)
    
    vx = rng.normal(0, sigma, size=size)
    vy = rng.normal(0, sigma, size=size)
    vz = rng.normal(0, sigma, size=size)
    velocities = np.sqrt(vx**2 + vy**2 + vz**2)
 
    if size == 1:
        return velocities[0]
    else:
        return velocities


def bounded_norm(loc: np.float64,scale: np.float64,size: Optional[Union[int, Tuple[int, ...]]]=1):
    """
    Sample from a truncated normal distribution that is bounded by 1-sigma stdev
    
    Parameters 
    ----------
    loc : float
       mean of the distribution
    scale : float
       standard deviation and bounds of the distribution
       
    Returns
    ----------
    float
       Truncated norm bounded by loc-scale, loc+scale
    """    
    
    lower_bound = loc - scale
    upper_bound = loc + scale
    truncated_normal = truncnorm(
          (lower_bound - loc) / scale,
            (upper_bound - loc) / scale,
            loc=loc, scale=scale
        )
    
    if size == 1:
        return truncated_normal.rvs(1)[0]
    else:
        return truncated_normal.rvs(size)