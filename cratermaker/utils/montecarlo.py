import numpy as np
from scipy.stats import truncnorm
def get_random_location(size=None, rng=None):
   """
   Computes random latitude and longitude values.
   
   Generates a set of latitude and longitude values that are uniformly distributed on the surface of a sphere.
   
   Parameters
   ----------
   size : int or tuple of ints, optional
      The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
   rng : numpy.random.Generator, optional
      An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.
   
   Returns
   -------
   ndarray[(lat,lon)]
      An array of pairs of latitude and longitude values.
   """

   if rng and not hasattr(rng, 'uniform'):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng()

   uv_shape = size if size is not None else (1,)
   u = rng.uniform(size=uv_shape)
   v = rng.uniform(size=uv_shape)
   
   # Compute the angles theta and phi
   theta = 2 * np.pi * u
   phi = np.arccos(2 * v - 1)
   
   # Convert to lat/lon
   lat = np.degrees(phi - np.pi / 2)
   lon = np.degrees(theta)
   
   # Reshape lat and lon to the original uv_shape if necessary
   lat = lat.reshape(uv_shape)
   lon = lon.reshape(uv_shape)
   
   # Combine lat and lon into a structured array
   latlon_arr = np.empty(uv_shape, dtype=[('lat', 'float'), ('lon', 'float')])
   latlon_arr['lat'] = lat
   latlon_arr['lon'] = lon
   
   return latlon_arr


def get_random_impact_angle(size=None, rng=None):
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
   ndarray 
      An array of impact angles (in radians).
   """   
   
   if rng and not hasattr(rng, 'uniform'):
      raise TypeError("The 'rng' argument must be a numpy.random.Generator instance or None")
   if rng is None:
      rng = np.random.default_rng()

   u_shape = size if size is not None else (1,)
   u = np.sqrt(rng.uniform(size=u_shape))
   impact_angle = np.arcsin(u)
   
   return impact_angle


def get_random_size(diameters, cdf=None, size=None, rng=None):
   """
   Sample diameters from a cumulative size-frequency distribution (SFD).
   
   Given an array of diameters and optionally a cumulative distribution function (CDF), this function generates new diameter values that follow the specified distribution. The SFD is treated as a continuous function that interpolates between the provided diameters, which are assumed to represent a power-law distribution.
   
   Parameters
   ----------
   diameters : array_like
      An array of diameters from which the SFD is constructed. Must be 1-dimensional.
   cdf : array_like, optional
      The cumulative distribution function corresponding to `diameters`. Must be the same size as `diameters`. If not provided, a CDF is calculated by cumulatively summing the sorted `diameters`.
   size : int or tuple of ints, optional
      The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
   rng : numpy.random.Generator, optional
      An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.
   
   Returns
   -------
   ndarray
      An array of sampled diameter values from the SFD. The size of the array is determined by the `size` parameter.
   
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
   if rng and hasattr(rng, 'uniform'):
      pass
   elif rng is None:  # Just use the basic uniform random number generator
      rng = np.random.default_rng()
   else:
      raise TypeError("The 'rng' argument must be a compatible with numpy random generator or None")

   # Sort the diameters in descending order and get the cumulative distribution if it was not supplied
   if cdf is None:
      sorted_indices = np.argsort(diameters)[::-1]
      sorted_diameters = diameters[sorted_indices]
      sorted_cdf = np.cumsum(sorted_diameters)
   else:
      sorted_indices = np.argsort(diameters)[::-1]
      sorted_diameters = diameters[sorted_indices]
      sorted_cdf = cdf[sorted_indices]

   # Normalize the cdf
   sorted_cdf /= sorted_cdf[-1]

   # Generate uniform random numbers for the entire sample size
   u_shape = size if size is not None else (1,)
   u = rng.uniform(low=sorted_cdf[0], high=sorted_cdf[-1], size=u_shape)

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
   noise = 1e-3 * rng.uniform(size=new_diameters.shape)
   new_diameters *= (1 + noise)

   return new_diameters


def get_random_velocity(vmean, target=None, size=None, rng=None):
   """
   Sample impact velocities from a Maxwell-Boltzmann distribution given a mean velocity.
   
   Parameters 
   ----------
   vmean : float
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
   
   u_shape = size if size is not None else (1,)
   sigma = vmean / np.sqrt(8/np.pi)
   
   vx = rng.normal(0, sigma, size=u_shape)
   vy = rng.normal(0, sigma, size=u_shape)
   vz = rng.normal(0, sigma, size=u_shape)
   velocities = np.sqrt(vx**2 + vy**2 + vz**2)
   return velocities

def bounded_norm(loc,scale):
   lower_bound = loc - scale
   upper_bound = loc + scale
   truncated_normal = truncnorm(
         (lower_bound - loc) / scale,
         (upper_bound - loc) / scale,
         loc=loc, scale=scale
      )

   # Generate a random number from this distribution
   return truncated_normal.rvs(1)