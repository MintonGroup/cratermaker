import numpy as np
from numpy.random import Generator
from numpy.typing import ArrayLike
from typing import Union, Optional, Tuple
from numpy.typing import NDArray
from scipy.stats import truncnorm
from scipy.stats import maxwell

def get_random_location(
                        size: int | Tuple[int, ...]=1, 
                        rng: Generator | None=None
                        ) -> Union[np.float64, Tuple[np.float64, np.float64], ArrayLike]:
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
        A pair or array of pairs of longitude and latitude values in degrees.
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
    lon = np.rad2deg(theta - np.pi) # Use the convention that longitude is in the range [-180, 180]
    lat = np.rad2deg(phi - np.pi / 2.0)
    
    if size == 1: 
        return np.float64(lon.item()),np.float64(lat.item())
    else:
        # Reshape lat and lon to the original size if necessary
        lon = lon.reshape(size)
        lat = lat.reshape(size)
  
        # Combine lat and lon into a structured array
        lonlat_arr = np.empty(size, dtype=[('lon', 'float'), ('lat', 'float')])
        lonlat_arr['lon'] = lon
        lonlat_arr['lat'] = lat
    
    return lonlat_arr


def get_random_impact_angle(
                            size: int | Tuple[int, ...]=1, 
                            rng: Generator | None=None
                            ) -> Union[np.float64,NDArray[np.float64]]:
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


def get_random_size(diameters: NDArray[np.float64], 
                    cdf: NDArray[np.float64], 
                    size: int | Tuple[int, ...] | None = None, 
                    mu: int | Tuple[int, ...] | None = None,
                    rng: Generator | None=None
                    ) -> Union[np.float64,NDArray[np.float64]]:
    """
    Sample diameters from a cumulative size-frequency distribution (SFD).
    
    Given an array of diameters and optionally a cumulative distribution function (CDF), this function generates new diameter values that follow the specified distribution. The SFD is treated as a continuous function that interpolates between the provided diameters, which are assumed to represent a power-law distribution.
    
    Parameters
    ----------
    diameters : array_like
        An array of diameters from which the SFD is constructed. Must be 1-dimensional.
    cdf : array_like
        The cumulative distribution function corresponding to `diameters`. Must be the same size as `diameters` and must be monotonically increasing with decreasing diameter value.
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None and mu is None then a single value is returned. Note: mu and size are mutually exclusive. 
    mu : int or tuple of ints, optional
        The expected number of samples to generate using a Poisson random number genertor. If the shape is (m, n, k), then m * n * k samples are drawn. Note: mu and size are mutually exclusive. 
    rng : numpy.random.Generator, optional
        An instance of a random number generator compatible with numpy's random generators. If not provided, `default_rng` is used to create a new instance.
    
    Returns
    -------
    np.float64 or ndarray of np.float 64 
        A scalar or array of sampled diameter values from the SFD. 
    
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
        
    # Check that the shapes and sizes of diameters and cdf are compatible
    if np.isscalar(diameters) or np.isscalar(cdf):
        raise ValueError("The 'diameters' and 'cdf' arguments must be arrays")
    if diameters.ndim != 1:
        raise ValueError("The 'diameters' argument must be a 1-dimensional array")
    if diameters.shape != cdf.shape:
        raise ValueError("The 'diameters' and 'cdf' arguments must have the same shape")
    if diameters.size != cdf.size:
        raise ValueError("The 'diameters' and 'cdf' arguments must have the same size")
    if diameters.size < 2:
        raise ValueError("The 'diameters' and 'cdf' arguments must have at least two elements")
    if np.any(diameters <= 0.0):
        raise ValueError("All values in the 'diameters' argument must be positive")
    if np.any(cdf <= 0.0):
        raise ValueError("All values in the 'cdf' argument must be positive")
    if size is None and mu is None:
        size = 1
    elif size is not None and mu is not None:
        raise ValueError("The 'size' and 'mu' arguments are mutually exclusive")
    elif size is None and mu is not None:
        size = rng.poisson(mu) 

    # Sort the diameters in descending order and get the cumulative distribution if it was not supplied
    sorted_indices = np.argsort(diameters)[::-1]
    sorted_diameters = diameters[sorted_indices]
    sorted_cdf = cdf[sorted_indices]
    
    # Check to make sure that the CDF is correctly specified so that as diameter is decreasing it is monotonically increasing
    is_monotonic_increasing = np.all(np.diff(sorted_cdf) >= 0)
    if not is_monotonic_increasing:
        raise ValueError("The CDF must be monotonically increasing with decreasing diameter.")
    
    # Normalize the cdf and put it in logspace
    sorted_cdf /= sorted_cdf[-1]
    log_sorted_cdf = np.log(sorted_cdf)

    # Generate uniform random numbers for the entire sample size
    u = rng.uniform(low=sorted_cdf[0], high=sorted_cdf[-1], size=size)
    u = np.log(u)

    # Handle the situation where u is a scalar (size is None and diameters is a scalar)
    if np.isscalar(u):
        u = np.array([u])
        
    # Flatten u to work with it as a 1D array
    original_shape = u.shape
    u = u.flatten()     

    # Find the indices where the random numbers would be inserted to maintain order
    # Use the right side of the interval to avoid edge effects for when u == sorted_cdf[0]
    # Because rng.uniform returns values in the half-open interval [sorted_cdf[0], sorted_cdf[-1]), u will never be exactly equal to sorted_cdf[-1]
    indices = np.searchsorted(log_sorted_cdf, u, side="right")

    # Initialize the new_diameters array
    new_diameters = np.empty(u.shape)

    log_diam_low = np.log(sorted_diameters[indices-1])
    log_diam_high = np.log(sorted_diameters[indices])
    log_nval_low = log_sorted_cdf[indices-1]
    log_nval_high = log_sorted_cdf[indices]
                             
    fractions = (u - log_nval_low) / (log_nval_high - log_nval_low)
    log_diam_interp = log_diam_low + fractions * (log_diam_high - log_diam_low)
    new_diameters = np.exp(log_diam_interp)
    
    # Reshape new_diameters to the original shape of u
    new_diameters = new_diameters.reshape(original_shape)

    # Add a small random noise to the diameters
    noise = 1e-8 * rng.uniform(size=new_diameters.shape)
    new_diameters *= (1 + noise)
    if size == 1:
        return new_diameters[0]
    else:
        return new_diameters


def get_random_velocity(
                        vmean: np.float64, 
                        size: int | Tuple[int, ...]=1, 
                        rng: Generator | None=None
                        ) -> Union[np.float64,NDArray[np.float64]]:
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
    
           
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    def plot_random_location():
        # Sample data generation
        size = 10000
        points = get_random_location(size=size)
        
        lons = points['lon']
        lats = points['lat']

        # Number of bins
        bins = 50 

        # Longitude plot
        observed_counts_lon, bins_lon_deg = np.histogram(lons, bins=bins, range=(-180, 180.0))
        expected_count_lon = size // bins

        # Latitude plot
        observed_counts_lat, bins_lat_deg = np.histogram(lats, bins=bins, range=(-90, 90))
        # Convert bins to degrees
        bins_lon = np.deg2rad(bins_lon_deg)
        bins_lat = np.deg2rad(bins_lat_deg)
        
        # For expected counts in latitude, adjust for area covered by each bin
        area_ratio = np.sin(bins_lat[1:]) - np.sin(bins_lat[:-1])
        total_area = np.sin(np.pi/2) - np.sin(-np.pi/2)  # Total area for the entire sphere
        expected_count_lat = size * area_ratio / total_area

        # Bar width in degrees
        bar_width_lon = np.diff(bins_lon_deg)
        bar_width_lat = np.diff(bins_lat_deg)

        # Plotting
        fig, axs = plt.subplots(2, 1, figsize=(8, 6))

        # Longitude plot
        axs[0].bar(bins_lon_deg[:-1], observed_counts_lon, width=bar_width_lon, align='edge', label='Observed')
        axs[0].plot(bins_lon_deg[:-1], [expected_count_lon] * len(bins_lon_deg[:-1]), color='red', label='Expected')
        axs[0].set_xlabel('Longitude (deg)')
        axs[0].set_ylabel('Number')
        axs[0].legend()
        axs[0].set_title('Number vs Longitude')
        
        # Latitude
        axs[1].bar(bins_lat_deg[:-1], observed_counts_lat, width=bar_width_lat, align='edge', alpha=0.5, label='Observed')
        axs[1].plot(bins_lat_deg[:-1], expected_count_lat, label='Expected', color='red')
        axs[1].set_xlabel('Latitude (deg)')
        axs[1].set_ylabel('Number')
        axs[1].legend()
        axs[1].set_title('Number vs Latitude')

        plt.tight_layout()
        plt.show()


    def plot_random_angle():
        # Sample data generation
        size = 10000
        angles = get_random_impact_angle(size=size)
        

        # Number of bins
        bins = 50 
        observed_counts, bins_ang = np.histogram(angles, bins=bins, range=(0.0, 90.0))

        # Calculate expected distribution
        uniform_dist = np.linspace(0, 1, size)
        transformed_angles = np.rad2deg(np.arcsin(np.sqrt(uniform_dist)))
        expected_counts, _ = np.histogram(transformed_angles, bins=bins, range=(0.0, 90.0))

        # Plotting
        fig, ax = plt.subplots(figsize=(8, 4))

        # Observed counts
        ax.bar(bins_ang[:-1], observed_counts, width=np.diff(bins_ang), align='edge', label='Observed', alpha=0.5)

        # Expected counts
        ax.plot(bins_ang[:-1], expected_counts, label='Expected', color='red')

        ax.set_xlabel('Impact Angle (deg)')
        ax.set_ylabel('Count')
        ax.legend()
        ax.set_title('Impact Angle Distribution')

        plt.show()


    def plot_random_velocity():
        vmean = 20e3
        size = 10000
        velocities = get_random_velocity(vmean, size)
        
        # Create histogram of the generated velocities
        bins = np.linspace(0, vmean*3, 50)
        histogram, bins = np.histogram(velocities, bins=bins, density=True)
        bin_centers = 0.5 * (bins[1:] + bins[:-1])

        # Generate theoretical Maxwell-Boltzmann distribution for the given vmean
        theoretical_dist = maxwell.pdf(bin_centers, scale=vmean/np.sqrt(8/np.pi))

        # Plotting
        plt.figure(figsize=(8, 4))
        plt.bar(bin_centers, histogram, width=bins[1]-bins[0], alpha=0.5, label='Generated Velocities')
        plt.plot(bin_centers, theoretical_dist, label='Theoretical Maxwell-Boltzmann', color='red')
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Probability Density')
        plt.title('Maxwell-Boltzmann Distribution of Velocities')
        plt.legend()
        plt.show()

           
    def plot_random_size():
        num_realizations = 100
        nbins = 10
        size = 100000
        Dhi = 100.0
        p = 2.0
        C =  Dhi**p 
        Dlo = (size/C)**(-1.0/p) 
        diameters = np.exp(np.linspace(np.log(Dlo), np.log(100*Dhi), nbins))
        # Test a simple power law SFD
        cdf = C * diameters**-p
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.set_title("Test SFD generation")
        ax.set_xlabel("$D$")
        ax.set_ylabel("$N_{>D}$")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim([1.0,1.1*size])

        for i in range(num_realizations):
            new_diameters = get_random_size(diameters, cdf, mu=size)
            new_diameters = np.sort(new_diameters)[::-1]
            nval = np.arange(1, new_diameters.size+1)
            if i == 0:
                label = 'Sampled SFD'
            else:
                label = None
                
            ax.plot(new_diameters, nval, label=label, alpha=0.25, color='orange')
            
        ax.plot(diameters, cdf, label='Model SFD')
        ax.legend(loc='upper right')
        plt.show()  
        
        
    plot_random_location()
    plot_random_angle()
    plot_random_velocity()
    plot_random_size()          