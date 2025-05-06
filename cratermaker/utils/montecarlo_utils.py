import numpy as np
from numpy.random import Generator
from numba import njit
from numpy.linalg import norm
from numpy import cross
from typing import Any
from numpy.typing import NDArray
from scipy.stats import truncnorm
from uxarray import Grid
from cratermaker.core.base import _rng_init
from cratermaker.utils.custom_types import FloatLike

def get_random_location(size: int=1, 
                        rng: Generator | None = None, 
                        rng_seed: int | None = None,
                        rng_state: dict | None = None,
                        **kwargs: Any) -> NDArray[np.float64]:
    """
    Computes random longitude and latitude values.
    
    Generates a set of latitude and longitude values that are uniformly distributed on the surface of a sphere.
    
    Parameters
    ----------
    size : int or tuple of ints, optional
        The number of samples to generate. If size is None (the default), a single tuple is returned. If size is greater than 1, 
        then a structured array with fields 'lon' and 'lat' is returned.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments.
    Returns
    -------
    A structured numpy array with the location data in the format [('lon', 'f8'), ('lat', 'f8')].
    """

    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    u = rng.uniform(size=size)
    v = rng.uniform(size=size)
    
    # Compute the angles theta and phi
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    
    # Convert to lon/lat
    lon = np.rad2deg(theta - np.pi) # Use the convention that longitude is in the range [-180, 180]
    lat = np.rad2deg(phi - np.pi / 2.0)
    
    lonlat_arr = np.empty(size, dtype=[('lon', 'f8'), ('lat', 'f8')])

    # Reshape lat and lon to the original size if necessary
    lon = lon.reshape(size)
    lat = lat.reshape(size)
  
    # Combine lat and lon into a structured array
    lonlat_arr['lon'] = lon
    lonlat_arr['lat'] = lat
    
    return lonlat_arr


@njit
def _get_one_random_location(face_nodes, node_x, node_y, node_z, rng_vals):
    valid_nodes = face_nodes[face_nodes > 0]
    n = len(valid_nodes)
    tris = [(0, j, j+1) for j in range(1, n - 1)]

    n_valid = len(valid_nodes)
    vertices = np.empty((n_valid, 3), dtype=np.float64)
    for i in range(n_valid):
        idx = valid_nodes[i]
        vertices[i, 0] = node_x[idx]
        vertices[i, 1] = node_y[idx]
        vertices[i, 2] = node_z[idx]

    areas = np.empty(len(tris))
    for i in range(len(tris)):
        j0, j1, j2 = tris[i]
        v0, v1, v2 = vertices[j0], vertices[j1], vertices[j2]
        areas[i] = 0.5 * norm(cross(v1 - v0, v2 - v0))
    areas /= areas.sum()
    cum_areas = np.cumsum(areas)

    # triangle selection
    tri_idx = np.searchsorted(cum_areas, rng_vals[0])
    j0, j1, j2 = tris[tri_idx]
    v0, v1, v2 = vertices[j0], vertices[j1], vertices[j2]

    r1, r2 = rng_vals[1], rng_vals[2]
    if r1 + r2 > 1.0:
        r1, r2 = 1.0 - r1, 1.0 - r2
    r0 = 1.0 - r1 - r2
    return r0 * v0 + r1 * v1 + r2 * v2


def get_random_location_on_face(grid: Grid, 
                                face_index: int | NDArray[np.int64], 
                                rng: Generator | None = None, 
                                rng_seed: int | None = None,
                                rng_state: dict | None = None,
                                **kwargs: Any) -> NDArray[np.float64]:
    """
    Generate a random coordinate within a given face of an unstructured mesh.

    Parameters
    ----------
    grid : uxarray.Grid
        The grid object containing the mesh information.
    face_index : int | NDArray[np.int64]
        The index or array of indices of the face within the grid to obtain the random sample.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 

    Returns
    -------
    A structured numpy array with the location data in the format of the same shape as face_index [('lon', 'f8'), ('lat', 'f8')].
    """
    from uxarray.grid import coordinates
    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    # Extract node indices for the given face(s)
    node_indices = grid.face_node_connectivity[face_index, :]
    face_index = np.atleast_1d(face_index)
    size = len(face_index)

    # Prepare index arrays
    locations = np.empty(size, dtype=[('lon', 'f8'), ('lat', 'f8')])

    node_x = grid.node_x.values
    node_y = grid.node_y.values
    node_z = grid.node_z.values
    node_indices = node_indices.values
    rng_vals = rng.random((size, 3))
    for i in range(size):
        p = _get_one_random_location(node_indices[i], node_x, node_y, node_z, rng_vals=rng_vals[i,:])
        lon, lat = coordinates._xyz_to_lonlat_deg(*p)
        locations['lon'][i] = np.float64(lon)
        locations['lat'][i] = np.float64(lat)

    return locations


def get_random_impact_angle(size: int | tuple[int, ...]=1, 
                            rng: Generator | None = None, 
                            rng_seed: int | None = None,
                            rng_state: dict | None = None,
                            **kwargs: Any) -> NDArray[np.float64]:
    """
    Sample impact angles from a distribution centered on 45deg.
    
    For the theory, see Shoemaker (1962) "Interpretation of lunar craters."
    
    Parameters 
    ----------
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 

    Returns
    ----------
    ndarray of np.float64 
        An array of impact angles (in degrees).
    """    
    
    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)

    u = np.sqrt(rng.uniform(size=size))
    impact_angle = np.arcsin(u)
    return np.rad2deg(impact_angle)


def get_random_impact_direction(size: int | tuple[int, ...]=1,
                                rng: Generator | None = None, 
                                rng_seed: int | None = None,
                                rng_state: dict | None = None,
                                **kwargs: Any) -> NDArray[np.float64]:
    """
    Sample impact direction from a uniform distribution.

    Parameters
    ----------
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single scalar value is returned.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 

    Returns
    -------
    ndarray of np.float64 
        An array of impact angles (in degrees).
    """
    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
    
    pdir = rng.uniform(0.0, 360.0, size=size)
    return pdir


def get_random_size(diameters: NDArray[np.float64], 
                    cdf: NDArray[np.float64], 
                    size: int | tuple[int, ...] | None = None, 
                    mu: int | tuple[int, ...] | None = None,
                    rng: Generator | None = None, 
                    rng_seed: int | None = None,
                    rng_state: dict | None = None,
                    **kwargs: Any) -> NDArray[np.float64]:
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
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 
    
    Returns
    -------
    ndarray of np.float64 
        An array of sampled diameter values from the SFD. 
    
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
    """

    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
        
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
    return new_diameters


def get_random_velocity(vmean: np.float64, 
                        size: int | tuple[int, ...]=1, 
                        rng: Generator | None = None, 
                        rng_seed: int | None = None,
                        rng_state: dict | None = None,
                        **kwargs: Any) -> NDArray[np.float64]:
    """
    Sample impact velocities from a Maxwell-Boltzmann distribution given a mean velocity.
    
    Parameters 
    ----------
    vmean : np.float64
        The mean velocity of the distribution.
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
        Additional keyword arguments. 

    Returns
    ----------
    ndarray 
        An array of impact velocities (in m/s).
    """
    
    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
    
    sigma = vmean / np.sqrt(8/np.pi)
    
    vx = rng.normal(0, sigma, size=size)
    vy = rng.normal(0, sigma, size=size)
    vz = rng.normal(0, sigma, size=size)
    velocities = np.sqrt(vx**2 + vy**2 + vz**2)
 
    return velocities


def bounded_norm(mean: FloatLike,
                 scale: FloatLike,
                 size: int | tuple[int, ...]=1,
                rng: Generator | None = None, 
                rng_seed: int | None = None,
                rng_state: dict | None = None,
                 **kwargs: Any) -> FloatLike: 
    """
    Sample from a truncated normal distribution that is bounded by 1-sigma stdev
    
    Parameters 
    ----------
    loc : FloatLike
       mean of the distribution
    scale : FloatLike
       standard deviation and bounds of the distribution
    size : int or tuple of ints, optional
        The number of samples to generate. If the shape is (m, n, k), then m * n * k samples are drawn. If size is None (the default), a single value is returned if `diameters` is a scalar, otherwise an array of samples is returned with the same size as `diameters`.
    rng : numpy.random.Generator | None
        A numpy random number generator. If None, a new generator is created using the rng_seed if it is provided.
    rng_seed : Any type allowed by the rng_seed argument of numpy.random.Generator, optional
        The rng_rng_seed for the RNG. If None, a new RNG is created.
    rng_state : dict, optional
        The state of the random number generator. If None, a new state is created.
    **kwargs : Any
       
    Returns
    ----------
    float
       Truncated norm bounded by loc-scale, loc+scale
    """    
    rng, _ = _rng_init(rng=rng, rng_seed=rng_seed, rng_state=rng_state, **kwargs)
    
    lower_bound = mean - scale
    upper_bound = mean + scale
    truncated_normal = truncnorm((lower_bound - mean) / scale, 
                                 (upper_bound - mean) / scale, 
                                 loc=mean, 
                                 scale=scale)
    
    return truncated_normal.rvs(size, random_state=rng)
    
           
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