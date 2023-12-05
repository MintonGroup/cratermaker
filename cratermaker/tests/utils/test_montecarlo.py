import unittest
import numpy as np
from cratermaker.utils.montecarlo import *
from scipy.stats import chisquare, ks_2samp, poisson

class TestMonteCarlo(unittest.TestCase):

    def test_get_random_location(self):
        """Test get_random_location with default parameters."""
        size = 10000
        points = get_random_location(size=size)

        lons = points['lon']
        lats = points['lat']
        self.assertIsInstance(lons[0], np.float64)
        self.assertIsInstance(lats[0], np.float64)
        self.assertTrue(np.all(0.0 <= lons) and np.all(lons <= 2 * np.pi))
        self.assertTrue(np.all(-np.pi/2 <= lats) and np.all(lats <= np.pi/2))
        
        # Test Longitude Uniformity
        
        # Set a significance level for your test, commonly 0.05
        alpha = 0.05
       
        bins = 100 
        observed_counts, bins_lon = np.histogram(lons, bins=bins, range=(0.0, 2 * np.pi))
        expected_count_lon = size // bins
        
        # Perform the chi-square test
        chi2_statistic, p_value = chisquare(f_obs=observed_counts, f_exp=[expected_count_lon]*bins)

        # Assert that the p-value is greater than the significance level
        self.assertTrue(p_value > alpha)        
        
        sin_lats = np.sin(lats)
        observed_counts, bins_lat = np.histogram(sin_lats, bins=bins, range=(-1, 1))
        expected_count_lat = size // bins
        
        # Perform the chi-square test
        chi2_statistic, p_value = chisquare(f_obs=observed_counts, f_exp=[expected_count_lat]*bins)

        # Assert that the p-value is greater than the significance level
        self.assertTrue(p_value > alpha)        
        
        return
    

    def test_get_random_impact_angle(self):
        """Test get_random_impact_angle function."""
        size = 10000
        angles = get_random_impact_angle(size=size)

        # Check type and shape
        self.assertIsInstance(angles, np.ndarray)
        self.assertEqual(angles.shape, (size,))

        # Check range (0 to np.pi/2 radians)
        self.assertTrue(np.all(angles >= 0) and np.all(angles <= 90.0))

        # Check distribution
        # Generate expected distribution centered around 45 degrees
        expected_angles = np.arcsin(np.sqrt(np.linspace(0, 1, size)))
        _, p_value = ks_2samp(angles, np.rad2deg(expected_angles))

        # Check if the p-value is greater than a chosen significance level (0.05)
        self.assertGreater(p_value, 0.05)

    def test_get_random_size(self):
        # Generate power law CDF and test that it results in a Poisson-like distribution of binned counts
        diameters = np.exp(np.linspace(np.log(1.0), np.log(1.0e2), 10))
        cdf = diameters**-3
        cdf /= cdf[-1]  # Normalize CDF

        # Generate random sizes
        size = 100000
        random_sizes = get_random_size(diameters, cdf, size=size)

        # Geometric binning
        bin_edges = np.geomspace(1.0, 1e5, 20)
        counts, _ = np.histogram(random_sizes, bins=bin_edges)

        # Expected counts based on Poisson distribution
        # Compute cumulative counts at bin edges
        cumulative_counts_at_edges = bin_edges**-3

        # Calculate expected counts in each bin as the difference in cumulative counts
        lambda_ = np.diff(cumulative_counts_at_edges)

        # Normalize lambda_ so that sum(lambda_) equals size
        lambda_ *= size / lambda_.sum()
        expected_counts = [poisson.pmf(k, lambda_[k]) * size for k in range(len(counts))]

        # Perform chi-square test
        chi2_statistic, p_value = chisquare(f_obs=counts, f_exp=expected_counts)

        # Check if p-value is significant
        self.assertGreater(p_value, 0.05)  # Common threshold for significance
        

if __name__ == '__main__':
    unittest.main()
