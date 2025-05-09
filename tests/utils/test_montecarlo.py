import unittest
import numpy as np
from cratermaker.utils.montecarlo_utils import *
from scipy.stats import chisquare, ks_2samp

class TestMonteCarlo(unittest.TestCase):

    def test_get_random_location(self):
        """Test get_random_location with default parameters."""
        # Do three trials. If all three fail, the test fails.
        ntrials = 3
        trial_results = []
        for _ in range(ntrials):
            size = 1000
            points = get_random_location(size=size)

            lons = points['lon']
            lats = points['lat']
            self.assertIsInstance(lons[0], np.float64)
            self.assertIsInstance(lats[0], np.float64)
            self.assertTrue(np.all(-180.0 <= lons) and np.all(lons <= 180.0), msg=f"Longitude range: {np.min(lons)} to {np.max(lons)}")
            self.assertTrue(np.all(-90.0 <= lats) and np.all(lats <= 90.0), msg=f"Latitude range: {np.min(lats)} to {np.max(lats)}")
            
            # Test Longitude Uniformity
            
            # Set a significance level for your test, commonly 0.05
            alpha = 0.05
        
            bins = 20 
            observed_counts, bins_lon = np.histogram(lons, bins=bins, range=(-180, 180))
            expected_count_lon = size // bins
            
            # Perform the chi-square test
            _, p_value = chisquare(f_obs=observed_counts, f_exp=[expected_count_lon]*bins)

            # Assert that the p-value is greater than the significance level
            result = p_value
            
            sin_lats = np.sin(np.deg2rad(lats))
            observed_counts, _ = np.histogram(sin_lats, bins=bins, range=(-1, 1))
            expected_count_lat = size // bins
            
            # Perform the chi-square test
            _, p_value = chisquare(f_obs=observed_counts, f_exp=[expected_count_lat]*bins)

            # Assert that the p-value is greater than the significance level
            result = result and p_value > alpha
            trial_results.append(result)
        self.assertTrue(any(trial_results)) 
        return
    

    def test_get_random_impact_angle(self):
        """Test get_random_impact_angle function."""
        size = 10000
        angles = get_random_impact_angle(size=size)

        # Check type and morphology
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
        nbins = 10
        size = 10000
        num_realizations = 1000
        
        Dhi = 100.0
        p = 3.0
        C =  Dhi**p 
        Dlo = (size/C)**(-1.0/p) 
        diameters = np.exp(np.linspace(np.log(Dlo), np.log(100*Dhi), nbins))
        model_cdf = C*diameters**-p

        # Generate multiple populations of sizes and bin the results from each population
        all_counts = []
        for _ in range(num_realizations):
            random_sizes = get_random_size(diameters, model_cdf, mu=size)
            counts, _ = np.histogram(random_sizes, bins=diameters)
            cumulative_counts = np.cumsum(counts[::-1])[::-1]
            all_counts.append(cumulative_counts)
       
        # Assuming 'all_counts' is a list of arrays, each array being the counts for one realization
        observed_counts = np.stack(all_counts)  # Stack all count arrays

        # Calculate the mean and standard deviation of observed counts for each bin
        observed_means = np.mean(observed_counts, axis=0)
        observed_std_devs = np.std(observed_counts, axis=0)

        # Calculate the expected standard deviation for each bin from the Poisson distribution
        expected_counts = model_cdf 
        expected_std_devs = np.sqrt(expected_counts)

        # Compare observed and expected standard deviations
        for i in range(len(observed_means)):
            if observed_means[i] == 0:
                continue
            self.assertAlmostEqual(observed_means[i], expected_counts[i], delta=4*expected_std_devs[i])
            
        # test with invalid diameter shape or size
        diameters = np.array([[100., 56.], [32., 18.]])  # 2D array
        cdf = np.array([1., 0.8, 0.6, 0.4])
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf)
        diameters = np.array([])  # Empty array
        cdf = np.array([])
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf) 
        diameters = np.array([1.0])  # Single element array
        cdf = np.array([1.0])
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf)
        diameters = 1.0  # Scalar
        cdf = 1.0
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf)                          
        diameters = np.array([100., 56., 32., 18.])
        cdf = np.array([0.4, 0.6, 0.8])
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf)
        
        # test the case when the cdf is not monotonically decreasing with increasing diameter 
        diameters = np.array([100., 56., 32., 18.])
        cdf = np.array([1., 0.8, 0.6, 0.4])
        with self.assertRaises(ValueError):
            get_random_size(diameters, cdf)
            
        # test the case when both size and mu are requested
        diameters = np.logspace(0,3,nbins)
        model_cdf = diameters**-3
        with self.assertRaises(ValueError):
            get_random_size(diameters, model_cdf, size=100, mu=100)
         
        # test the case when negative values are requested
        diameters = np.logspace(0,3,nbins)
        model_cdf = diameters**-3
        with self.assertRaises(ValueError):
            get_random_size(diameters, model_cdf, size=-100)
        with self.assertRaises(ValueError):
            get_random_size(diameters, model_cdf, mu=-100)
        with self.assertRaises(ValueError):
            get_random_size(diameters, -model_cdf)
        with self.assertRaises(ValueError):
            get_random_size(-diameters, model_cdf)

        return
        

if __name__ == '__main__':
    unittest.main()
