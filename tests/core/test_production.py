import unittest
from cratermaker import NeukumProduction, Production
import numpy as np

class TestProduction(unittest.TestCase):
    
    def test_production_N_to_time(self):

        plaw = Production()
        neukum = NeukumProduction(model="Moon")

        for age_orig in np.linspace(0,4500,num=10):
            D = np.logspace(-1,3,num=1000)
            N = neukum.function(diameter=D,age=age_orig)
            age_new = neukum.function_inverse(cumulative_number=N,diameter=D)
            np.testing.assert_array_almost_equal(age_orig,age_new,decimal=2)
            
            N = plaw.function(diameter=D,age=age_orig)
            age_new = plaw.function_inverse(cumulative_number=N,diameter=D)
            np.testing.assert_array_almost_equal(age_orig,age_new,decimal=2)
            
    def test_sample_arguments(self):
        production = [Production(), NeukumProduction()]
        for prod in production:

            # Test valid arguments
            try:
                prod.sample(age=500, diameter_range=(10, 100), area=1000)
            except ValueError:
                self.fail("sample() raised ValueError unexpectedly with valid arguments!")

            # Test providing both age and cumulative_number_at_diameter
            with self.assertRaises(ValueError):
                prod.sample(age=500, cumulative_number_at_diameter=(100, 10), diameter_range=(10, 100))
               
            # Test providing both reference_age and reference_cumulative_number_at_diameter
            with self.assertRaises(ValueError):
                prod.sample(age=500, reference_age=1000, reference_cumulative_number_at_diameter=(100, 10), diameter_range=(10, 100))

            # Test missing diameter_range
            with self.assertRaises(ValueError):
                prod.sample(age=500)

            # Test diameter_range with invalid length
            with self.assertRaises(ValueError):
                prod.sample(age=500, diameter_range=(10,), area=1000)

            # Test diameter_range with non-positive minimum diameter
            with self.assertRaises(ValueError):
                prod.sample(age=500, diameter_range=(-1, 100), area=1000)

            # Test diameter_range with max diameter less than min diameter
            with self.assertRaises(ValueError):
                prod.sample(age=500, diameter_range=(100, 10), area=1000)

            # Test invalid cumulative_number_at_diameter length
            with self.assertRaises(ValueError):
                prod.sample(cumulative_number_at_diameter=(100,), diameter_range=(10, 100))

            # Test negative value in cumulative_number_at_diameter
            with self.assertRaises(ValueError):
                prod.sample(cumulative_number_at_diameter=(-100, 10), diameter_range=(10, 100))

            # Test negative area
            with self.assertRaises(ValueError):
                prod.sample(age=500, diameter_range=(10, 100), area=-1000)
                
            # Test non-scalar age
            with self.assertRaises(ValueError):
                prod.sample(age=[500, 1000], diameter_range=(10, 100))
            
            # This checks if a ValueError is raised for an invalid length of reference_cumulative_number_at_diameter.    
            with self.assertRaises(ValueError):
                prod.sample(reference_cumulative_number_at_diameter=(100,), diameter_range=(10, 100))
           
            #Testing when a negative value is provided in reference_cumulative_number_at_diameter. 
            with self.assertRaises(ValueError):
                prod.sample(reference_cumulative_number_at_diameter=(-100, 10), diameter_range=(10, 100))

        
if __name__ == '__main__':
    unittest.main()