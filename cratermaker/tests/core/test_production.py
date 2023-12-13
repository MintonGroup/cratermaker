import unittest
import cratermaker
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
         
        
if __name__ == '__main__':
    unittest.main()