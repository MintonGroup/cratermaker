import unittest
import cratermaker
from cratermaker import NeukumProduction, Production
import numpy as np

class TestProduction(unittest.TestCase):
    
    def test_production_N_to_time(self):

        plaw = Production()
        neukum = NeukumProduction(model="Moon")

        for t in np.linspace(0,4500,num=10):
            D = np.logspace(-1,3,num=1000)
            N = neukum.function(diameter=D,time=t)
            tnew = neukum.function_inverse(cumulative_number=N,diameter=D)
            np.testing.assert_array_almost_equal(t,tnew,decimal=2)
            
            N = plaw.function(diameter=D,time=t)
            tnew = plaw.function_inverse(cumulative_number=N,diameter=D)
            np.testing.assert_array_almost_equal(t,tnew,decimal=2)
         
        
if __name__ == '__main__':
    unittest.main()