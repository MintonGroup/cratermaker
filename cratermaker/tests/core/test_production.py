import unittest
import cratermaker
from cratermaker import NeukumProductionFunction
import numpy as np

class TestNPF(unittest.TestCase):
    
    def test_npf_N_to_time(self):
        npf = NeukumProductionFunction()
        time_input = np.logspace(-4,np.log10(4.4),num=100)
        time_output = npf._Nrel_to_time(npf._time_to_Nrel(time_input))
        np.testing.assert_array_almost_equal(time_input, time_output, decimal=6)
        
        bad_N = np.array([-1.0,1.0,npf._time_to_Nrel(4.6,check_valid_time=False)])
        expected_time = np.array([-1.0,1.0,np.nan])
        time_output = npf._Nrel_to_time(bad_N)
        np.testing.assert_array_almost_equal(expected_time, time_output, decimal=1)
        
    def test_arr_vs_scalar_args(self):
        npf = NeukumProductionFunction()
        time_array = np.logspace(-4,np.log10(4.4))
        time_scalar = 2.0
        
        N1_scalar = npf._N1(time_scalar)
        N1_array = npf._N1(time_array)
        self.assertTrue(np.isscalar(N1_scalar)) 
        self.assertTrue(N1_array.shape == time_array.shape)
        
        D_scalar = 50.0
        D_array = np.logspace(0,2)
        N_scalar = npf._CSFD(D_scalar)
        N_array = npf._CSFD(D_array)
        self.assertTrue(np.isscalar(N_scalar))
        self.assertTrue(N_array.shape == D_array.shape)
        
        N_scalar = npf._dNdD(D_scalar)
        N_array = npf._dNdD(D_array)
        self.assertTrue(np.isscalar(N_scalar))
        self.assertTrue(N_array.shape == D_array.shape)
        
        N_scalar = npf._DSFD(D_scalar)
        N_array = npf._DSFD(D_array)
        self.assertTrue(np.isscalar(N_scalar))
        self.assertTrue(N_array.shape == D_array.shape) 
        
        N_scalar = npf._time_to_Nrel(time_scalar)
        N_array = npf._time_to_Nrel(time_array)
        self.assertTrue(np.isscalar(N_scalar)) 
        self.assertTrue(N_array.shape == time_array.shape)
        
         
        
if __name__ == '__main__':
    unittest.main()